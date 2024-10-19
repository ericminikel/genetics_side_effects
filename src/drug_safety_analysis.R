##############
# STARTUP
##############

overall_start_time = Sys.time()
tell_user = function(x) { cat(file=stderr(), x); flush.console() }

tell_user('Loading dependencies...')

options(stringsAsFactors=F)
if(interactive()) setwd('~/src/genetics_side_effects')
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(binom))
suppressMessages(library(weights))
suppressMessages(library(openxlsx))
suppressMessages(library(Hmisc))
suppressMessages(library(lawstat))
suppressMessages(library(MASS)); select = dplyr::select; summarize = dplyr::summarize


##############
# OUTPUT STREAMS
##############

tell_user('done.\nCreating output streams...')

text_stats_path = 'display_items/stats_for_text.txt'
write(paste('Last updated: ',Sys.Date(),'\n',sep=''),text_stats_path,append=F) # start anew - but all subsequent writings will be append=T
write_stats = function(...) {
  write(paste(list(...),collapse='',sep=''),text_stats_path,append=T)
  write('\n',text_stats_path,append=T)
}

supplement_path = 'display_items/supplement.xlsx'
supplement = createWorkbook()
# options("openxlsx.numFmt" = "0.00") # this looks better for residuals but terrible for p values and weeks post-dose
supplement_directory = tibble(name=character(0), title=character(0))
write_supp_table = function(tbl, title='') {
  # write Excel sheet for supplement
  table_number = length(names(supplement)) + 1
  table_name = paste0('s',formatC(table_number,'d',digits=0,width=2,flag='0'))
  addWorksheet(supplement,table_name)
  bold_style = createStyle(textDecoration = "Bold")
  writeData(supplement,table_name,tbl,headerStyle=bold_style,withFilter=T)
  freezePane(supplement,table_name,firstRow=T)
  saveWorkbook(supplement,supplement_path,overwrite = TRUE)
  
  # also write tab-sep version for GitHub repo
  write_tsv(tbl,paste0('display_items/table-',table_name,'.tsv'), na='')
  
  # and save the title in the directory tibble for later
  assign('supplement_directory',
         supplement_directory %>% add_row(name=table_name, title=title),
         envir = .GlobalEnv)
}




##############
# FUNCTIONS & CONSTANTS
##############

tell_user('done.\nSetting constants and functions...')

whisker_factor = 1.5

percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x <= 0, '', ifelse(signed, '+', '')),formatC(100*x,format='f',digits=digits),'%'))

upper = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) + sds*sd(x)/sqrt(sum(!is.na(x)))
}
lower = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) - sds*sd(x)/sqrt(sum(!is.na(x)))
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot


clipcopy = function(tbl) {
  clip = pipe("pbcopy", "w")  
  write.table(tbl, file=clip, sep = '\t', quote=F, row.names = F, na='')
  close(clip)
}

clipdist = function(x, minx, maxx) {
  return (pmin(maxx,pmax(minx,x)))
}

stderrmsg = function(msg, verbose=T) {
  if (verbose) {
    cat(file=stderr(),msg)
    flush.console()
  }
}

rbind_files = function(path, grepstring) {
  all_files = list.files(path, full.names=T)
  these_files = all_files[grepl(grepstring,all_files)]
  if (exists('rbound_table')) rm('rbound_table')
  for (this_file in these_files) {
    this_tbl = read_delim(this_file, col_types=cols()) %>% clean_names()
    this_tbl$file = gsub('.*\\/','',gsub('\\.[tc]sv','',this_file))
    if (exists('rbound_table')) {
      rbound_table = rbind(rbound_table, this_tbl)
    } else {
      rbound_table = this_tbl
    }
  }
  return (rbound_table)
}

make_dse = function(drugtable, setable, drugsematch, simmat, assocs, se_metadata, verbose=T, skip_enrich=F, presimindic=NULL) {
  
  ####
  # 0. prep the tables
  drugtable$drug_name = tolower(drugtable$drug_name)
  drugsematch$drug_name = tolower(drugsematch$drug_name)
  drugsematch$dse_uid = paste0(drugsematch$drug_name,'-',drugsematch$mesh_id)
  
  ####
  # 1. get pharmaprojects info - targets and approved indications
  if (verbose) {
    cat(file=stderr(),paste0('1. Mapping Pharmaprojects targets and indications.\n'))
    flush.console()
  }
  if (verbose) {
    cat(file=stderr(),paste0('- Mapping drug targets.\n'))
    flush.console()
  }
  drugtable %>%
    left_join(pp_targets, by='ppid') %>%
    distinct(drug_name, gene) %>%
    filter(!is.na(gene) & !is.na(drug_name)) -> drug_target_match
  
  
  ####
  # 2. subset to known targets and known approved indications
  if (verbose) {
    cat(file=stderr(),paste0('2. Subsetting to known targets, approved indication, available similaritys.\n'))
    flush.console()
  }
  drugtable %>%
    inner_join(pp_launched_indic, by='ppid') %>%
    filter(!is.na(indication_mesh_id) & indication_mesh_id!='') %>%
    distinct(drug_name, indication_mesh_id) -> approved_indics
  orig_nrow = nrow(drugtable)
  drugtable %>%
    filter(drug_name %in% drug_target_match$drug_name) %>%
    filter(drug_name %in% approved_indics$drug_name) -> drugtable
  new_nrow = nrow(drugtable)
  if (verbose) {
    cat(file=stderr(),paste0('- ',new_nrow,'/',orig_nrow,' drugs remain.\n'))
    flush.console()
  }
  orig_nrow = nrow(setable)
  setable %>%
    filter(mesh_id %in% simmat$meshcode_a) -> setable
  new_nrow = nrow(setable)
  if (verbose) {
    cat(file=stderr(),paste0('- ',new_nrow,'/',orig_nrow,' SEs remain.\n'))
    flush.console()
  }
  
  
  ####
  # 3. make cartesian product
  if (verbose) {
    cat(file=stderr(),paste0('3. Creating the cartesian product of ',nrow(drugtable),' drugs x ',nrow(setable),' SEs.\n'))
    flush.console()
  }
  expand.grid(drug_name=drugtable$drug_name, se_mesh_id=setable$mesh_id) %>%
    mutate(dse_uid = paste0(drug_name, '-', se_mesh_id)) %>%
    inner_join(mesh_best_names, by=c('se_mesh_id'='id')) %>%
    rename(se_mesh_term = labeltext) %>%
    filter(se_mesh_id %in% simmat$meshcode_a) %>%
    as_tibble() -> cp
  
  ####
  # 4. annotate in most similar indication by sort and de-dup
  if (is.null(presimindic)) {
    if (verbose) {
      cat(file=stderr(),paste0('4. Selecting most similar approved indications.\n'))
      flush.console()
    }
    simmat %>%
      filter(meshcode_a %in% cp$se_mesh_id) %>%
      filter(meshcode_b %in% approved_indics$indication_mesh_id) -> sims
    if (verbose) {
      n_missing = sum(!(approved_indics$indication_mesh_id %in% simmat$meshcode_b))
      n_unique_missing = length(unique(approved_indics$indication_mesh_id[!approved_indics$indication_mesh_id %in% simmat$meshcode_b])) 
      cat(file=stderr(),paste0('- ',n_missing,' indications (',n_unique_missing,' unique) are missing from similarity matrix.\n'))
      missing_from_matrix = tibble(mesh_id=unique(approved_indics$indication_mesh_id[!approved_indics$indication_mesh_id %in% simmat$meshcode_b]))
      missing_from_matrix$mesh_term = mesh_best_names$labeltext[match(missing_from_matrix$mesh_id, mesh_best_names$id)]
      write_tsv(missing_from_matrix, 'ignore/pp_indics_missing_from_sim_matrix.tsv')
      flush.console()
    }
    
    cp %>%
      inner_join(approved_indics, by=c('drug_name'='drug_name')) %>%
      left_join(sims, by=c('se_mesh_id'='meshcode_a', 'indication_mesh_id'='meshcode_b')) %>%
      mutate(comb_norm = replace_na(comb_norm, 0)) %>%
      mutate(dse_uid = paste0(drug_name,'-',se_mesh_id)) %>%
      select(dse_uid, drug_name, se_mesh_id, se_mesh_term, indication_mesh_id, se_indic_similarity=comb_norm) %>%
      arrange(dse_uid, desc(se_indic_similarity)) %>%
      group_by(dse_uid) %>%
      slice(1) %>%
      ungroup() %>%
      inner_join(mesh_best_names, by=c('indication_mesh_id'='id')) %>%
      rename(indication_mesh_term = labeltext) -> cp_indic
    
    if (verbose) {
      cat(file=stderr(),paste0('- ',sum(is.na(cp_indic$indication_mesh_id)),'/',nrow(cp_indic),' rows lack a most similar indication.\n'))
      flush.console()
    }
  } else {
    if (verbose) {
      cat(file=stderr(),paste0('4. Using most similar approved indications provided.\n'))
      flush.console()
    }
    presimindic %>%
      select(dse_uid, indication_mesh_id, indication_mesh_term, se_indic_similarity) -> presimindic_skinny
    cp %>%
      inner_join(presimindic_skinny, by='dse_uid') -> cp_indic
    
  }
  
  
  if (verbose) {
    cat(file=stderr(),paste0('5. Annotating most similar genetic associations.\n'))
    flush.console()
  }
  if (verbose) {
    cat(file=stderr(),paste0('- Beginning from cartesian product table, ',nrow(cp),' rows.\n'))
    flush.console()
  }
  cp_indic %>%
    inner_join(drug_target_match, by='drug_name') %>%
    select(dse_uid, drug_name, se_mesh_id, se_mesh_term, gene) -> cp_with_targets
  if (verbose) {
    cat(file=stderr(),paste0('- Added targets to cartesian product table, ',nrow(cp_with_targets),' rows total.\n'))
    flush.console()
  }
  simmat %>%
    filter(meshcode_a %in% cp$se_mesh_id) %>%
    filter(meshcode_b %in% assocs$mesh_id) -> sims
  assocs %>%
    mutate(gmid = paste0(gene, '-', mesh_id)) %>%
    group_by(gmid) %>%
    slice(1) %>%
    ungroup() -> assoc1
  cp_with_targets %>%
    inner_join(assoc1, by='gene') %>%
    select(dse_uid, drug_name, se_mesh_id, se_mesh_term, gene,
           arow, mesh_id, mesh_term, source, extra_info, original_trait) %>%
    rename(assoc_mesh_id = mesh_id,
           assoc_mesh_term = mesh_term,
           assoc_source = source,
           assoc_info = extra_info,
           assoc_orig = original_trait) -> merge1_present
  cp_with_targets %>%
    mutate(arow = as.integer(NA),
           assoc_mesh_id = as.character(NA),
           assoc_mesh_term = as.character(NA),
           assoc_source = as.character(NA),
           assoc_info = as.character(NA),
           assoc_orig = as.character(NA)) -> merge1_blank
  merge1 = rbind(merge1_present, merge1_blank)
  
  if (verbose) {
    cat(file=stderr(),paste0('- Annotated all unique genetic associations, ',nrow(merge1),' rows.\n'))
    flush.console()
  }
  
  # running line-by-line this works, and running it as a command line script,
  # this works. 
  # but when running make_dse as a function call in interactive mode,
  # this step gives:
  # Error: vector memory exhausted (limit reached?)
  merge1 %>%
    left_join(sims, by=c('se_mesh_id'='meshcode_a', 'assoc_mesh_id'='meshcode_b')) %>%
    mutate(comb_norm = case_when(se_mesh_id==assoc_mesh_id ~ 1,
                                 is.na(comb_norm) ~ 0,
                                 TRUE ~ comb_norm)) %>%
    rename(se_assoc_similarity = comb_norm) -> merge2
  
  if (verbose) {
    cat(file=stderr(),paste0('- Added similarities for genetic associations, ',nrow(merge2),' rows.\n'))
    flush.console()
  }
  
  merge2 %>%
    group_by(dse_uid) %>%
    arrange(desc(se_assoc_similarity)) %>%
    slice(1) %>%
    ungroup() -> mtable
  # note you could make an argument for this implementation:
  # slice_max(similarity, n=1, with_ties=FALSE) 
  # instead of arrange followed by slice(1) 
  # however i found that it for some reason was far slower
  
  # note: from here on out there should be no more one-to-many drug-target matches.
  # we have preserved for each drug only the one target that yields the most similar genetic association.
  
  if (verbose) {
    cat(file=stderr(),paste0('- Removed duplicates, preserving only most similar association, ',nrow(mtable),' rows remain.\n'))
    flush.console()
  }
  
  ####
  # 6. test for target enrichment with nominal fisher exact test
  if (!skip_enrich) {
    
    if (verbose) {
      cat(file=stderr(),paste0('6. Testing for target enrichment.\n'))
      flush.console()
    }
    # only bother with:
    # SEs observed for at least 2 drugs, 
    # targets observed for at least 2 drugs, and 
    # combination observed ever
    drugsematch %>% 
      filter(!is.na(mesh_id)) %>%
      group_by(mesh_id) %>%
      summarize(.groups='keep', n_drugs = length(unique(drug_name))) %>%
      ungroup() %>%
      filter(n_drugs > 1) -> se_use
    drug_target_match %>%
      filter(!is.na(gene)) %>%
      group_by(gene) %>%
      summarize(.groups='keep', n_drugs = length(unique(drug_name))) %>%
      ungroup() %>%
      filter(n_drugs > 1) -> gene_use
    drugsematch %>%
      inner_join(se_use, by='mesh_id') %>%
      inner_join(drug_target_match, by='drug_name') %>%
      inner_join(gene_use, by='gene') %>%
      distinct(gene, se_mesh_id = mesh_id) %>%
      mutate(target_or = as.numeric(NA),
             target_p = as.numeric(NA)) -> enrich
    
    if (verbose) {
      cat(file=stderr(),paste0('- Identified ',nrow(enrich),' gene-SE combinations that qualify for enrichment testing.\n'))
      flush.console()
    }
    cp %>%
      distinct(drug_name) -> dt
    for (i in 1:nrow(enrich)) {
      
      if (verbose) {
        cat(file=stderr(),paste0('\r- Enrichment testing row ',i,'/',nrow(enrich),'...'))
        flush.console()
      }
      dt$target = dt$drug_name %in% drug_target_match$drug_name[drug_target_match$gene==enrich$gene[i]]
      dt$se = dt$drug_name %in% drugsematch$drug_name[drugsematch$mesh_id==enrich$se_mesh_id[i]]
      ctable = table(dt[,c('target','se')])
      if (prod(dim(ctable)) < 4) {
        enrich$target_or[i] = 1
        enrich$target_p[i] = 1
      } else {
        fobj = fisher.test(ctable)
        enrich$target_or[i] = fobj$estimate
        enrich$target_p[i] = fobj$p.value
      }
    }
    enrich$fdr = p.adjust(enrich$target_p, method='BH')
    enrich$target_enriched = enrich$target_or >=2 & enrich$fdr < 0.05
    if (verbose) {
      cat(file=stderr(),paste0('done!\n- ',sum(enrich$target_enriched),' gene-SE combinations show enrichment at OR >= 2 and nominal P < 0.01.\n'))
      flush.console()
    }
    
  } else {
    cp_indic %>%
      distinct(se_mesh_id) %>%
      inner_join(drugsematch, by=c('se_mesh_id'='mesh_id')) %>%
      inner_join(drug_target_match, by=c('drug_name')) %>%
      distinct(gene, se_mesh_id) %>%
      mutate(target_or=NA, target_p=NA, target_enriched=NA) -> enrich
  }
  
  ###
  # 7. merge & bring back any extra columns from drugsematch
  if (verbose) {
    cat(file=stderr(),paste0('7. Merging all data back together\n'))
    flush.console()
  }
  cp_indic %>%
    inner_join(mtable, by=c('dse_uid','drug_name','se_mesh_id','se_mesh_term')) %>%
    select(dse_uid,
           drug_name,
           se_mesh_id,
           se_mesh_term,
           indication_mesh_id,
           indication_mesh_term,
           se_indic_similarity,
           gene,
           assoc_mesh_id,
           assoc_mesh_term,
           se_assoc_similarity,
           arow,
           assoc_source,
           assoc_info,
           assoc_orig
           ) %>%
    left_join(enrich, by=c('gene','se_mesh_id')) %>%
    mutate(target_or = replace_na(target_or, 1)) %>%
    mutate(target_p = replace_na(target_p, 1)) %>%
    mutate(target_enriched = replace_na(target_enriched, F)) %>%
    mutate(observed = dse_uid %in% drugsematch$dse_uid) -> dse_info
  
  dse_info$genetic_insight = se_metadata$genetic_insight[match(dse_info$se_mesh_id, se_metadata$mesh_id)]
  dse_info$genetic_insight[is.na(dse_info$genetic_insight)] = 'none'
  
  # more_dse_cols = setdiff(colnames(drugsematch), colnames(dse_info))
  more_dse_cols = c('numfreq','wordfreq','wordrank','placebo')
  for (col in more_dse_cols) {
    dse_info[,col] = drugsematch[match(dse_info$dse_uid, drugsematch$dse_uid),col]
  }
  # more_drug_cols = setdiff(colnames(drugtable), colnames(dse_info))
  more_drug_cols = c('n_se')
  for (col in more_drug_cols) {
    dse_info[,col] = drugtable[match(dse_info$drug_name, drugtable$drug_name),col]
  }
  # more_se_cols = setdiff(colnames(setable), colnames(dse_info))
  more_se_cols = c('n_drugs','max_severity')
  for (col in more_se_cols) {
    dse_info[,col] = setable[match(dse_info$se_mesh_id, setable$mesh_id),col]
  }
  if (verbose) {
    cat(file=stderr(),paste0('- Final row count: ',nrow(dse_info),'.\n'))
    cat(file=stderr(),paste0('- *Of which ',sum(dse_info$genetic_insight!='none'),' are genetic insight SEs.\n'))
    cat(file=stderr(),paste0('- *And ',sum(dse_info$target_enriched),' are target-enriched.\n'))
    cat(file=stderr(),paste0('- Final columns: ',paste0(colnames(dse_info),collapse=','),'.\n'))
    flush.console()
  }
  
  return(dse_info)
}


careful_fisher_test = function(ctable) {
  if (prod(dim(ctable)) < 4) {
    fobj = list(estimate=1, conf.int=c(0,Inf), p.value=1)
    if ('TRUE' %in% rownames(ctable)) {
      fobj$trtot = sum(ctable['TRUE',])
      if ('TRUE' %in% colnames(ctable)) {
        fobj$yy = ctable['TRUE','TRUE']
        fobj$tctot = sum(ctable[,'TRUE'])
      } else {
        fobj$tctot = 0
        fobj$yy = 0
      }
    } else {
      fobj$trtot = 0
      fobj$yy = 0
    }
  } else {
    fobj = fisher.test(ctable)
    fobj$yy = ctable['TRUE','TRUE']
    fobj$trtot = sum(ctable['TRUE',])
    fobj$tctot = sum(ctable[,'TRUE'])
  }
  return (fobj)
}




##############
# DATA
##############

tell_user('done.\nReading in data...')

# need guess_max
assoc = read_tsv('data/assocs/assoc.tsv.gz', col_types=cols(), guess_max=1e5) %>%
  filter(!(source %in% 'OTG') | l2g_share >= 0.5)

all_se = read_tsv('data/sider/meddra_all_se.tsv.gz', col_types=cols(),
                  col_names = c('drug_id','id2','umls_label','meddra_type','umls_meddra','se_name'))

drug_names = read_tsv('data/synthesis/sider_drug_names.tsv', col_types=cols())

freq = read_tsv('data/sider/meddra_freq.tsv.gz', col_types=cols(),
                col_names = c('drug_id','id2','umls_label','placebo','freq','lb','ub','meddra_type','umls_meddra','se_name'))

pp_targets = read_tsv('data/synthesis/pp_targets.tsv', col_types=cols())

pp_launched_indic = read_tsv('data/synthesis/pp_launched_indic.tsv', col_types=cols())

mesh_best_names = read_tsv('data/mesh/mesh_best_names.tsv.gz', col_types=cols())

sider_mesh_map = read_tsv('data/mesh/sider_mesh_map.tsv', col_types=cols()) %>%
  left_join(mesh_best_names, by=c('mesh_id'='id')) %>%
  mutate(mesh_term = labeltext) %>%
  select(se_name, mesh_id, mesh_term, map_type)
# check for missing types # sider_mesh_map[!is.na(sider_mesh_map$mesh_id) & is.na(sider_mesh_map$map_type),] %>% View()

se_topl_match = read_tsv('data/synthesis/se_topl_match.tsv', col_types=cols())

se_insight = read_tsv('data/synthesis/se_insight.tsv', col_types=cols())

sim = read_tsv('data/mesh/mesh_sim.tsv.gz', col_types=cols())

hubal_day = read_tsv('data/other/hubal-day-2006-table-3.tsv', col_types=cols()) %>%
  arrange(task_value) %>%
  mutate(wordrank = rank(task_value))
gottlieb = read_tsv('data/other/gottlieb-2015-table-s2.tsv', col_types=cols()) %>% 
  clean_names()

sider_mesh_map$severity = gottlieb$rank_score[match(tolower(sider_mesh_map$se_name), tolower(gottlieb$name))]

# check none are invalid - 0, great
# sider_mesh_map %>%
#   distinct(mesh_id, mesh_term) %>%
#   filter(!(mesh_id %in% mesh_best_names$id))

# sider_mesh_map %>%
#   distinct(mesh_id, mesh_term) %>%
#   filter(!(mesh_id %in% sim$meshcode_a))
# 
# sider_mesh_map %>%
#   group_by(map_type) %>%
#   summarize(.groups='keep', n=n()) %>%
#   ungroup() %>%
#   arrange(desc(n)) %>%
#   mutate(proportion = n/sum(n)) %>% View()
# 
# sample(sider_mesh_map$se_name[is.na(sider_mesh_map$mesh_id)], size=10)


# join together the datasets to feed into make_dse

drug_names %>%
  inner_join(all_se, by='drug_id') %>%
  inner_join(sider_mesh_map, by='se_name') %>%
  filter(!is.na(mesh_id) & !is.na(ppid)) %>% # keep only SEs mapped to MeSH & drugs mapped to PP
  distinct(drug_name, ppid, mesh_id, mesh_term, severity) -> drug_mesh_match
  
drug_mesh_match %>%
  group_by(drug_name) %>%
  summarize(.groups='keep', ppid=min(ppid), n_se=length(unique(mesh_id))) -> sider_drugs

drug_mesh_match %>%
  group_by(mesh_id) %>%
  summarize(.groups='keep', 
            n_drugs=length(unique(drug_name)),
            max_severity = suppressWarnings(max(severity, na.rm=T))) %>% # many have no values, will get -Inf, don't need to see warning message
  ungroup() %>%
  left_join(mesh_best_names, by=c('mesh_id'='id')) %>%
  mutate(mesh_term=labeltext) %>%
  mutate(max_severity = replace(max_severity, is.infinite(max_severity), NA)) %>%
  select(mesh_id, mesh_term, n_drugs, max_severity) -> sider_ses

# note the group-bys above serve to de-dup the drug names table, see example of adefovir 13932, it has two distinct CIDs

freq %>%
  select(drug_id, se_name, placebo, chrfreq=freq) %>%
  mutate(numfreq = suppressWarnings(as.numeric(gsub('%','',chrfreq))/100)) %>%
  mutate(wordfreq = case_when(chrfreq %in% hubal_day$term ~ chrfreq,
                              TRUE ~ as.character(NA))) %>%
  left_join(hubal_day, by=c('wordfreq'='term')) %>%
  select(drug_id, se_name, numfreq, wordfreq, wordrank, placebo) %>%
  inner_join(sider_mesh_map, by='se_name') %>%
  inner_join(drug_names, by='drug_id') %>%
  filter(!is.na(mesh_id) & !is.na(ppid)) %>%
  group_by(drug_name, mesh_id) %>%
  summarize(.groups='keep',
            minppid = min(ppid),
            maxnumfreq = max(numfreq),
            maxwordrank=max(wordrank),
            anyplacebo=any(placebo %in% 'placebo')) %>%
  ungroup() %>%
  left_join(hubal_day, by=c('maxwordrank'='wordrank')) %>%
  mutate(dse_uid = paste0(drug_name,'-',mesh_id)) %>%
  select(dse_uid, drug_name, ppid=minppid, mesh_id, numfreq=maxnumfreq, wordfreq=term, wordrank=maxwordrank, placebo=anyplacebo) -> sider_freq
  
drug_mesh_match %>%
  group_by(drug_name, mesh_id) %>%
  summarize(.groups='keep', ppid=min(ppid), mesh_term=min(mesh_term)) %>%
  ungroup() %>%
  mutate(observed=TRUE) %>%
  mutate(dse_uid = paste0(drug_name, '-', mesh_id)) %>%
  left_join(select(sider_freq, dse_uid, numfreq, wordfreq, wordrank, placebo), by='dse_uid') %>%
  select(dse_uid, drug_name, ppid, mesh_id, mesh_term, numfreq, wordfreq, wordrank, placebo) -> sider_drugse_match

tell_user('done.\nCreating and annotating drug-SE matrix...')

dse_info = make_dse(sider_drugs, 
         sider_ses,
         sider_drugse_match,
         sim,
         assoc,
         se_insight,
         verbose=F, skip_enrich=F)

sim_threshold = 0.9
dse_info$sim_assoc = dse_info$se_assoc_similarity >= sim_threshold
dse_info$sim_indic = dse_info$se_indic_similarity >= sim_threshold

write_tsv(dse_info, 'intermediate/dse_info.tsv.gz')

# if debugging in interactive mode, you can just run this to save several minutes of creating the drug-SE matrix:
# dse_info = read_tsv('intermediate/dse_info.tsv.gz', col_types=cols())

###
# SUPPLEMENTARY TABLES
###

tell_user('done.\nCreating supplementary tables...')


tell_user('done.\nCreating Table S1...')

dse_info %>%
  summarize( n_drugs = length(unique(drug_name)),
             n_ses = length(unique(se_mesh_id)),
             n_rows = n(),
             n_observed = sum(observed),
             n_simassoc = sum(sim_assoc),
             n_simindic = sum(sim_indic),
             n_geninsight = sum(genetic_insight != 'none'),
             n_targenrich = sum(target_enriched)
  ) %>%
  ungroup() %>%
  t() -> temp
colnames(temp) = 'count'
temp %>%
  as_tibble(rownames='attribute', .name_repair = 'universal') %>%
  mutate(proportion = count/count[attribute=='n_rows']) -> table_1

write_supp_table(table_1, 'Properties of drug-side effect matrix.')


# types of matches to PP
drug_names %>%
  group_by(pp_match_type) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  mutate(proportion = n/sum(n)) -> pp_match_types

write_supp_table(pp_match_types, 'Methods of matching SIDER drug names to Pharmaprojects.')

sider_mesh_map %>%
  group_by(map_type) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  mutate(proportion = n/sum(n)) %>%
  mutate(map_type=replace_na(map_type,'UNMAPPED')) -> sider_mesh_sources

write_supp_table(sider_mesh_sources, 'Methods of matching SIDER side effect names to MeSH terms.')

# reasons drugs drop out of analysis
drug_names %>%
  mutate(duplicate_drug_name = duplicated(tolower(drug_name))) %>%
  mutate(unmapped_to_pharmaprojects = is.na(pp_match_type) | pp_match_type == 'missing') %>%
  mutate(lacks_approved_indications = !(ppid %in% pp_launched_indic$ppid)) %>%
  mutate(lacks_annotated_human_target = !(ppid %in% pp_targets$ppid)) %>%
  pivot_longer(cols=c(duplicate_drug_name,unmapped_to_pharmaprojects,lacks_approved_indications,lacks_annotated_human_target)) %>%
  group_by(name) %>%
  summarize(.groups='keep', n=sum(value)) %>%
  ungroup() -> dropout_stats

write_supp_table(dropout_stats, 'Reasons SIDER drugs drop out of analysis.')


######
# FIGURE S1
######

tell_user('done.\nCreating Figure S1...')

make_forest = function(forest_tibble, verbose=T) {
  
  for (i in 1:nrow(forest_tibble)) {
    
    if (verbose) {
      write(file=stderr(), paste0('\rNow on row ',i,'...'))
      flush.console()
    }
    
    # this just redoes everything - very inefficient:
    # dsetable = make_dse(sider_drugs, get(forest_tibble$setable[i]), sider_drugse_match, sim, get(forest_tibble$assoctable[i]), se_insight, verbose=T, skip_enrich=T)
    
    get(forest_tibble$setable[i]) %>%
      filter(mesh_id %in% se_insight$mesh_id[se_insight$genetic_insight != 'none']) -> this_ses
    get(forest_tibble$assoctable[i]) %>%
      group_by(gene, mesh_id) %>%
      slice(1) %>%
      ungroup() -> this_assocs
    
    # use the precast SE-indic similarities to avoid redoing the hardest step
    dse_info %>%
      select(dse_uid, indication_mesh_id, indication_mesh_term, se_indic_similarity) -> presimindic
    dsetable = make_dse(sider_drugs, this_ses, sider_drugse_match, sim, this_assocs, se_insight, verbose=verbose, skip_enrich=T, presimindic=presimindic)
    
    dsetable$sim_assoc = dsetable$se_assoc_similarity >= sim_threshold
    dsetable$sim_indic = dsetable$se_indic_similarity >= sim_threshold
    dsetable$geninsight = dsetable$genetic_insight != 'none'
    
    ctable = table(dsetable[dsetable$genetic_insight !='none' & !dsetable$sim_indic,c('sim_assoc','observed')])
    fobj = careful_fisher_test(ctable)
    forest_tibble$or[i] = as.numeric(fobj$estimate)
    forest_tibble[i,c('or_l95','or_u95')] = as.list(as.numeric(fobj$conf.int))
    forest_tibble$p[i] = as.numeric(fobj$p.value)
    forest_tibble$sim_obs[i] = ctable['TRUE','TRUE']
    forest_tibble$sim_not[i] = ctable['TRUE','FALSE']
    forest_tibble$dis_obs[i] = ctable['FALSE','TRUE']
    forest_tibble$dis_not[i] = ctable['FALSE','FALSE']
    forest_tibble$n_rows[i] = sum(ctable)
  }
  
  return(forest_tibble)
  
}


sider_ses %>%
  filter(mesh_id %in% se_topl_match$mesh_id[se_topl_match$topl=='C04']) -> sider_onco_ses

sider_ses %>%
  filter(!(mesh_id %in% se_topl_match$mesh_id[se_topl_match$topl=='C04'])) -> sider_nononco_ses

assoc_omim = assoc[assoc$source=='OMIM',]
assoc_intogen = assoc[assoc$source=='intOGen',] 
assoc_non_intogen = assoc[assoc$source!='intOGen',] 
assoc_gwas = assoc[assoc$source %in% c('PICCOLO','OTG','Genebass'),]
assoc_otg = assoc[assoc$source %in% c('OTG'),]
assoc_piccolo = assoc[assoc$source %in% c('PICCOLO'),]
assoc_genebass = assoc[assoc$source %in% c('Genebass'),]

assoc_source_forest_setup = tibble(label=c('all','OMIM','IntOGen oncology','GWAS','OTG','PICCOLO','Genebass'), 
                                   assoctable=c('assoc','assoc_omim','assoc_intogen','assoc_gwas','assoc_otg','assoc_piccolo','assoc_genebass'), 
                                   setable = c('sider_ses', 'sider_ses', 'sider_onco_ses', 'sider_ses', 'sider_ses', 'sider_ses', 'sider_ses'),
                                   or      = as.numeric(NA), 
                                   or_l95  = as.numeric(NA), 
                                   or_u95  = as.numeric(NA), 
                                   p       = as.numeric(NA), 
                                   sim_obs = as.integer(NA),
                                   sim_not = as.integer(NA),
                                   dis_obs = as.integer(NA),
                                   dis_not = as.integer(NA),
                                   n_rows  = as.integer(NA)) %>%
  mutate(y = max(row_number()) - row_number() + 1)

assoc_source_forest = make_forest(assoc_source_forest_setup, verbose=F)


make_figure_s1 = TRUE
if (make_figure_s1) {
  resx=300
  png(paste0('display_items/figure-s1.png'),width=6.5*resx,height=2.5*resx,res=resx)
  
  layout_matrix = matrix(c(1,1,2), nrow=1, byrow=T)
  layout(layout_matrix)
  panel = 1
  
  
  dse_info$geninsight = dse_info$genetic_insight != 'none'
  test_columns = tibble(var1=c('observed','sim_assoc','sim_indic','geninsight','target_enriched'),
                        var2=c('observed','sim_assoc','sim_indic','geninsight','target_enriched'),
                        longdisp = c('SE observed','similar association','similar indication','studied genetically','target enriched')) %>%
    mutate(y = max(row_number()) - row_number() + 1,
           x = row_number())
  crossing(select(test_columns, var1, x), select(test_columns, var2, y)) %>%
    filter(y + x <= max(y+x)/2+1) %>%
    mutate(or = as.numeric(NA), p = as.numeric(NA), yy = as.integer(NA)) -> crosstabs
  for (i in 1:nrow(crosstabs)) {
    crosstabs$yy[i] = sum(dse_info[,crosstabs$var1[i]] & dse_info[,crosstabs$var2[i]])
    fobj = careful_fisher_test(table(dse_info[,c(crosstabs$var1[i],crosstabs$var2[i])]))
    crosstabs$or[i] = as.numeric(fobj$estimate)
    crosstabs$p[i] = fobj$p.value
  }
  
  color_ramp = c('#e5f5e0','#c7e9c0','#a1d99b','#74c476','#41ab5d')
  color_ramper = colorRamp(color_ramp)
  color_scale = tibble(floor = 100*(10^seq(-2,0,length.out=20)))
  color_scale$color = rgb(color_ramper((2 + log10(color_scale$floor/100))/2), maxColorValue=255)
  # orparams = data.frame(floor=c(1,3,10,30,100), color=color_ramp)
  crosstabs$idx = floor(log10(crosstabs$or)*2)+1
  # crosstabs$floor = orparams$floor[crosstabs$idx]
  crosstabs$color = as.character(NA)
  crosstabs$color[!is.infinite(crosstabs$or)] = rgb(color_ramper((2 + log10(crosstabs$or[!is.infinite(crosstabs$or)]/100))/2), maxColorValue=255)
  # crosstabs$color = orparams$color[match(crosstabs$floor, orparams$floor)]
  crosstabs$color[is.infinite(crosstabs$or)] = '#000000'
  
  
  par(mar=c(5,7,3,2))
  xlims = range(crosstabs$x) + c(-0.5, 0.5)
  ylims = range(crosstabs$y) + c(-0.5, 0.5)
  radius = 0.4
  plot(x=NA, y=NA, xlim = xlims, ylim=ylims, axes=F, ann=F)
  rect(xleft=crosstabs$x-radius, xright=crosstabs$x+radius, ybottom=crosstabs$y-radius, ytop=crosstabs$y+radius,
       col=crosstabs$color, border=NA)
  rect(xleft=crosstabs$x-radius, xright=crosstabs$x+radius, ybottom=crosstabs$y-radius, ytop=crosstabs$y+radius,
       col=NA, border='#000000')
  crosstabs %>%
    filter(!is.infinite(or)) -> subs
  text(x=subs$x, y=subs$y, labels=formatC(subs$or, format='f',digits=1))
  par(xpd=T)
  text(x=rep(0.4, nrow(test_columns)), y=test_columns$y, labels=test_columns$longdisp, adj=1, cex=1)
  text(x=test_columns$x, y=rep(0.4, nrow(test_columns)), labels=test_columns$longdisp, srt=30, adj=1, cex=1)
  par(xpd=F)
  
  orlabels = tibble(floor=c(1,3,10,30,100))
  orlabels$color = rgb(color_ramper((2 + log10(orlabels$floor/100))/2), maxColorValue=255)
  scale_points = nrow(color_scale)
  scale_xleft = 3.3
  scale_xright = 4.8
  scale_factor = (scale_xright - scale_xleft) / scale_points
  scale_ybot = 5.1
  scale_ytop = 5.4
  y_offset = 0.1
  par(xpd=T)
  rect(xleft=scale_xleft,xright=scale_xright,ybottom=scale_ybot,ytop=scale_ytop,lwd=2.5)
  rect(xleft=scale_xleft+(0:(scale_points-1))*scale_factor, xright=scale_xleft+(1:scale_points)*scale_factor, ybottom=rep(scale_ybot,scale_points), ytop=rep(scale_ytop,scale_points), border=NA, col=color_scale$color)
  text(x=(scale_xright + scale_xleft)/2, y=scale_ytop-y_offset, pos=3, labels="odds ratio", cex=0.9)
  text(x=seq(scale_xleft, scale_xright, length.out=nrow(orlabels)),y=rep(scale_ybot,nrow(orlabels))+y_offset,pos=1,labels=orlabels$floor,cex=0.9,font=3)
  par(xpd=F)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = -0.4)
  panel = panel + 1
  
  write_supp_table(crosstabs, 'Cross-tabulation of drug-SE properties.')
  
  
  metrics_tbl = tibble(disp=c('no filters','require studied genetically','remove similar indications','both'),
                       require_insight = c(F,T,F,T),
                       remove_simindic = c(F,F,T,T),
                       y = 4:1,
                       or = as.numeric(NA),
                       or_l95 = as.numeric(NA),
                       or_u95 = as.numeric(NA),
                       p = as.numeric(NA),
                       sim_obs = as.integer(NA),
                       sim_not = as.integer(NA),
                       dis_obs = as.integer(NA),
                       dis_not = as.integer(NA),
                       n_rows = as.integer(NA))
  for (i in 1:nrow(metrics_tbl)) {
    dse_info %>%
      filter(geninsight | !metrics_tbl$require_insight[i]) %>%
      filter(!sim_indic | !metrics_tbl$remove_simindic[i]) -> subs
    ctable = table(subs[,c('sim_assoc','observed')])
    fobj = fisher.test(ctable)
    metrics_tbl$sim_obs[i] = ctable['TRUE','TRUE']
    metrics_tbl$sim_not[i] = ctable['TRUE','FALSE']
    metrics_tbl$dis_obs[i] = ctable['FALSE','TRUE']
    metrics_tbl$dis_not[i] = ctable['FALSE','FALSE']
    metrics_tbl$n_rows[i] = sum(ctable)
    metrics_tbl$or[i] = fobj$estimate
    metrics_tbl$p[i] = fobj$p.value
    metrics_tbl[i,c('or_l95','or_u95')] = as.list(as.numeric(fobj$conf.int))
  }
  metrics_tbl %>%
    mutate(denominator = sim_obs + sim_not) %>%
    rename(mean=or, l95=or_l95, u95=or_u95, label=disp, numerator=sim_obs) -> metrics_frst
  
  xlims = c(0,3)
  ylims = range(metrics_frst$y) + c(-0.5, 0.5)
  xats = 0:8/2
  xbigs = 0:4
  par(mar=c(3,4,3,1))
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, lwd=0, line=-0.25, cex.axis=1)
  abline(v=1, lty=3)
  mtext(side=1, line=1.75, text='OR')
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  par(lheight=0.8)
  mtext(side=2, at=metrics_tbl$y, text=gsub(' ','\n',metrics_tbl$disp), cex=0.7, line=0.25, las=2)
  par(lheight=1.0)
  points(metrics_tbl$or, metrics_tbl$y, pch=19)
  segments(x0=metrics_tbl$or_l95, x1=metrics_tbl$or_u95, y0=metrics_tbl$y, lwd=1.5)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5)
  panel = panel + 1
  
  write_supp_table(metrics_tbl, 'Effects of requiring genetic insight and removing similar indications.')
  
  end_of_figure_s1 = dev.off()
}


######
# FIGURE 1
######

tell_user('done.\nCreating Figure 1...')


make_figure_1 = TRUE
if (make_figure_1) {
  
  resx=300
  png(paste0('display_items/figure-1.png'),width=6.5*resx,height=2.5*resx,res=resx)
  
  layout_matrix = matrix(c(1,2,3), nrow=1, byrow=T)
  layout(layout_matrix)
  panel = 1
  
  rocassoc_tbl = tibble(threshold = 1:20/20,
                        or = as.numeric(NA),
                        or_l95 = as.numeric(NA),
                        or_u95 = as.numeric(NA),
                        p = as.numeric(NA),
                        sim_obs = as.integer(NA),
                        sim_not = as.integer(NA),
                        dis_obs = as.integer(NA),
                        dis_not = as.integer(NA),
                        n_rows = as.integer(NA))
  dse_info %>%
    filter(geninsight & !sim_indic) -> dse_strip
  for (i in 1:nrow(rocassoc_tbl)) {
    dse_strip$temp_sim = dse_strip$se_assoc_similarity >= rocassoc_tbl$threshold[i]
    ctable = table(dse_strip[,c('temp_sim','observed')])
    fobj = fisher.test(ctable)
    rocassoc_tbl$sim_obs[i] = ctable['TRUE','TRUE']
    rocassoc_tbl$sim_not[i] = ctable['TRUE','FALSE']
    rocassoc_tbl$dis_obs[i] = ctable['FALSE','TRUE']
    rocassoc_tbl$dis_not[i] = ctable['FALSE','FALSE']
    rocassoc_tbl$n_rows[i] = sum(ctable)
    rocassoc_tbl$or[i] = fobj$estimate
    rocassoc_tbl$p[i] = fobj$p.value
    rocassoc_tbl[i,c('or_l95','or_u95')] = as.list(as.numeric(fobj$conf.int))
  }
  roc_col = '#FF9912'
  xlims = c(0, 45000)
  ylims = c(0, 2.5)
  xbigs = 0:5*10000
  xats = 0:50*1000
  ybigs = 0:3
  yats = 0:6/2
  
  par(mar=c(5,3,3,1))
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, labels=gsub('000$','K',xbigs), lwd=0, line=-0.25, cex.axis=1)
  mtext(side=1, line=3.0, text='observed SEs with\nsimilar association', cex=0.8, padj=0)
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, lwd=0, line=-0.25, cex.axis=1, las=2)
  mtext(side=2, line=1.25, text='OR', cex=0.8)
  abline(h=1, lty=3)
  points(rocassoc_tbl$sim_obs, rocassoc_tbl$or, type='l', lwd=1.5, col=roc_col)
  points(rocassoc_tbl$sim_obs, rocassoc_tbl$or, pch=20, col=roc_col)
  polygon(x=c(rocassoc_tbl$sim_obs, rev(rocassoc_tbl$sim_obs)), y=c(rocassoc_tbl$or_l95, rev(rocassoc_tbl$or_u95)), border=NA, col=alpha(roc_col,ci_alpha))
  rocassoc_tbl %>%
    filter(threshold %in% c(1, .9, .8, .5, .2)) -> to_label
  par(xpd=T)
  text(x=to_label$sim_obs, y=to_label$or, pos=4, labels=formatC(to_label$threshold, digits=2, format='f'), col=roc_col, cex=0.8)
  par(xpd=F)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5)
  panel = panel + 1
  
  write_supp_table(rocassoc_tbl, 'Sensitivity to similarity threshold for inclusion of similar genetic associations.')
  
  rocindic_tbl = tibble(threshold = 1:20/20,
                        or = as.numeric(NA),
                        or_l95 = as.numeric(NA),
                        or_u95 = as.numeric(NA),
                        p = as.numeric(NA),
                        sim_obs = as.integer(NA),
                        sim_not = as.integer(NA),
                        dis_obs = as.integer(NA),
                        dis_not = as.integer(NA),
                        n_rows = as.integer(NA))
  for (i in 1:nrow(rocindic_tbl)) {
    dse_info %>%
      filter(geninsight & se_indic_similarity < rocindic_tbl$threshold[i]) -> dse_strip
    ctable = table(dse_strip[,c('sim_assoc','observed')])
    fobj = fisher.test(ctable)
    rocindic_tbl$sim_obs[i] = ctable['TRUE','TRUE']
    rocindic_tbl$sim_not[i] = ctable['TRUE','FALSE']
    rocindic_tbl$dis_obs[i] = ctable['FALSE','TRUE']
    rocindic_tbl$dis_not[i] = ctable['FALSE','FALSE']
    rocindic_tbl$n_rows[i] = sum(ctable)
    rocindic_tbl$or[i] = fobj$estimate
    rocindic_tbl$p[i] = fobj$p.value
    rocindic_tbl[i,c('or_l95','or_u95')] = as.list(as.numeric(fobj$conf.int))
  }
  
  roc_col = '#FF9912'
  xlims = c(0, 800)
  ylims = c(0, 2.5)
  xbigs = 0:3*500
  xats = 0:15*100
  ybigs = 0:3
  yats = 0:6/2
  par(mar=c(5,3,3,1))
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, labels=formatC(xbigs, big.mark=',', format='d'), lwd=0, line=-0.25, cex.axis=1)
  mtext(side=1, line=3.0, text='observed SEs with\nsimilar association', padj=0, cex=0.8)
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, lwd=0, line=-0.25, cex.axis=1, las=2)
  mtext(side=2, line=1.25, text='OR', cex=0.8)
  abline(h=1, lty=3)
  points(rocindic_tbl$sim_obs, rocindic_tbl$or, type='l', lwd=1.5, col=roc_col)
  points(rocindic_tbl$sim_obs, rocindic_tbl$or, pch=20, col=roc_col)
  polygon(x=c(rocindic_tbl$sim_obs, rev(rocindic_tbl$sim_obs)), y=c(rocindic_tbl$or_l95, rev(rocindic_tbl$or_u95)), border=NA, col=alpha(roc_col,ci_alpha))
  rocindic_tbl %>%
    filter(threshold %in% c(.9, .5, .2)) -> to_label
  par(xpd=T)
  text(x=to_label$sim_obs, y=to_label$or, pos=3, labels=formatC(to_label$threshold, digits=2, format='f'), col=roc_col, cex=0.8)
  par(xpd=F)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5)
  panel = panel + 1
  
  write_supp_table(rocindic_tbl, 'Sensitivity to similarity threshold for removal of similar indications.')
  
  

  
  xlims = c(0,4)
  ylims = range(assoc_source_forest$y) + c(-0.5, 0.5)
  xats = 0:8/2
  xbigs = 0:4
  par(mar=c(3,6,3,1))
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, lwd=0, line=-0.25, cex.axis=1)
  abline(v=1, lty=3)
  mtext(side=1, line=1.75, text='OR')
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  par(lheight=0.8)
  mtext(side=2, at=assoc_source_forest$y, text=gsub(' ','\n',assoc_source_forest$label), cex=0.7, line=0.25, las=2)
  par(lheight=1.0)
  points(assoc_source_forest$or, assoc_source_forest$y, pch=19)
  segments(x0=assoc_source_forest$or_l95, x1=assoc_source_forest$or_u95, y0=assoc_source_forest$y, lwd=1.5)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5)
  panel = panel + 1
  
  write_supp_table(assoc_source_forest, 'Breakdown by source of genetic evidence.')
  
  end_of_figure_1 = dev.off()

}





















######
# FIGURE S2
######

tell_user('done.\nCreating Figure S2...')

onco_breakdown_forest_setup = tibble(label=c('IntOGen all', 'IntOGen oncology','IntOGen non-oncology', 'germline oncology', 'germline non-oncology'), 
                             assoctable=c('assoc_intogen','assoc_intogen','assoc_intogen','assoc_non_intogen','assoc_non_intogen'), 
                             setable = c('sider_ses', 'sider_onco_ses', 'sider_nononco_ses', 'sider_onco_ses', 'sider_nononco_ses'),
                             or      = as.numeric(NA), 
                             or_l95  = as.numeric(NA), 
                             or_u95  = as.numeric(NA), 
                             p       = as.numeric(NA), 
                             sim_obs = as.integer(NA),
                             sim_not = as.integer(NA),
                             dis_obs = as.integer(NA),
                             dis_not = as.integer(NA),
                             n_rows  = as.integer(NA)) %>%
  mutate(y = max(row_number()) - row_number() + 1)

onco_breakdown_forest = make_forest(onco_breakdown_forest_setup, verbose=F)
write_supp_table(onco_breakdown_forest, 'Breakdown by somatic vs. germline and oncology vs. non-oncology.')


make_figure_s2 = TRUE
if (make_figure_s2) {
  
resx=300
png('display_items/figure-s2.png',width=6.5*resx,height=4*resx,res=resx)

layout_matrix = matrix(1:2, nrow=2, byrow=T)
layout(layout_matrix, heights=c(1,.5))

xlims = c(0,10)
ylims = range(onco_breakdown_forest$y) + c(-0.5, 0.5)
xats = 0:20/2
xbigs = 0:10
par(mar=c(3,8,1,1))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=xats, tck=-0.025, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xbigs, lwd=0, line=-0.5, cex.axis=0.8)
abline(v=1, lty=3)
mtext(side=1, line=1.75, text='OR')
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
mtext(side=2, at=onco_breakdown_forest$y, text=onco_breakdown_forest$label, cex=0.8, line=0.25, las=2)
points(onco_breakdown_forest$or, onco_breakdown_forest$y, pch=19)
segments(x0=onco_breakdown_forest$or_l95, x1=onco_breakdown_forest$or_u95, y0=onco_breakdown_forest$y, lwd=1.5)

onco_color = '#66A61E'
nononco_color = '#9959E1'

sider_ses %>%
  mutate(onco = mesh_id %in% se_topl_match$mesh_id[se_topl_match$topl=='C04']) %>%
  mutate(y=as.integer(onco) + 1) %>%
  mutate(color = case_when(onco ~ onco_color, !onco ~ nononco_color)) -> sider_ses_by_onco

sider_ses_by_onco %>%
  group_by(onco, y, color) %>%
  summarize(.groups='keep',
            n    =     n(),
            mean =  mean(n_drugs),
            sd   =    sd(n_drugs),
            l95  = lower(n_drugs),
            u95  = upper(n_drugs)) %>%
  ungroup() -> sider_ses_by_onco_smry

write_supp_table(sider_ses_by_onco_smry, 'Properties of oncology vs. non-oncology SEs')
  
par(mar=c(3,8,1,1))
xlims = c(0.7,1100)
ylims = c(0.5, 2.5) 
xbigs = c(1,10,100,1000,10000)
xats = rep(1:9,3) * 10^rep(0:2,each=9)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='x')
axis(side=1, at=xats, labels=NA, tck=-0.025)
axis(side=1, at=xbigs, labels=NA, tck=-0.05)
axis(side=1, line=-0.5, at=xbigs, lwd=0)
mtext(side=1, line=1, text='N drugs / SE')
# axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
mtext(side=2, line=0.25, las=2, at=sider_ses_by_onco_smry$y, text=paste0(ifelse(sider_ses_by_onco_smry$onco,'','non-'),'oncology'))
set.seed(1)
points(sider_ses_by_onco$n_drugs, jitter(sider_ses_by_onco$y,amount=.25), pch=20, cex=0.5, col=alpha(sider_ses_by_onco$color,0.1))
barwidth = 0.3
segments(x0=sider_ses_by_onco_smry$mean, y0=sider_ses_by_onco_smry$y-barwidth, y1=sider_ses_by_onco_smry$y+barwidth, col=sider_ses_by_onco_smry$color, lwd=2)
arrows(x0=sider_ses_by_onco_smry$l95, x1=sider_ses_by_onco_smry$u95, y0=sider_ses_by_onco_smry$y, code=3, angle=90, length=0.03, col=sider_ses_by_onco_smry$color, lwd=2)
end_of_figure_s2 = dev.off()

}



######
# FIGURE 2
######

tell_user('done.\nCreating Figure 2...')

binned_analysis = function(dsetable, bincolumn, binmins, binmaxs, binlabels, flip) {
  
  stopifnot(bincolumn %in% colnames(dsetable))
  
  forest = tibble(labels=binlabels,
                  min=binmins,
                  max=binmaxs,
                  mean=as.numeric(NA),
                  l95=as.numeric(NA),
                  u95=as.numeric(NA),
                  ntot=as.numeric(NA),
                  numerator=as.numeric(NA),
                  denominator=as.numeric(NA),
                  fracdisp = as.character(NA)) %>%
    mutate(y = max(row_number()) - row_number() + 1)
  
  # "flip" is for variables only defined for observed SEs. for instance, SIDER frequency
  # it flips the analysis so we look at assoc/obs rather than obs/assoc
  # and to create bins, instead of subsetting the denominator, we create conditional "observed" variables
  if (flip) {
    
    for (i in 1:nrow(forest)) {
      
      if (is.na(forest$min[i])) {
        dsetable$conditional_observed = dsetable$observed & is.na(dsetable[,bincolumn])
      } else {
        dsetable$conditional_observed = dsetable$observed & !is.na(dsetable[,bincolumn]) & dsetable[,bincolumn] > forest$min[i] & dsetable[,bincolumn] <= forest$max[i]
      }
      
      dse_subs = dsetable # do not subset - we use conditional_observed instead
      
      ctable = table(dse_subs[dse_subs$genetic_insight !='none' & !dse_subs$sim_indic,c('sim_assoc','conditional_observed')])
      forest$ntot[i] = sum(ctable)
      fobj = careful_fisher_test(ctable)
      forest$mean[i] = as.numeric(fobj$estimate)
      forest[i,c('l95','u95')] = as.list(fobj$conf.int)
      forest$numerator[i] = fobj$yy
      forest$denominator[i] = fobj$tctot
      
    }
    
    logit_model = glm(sim_assoc ~ get(bincolumn), data=dsetable, family='binomial')
    logit_model$call = paste0('sim_assoc ~ ',bincolumn)
    names(logit_model$coefficients)[2] = bincolumn
    
  } else {
    
    for (i in 1:nrow(forest)) {
      
      if (is.na(forest$min[i])) {
        dse_subs = dsetable[is.na(dsetable[,bincolumn]),]
      } else {
        dse_subs = dsetable[dsetable[,bincolumn] > forest$min[i] & dsetable[,bincolumn] <= forest$max[i],]
      }
      
      ctable = table(dse_subs[dse_subs$genetic_insight !='none' & !dse_subs$sim_indic,c('sim_assoc','observed')])
      forest$ntot[i] = sum(ctable)
      fobj = careful_fisher_test(ctable)
      forest$mean[i] = as.numeric(fobj$estimate)
      forest[i,c('l95','u95')] = fobj$conf.int
      forest$numerator[i] = fobj$yy
      forest$denominator[i] = fobj$trtot
      
      # set up logistic model
      logit_data = dsetable[dsetable$genetic_insight != 'none' & !dsetable$sim_indic,]
      logit_model = glm(conditional_observed ~ get(bincolumn) * sim_assoc, data=logit_data, family='binomial')
    }
    
  }
  
  forest$fracdisp = paste0(formatC(forest$numerator,format='d',big.mark=','), '/', formatC(forest$denominator,format='d',big.mark=','))
  
  return( list(forest=forest, model=logit_model))
  
}

make_figure_2 = TRUE
if (make_figure_2) {
  resx=300
  png('display_items/figure-2.png',width=3.25*resx,height=3.0*resx,res=resx)
  
  layout_matrix = matrix(1:3, nrow=3, byrow=T)
  layout(layout_matrix, heights=c(1,1,1))
  
  par(mar = c(0.5, 10, 1.5, 6))
  
  numfreq_list = binned_analysis(dse_info,
                                 bincolumn='numfreq', 
                                 binmins=c(1e-6,.001,.01,.1), 
                                 binmaxs=c(.001,.01,.1,1), 
                                 binlabels = c('<0.1%','0.1%-1%','1-10%','10%+'),
                                 flip=T)

  write_supp_table(numfreq_list[['forest']], 'Binned analysis of numerical SE frequency.')
  summary(numfreq_list[['model']])$coefficients %>%
    as_tibble(rownames='variable') %>%
    clean_names() -> numfreq_model
  write_supp_table(numfreq_model, 'Logit model coefficients for numerical SE frequency.')
    
  xlims = c(0,5)
  xats = 0:10/2
  xbigs = 0:5
  ylims = range(numfreq_list[['forest']]$y) + c(-0.5,0.5)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  abline(v=1, lwd=0.5, lty=3)
  # axis(side=1, at=xbigs, lwd=0, line=-0.5, cex.axis=0.8)
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  mtext(side=2, line=0.25, las=2, at=numfreq_list[['forest']]$y, text=numfreq_list[['forest']]$labels, cex=0.8)
  mtext(side=2, line=6, text='numerical\nfrequency')
  points(x=numfreq_list[['forest']]$mean, y=numfreq_list[['forest']]$y, pch=19)
  segments(x0=numfreq_list[['forest']]$l95, x1=numfreq_list[['forest']]$u95, y0=numfreq_list[['forest']]$y, lwd=1.5)
  mtext(side=4, line=0.25, at=numfreq_list[['forest']]$y, text=numfreq_list[['forest']]$fracdisp, las=2, cex=0.8)
  mtext(side=4, line=0.25, at=max(ylims)+0.5, text='assoc/obs', font=2, cex=0.8, las=2)
  
  summary(numfreq_list[['model']])
  
  
  wordfreq_list = binned_analysis(dse_info,
                                  bincolumn='wordrank', 
                                  binmins=hubal_day$wordrank-0.1, 
                                  binmaxs=hubal_day$wordrank+0.1, 
                                  binlabels = hubal_day$term,
                                  flip=T)
  par(mar = c(0.5, 10, 0.5, 6))
  xlims = c(0,5)
  xats = 0:10/2
  xbigs = 0:5
  ylims = range(wordfreq_list[['forest']]$y) + c(-0.5,0.5)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  abline(v=1, lwd=0.5, lty=3)
  # axis(side=1, at=xbigs, lwd=0, line=-0.5, cex.axis=0.8)
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  mtext(side=2, line=0.25, las=2, at=wordfreq_list[['forest']]$y, text=wordfreq_list[['forest']]$labels, cex=0.8)
  mtext(side=2, line=6, text='frequency\nword')
  points(x=wordfreq_list[['forest']]$mean, y=wordfreq_list[['forest']]$y, pch=19)
  segments(x0=wordfreq_list[['forest']]$l95, x1=wordfreq_list[['forest']]$u95, y0=wordfreq_list[['forest']]$y, lwd=1.5)
  mtext(side=4, line=0.25, at=wordfreq_list[['forest']]$y, text=wordfreq_list[['forest']]$fracdisp, las=2, cex=0.8)
  #mtext(side=4, line=0.25, at=max(ylims)+0.5, text='assoc/obs', font=2, cex=0.8, las=2)
  
  
  # see https://stackoverflow.com/questions/57297771/interpretation-of-l-q-c-4-for-logistic-regression
  # for interpretation of the binomial logit with an ordered (ordinal) independent variable
  # remember that polr(Hess=T) ordinal logits are for ordinal dependent variables, that's why they are not used here
  dse_info$wordrank_ordinal = factor(dse_info$wordrank, ordered=T)
  wordfreq_logit_model = glm(sim_assoc ~ wordrank_ordinal, data=dse_info, family='binomial')
  summary(wordfreq_logit_model) 
  # the linear term .L trends in the above. try also modeling it with only the linear term
  # which we can do by just feeding in rank as an integer instead of an ordered factor:
  wordfreq_logit_model_linearonly = glm(sim_assoc ~ wordrank, data=dse_info, family='binomial')
  summary(wordfreq_logit_model_linearonly) 

  
  write_supp_table(wordfreq_list[['forest']], 'Binned analysis of SE frequency terms.')
  summary(wordfreq_logit_model)$coefficients %>%
    as_tibble(rownames='variable') %>%
    clean_names() -> wordfreq_model_all
  write_supp_table(wordfreq_model_all, 'Logit model coefficients for SE frequency terms, ordinal model.')
  summary(wordfreq_logit_model_linearonly)$coefficients %>%
    as_tibble(rownames='variable') %>%
    clean_names() -> wordfreq_model_linearonly
  write_supp_table(wordfreq_model_linearonly, 'Logit model coefficients for SE frequency terms, linear term only.')
  
  
  
  
  placebo_list = binned_analysis(dse_info,
                                 bincolumn='placebo', 
                                 binmins=c(0,1)-0.1, 
                                 binmaxs=c(0,1)+0.1, 
                                 binlabels = c('non-placebo','placebo'),
                                 flip=T)
  
  # redo logit
  placebo_logit_model = glm(sim_assoc ~ placebo, data=dse_info, family='binomial')
  summary(placebo_logit_model) 
  
  write_supp_table(placebo_list[['forest']], 'Binned analysis of placebo status.')
  summary(placebo_logit_model)$coefficients %>%
    as_tibble(rownames='variable') %>%
    clean_names() -> placebo_model
  write_supp_table(placebo_model, 'Logit model coefficients for placebo status.')
  
  
  
  par(mar=c(3, 10, 0.5, 6))
  xlims = c(0,5)
  xats = 0:10/2
  xbigs = 0:5
  ylims = range(placebo_list[['forest']]$y) + c(-0.5,0.5)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, lwd=0, line=-0.75, cex.axis=0.8)
  abline(v=1, lwd=0.5, lty=3)
  mtext(side=1, line=1.75, text='OR')
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  mtext(side=2, line=0.25, las=2, at=placebo_list[['forest']]$y, text=placebo_list[['forest']]$labels, cex=0.8)
  mtext(side=2, line=6, text='evidence\nbasis')
  points(x=placebo_list[['forest']]$mean, y=placebo_list[['forest']]$y, pch=19)
  segments(x0=placebo_list[['forest']]$l95, x1=placebo_list[['forest']]$u95, y0=placebo_list[['forest']]$y, lwd=1.5)
  mtext(side=4, line=0.25, at=placebo_list[['forest']]$y, text=placebo_list[['forest']]$fracdisp, las=2, cex=0.8)
  #mtext(side=4, line=0.25, at=max(ylims)+0.5, text='assoc/obs', font=2, cex=0.8, las=2)
  
  end_of_figure_2 = dev.off()
}











######
# FIGURE 3
######

tell_user('done.\nCreating Figure 3...')


make_figure_3 = TRUE

if (make_figure_3) {
  resx=300
  png('display_items/figure-3.png',width=6.5*resx,height=3.0*resx,res=resx)
  
  layout_matrix = matrix(c(1:3,7,4:6,7), nrow=2, byrow=T)
  layout(layout_matrix, widths=c(0.6, 1, 1, 1.5))
  panel = 1
  
  ### PPV and OR by specificity (nee arcanity)
  dse_info %>%
    filter(observed) %>%
    group_by(se_mesh_id, se_mesh_term) %>%
    summarize(.groups='keep', n_drugs=length(unique(drug_name))) %>%
    ungroup() %>%
    mutate(log10bin = case_when(n_drugs==1 ~ '1',
                                n_drugs > 1 & n_drugs < 10 ~ '2-9',
                                n_drugs >= 10 & n_drugs < 100 ~ '10-99',
                                n_drugs >= 100 ~ '100+')) -> se_drug_counts
  
  specificity_bins = tibble(bin_name = c('1','2-9','10-99','100+'),
                            bin_order = c(1,2,3,4),
                            y=c(4,3,2,1),
                            or=numeric(4),
                            or_l95=numeric(4),
                            or_u95=numeric(4),
                            ntot=integer(4),
                            nse=integer(4),
                            sim_obs=integer(4),
                            sim_not=integer(4),
                            dis_obs=integer(4),
                            dis_not=integer(4),
                            baserate = numeric(4),
                            baserate_l95 = numeric(4),
                            baserate_u95 = numeric(4),
                            ppv = numeric(4),
                            ppv_l95 = numeric(4),
                            ppv_u95 = numeric(4),
                            npv = numeric(4),
                            npv_l95 = numeric(4),
                            npv_u95 = numeric(4))
  
  
  for (i in 1:nrow(specificity_bins)) {
    dse_info %>%
      filter(se_mesh_id %in% se_drug_counts$se_mesh_id[se_drug_counts$log10bin == specificity_bins$bin_name[i]]) -> dse_subs
    ctable = table(dse_subs[dse_subs$genetic_insight !='none' & !dse_subs$sim_indic,c('sim_assoc','observed')])
    specificity_bins$ntot[i] = sum(ctable)
    specificity_bins$nse[i] = length(unique(dse_subs$se_mesh_id))
    fobj = careful_fisher_test(ctable)
    specificity_bins$or[i] = as.numeric(fobj$estimate)
    specificity_bins[i,c('or_l95','or_u95')] = as.list(fobj$conf.int)
    specificity_bins$sim_obs[i] = ctable['TRUE','TRUE']
    specificity_bins$sim_not[i] = ctable['TRUE','FALSE']
    specificity_bins$dis_obs[i] = ctable['FALSE','TRUE']
    specificity_bins$dis_not[i] = ctable['FALSE','FALSE']
    baserate_binom = binom.confint(x=specificity_bins$sim_obs[i] + specificity_bins$dis_obs[i], n=sum(ctable), method='wilson')
    specificity_bins[i,c('baserate','baserate_l95','baserate_u95')] = as.list(baserate_binom[,c('mean','lower','upper')])
    ppv_binom = binom.confint(x=specificity_bins$sim_obs[i], n=specificity_bins$sim_obs[i] + specificity_bins$sim_not[i], method='wilson')
    specificity_bins[i,c('ppv','ppv_l95','ppv_u95')] = as.list(ppv_binom[,c('mean','lower','upper')])
    npv_binom = binom.confint(x=specificity_bins$dis_not[i], n=specificity_bins$dis_not[i] + specificity_bins$dis_obs[i], method='wilson')
    specificity_bins[i,c('npv','npv_l95','npv_u95')] = as.list(npv_binom[,c('mean','lower','upper')])
    
  }
  specificity_bins$label = specificity_bins$bin_name
  
  # layout_matrix = c(1:4, nrow=1, byrow=T)
  # layout(layout_matrix, widths=c(0.6, 1, 1, 1))
  ylims = range(specificity_bins$y) + c(-0.5, 0.5)
  par(mar=c(3,0,2,0.25))
  plot(NA, NA, xlim=c(0,1), ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  mtext(side=4, adj=1, las=2, at=specificity_bins$y, text=specificity_bins$bin_name, cex=0.7)
  mtext(side=2, line=-2, text='N drugs')
  xlims = c(0, 3.75)
  par(mar=c(3,0,2,1))
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  xats = 0:8/2
  xbigs = 0:4
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, labels=xbigs, lwd=0, line=-0.5)
  mtext(side=1, line=1.6, text='OR', cex=0.75)
  abline(v=1, lty=3, lwd=0.5)
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  points(specificity_bins$or, specificity_bins$y, pch=19)
  segments(x0=specificity_bins$or_l95, x1=specificity_bins$or_u95, y0=specificity_bins$y, lwd=1.5)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
  panel = panel + 1
  
  par(mar=c(3,0,2,1.5))
  xlims = c(0.001, 1)
  xats = rep(1:9,4) * rep(10^(-3:0),each=9)
  xbigs = 10^(-3:0)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='x')
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, labels=percent(xbigs), lwd=0, line=-0.5)
  par(xpd=T)
  legend(x=min(xlims), y=min(ylims)-0.75, c('base rate','PPV'), pch=c(1,19), horiz = T, bty='n', cex=1)
  par(xpd=F)
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  points(specificity_bins$baserate, specificity_bins$y, pch=1)
  points(specificity_bins$ppv, specificity_bins$y, pch=19)
  segments(x0=specificity_bins$ppv_l95, x1=specificity_bins$ppv_u95, y0=specificity_bins$y, lwd=1.5)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
  panel = panel + 1
  
  write_supp_table(specificity_bins, 'Binned analysis of SEs by drug specificity (number of drugs where the SE is observed).')
  
  quartiles = quantile(sider_ses$max_severity,0:4/4, na.rm=T)
  severity_quartiles = tibble(labels=c('null','0-25th','26-50th','50-75th','76-100th'),
                              min=c(NA,0.00,quartiles[2:4]),
                              max=c(NA,quartiles[2:5]),
                              y=5:1,
                              or=numeric(5),
                              or_l95=numeric(5),
                              or_u95=numeric(5),
                              ntot=integer(5),
                              nse=integer(5),
                              sim_obs=integer(5),
                              sim_not=integer(5),
                              dis_obs=integer(5),
                              dis_not=integer(5),
                              baserate = numeric(5),
                              baserate_l95 = numeric(5),
                              baserate_u95 = numeric(5),
                              ppv = numeric(5),
                              ppv_l95 = numeric(5),
                              ppv_u95 = numeric(5),
                              npv = numeric(5),
                              npv_l95 = numeric(5),
                              npv_u95 = numeric(5))
  
  for (i in 1:nrow(severity_quartiles)) {
    if (is.na(severity_quartiles$min[i])) {
      dse_subs = dse_info %>% filter(is.na(max_severity))
    } else {
      dse_subs = dse_info %>% filter(max_severity > severity_quartiles$min[i] & max_severity <= severity_quartiles$max[i])
    }
    ctable = table(dse_subs[dse_subs$genetic_insight !='none' & !dse_subs$sim_indic,c('sim_assoc','observed')])
    severity_quartiles$ntot[i] = sum(ctable)
    fobj = careful_fisher_test(ctable)
    severity_quartiles$or[i] = as.numeric(fobj$estimate)
    severity_quartiles[i,c('or_l95','or_u95')] = as.list(fobj$conf.int)
    severity_quartiles$sim_obs[i] = ctable['TRUE','TRUE']
    severity_quartiles$sim_not[i] = ctable['TRUE','FALSE']
    severity_quartiles$dis_obs[i] = ctable['FALSE','TRUE']
    severity_quartiles$dis_not[i] = ctable['FALSE','FALSE']
    baserate_binom = binom.confint(x=severity_quartiles$sim_obs[i] + severity_quartiles$dis_obs[i], n=sum(ctable), method='wilson')
    severity_quartiles[i,c('baserate','baserate_l95','baserate_u95')] = as.list(baserate_binom[,c('mean','lower','upper')])
    ppv_binom = binom.confint(x=severity_quartiles$sim_obs[i], n=severity_quartiles$sim_obs[i] + severity_quartiles$sim_not[i], method='wilson')
    severity_quartiles[i,c('ppv','ppv_l95','ppv_u95')] = as.list(ppv_binom[,c('mean','lower','upper')])
    npv_binom = binom.confint(x=severity_quartiles$dis_not[i], n=severity_quartiles$dis_not[i] + severity_quartiles$dis_obs[i], method='wilson')
    severity_quartiles[i,c('npv','npv_l95','npv_u95')] = as.list(npv_binom[,c('mean','lower','upper')])
  }
  
  dse_info %>%
    filter(!is.na(max_severity) & genetic_insight != 'none' & !sim_indic) -> severity_logit_data
  severity_logit = glm(observed ~ max_severity * sim_assoc, data=severity_logit_data, family='binomial')
  
  
  ylims = range(severity_quartiles$y) + c(-0.5, 0.5)
  par(mar=c(3,0,2,0.25))
  plot(NA, NA, xlim=c(0,1), ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  mtext(side=4, adj=1, las=2, at=severity_quartiles$y, text=severity_quartiles$labels, cex=0.7)
  mtext(side=2, line=-2, text='severity quartile')
  xlims = c(0, 3.75)
  par(mar=c(3,0,2,1))
  xats = 0:8/2
  xbigs = 0:4
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, labels=xbigs, lwd=0, line=-0.5)
  mtext(side=1, line=1.6, text='OR', cex=0.75)
  abline(v=1, lty=3, lwd=0.5)
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  points(severity_quartiles$or, severity_quartiles$y, pch=19)
  segments(x0=severity_quartiles$or_l95, x1=severity_quartiles$or_u95, y0=severity_quartiles$y, lwd=1.5)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
  panel = panel + 1
  
  par(mar=c(3,0,2,1.5))
  xlims = c(0.001, 1)
  xats = rep(1:9,4) * rep(10^(-3:0),each=9)
  xbigs = 10^(-3:0)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='x')
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, labels=percent(xbigs), lwd=0, line=-0.5)
  par(xpd=T)
  legend(x=min(xlims), y=min(ylims)-0.75, c('base rate','PPV'), pch=c(1,19), horiz = T, bty='n', cex=1)
  par(xpd=F)
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  points(severity_quartiles$baserate, severity_quartiles$y, pch=1)
  points(severity_quartiles$ppv, severity_quartiles$y, pch=19)
  segments(x0=severity_quartiles$ppv_l95, x1=severity_quartiles$ppv_u95, y0=severity_quartiles$y, lwd=1.5)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
  panel = panel + 1
  
  write_supp_table(severity_quartiles, 'Binned analysis of SEs by severity quartile.')
  summary(severity_logit)$coefficients %>%
    as_tibble(rownames='variable') %>%
    clean_names() -> severity_logit_coefs
  write_supp_table(severity_logit_coefs, 'Logit model coefficients for severity analysis')
  
  dse_info %>%
    filter(genetic_insight != 'none') %>%
    filter(!sim_indic) %>%
    filter(observed) %>%
    group_by(se_mesh_id, se_mesh_term, max_severity) %>%
    summarize(.groups='keep', n_drugs = length(unique(drug_name))) %>%
    ungroup() -> severity_vs_specificity
  
  n_bins = 5
  severity_vs_specificity %>%
    filter(!is.na(max_severity)) %>%
    mutate(bin_max = ceiling(max_severity*n_bins)/n_bins) %>%
    mutate(bin_disp = paste0(formatC(bin_max-1/n_bins,format='f',digits=1), '-', formatC(bin_max,format='f',digits=1))) %>%
    group_by(bin_max, bin_disp) %>%
    summarize(.groups='keep',
              n_se = n(),
              min_n = min(n_drugs),
              q25_n = quantile(n_drugs,.25),
              q50_n = median(n_drugs),
              q75_n = quantile(n_drugs,.75),
              max_n = max(n_drugs)) %>%
    ungroup() %>%
    mutate(iqr = q75_n - q25_n) %>%
    mutate(minwhisk_n = as.numeric(NA),
           maxwhisk_n = as.numeric(NA)) -> severity_specificity_binned

  severity_vs_specificity$bin_max = ceiling(severity_vs_specificity$max_severity*n_bins)/n_bins
  for (i in 1:nrow(severity_specificity_binned)) {
    lowest_whisker = severity_specificity_binned$q25_n[i] - whisker_factor * severity_specificity_binned$iqr[i]
    highest_whisker = severity_specificity_binned$q75_n[i] + whisker_factor * severity_specificity_binned$iqr[i]
    severity_specificity_binned$minwhisk_n[i] = min(severity_vs_specificity$n_drugs[severity_vs_specificity$bin_max==severity_specificity_binned$bin_max[i] & severity_vs_specificity$n_drugs >= lowest_whisker],na.rm=T)
    severity_specificity_binned$maxwhisk_n[i] = max(severity_vs_specificity$n_drugs[severity_vs_specificity$bin_max==severity_specificity_binned$bin_max[i] & severity_vs_specificity$n_drugs <= highest_whisker],na.rm=T)
  }
    
  par(mar=c(3,3,2,1))
  xlims = c(0, 1)
  xbigs = 0:5/5
  ylims = c(0.7,550)
  yats = rep(1:9,4) * rep(10^(0:3),each=9)
  ybigs = 10^(0:3)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='y')
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  mtext(side=1, at=severity_specificity_binned$bin_max-1/n_bins/2, text=severity_specificity_binned$bin_disp, line=0.25, cex=0.6)
  mtext(side=1, line=1.75, text='severity quantile', cex=0.8)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=-0.25)
  boxwidth = 0.15
  rect(xleft=severity_specificity_binned$bin_max - 1/n_bins/2 - boxwidth/2, xright=severity_specificity_binned$bin_max - 1/n_bins/2 + boxwidth/2, ybottom=severity_specificity_binned$q25_n, ytop=severity_specificity_binned$q75_n, lwd=1.5, border='#000000', col=NA)
  segments(x0=severity_specificity_binned$bin_max - 1/n_bins/2 - boxwidth/2, x1=severity_specificity_binned$bin_max - 1/n_bins/2 + boxwidth/2, y0=severity_specificity_binned$q50_n, lwd=1.5, col='#000000')
  arrows(x0=severity_specificity_binned$bin_max - 1/n_bins/2, y0=severity_specificity_binned$q25_n, y1=severity_specificity_binned$min_n, code=2, length=0.05, angle=90, col='#000000')
  arrows(x0=severity_specificity_binned$bin_max - 1/n_bins/2, y0=severity_specificity_binned$q75_n, y1=severity_specificity_binned$max_n, code=2, length=0.05, angle=90, col='#000000')
  mtext(side=2, line=2.25, text='N drugs', cex=0.8)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
  panel = panel + 1
  
  write_supp_table(severity_specificity_binned, 'SE drug specificity versus severity bin.')
  
  m = lm(log(n_drugs) ~ max_severity, data=severity_vs_specificity)
  summary(m)$coefficients %>%
    as_tibble(rownames='variable') %>%
    clean_names() -> sev_spec_linear_coefs
  write_supp_table(sev_spec_linear_coefs, 'Linear model coefficients for SE severity vs. drug specificity.')
  
  end_of_figure_3 = dev.off()
}















######
# FIGURE 4
######

tell_user('done.\nCreating Figure 4...')

areas = read_tsv('data/mesh/areas.tsv', col_types=cols()) 
se_topl_match = read_tsv('data/synthesis/se_topl_match.tsv', col_types=cols())
topl_maps = read_tsv('data/mesh/mesh_2023_topl_maps.tsv', col_types=cols())

# mean(se_topl_match$topl %in% topl_maps$topl) # 100%
# which(!(se_topl_match$topl %in% topl_maps$topl)) # none
# mean(dse_info$se_mesh_id %in% se_topl_match$mesh_id) # 100% -OK 
# mean(dse_info$assoc_mesh_id %in% se_topl_match$mesh_id) # 95.6% - rest are NA - best we can do
# head(which(!(dse_info$assoc_mesh_id %in% se_topl_match$mesh_id) & !is.na(dse_info$assoc_mesh_id)))

area_modes = tibble(filter_by=c('SE','SE'),
                    filter_mode = c('soft','hard'))

areas %>%
  add_row(.before=1, topl='ALL',area='all',color='#000000',filter='none') %>%
  mutate(y = max(row_number()) + 1 - row_number()) %>%
  mutate(sim_obs = as.integer(NA),
         sim_not = as.integer(NA),
         dis_obs = as.integer(NA),
         dis_not = as.integer(NA),
         or = as.numeric(NA),
         p = as.numeric(NA),
         or_l95 = as.numeric(NA),
         or_u95 = as.numeric(NA),
         baserate     = as.numeric(NA),
         baserate_l95 = as.numeric(NA),
         baserate_u95 = as.numeric(NA),
         median_severity = as.numeric(NA),
         q25_severity = as.numeric(NA),
         q75_severity = as.numeric(NA),
         minwhisk_severity = as.numeric(NA),
         maxwhisk_severity = as.numeric(NA),
         ppv = as.numeric(NA),
         ppv_l95 = as.numeric(NA),
         ppv_u95 = as.numeric(NA),
         npv = as.numeric(NA),
         npv_l95 = as.numeric(NA),
         npv_u95 = as.numeric(NA)) -> area_stats_one
  
crossing(area_stats_one, area_modes) %>%
  relocate(c(filter_by, filter_mode), .after=filter) -> area_stats

skip_intensive = F
for (i in 1:nrow(area_stats)) {

  if (area_stats$topl[i] == 'ALL') {
    this_area_ses = unique(dse_info$se_mesh_id)
    this_area_assocs = unique(dse_info$assoc_mesh_id)
    this_area_indics = unique(pp_launched_indic$indication_mesh_id)
  } else {
    this_area_topls = topl_maps$topl[topl_maps$map_to == area_stats$topl[i]]
    this_area_ses = se_topl_match$mesh_id[se_topl_match$topl %in% this_area_topls]
    this_area_assocs = se_topl_match$mesh_id[se_topl_match$topl %in% this_area_topls]
    this_area_indics = intersect(pp_launched_indic$indication_mesh_id, se_topl_match$mesh_id[se_topl_match$topl %in% this_area_topls])
  }
  
  if (area_stats$filter_by[i]=='SE') {
    dse_info %>%
      filter(se_mesh_id %in% this_area_ses) -> current_dse
  }

  if (area_stats$filter_mode[i]=='hard') {
    
    # "hard filter mode" - remove any drug with any indication in the same MeSH tree branch
    pp_drugs_this_area = pp_launched_indic$ppid[pp_launched_indic$indication_mesh_id %in% this_area_indics]
    sider_drugs_this_area = sider_drugs$drug_name[sider_drugs$ppid %in% pp_drugs_this_area]
    current_dse %>%
      filter(!(drug_name %in% sider_drugs_this_area)) -> current_dse
  }
  
  ctable = table(current_dse[current_dse$genetic_insight != 'none' &
                                 !current_dse$sim_indic, c('sim_assoc', 'observed')])
  
  if (prod(dim(ctable)) < 4) {
    next
  }
  area_stats$sim_obs[i] = ctable['TRUE','TRUE']
  area_stats$sim_not[i] = ctable['TRUE','FALSE']
  area_stats$dis_obs[i] = ctable['FALSE','TRUE']
  area_stats$dis_not[i] = ctable['FALSE','FALSE']
  baserate_binom = binom.confint(x=area_stats$sim_obs[i] + area_stats$dis_obs[i], n=sum(ctable), method='wilson')
  area_stats[i,c('baserate','baserate_l95','baserate_u95')] = as.list(baserate_binom[,c('mean','lower','upper')])
  area_stats$median_severity[i] = median(current_dse$max_severity[current_dse$genetic_insight != 'none' & !current_dse$sim_indic], na.rm=T)
  area_stats$q25_severity[i] = quantile(current_dse$max_severity[current_dse$genetic_insight != 'none' & !current_dse$sim_indic], 0.25, na.rm=T)
  area_stats$q75_severity[i] = quantile(current_dse$max_severity[current_dse$genetic_insight != 'none' & !current_dse$sim_indic], 0.75, na.rm=T)
  this_area_iqr = area_stats$q75_severity[i] - area_stats$q25_severity[i]
  lowest_whisker = area_stats$q25_severity[i] - this_area_iqr*whisker_factor
  highest_whisker = area_stats$q75_severity[i] + this_area_iqr*whisker_factor
  area_stats$minwhisk_severity[i] = min(current_dse$max_severity[current_dse$genetic_insight != 'none' & !current_dse$sim_indic & current_dse$max_severity >= lowest_whisker], na.rm=T)
  area_stats$maxwhisk_severity[i] = max(current_dse$max_severity[current_dse$genetic_insight != 'none' & !current_dse$sim_indic & current_dse$max_severity <= highest_whisker], na.rm=T)
  fobj = fisher.test(ctable)
  area_stats$or[i] = fobj$estimate
  area_stats$p[i] = fobj$p.value
  area_stats[i,c('or_l95','or_u95')] = as.list(as.numeric(fobj$conf.int))
  ppv_binom = binom.confint(x=area_stats$sim_obs[i], n=area_stats$sim_obs[i] + area_stats$sim_not[i], method='wilson')
  area_stats[i,c('ppv','ppv_l95','ppv_u95')] = as.list(ppv_binom[c('mean','lower','upper')])
  npv_binom = binom.confint(x=area_stats$dis_not[i], n=area_stats$dis_obs[i] + area_stats$dis_not[i], method='wilson')
  area_stats[i,c('npv','npv_l95','npv_u95')] = as.list(npv_binom[c('mean','lower','upper')])
}
if (!skip_intensive) {
  write_supp_table(area_stats, 'Breakdown by MeSH area.')
  write_tsv(area_stats, 'ignore/area_stats.tsv')
} else {
  area_stats = read_tsv('ignore/area_stats.tsv', col_types=cols())
}

area_stats %>%
  filter(filter_by=='SE' & filter_mode=='soft') %>%
  mutate(isnt_all = area != 'all') %>%
  arrange(isnt_all, desc(or)) %>%
  mutate(y = max(row_number()) - row_number() + 1) -> area_ys

area_stats_one$y = area_ys$y[match(area_stats_one$area, area_ys$area)]
area_stats$y = area_ys$y[match(area_stats$area, area_ys$area)]

# note that if anyone complains the plot doesn't show individual points, you can check the number of points in each bin:
# sum(!is.na(current_dse$max_severity))
# for 'all' it is 849,366. remember that in order to remove drug-SE tuples with similar indication, parallel with other
# analyses in this figure, we do the quantiling of severity in drug-SE mode rather than in unique SE mode, so there are
# hundreds of thousands of points rather than thousands.





area_stats %>%
  filter(filter_mode=='soft' & area != 'all') -> area_data
cmh_array = array(data=rep(0,2*2*(nrow(area_data))), dim = c(2,2,(nrow(area_data))))
for (i in 1:nrow(areas)) {
  cmh_array[1,1,i] = as.numeric(area_data[i,'sim_obs'])
  cmh_array[1,2,i] = as.numeric(area_data[i,'sim_not'])
  cmh_array[2,1,i] = as.numeric(area_data[i,'dis_obs'])
  cmh_array[2,2,i] = as.numeric(area_data[i,'dis_not'])
}
cmh_result = suppressMessages(cmh.test(cmh_array))
extract_cmh_p = function(cmh_obj) {
  cmh_minp = 1e-15 # 1 - pchisq(q=70, df=1) = 1.1e-16, this is the lowest non-zero p value it ever returns, after that all 0. so we can just say P < 1e-15
  chisq_mh = cmh_obj$parameter['CMH statistic']
  df = cmh_obj$parameter['df']
  p = 1 - pchisq(q=chisq_mh, df=df)
  return (max(p,cmh_minp))
}
cmh_p = extract_cmh_p(cmh_result)
write_stats('CMH test for heterogeneity between MeSH areas: P = ', formatC(cmh_p, format='e', digits=1))

make_figure_4 = TRUE
if (make_figure_4) {
  
  resx=300
  png(paste0('display_items/figure-4.png'),width=8.5*resx,height=4*resx,res=resx)
  
  this_by = 'SE'
  this_mode = 'soft'
  
  layout_matrix = matrix(c(1:4,5, 1:4,6),nrow=2, byrow=T)
  layout(layout_matrix, widths=c(.5, 1, .75, .75, 1))
  panel = 1
  
  xlims=c(0,8)
  ylims=range(area_stats_one$y) + c(-0.5, 0.5)
  
  par(mar=c(4,0,2.5,0))
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  mtext(side=4, at=area_stats_one$y, text=area_stats_one$area, col=area_stats_one$color, las=2, line=0.25, adj=1, cex=0.9)
  
  area_stats %>%
    filter(filter_by==this_by & filter_mode==this_mode) -> subs
  
  par(mar=c(4,1,2.5,7))
  xats = 0:16/2
  xbigs = 0:8
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, labels=xbigs, lwd=0, line=-0.5, cex.axis=0.8)
  mtext(side=1, line=1.6, text='OR', cex=0.75)
  #par(xpd=T)
  #legend(x=min(xlims), y=min(ylims)-0.75, c('OR'), pch=c(19), horiz = T, bty='n', cex=1.25)
  #par(xpd=F)
  abline(v=1, lty=3)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  points(x=subs$or, y=subs$y, col=subs$color, pch=19)
  segments(x0=subs$or_l95, x1=pmin(subs$or_u95,max(xlims)), y0=subs$y, lwd=2, col=subs$color)
  par(xpd=T)
  mtext(side=4, at=max(ylims)+1, text='obs/assoc', las=2, font=2, cex=0.8)
  par(xpd=F)
  mtext(side=4, at=subs$y, text=paste0(subs$sim_obs,'/',formatC(subs$sim_obs + subs$sim_not, format='d', big.mark=',')), col=subs$color, line=0.5, las=2, cex=0.8)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
  panel = panel + 1
  
  xlims = c(0, 0.30)
  xats = 0:30/100
  xbigs = 0:20/20
  par(mar=c(4,1,2.5,1))
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, labels=percent(xbigs), lwd=0, line=-0.5)
  par(xpd=T)
  legend(x=min(xlims), y=min(ylims)-0.75, c('base rate','PPV'), pch=c(1,19), horiz = T, bty='n', cex=1.25)
  par(xpd=F)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  points(x=subs$ppv, y=subs$y, col=subs$color, pch=19)
  segments(x0=subs$ppv_l95, x1=pmin(subs$ppv_u95,max(xlims)), y0=subs$y, lwd=2, col=subs$color)
  points(x=subs$baserate, y=subs$y, pch=21, bg='#FFFFFF', col=subs$color)
  barwidth = 0.3
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
  panel = panel + 1

  ylims = range(subs$y) + c(-0.5, 0.5)
  xlims = c(0, 1.05)
  xats = 0:20/20
  xbigs = 0:10/10
  par(mar=c(4,1,2.5,1))
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, labels=xbigs, lwd=0, line=-0.5)
  mtext(side=1, line=1.75, text='med/IQR severity', cex=0.8)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  boxwidth = 0.3
  rect(xleft=subs$q25_severity, xright=subs$q75_severity, ybottom=subs$y-boxwidth, ytop=subs$y+boxwidth, lwd=1.5, border=subs$color, col=NA)
  arrows(x0=subs$q25_severity, x1=subs$minwhisk_severity, y0=subs$y, code=2, length=0.05, angle=90, col=subs$color)
  arrows(x0=subs$q75_severity, x1=subs$maxwhisk_severity, y0=subs$y, code=2, length=0.05, angle=90, col=subs$color)
  segments(x0=subs$median_severity, y0=subs$y-boxwidth, y1=subs$y+boxwidth, lwd=1.5, col=subs$color)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
  panel = panel + 1
  
  area_stats %>%
    filter(filter_by==this_by & filter_mode==this_mode & area != 'all') -> area_scatter
  spearman_obj = suppressWarnings(cor.test(area_scatter$or, area_scatter$ppv, method='spearman'))
  write_stats('Spearman test OR vs. PPV across areas (soft filter mode), rho = ',as.numeric(spearman_obj$estimate),'P = ',spearman_obj$p.value)
  subs = area_scatter
  par(mar=c(4,4,2.5,1))
  xlims = c(0, 8)
  ylims = c(0, 0.30)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=0:max(xlims))
  mtext(side=1, line=2.5, text='OR')
  abline(v=1, lty=3)
  axis(side=2, at=0:6/20, labels=percent(0:6/20), las=2)
  mtext(side=2, line=2.75, text='PPV')
  points(x=subs$or, y=subs$ppv, col=subs$color, pch=19, cex=1.5)
  segments(x0=subs$or_l95, x1=pmin(subs$or_u95,max(xlims)), y0=subs$ppv, lwd=1, col=subs$color)
  segments(x0=subs$or, y0=subs$ppv_l95, y1=subs$ppv_u95, lwd=1, col=subs$color)
  subs %>% 
    filter(ppv > .15 | or > 4) -> to_label
  text(x=to_label$or, y=to_label$ppv + 0.01, pos=ifelse(to_label$or > 6, 2, 4), labels=to_label$area, col=to_label$color, cex=0.7)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
  panel = panel + 1
  
  spearman_obj = suppressWarnings(cor.test(area_scatter$baserate, area_scatter$median_severity, method='spearman'))
  write_stats('Spearman test base rate vs. median severity across areas (soft filter mode), rho = ',as.numeric(spearman_obj$estimate),'P = ',spearman_obj$p.value)
  subs = area_scatter
  par(mar=c(4,4,2.5,1))
  xlims = c(0, 0.12)
  ylims = c(0.3, 0.8)
  xats = 0:100/100
  xbigs = 0:20/20
  yats = 0:20/20
  ybigs = 0:10/10
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, labels=percent(xbigs), lwd=0, line=-0.5)
  mtext(side=1, line=1.5, text='base rate')
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, labels=ybigs, lwd=0, line=-0.5, las=2)
  mtext(side=2, line=2.75, text='median severity')
  points(x=subs$baserate, y=subs$median_severity, col=subs$color, pch=19, cex=1.5)
  subs %>% 
    filter(median_severity > 0.6 | baserate > 0.1) -> to_label
  par(xpd=T)
  text(x=to_label$baserate, y=to_label$median_severity + 0.01, pos=ifelse(to_label$baserate > 0.1, 2, 4), labels=to_label$area, col=to_label$color, cex=0.7)
  par(xpd=F)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
  panel = panel + 1
  
  end_of_figure_4 = dev.off()
}




######
# FIGURE S3
######

tell_user('done.\nCreating Figure S3...')

make_figure_s3 = TRUE
if (make_figure_s3) {
  
  resx=300
  png(paste0('display_items/figure-s3.png'),width=8.5*resx,height=4*resx,res=resx)
  
  this_by = 'SE'
  this_mode = 'hard'
  
  area_stats %>%
    filter(area != 'all') %>%
    filter(filter_by==this_by & filter_mode==this_mode) -> subs
  
  layout_matrix = matrix(c(1:4,5, 1:4,6),nrow=2, byrow=T)
  layout(layout_matrix, widths=c(.5, 1, .75, .75, 1))
  panel = 1
  
  xlims=c(0,14)
  ylims=range(subs$y) + c(-0.5, 0.5)
  
  par(mar=c(4,0,2.5,0))
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  mtext(side=4, at=subs$y, text=subs$area, col=subs$color, las=2, line=0.25, adj=1, cex=0.8)
  
  par(mar=c(4,1,2.5,7))
  xats = 0:28/2
  xbigs = 0:14
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, labels=xbigs, lwd=0, line=-0.5, cex.axis=0.8)
  par(xpd=T)
  legend(x=min(xlims), y=min(ylims)-0.75, c('OR'), pch=c(19), horiz = T, bty='n', cex=1.25)
  par(xpd=F)
  abline(v=1, lty=3)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  points(x=subs$or, y=subs$y, col=subs$color, pch=19)
  segments(x0=subs$or_l95, x1=pmin(subs$or_u95,max(xlims)), y0=subs$y, lwd=2, col=subs$color)
  par(xpd=T)
  mtext(side=4, at=max(ylims)+1, text='obs/assoc', las=2, font=2, cex=0.8)
  par(xpd=F)
  mtext(side=4, at=subs$y, text=paste0(subs$sim_obs,'/',formatC(subs$sim_obs + subs$sim_not, format='d', big.mark=',')), col=subs$color, line=0.5, las=2, cex=0.8)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
  panel = panel + 1
  
  
  xlims = c(0, 0.30)
  xats = 0:30/100
  xbigs = 0:6/20
  par(mar=c(4,1,2.5,1))
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, labels=percent(xbigs), lwd=0, line=-0.5)
  par(xpd=T)
  legend(x=min(xlims), y=min(ylims)-0.75, c('base rate','PPV'), pch=c(1,19), horiz = T, bty='n', cex=1.25)
  par(xpd=F)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  points(x=subs$ppv, y=subs$y, col=subs$color, pch=19)
  segments(x0=subs$ppv_l95, x1=pmin(subs$ppv_u95,max(xlims)), y0=subs$y, lwd=2, col=subs$color)
  points(x=subs$baserate, y=subs$y, pch=21, bg='#FFFFFF', col=subs$color)
  barwidth = 0.3
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
  panel = panel + 1
  
  
  ylims = range(subs$y) + c(-0.5, 0.5)
  xlims = c(0.28, 0.88)
  xats = 0:20/20
  xbigs = 0:10/10
  par(mar=c(4,1,2.5,1))
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, labels=xbigs, lwd=0, line=-0.5)
  mtext(side=1, line=1.75, text='med/IQR severity', cex=0.8)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  boxwidth = 0.3
  rect(xleft=subs$q25_severity, xright=subs$q75_severity, ybottom=subs$y-boxwidth, ytop=subs$y+boxwidth, lwd=1.5, border=subs$color, col=NA)
  arrows(x0=subs$q25_severity, x1=subs$minwhisk_severity, y0=subs$y, code=2, length=0.05, angle=90, col=subs$color)
  arrows(x0=subs$q75_severity, x1=subs$maxwhisk_severity, y0=subs$y, code=2, length=0.05, angle=90, col=subs$color)
  segments(x0=subs$median_severity, y0=subs$y-boxwidth, y1=subs$y+boxwidth, lwd=1.5, col=subs$color)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
  panel = panel + 1
  
  
  area_stats %>%
    filter(filter_by==this_by & filter_mode==this_mode & area != 'all') -> area_scatter
  spearman_obj = suppressWarnings(cor.test(area_scatter$or, area_scatter$ppv, method='spearman'))
  write_stats('Spearman test OR vs. PPV across areas (hard filter mode), rho = ',as.numeric(spearman_obj$estimate),'P = ',spearman_obj$p.value)
  subs = area_scatter
  par(mar=c(4,4,2.5,1))
  xlims = c(0, 14)
  ylims = c(0, 0.30)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=0:max(xlims))
  mtext(side=1, line=2.5, text='OR')
  abline(v=1, lty=3)
  axis(side=2, at=0:6/20, labels=percent(0:6/20), las=2)
  mtext(side=2, line=2.75, text='PPV')
  points(x=subs$or, y=subs$ppv, col=subs$color, pch=19, cex=1.5)
  segments(x0=subs$or_l95, x1=pmin(subs$or_u95,max(xlims)), y0=subs$ppv, lwd=1, col=subs$color)
  segments(x0=subs$or, y0=subs$ppv_l95, y1=subs$ppv_u95, lwd=1, col=subs$color)
  subs %>% 
    filter(ppv > .15 | or > 5) -> to_label
  text(x=to_label$or, y=to_label$ppv + 0.01, pos=ifelse(to_label$or > 5 & to_label$ppv > .10, 2, 4), labels=to_label$area, col=to_label$color, cex=0.7)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
  panel = panel + 1
  
  
  
  spearman_obj = suppressWarnings(cor.test(area_scatter$baserate, area_scatter$median_severity, method='spearman'))
  write_stats('Spearman test base rate vs. median severity across areas (hard filter mode), rho = ',as.numeric(spearman_obj$estimate),'P = ',spearman_obj$p.value)
  subs = area_scatter
  par(mar=c(4,4,2.5,1))
  xlims = c(0, 0.12)
  ylims = c(0.3, 0.8)
  xats = 0:100/100
  xbigs = 0:20/20
  yats = 0:20/20
  ybigs = 0:10/10
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, labels=percent(xbigs), lwd=0, line=-0.5)
  mtext(side=1, line=1.5, text='base rate')
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, labels=ybigs, lwd=0, line=-0.5, las=2)
  mtext(side=2, line=2.75, text='median severity')
  points(x=subs$baserate, y=subs$median_severity, col=subs$color, pch=19, cex=1.5)
  subs %>% 
    filter(median_severity > 0.6 | baserate > 0.1) -> to_label
  text(x=to_label$baserate, y=to_label$median_severity + 0.01, pos=ifelse(to_label$baserate > 0.1, 2, 4), labels=to_label$area, col=to_label$color, cex=0.7)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.2)
  panel = panel + 1
  
  
  end_of_figure_s3 = dev.off()
  
}


####
# Supplementary tables grouping by SE and association
####

dse_info %>%
  filter(genetic_insight != 'none') %>%
  filter(!sim_indic) %>%
  filter(sim_assoc) %>%
  filter(observed) %>%
  group_by(assoc_mesh_id, assoc_mesh_term) %>%
  summarize(.groups='keep',
            genes = toString(unique(gene)),
            drugs = toString(unique(drug_name)),
            side_effects = toString(unique(se_mesh_term))) %>%
  ungroup() %>%
  mutate(or = as.numeric(NA),
         or_l95 = as.numeric(NA),
         or_u95 = as.numeric(NA),
         fisher_p = as.numeric(NA),
         sim_obs = as.integer(NA),
         sim_not = as.integer(NA),
         dis_obs = as.integer(NA),
         dis_not = as.integer(NA)) -> assocs

get_ctable_cell = function(ctable, row, col) {
  if (row %in% rownames(ctable) & col %in% colnames(ctable)) {
    return (ctable[row,col])
  } else {
    return (0)
  }
}

for (i in 1:nrow(assocs)) {
  dse_subs = dse_info %>% 
    filter(assoc_mesh_id == assocs$assoc_mesh_id[i]) %>%
    filter(genetic_insight != 'none') %>%
    filter(!sim_indic)
  ctable = table(dse_subs[,c('sim_assoc','observed')])
  if (!('TRUE' %in% colnames(ctable))) {
    
  }
  fobj = careful_fisher_test(ctable)
  assocs$or[i] = as.numeric(fobj$estimate)
  assocs$or_l95[i] = as.numeric(fobj$conf.int[1])
  assocs$or_u95[i] = as.numeric(fobj$conf.int[2])
  assocs$fisher_p[i] = as.numeric(fobj$p.value)
  assocs$sim_obs[i] = get_ctable_cell(ctable, 'TRUE', 'TRUE')
  assocs$dis_obs[i] = get_ctable_cell(ctable, 'FALSE','TRUE')
  assocs$sim_not[i] = get_ctable_cell(ctable, 'TRUE', 'FALSE')
  assocs$dis_not[i] = get_ctable_cell(ctable, 'FALSE','FALSE')
}

assocs %>% 
  arrange(fisher_p) %>%
  relocate(genes, drugs, side_effects, .after=dis_not) -> assocs


write_supp_table(assocs, 'Enrichment statistics by GWAS association MeSH term.')








dse_info %>%
  filter(genetic_insight != 'none') %>%
  filter(!sim_indic) %>%
  filter(sim_assoc) %>%
  filter(observed) %>%
  group_by(se_mesh_id, se_mesh_term) %>%
  summarize(.groups='keep',
            genes = toString(unique(gene)),
            drugs = toString(unique(drug_name)),
            assocs = toString(unique(assoc_mesh_term))) %>%
  ungroup() %>%
  mutate(or = as.numeric(NA),
         or_l95 = as.numeric(NA),
         or_u95 = as.numeric(NA),
         fisher_p = as.numeric(NA),
         sim_obs = as.integer(NA),
         sim_not = as.integer(NA),
         dis_obs = as.integer(NA),
         dis_not = as.integer(NA)) -> se_stats

get_ctable_cell = function(ctable, row, col) {
  if (row %in% rownames(ctable) & col %in% colnames(ctable)) {
    return (ctable[row,col])
  } else {
    return (0)
  }
}

for (i in 1:nrow(se_stats)) {
  dse_subs = dse_info %>% 
    filter(se_mesh_id == se_stats$se_mesh_id[i]) %>%
    filter(genetic_insight != 'none') %>%
    filter(!sim_indic)
  ctable = table(dse_subs[,c('sim_assoc','observed')])
  fobj = careful_fisher_test(ctable)
  se_stats$or[i] = as.numeric(fobj$estimate)
  se_stats$or_l95[i] = as.numeric(fobj$conf.int[1])
  se_stats$or_u95[i] = as.numeric(fobj$conf.int[2])
  se_stats$fisher_p[i] = as.numeric(fobj$p.value)
  se_stats$sim_obs[i] = get_ctable_cell(ctable, 'TRUE', 'TRUE')
  se_stats$dis_obs[i] = get_ctable_cell(ctable, 'FALSE','TRUE')
  se_stats$sim_not[i] = get_ctable_cell(ctable, 'TRUE', 'FALSE')
  se_stats$dis_not[i] = get_ctable_cell(ctable, 'FALSE','FALSE')
}

se_stats %>% 
  arrange(fisher_p) %>%
  relocate(genes, drugs, assocs, .after=dis_not) -> se_stats


write_supp_table(se_stats, 'Enrichment statistics by side effect MeSH term.')


dse_info %>% 
  filter(genetic_insight=='none') %>% 
  group_by(se_mesh_id, se_mesh_term) %>% 
  summarize(.groups='keep', n_drugs_with_this_se=sum(observed)) %>% 
  arrange(desc(n_drugs_with_this_se)) -> lacking_genetic_insight

write_supp_table(lacking_genetic_insight, 'Side effects lacking genetic insight, by number of drugs ')

assoc %>%
  group_by(gene) %>%
  summarize(.groups='keep',
            n_associated_mesh_ids = length(unique(mesh_id))) %>%
  ungroup() -> assoc_traits_per_gene

dse_info %>%
  filter(observed & genetic_insight != 'none' & !sim_indic) %>%
  group_by(drug_name, gene) %>%
  summarize(.groups='keep',
            n_observed_ses = length(unique(se_mesh_id))) %>%
  ungroup() -> obs_ses_per_drug

assoc_traits_per_gene %>%
  right_join(obs_ses_per_drug, by='gene') %>%
  mutate(n_associated_mesh_ids = replace_na(n_associated_mesh_ids, 0)) -> assocs_vs_ses

cor_obj = cor.test(assocs_vs_ses$n_associated_mesh_ids, assocs_vs_ses$n_observed_ses, method='pearson')
write_stats_text("Pearson's correlation between a drug's number of side effects and the number of unique MeSH IDs to which it is genetically associated: ",
                 "rho = ",formatC(cor_obj$estimate,format='f',digits=2),
                 "P = ",formatC(cor_obj$p.value,format='f',digits=2))



dse_info %>%
  filter(genetic_insight!='none' & observed) %>%
  distinct(drug_name, se_mesh_id, indication_mesh_id, se_indic_similarity, sim_indic) %>%
  summarize(.groups='keep',
            without_sim_indic = sum(!sim_indic),
            with_sim_indic = sum(sim_indic),
            total_drugse_pairs = n(),
            proportion_sim_indic = mean(sim_indic)) %>%
  ungroup() -> sim_indic_stats

write_stats_text('At our 0.9 similarity threshold, ',
                 sim_indic_stats$with_sim_indic,'/',
                 sim_indic_stats$total_drugse_pairs,
                 ' (',percent(sim_indic_stats$proportion_sim_indic),
                 ') were removed.')


dse_info %>%
  left_join(sim, by=c('indication_mesh_id'='meshcode_a', 'assoc_mesh_id'='meshcode_b')) %>%
  rename(indic_assoc_sim = comb_norm) %>%
  mutate(indic_assoc_sim = replace_na(indic_assoc_sim, 0)) %>%
  mutate(indic_supported = indic_assoc_sim >= 0.8) -> dse_info_with_indic_assoc_sim

dse_info_with_indic_assoc_sim %>%
  distinct(drug_name, gene, assoc_mesh_id, indication_mesh_id, indic_supported) %>%
  group_by(indic_supported) %>%
  summarize(.groups='keep',
            n=n()) %>%
  ungroup() %>%
  mutate(proportion = n/sum(n)) -> prevalence_of_gensup

write_stats_text(prevalence_of_gensup$n[prevalence_of_gensup$indic_supported],"(",
                 percent(prevalence_of_gensup$n[prevalence_of_gensup$proportion]),")",
                 " of drug-gene-association-indication tuples are genetically supported")

dse_info_with_indic_assoc_sim %>%
  filter(genetic_insight != 'none' & observed) %>%
  summarize(.groups='keep',
            total_n = n(),
            n_indic_supported = sum(indic_supported),
            proportion_indic_supported = mean(indic_supported)) %>%
  ungroup() -> drugse_observed_supported


write_stats_text(drugse_observed_supported$n_indic_supported,"(",
                 percent(drugse_observed_supported$proportion_indic_supported),")",
                 " of observed drug-SE pairs have a supported indication.")

dse_info_with_indic_assoc_sim %>%
  filter(genetic_insight != 'none') %>%
  group_by(indic_supported) %>%
  summarize(.groups='keep',
            n_rows = n(),
            n_observed = sum(observed),
            proportion_observed = mean(observed)) %>%
  ungroup() -> baserates_by_gensup

write_stats_text("Base rate (overall proportion of possible SEs observed) is ",
                 baserates_by_gensup$n_observed[baserates_by_gensup$indic_supported],"/", baserates_by_gensup$n_rows[baserates_by_gensup$indic_supported]," (",
                 percent(baserates_by_gensup$proportion_observed[baserates_by_gensup$indic_supported],digits=1),") with genetic support or ",
                 baserates_by_gensup$n_observed[!baserates_by_gensup$indic_supported],"/", baserates_by_gensup$n_rows[!baserates_by_gensup$indic_supported]," (",
                 percent(baserates_by_gensup$proportion_observed[!baserates_by_gensup$indic_supported],digits=1),") without.")


dse_info_with_indic_assoc_sim %>%
  filter(genetic_insight != 'none'  & observed) -> temp
ctable = table(temp[,c('sim_assoc','indic_supported')])
fobj = fisher.test(ctable)

write_stats_text("Among observed drug-SE pairs, genetic evidence for the SE and genetic evidence for the indication are enriched for each other, OR = ",
                 formatC(fobj$estimate,format='f',digits=1),", P = ",
                 formatC(fobj$p.value,format='e',digits=1))


# this test asks whether genetic support for indications predicts
# that a side effect will be osberved. but it's not a good test
# because each row is aligning the SE to the drug's indication (out of 
# many possible indications) that is most similar, so it's 
# kind of just a roundabout way of asking about genetic evidence
# for the SE itself.
# dse_info_with_indic_assoc_sim %>%
#   filter(genetic_insight != 'none') -> temp
# ctable = table(temp[,c('observed','indic_supported')])
# fisher.test(ctable) 

dse_info_with_indic_assoc_sim %>%
  filter(genetic_insight != 'none' & !sim_indic) -> temp
ctable = table(temp[,c('observed','indic_supported')])
fobj = fisher.test(ctable)
write_stats_text("Among observed drug-SE pairs, genetic evidence for the SE and ",
                 "genetic evidence for the indication are not enriched after the ",
                 "similar indication filted is applied, OR = ",
                 formatC(fobj$estimate,format='f',digits=1),", P = ",
                 formatC(fobj$p.value,format='e',digits=1))

#####
# SUPPLEMENT
#####

tell_user('done.\nFinalizing supplementary tables...')

# write the supplement directory / table of contents
supplement_directory %>% rename(table_number = name, description=title) -> contents
addWorksheet(supplement,'contents')
bold_style = createStyle(textDecoration = "Bold")
writeData(supplement,'contents',contents,headerStyle=bold_style,withFilter=T)
freezePane(supplement,'contents',firstRow=T)
# move directory to the front
original_order = worksheetOrder(supplement)
n_sheets = length(original_order)
new_order = c(n_sheets, 1:(n_sheets-1))
worksheetOrder(supplement) = new_order
activeSheet(supplement) = 'contents'
# now save
saveWorkbook(supplement,supplement_path,overwrite = TRUE)

elapsed_time = Sys.time() - overall_start_time
cat(file=stderr(), paste0('done.\nAll tasks complete in ',round(as.numeric(elapsed_time),1),' ',units(elapsed_time),'.\n'))

