library(tidyverse)
library(janitor)
library(openxlsx)
setwd('~/d/sci/src/genetic_safety')


clipcopy = function(tbl) {
  clip = pipe("pbcopy", "w")  
  write.table(tbl, file=clip, sep = '\t', quote=F, row.names = F, na='')
  close(clip)
}


map_pp_ids = function(drugtable, verbose=T) {
  drugtable$drug_name = tolower(drugtable$drug_name)
  drugtable$db_name = tolower(drugtable$db_name)
  drugtable$ppid = as.integer(NA)
  drugtable$pp_match_type = 'missing'
  pp$drugnamelower = tolower(iconv(pp$drug_primary_name, 'utf-8', 'utf-8', sub=' '))
  pp$synonymlower = tolower(iconv(pp$drug_name_synonyms, 'utf-8', 'utf-8', sub=' '))
  
  primary_ppid = pp$drug_id[match(drugtable$drug_name, pp$drugnamelower)]
  to_update = !is.na(primary_ppid)
  drugtable$pp_match_type[to_update] = 'exact primary'
  drugtable$ppid[to_update] = primary_ppid[to_update]
  
  primary_db_ppid = pp$drug_id[match(drugtable$db_name, pp$drugnamelower)]
  to_update = is.na(drugtable$ppid) & !is.na(primary_db_ppid)
  drugtable$pp_match_type[to_update] = 'exact db primary'
  drugtable$ppid[to_update] = primary_db_ppid[to_update]
  
  synonym_ppid = pp$drug_id[match(drugtable$drug_name, pp$synonymlower)]
  to_update = is.na(drugtable$ppid) & !is.na(synonym_ppid)
  drugtable$pp_match_type[to_update] = 'exact synonym'
  drugtable$ppid[to_update] = synonym_ppid[to_update]

  synonym_db_ppid = pp$drug_id[match(drugtable$db_name, pp$synonymlower)]
  to_update = is.na(drugtable$ppid) & !is.na(synonym_db_ppid)
  drugtable$pp_match_type[to_update] = 'exact db synonym'
  drugtable$ppid[to_update] = synonym_db_ppid[to_update]
  
  cas_ppid = pp$drug_id[match(drugtable$db_cas, pp$cas_number)]
  to_update = is.na(drugtable$ppid) & !is.na(cas_ppid)
  drugtable$pp_match_type[to_update] = 'db cas match'
  drugtable$ppid[to_update] = cas_ppid[to_update]
    
  for (i in 1:nrow(drugtable)) {
    if (!is.na(drugtable$ppid[i])) {
      next
    }
    if (verbose) {
      cat(file=stderr(),paste0('\r- Matching drug names, row ',i,'/',nrow(drugtable),'...'))
      flush.console()
    }
    primary_matches = grep(drugtable$drug_name[i], pp$drugnamelower, fixed=T) 
    if (length(primary_matches) > 0) {
      drugtable$ppid[i] = pp$drug_id[primary_matches[1]]
      drugtable$pp_match_type[i] = 'substring primary'
      next
    } 
    synonym_matches = grep(drugtable$drug_name[i], pp$synonymlower, fixed=T)
    if (length(synonym_matches) > 0) {
      drugtable$ppid[i] = pp$drug_id[synonym_matches[1]]
      drugtable$pp_match_type[i] = 'substring synonym'
      next
    }
    primary_matches = grep(drugtable$db_name[i], pp$drugnamelower, fixed=T) 
    if (length(primary_matches) > 0) {
      drugtable$ppid[i] = pp$drug_id[primary_matches[1]]
      drugtable$pp_match_type[i] = 'substring db primary'
      next
    } 
    synonym_matches = grep(drugtable$db_name[i], pp$synonymlower, fixed=T)
    if (length(synonym_matches) > 0) {
      drugtable$ppid[i] = pp$drug_id[synonym_matches[1]]
      drugtable$pp_match_type[i] = 'substring db synonym'
      next
    }
  }
  if (verbose) {
    cat(file=stderr(),paste0('done!\n'))
    flush.console()
  }
  
  return (drugtable)
}


assoc = read_tsv('data/assocs/assoc.tsv.gz', col_types=cols())
scr_best = read_tsv('data/mesh/mesh_scr_to_best.tsv.gz', col_types=cols())
mesh_all_vocab = read_tsv('data/mesh/mesh_all_vocab.tsv.gz', col_types=cols())
mesh_best_names = read_tsv('data/mesh/mesh_best_names.tsv.gz', col_types=cols())

sim = read_tsv('data/mesh/mesh_sim.tsv.gz', col_types=cols())

drug_names = read_tsv('data/sider/drug_names.tsv', 
                      col_types=cols(),
                      col_names=c('drug_id','drug_name'))

db = read_tsv('data/other/drugbank.tsv', col_types=cols())
drug_atc = read_tsv('data/sider/drug_atc.tsv', col_types=cols(), col_names=c('drug_id','atc'))
drug_atc$db_name = db$drug[match(drug_atc$atc, db$atc)]
drug_atc$db_cas = db$cas[match(drug_atc$atc, db$atc)]
drug_names$db_name = drug_atc$db_name[match(drug_names$drug_id, drug_atc$drug_id)]
drug_names$db_cas = drug_atc$db_cas[match(drug_names$drug_id, drug_atc$drug_id)]

pp_full = read.xlsx('ignore/launched_drug.xlsx') %>% 
  as_tibble() %>%
  mutate(drug_id = as.integer(drug_id))

pp_full %>%
  rename(orig_mesh = mesh_id) %>%
  left_join(scr_best, by=c('orig_mesh'='scr')) %>% # join to mapping of SCRs to main headings. note that sum(duplicated(scr_best$scr))==0 so there are no many-to-one joins here
  mutate(mesh_id = coalesce(main, orig_mesh)) %>% # choose the main that the SCR is mapped to, but if not an SCR, just use the original
  left_join(mesh_best_names, by=c('mesh_id'='id')) %>%
  rename(indication_mesh_id = mesh_id) %>%
  rename(indication_mesh_term=labeltext) %>%
  filter(!grepl('iagnos',indication_mesh_term)) %>%
  filter(!grepl('iagnos',disease_name)) %>%
  filter(!is.na(indication_mesh_id) & !(indication_mesh_id %in% c('Unknown Mesh ID','Not applicable'))) %>%
  distinct(drug_id, indication_mesh_id) -> launched_indic

launched_indic %>%
  group_by(drug_id) %>%
  summarize(.groups='keep', n_launched_indic = length(unique(indication_mesh_id))) %>%
  ungroup() %>%
  select(drug_id, n_launched_indic) -> launched_status

pp_full %>%
  left_join(mesh_best_names, by=c('mesh_id'='id')) %>%
  rename(indication_mesh_term=labeltext) %>%
  filter(!grepl('iagnos',indication_mesh_term)) %>%
  filter(!grepl('iagnos',disease_name)) %>%
  filter(!is.na(gene)) %>%
  distinct(drug_id, gene) -> all_targets

all_targets %>%
  group_by(drug_id) %>%
  summarize(.groups='keep', n_targets = length(unique(gene))) -> target_count

pp = read_tsv('ignore/pp_drug_synonyms.tsv', col_types=cols()) %>%
  mutate(drug_id = as.integer(drug_id)) %>%
  left_join(launched_status, by='drug_id') %>%
  mutate(n_launched_indic = replace_na(n_launched_indic, 0)) %>%
  left_join(target_count, by='drug_id') %>%
  mutate(n_targets = replace_na(n_targets, 0)) %>%
  rename(drug_name_synonyms = drug_name_synomyns) %>%
  select(drug_id, drug_primary_name, drug_name_synonyms, cas_number, n_launched_indic, n_targets) %>%
  mutate(priority = case_when(n_launched_indic > 0 & n_targets > 0 ~ 1,
                              n_launched_indic == 0 & n_targets > 0 ~ 2,
                              n_launched_indic > 0 & n_targets==0 ~ 3,
                              n_launched_indic == 0 & n_targets==0 ~ 4)) %>%
  arrange(priority, desc(n_launched_indic), desc(n_targets), drug_id)
  
drug_names = map_pp_ids(drug_names)

# check how many are launched
# sum(drug_names$ppid %in% launched_status$drug_id)

# set.seed(1)
# sample(drug_names$drug_name[drug_names$pp_match_type=='missing'], size=10)

# clipcopy(drug_names)
# clipcopy(pp)

# drug_names %>%
#   group_by(pp_match_type) %>%
#   summarize(.groups='keep', n=n()) %>%
#   ungroup() %>%
#   arrange(desc(n)) %>%
#   mutate(proportion = n/sum(n))

write_tsv(drug_names, 'data/synthesis/sider_drug_names.tsv')

drug_names = read_tsv('data/synthesis/sider_drug_names.tsv', col_types=cols())

all_targets %>%
  filter(drug_id %in% drug_names$ppid) %>%
  rename(ppid = drug_id) %>%
  distinct(ppid, gene) -> pp_targets

write_tsv(pp_targets, 'data/synthesis/pp_targets.tsv')


launched_indic %>%
  filter(drug_id %in% drug_names$ppid) %>%
  rename(ppid = drug_id) %>%
  distinct(ppid, indication_mesh_id) -> pp_launched_indic

write_tsv(pp_launched_indic, 'data/synthesis/pp_launched_indic.tsv')

#####
# SE - MeSH mapping
#####

sider_mesh_map = read_tsv('data/mesh/sider_mesh_map.tsv', col_types=cols())
# sider_mesh_map %>%
#   group_by(map_type) %>%
#   summarize(.groups='keep', n=n()) %>%
#   ungroup() %>%
#   arrange(desc(n)) %>%
#   mutate(proportion = n/sum(n))

all_se = read_tsv('data/sider/meddra_all_se.tsv.gz', col_types=cols(),
                  col_names=c('drug_id','id2','umls_label','meddra_type','umls_meddra','se_name'))

# set.seed(1)
# sample(sider_mesh_map$se_name[is.na(sider_mesh_map$mesh_id)], size=10)

mesh_tree = read_delim('data/mesh/mesh_tree.tsv.gz', col_names = c("id", "tree"), col_types=cols())
mesh_tree$topl = substr(mesh_tree$tree,1,3)
mesh_tree$letter = substr(mesh_tree$tree,1,1)

# first create the general match table without filters
sider_mesh_map %>%
  distinct(mesh_id) -> se_mesh_ids
assoc %>%
  distinct(mesh_id) -> assoc_mesh_ids
launched_indic %>%
  rename(mesh_id = indication_mesh_id) %>%
  distinct(mesh_id) -> indic_mesh_ids
rbind(se_mesh_ids, assoc_mesh_ids, indic_mesh_ids) %>%
  distinct(mesh_id) %>%
  left_join(mesh_tree, by=c('mesh_id' = 'id')) %>%
  group_by(mesh_id, topl) %>%
  summarize(.groups='keep') %>%
  arrange(mesh_id, topl) -> se_topl_match_raw
se_topl_match_raw$topl[is.na(se_topl_match_raw$topl)] = 'OTH'
topl_maps = read_tsv('data/mesh/mesh_2023_topl_maps.tsv',col_types=cols())
se_topl_match = se_topl_match_raw
keeps = se_topl_match$topl %in% topl_maps$topl[topl_maps$topl == topl_maps$map_to]
potential_remaps = se_topl_match$topl %in% topl_maps$topl[topl_maps$map_to != topl_maps$topl]
# before remapping, if the SE maps to a kept area, simply use that - delete other mappings
to_delete = potential_remaps & se_topl_match$mesh_id %in% se_topl_match$mesh_id[keeps]
if (sum(to_delete) > 0) {
  se_topl_match = se_topl_match[-which(to_delete),]
}
# now apply remapping
remaps = se_topl_match$topl %in% topl_maps$topl[topl_maps$map_to != topl_maps$topl]
se_topl_match$topl[remaps] = topl_maps$map_to[match(se_topl_match$topl[remaps], topl_maps$topl)]
# bake in the non-cancer filter: delete all rows that are non-C04 but DO have a match to a thing in C04
delete_cancer = se_topl_match$topl != 'C04' & se_topl_match$mesh_id %in% se_topl_match$mesh_id[se_topl_match$topl == 'C04']
if (sum(delete_cancer) > 0) {
  se_topl_match = se_topl_match[-which(delete_cancer),]
}
# bake in the other filter: delete all rows that are OTH where the SE does match a non-OTH area
delete_other = se_topl_match$topl == 'OTH' & se_topl_match$mesh_id %in% se_topl_match$mesh_id[se_topl_match$topl != 'OTH']
if (sum(delete_other) > 0) {
  se_topl_match = se_topl_match[-which(delete_other),]
}
# remove dups (some that mapped to >1 other area have multiple OTH rows)
se_topl_match %>%
  filter(!is.na(mesh_id)) %>%
  group_by(mesh_id, topl) %>%
  slice(1) %>%
  ungroup() -> se_topl_match
write_tsv(se_topl_match, 'data/synthesis/se_topl_match.tsv')


# sider_mesh_map %>%
#   inner_join(all_se, by='se_name') %>%
#   group_by(se_name) %>%
#   summarize(.groups='keep', drugs=length(unique(drug_id))) 




######
# GENETIC INSIGHT FOR SEs
######

# indications table
v2d_uniq_assoc = read_tsv('../digap/output/otg2209/v2d_uniq_assoc.tsv.gz',col_types=cols())
v2d_uniq_assoc$k1m = formatC(round(v2d_uniq_assoc$lead_pos/1e6)*1e6,format='d',width=9,flag='0')
v2d_uniq_assoc$locus = paste0(v2d_uniq_assoc$lead_chrom,':',v2d_uniq_assoc$k1m)
assoc %>%
  group_by(original_trait, mesh_id) %>%
  summarize(.groups='keep') %>%
  arrange(original_trait, mesh_id) -> assoc_map
v2d_uniq_assoc$mesh_id = assoc_map$mesh_id[match(v2d_uniq_assoc$trait_reported,assoc_map$original_trait)]

sider_mesh_map %>%
  filter(!is.na(mesh_id)) %>%
  distinct(mesh_id, mesh_term) -> se_mesh

sim %>%
  filter(comb_norm >= 0.8,
         meshcode_a %in% se_mesh$mesh_id) -> mesh_sim_to_join

assoc %>%
  filter(source %in% c('OMIM','intOGen')) %>% 
  inner_join(mesh_sim_to_join, by=c('mesh_id' = 'meshcode_b')) %>%
  inner_join(se_mesh, by=c('meshcode_a' = 'mesh_id'), keep=T, suffix=c('_a','_s')) %>%
  group_by(meshcode_a) %>%
  summarize(.groups='keep', n_omim_genes=n_distinct(gene)) %>%
  ungroup() %>%
  filter(n_omim_genes >= 1)  -> omimintogen_insight

v2d_uniq_assoc %>%
  filter(!is.na(locus)) %>%
  group_by(mesh_id, locus) %>%
  summarize(.groups='keep') -> gwas_loci

se_mesh %>%
  filter(mesh_id %in% mesh_sim_to_join$meshcode_a) -> se_to_join

gwas_loci %>%
  filter(mesh_id %in% mesh_sim_to_join$meshcode_b) %>%
  inner_join(mesh_sim_to_join, by=c('mesh_id' = 'meshcode_b')) %>%
  inner_join(se_to_join, by=c('meshcode_a' = 'mesh_id'), keep=T) %>%
  group_by(meshcode_a) %>%
  summarize(.groups='keep', n_gwas_loci=n_distinct(locus)) %>%
  ungroup() %>%
  filter(n_gwas_loci >= 1) -> gwas_insight

se_mesh$n_omimintogen_genes = omimintogen_insight$n_omim_genes[match(se_mesh$mesh_id, omimintogen_insight$meshcode_a)]
se_mesh$n_omimintogen_genes[is.na(se_mesh$n_omimintogen_genes)] = 0
se_mesh$n_gwas_loci = gwas_insight$n_gwas_loci[match(se_mesh$mesh_id, gwas_insight$meshcode_a)]
se_mesh$n_gwas_loci[is.na(se_mesh$n_gwas_loci)] = 0
se_mesh$genetic_insight = 'none'
se_mesh$genetic_insight[se_mesh$n_omimintogen_genes >= 1] = 'omim/intogen'
se_mesh$genetic_insight[se_mesh$n_gwas_loci >= 3] = 'gwas'
se_mesh$genetic_insight[se_mesh$n_omimintogen_genes >= 1 & se_mesh$n_gwas_loci >= 3] = 'both'
se_mesh$genetic_insight[is.na(se_mesh$genetic_insight)] = 'none'

# require presence in sim matrix
se_mesh$genetic_insight[!(se_mesh$mesh_id %in% sim$meshcode_a)] = 'none'

write_tsv(se_mesh, 'data/synthesis/se_insight.tsv')

