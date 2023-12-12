#!/usr/bin/python3

import sys
import gzip
from xml.etree import cElementTree

# adapted from https://github.com/ericminikel/nd_trials/blob/sio/src/drugbank.py

# this script takes several minutes, mostly due to a small handful (<10) drugs where it hangs
# I decided to just live with the slowness because the top solution online is pretty non-trivial: https://stackoverflow.com/a/56090601/3806692

hgnc_path = '/Users/eminikel/d/sci/src/nd_trials/other_data/hgnc.txt'
with open(hgnc_path) as h:
	lines = h.readlines()
	hgnc_dict = {}
	for line in lines:
		cells = line.split('\t')
		hgnc_dict[cells[0]] = cells[1]

in_path = '/Users/eminikel/d/sci/src/nd_trials/bigdata/drugbank.xml.gz'
#in_path = 'bigdata/drugbank_sample.xml'

if in_path[-3:] == '.gz':
    f = gzip.open(in_path)
else:
    f = open(in_path)

sys.stdout.write('drug'+'\t'+'status'+'\t'+'target'+'\t'+'atc'+'\t'+'cas'+'\n')
counter = 0
iterparser = iter(cElementTree.iterparse(f, events=("start", "end")))
for event, elem in iterparser:
    if elem.tag == 'drug' and event == 'end':
        drug = elem
        counter += 1
        sys.stderr.write('\rProcessing drug '+str(counter))
        name = ''
        status = 'non'
        target = ''
        atc = ''
        #### NAME
        name_node = drug.find('name')
        if name_node is not None:
        	name = name_node.text.lower()
        #### APPROVAL STATUS
        groups_node = drug.find("groups")
        if groups_node is not None:
            for group_node in groups_node.findall("group"):
                if group_node.text == "approved":
                    status = 'approved'
        #### CAS
        cas_node = drug.find("cas-number")
        if cas_node is not None:
            if cas_node.text is not None:
                cas = cas_node.text
            else:
                cas = ''
        else:
            cas = ''
        #### TARGET
        targets_node = drug.find("targets")
        if targets_node is not None:
            target_nodes = targets_node.findall("target")
            for target_node in target_nodes: # loop through targets listed
                # remove if not ranked #1
                #if 'position' in target_node.attrib:
                #    if target_node.attrib["position"] != '1':
                #        continue # on to next target in for loop
                # remove if not "known action"
                known_action = target_node.find("known-action")
                if known_action.text != "yes":
                    continue
                # remove if not human polypeptide target
                polypeptide_node = target_node.find("polypeptide")
                if polypeptide_node is None:
                    continue
                organism_node = polypeptide_node.find("organism")
                if organism_node.text != "Humans":
                    continue
                extid_nodes = polypeptide_node.find("external-identifiers")
                for extid_node in extid_nodes.findall("external-identifier"):
                    if extid_node.find("resource").text == "HUGO Gene Nomenclature Committee (HGNC)":
                        hgnc_id = extid_node.find("identifier").text
                        gene_symbol = hgnc_dict[hgnc_id]
                        if target == '':
                            target = gene_symbol
                        else:
                        	target = target + ',' + gene_symbol
        #### ATC
        atc_codes_parent = drug.findall("atc-codes")
        if atc_codes_parent is not None and len(atc_codes_parent) > 0:
            atc_codes = atc_codes_parent[0].findall("atc-code")
            if atc_codes is not None and len(atc_codes) > 0:
                atc_code = atc_codes[0]
                atc = atc_code.attrib['code']
        sys.stdout.write(name+'\t'+status+'\t'+target+'\t'+atc+'\t'+cas+'\n')
        elem.clear()

