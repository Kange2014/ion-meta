#!/usr/bin/env python

from bioservices import KEGG
import json
import time
import re

k = KEGG()

organisms = k.list("genome")
organisms = organisms.rstrip().split("\n")

# get all organisms codes
org_codes = []
for org in organisms:
    org_code = org.rstrip().split("\t")
    org_codes.append(org_code[0])

records = {}
fout = open("kegg_genomeinfo.tsv","w")
fout.write("Taxonomy ID\tDefinition\tTaxonomy LINEAGE\tKeywords\tDisease\tComment\tReference\n")

length = len(org_codes)
for i in range(length):
    record = k.get(org_codes[i])
    if record == 403: 
        time.sleep(1)
        record = k.get(org_codes[i])
     
    # exclude viruses since parse function cannot deal with them
    m = re.search(r"Viruses;", str(record))
    if m: continue
    
    dict = k.parse(record)
    if not 'KEYWORDS' in dict: dict['KEYWORDS'] = ""
    if not 'DISEASE' in dict: dict['DISEASE'] = {}
    if not 'COMMENT' in dict: dict['COMMENT'] = ""
    ref = ""
    if 'REFERENCE' in dict:
        for key, value in dict["REFERENCE"][0].items():
            if ref: ref = ref + ", " + key + ":" + value
            else: ref = key + ":" + value
            
    fout.write(dict['TAXONOMY'][0]['TAXONOMY'].split(":")[1]+"\t"+dict['DEFINITION']+"\t"+dict['TAXONOMY'][0]['LINEAGE']+"\t"+dict['KEYWORDS']+"\t"+", ".join(list(dict['DISEASE'].values()))+"\t"+"; ".join(dict['COMMENT'])+"\t"+ref+"\n")
    records[dict['TAXONOMY'][0]['TAXONOMY'].split(":")[1]] = dict
    time.sleep(0.1)

fout.close()

json_str = json.dumps(records, indent=4)
with open('kegg_genomeinfo.json', 'w+') as json_file:
    json_file.write(json_str)
