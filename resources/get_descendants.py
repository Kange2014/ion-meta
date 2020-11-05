#!/usr/bin/python3

import sys
from ete3 import NCBITaxa
import pandas as pd

input_file = sys.argv[1]
output_file = sys.argv[2]
taxdb = sys.argv[3] # "resources/taxonomy/taxdb.sqlite"

def get_descendant(taxid, ncbitaxa):
    descendants = ncbitaxa.get_descendant_taxa(taxid,intermediate_nodes=True)
    if taxid not in descendants: descendants = descendants + [taxid]
    return descendants

ncbi = NCBITaxa(taxdb)
fw = open(output_file,"w")

df = pd.read_csv(input_file, header=0, sep="\t")
names = df.columns.values.tolist()
fw.write("\t".join(names) + "\t" + "Descendants"+ "\n")

for index, row in df.iterrows():
    try:
        name2taxid = ncbi.get_name_translator([row['TaxonomyName']])
    except ValueError:
        sys.stderr.write("WARNING: Name {} missing from db\n".format(row['TaxonomyName']))
        continue
    for name, taxid in name2taxid.items():
        if(int(taxid[0]) != int(row['TaxonomyID'])):
            print("WARNING: Name {} has a different ID {} in db from {}.\n".format(name, str(taxid[0]), str(row['TaxonomyID'])))
        descendants = get_descendant(taxid[0],ncbi)
        fw.write("\t".join(map(str,row)) + "\t" + ",".join(map(str,descendants))+ "\n")

fw.close()
