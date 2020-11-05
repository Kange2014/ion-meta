#!/usr/bin/env python

from ete3 import NCBITaxa
import pandas as pd
from argparse import ArgumentParser

def get_seqids(df, seqid2taxid, ncbitaxa):
    filtered = pd.DataFrame()
    for i in df.index:
        taxid = df.loc[i,"taxID"]
        # Intersect descendants with those we have genomes for
        descendants = seqid2taxid.loc[seqid2taxid["taxID"].isin(ncbitaxa.get_descendant_taxa(taxid)+[taxid])].copy()
        if len(descendants) == 0:
            continue
        # we still use the ancestry taxid for those descendant genomes
        descendants["taxID"] = taxid
        filtered = pd.concat([filtered,descendants])
    return filtered

def filter_genomes(input,db,mapfile,min_read_count,min_abundance,host_taxid,output):
    df = pd.read_csv(input, sep="\t")
    # Filter to number of reads >= 3
    df = df.loc[df.numReads >= 3] # for each classification, at least 3 different reads
    # Filter to ranks genus, species, subspecies, and leaf
    df = df.loc[df.taxRank.isin(["leaf", "subspecies", "species", "genus"])]
    # discard host: human genome
    df = df.loc[df.taxID != int(host_taxid)]
    
    # Read the centrifuge sqlite database
    ncbitaxa = NCBITaxa(db)
    
    df["genus"] = ""
    df["family"] = ""
    df["kingdom"] = ""
    genus = []
    family = []
    kingdom = []
    HERV_taxids = []
    for taxid in df.taxID:
        lineage = ncbitaxa.get_lineage(taxid)
        ranks = ncbitaxa.get_rank(lineage)
        genus_name = ""
        family_name = ""
        kingdom_name = ""
        for key, value in ranks.items():
            name = ncbitaxa.get_taxid_translator([key]).get(key)
            if(name == "Human endogenous retroviruses"): HERV_taxids.append(taxid)
            if(value == "genus"): genus_name = name
            if(value == "family"): family_name = name
            if(value == "superkingdom"): kingdom_name = name
        genus.append(genus_name)
        family.append(family_name)
        kingdom.append(kingdom_name)
    
    df["genus"] = genus
    df["family"] = family
    df["kingdom"] = kingdom
    # for virus, keep all records, and 
    # for non-virus, filter based on defined cutoffs
    df_virus = df.loc[df.kingdom == "Viruses"]
    # exclude human endogenous retroviruses
    df_virus = df_virus.loc[ -df_virus['taxID'].isin(HERV_taxids) ]
    df_nonvirus = df.loc[df.kingdom != "Viruses"]
    df_filtered = df_nonvirus.loc[(df_nonvirus.numReads >= min_read_count) | (df_nonvirus.abundance >= min_abundance)]
    df_filtered = df_filtered.append(df_virus)
    
    # if empty, then keep both the most abundant one and one with most unique mapped reads
    if df_filtered.empty:
        loc1 = df.numReads.idxmax()
        loc2 = df.abundance.idxmax()
        df_filtered = df.loc[list(set([loc1,loc2]))]
    
    # Read the seqid -> taxid map file
    seqid2taxid = pd.read_csv(mapfile, header=None, names=["seq", "taxID"], sep="\s+")
    filtered = get_seqids(df_filtered, seqid2taxid, ncbitaxa)
    if len(filtered) == 0:
        filtered = pd.DataFrame(columns=["seq","taxID"])
    # add taxonomy information based on taxID by merging df_filtered
    ids2tax = pd.merge(filtered,df_filtered[["taxID","name","genus","family","kingdom"]],on="taxID")
    ids2tax.to_csv(output, sep="\t",index=False)

#def main(args):
    #filter_genomes(args.input, args.db, args.mapfile, args.min_read_count, args.min_abundance, args.host_taxid, args.output)

## -i snakemake.input[0] -d snakemake.input[1] -M snakemake.input[2] -r snakemake.params[0] 
## -a snakemake.params[1] -t snakemake.params[2] -o snakemake.output[0]

def main():
    filter_genomes(snakemake.input[0],snakemake.input[1],snakemake.input[2],snakemake.params[0],\
                   snakemake.params[1],snakemake.params[2],snakemake.output[0])

if __name__ == '__main__':
    #parser = ArgumentParser()
    #parser.add_argument("-i", "--input", type=str, required=True)
    #parser.add_argument("-d", "--db", type=str, required=True)
    #parser.add_argument("-r", "--min_read_count", type=int, required=True)
    #parser.add_argument("-a", "--min_abundance", type=float, required=True)
    #parser.add_argument("-t", "--host_taxid", type=str, required=True)
    #parser.add_argument("-M", "--mapfile", type=str, required=True)
    #parser.add_argument("-o", "--output", type=str, required=True)
    #args = parser.parse_args()
    #main(args)
    main()
