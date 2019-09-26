#!/usr/bin/env python

from ete3 import NCBITaxa
import pandas as pd
from argparse import ArgumentParser


def get_seqids(df, seqid2taxid, ncbitaxa, host_taxid):
    filtered = pd.DataFrame()
    for i in df.index:
        taxid = df.loc[i,"taxID"]
        # discard host: human genome
        if str(taxid) == host_taxid:
            continue
        # Intersect descendants with those we have genomes for
        descendants = seqid2taxid.loc[seqid2taxid["taxID"].isin(ncbitaxa.get_descendant_taxa(taxid)+[taxid])]
        if len(descendants) == 0:
            continue
        filtered = pd.concat([filtered,descendants])
    return filtered


def filter_genomes(input,db,mapfile,min_read_count,min_abundance,host_taxid,output):
    df = pd.read_csv(input, sep="\t")
    # Filter to number of reads
    df = df.loc[df.numReads >= 3] # for each classification, at least 3 different reads
    df = df.loc[(df.numReads >= min_read_count) | (df.abundance >= min_abundance)]
    # Filter to ranks genus, species and leaf
    df = df.loc[df.taxRank.isin(["leaf", "species", "genus"])]
    # Read the centrifuge sqlite database
    ncbitaxa = NCBITaxa(db)
    # Read the seqid -> taxid map file
    seqid2taxid = pd.read_csv(mapfile, header=None, names=["seq", "taxID"], index_col=0, sep="\s+")
    filtered = get_seqids(df, seqid2taxid, ncbitaxa, host_taxid)
    if len(filtered) == 0:
        filtered = pd.DataFrame(columns=["taxID"], index=["seq"])
    filtered.to_csv(output, sep="\t", index=True)


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
