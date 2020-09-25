#!/usr/bin/env python

from argparse import ArgumentParser
from ete3 import NCBITaxa
import sys

def summarize_to_ranks(mapfile, summary_rank, mapTaxafile, taxdbfile):
    ncbitaxa = NCBITaxa(taxdbfile)
    summary1 = {}
    summary2 = {}
    lineage_dict = {}
    sys.stderr.write("Reading seq2taxid mapfile and summarizing at rank {}\n".format(summary_rank))
    fhout = open(mapTaxafile, 'w')
    with open(mapfile, 'r') as fh:
        for i, line in enumerate(fh, start=1): # start=1: index starts from 1 
            if i%1000 == 0:
                sys.stderr.write("Read {} lines...\n".format(i))
            line = line.rstrip()
            try:
                seqid, taxid = line.split("\t")
            except ValueError:
                seqid, taxid = line.split()
            taxid = int(taxid)
            seq_name = ""
            try:
                seq_name = ncbitaxa.get_taxid_translator([taxid])[taxid]
            except ValueError:
                sys.stderr.write("WARNING: Taxid {} missing from db\n".format(taxid))
                continue
            
            ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
            names = {}
            for rank in ranks:
                names[rank] = "unclassified"
            
            lineage =  ncbitaxa.get_rank(ncbitaxa.get_lineage(taxid))
            for id, rank in lineage.items():
                if rank not in ranks:
                    continue
                name = ncbitaxa.translate_to_names([id])[0]
                if rank == "superkingdom": names[rank] = "k__" + name
                if rank == "phylum": names[rank] = "p__" + name
                if rank == "class": names[rank] = "c__" + name
                if rank == "order": names[rank] = "o__" + name
                if rank == "family": names[rank] = "f__" + name
                if rank == "genus": names[rank] = "g__" + name
                if rank == "species": names[rank] = "s__" + name
            
            if summary_rank not in names.keys():
                rank_name = "Unknown"
            else:
                rank_name = names[summary_rank]
            
            kingdom = names["superkingdom"]
            
            if not kingdom in summary1:
                summary1[kingdom] = {taxid: 0}
            else:
                if not taxid in summary1[kingdom]:
                    summary1[kingdom][taxid] = 0
            
            if not kingdom in summary2:
                summary2[kingdom] = {rank_name: 0}
            else:
                if not rank_name in summary2[kingdom]:
                    summary2[kingdom][rank_name] = 0
            
            summary1[kingdom][taxid] += 1
            summary2[kingdom][rank_name] += 1
            
            name_string = ""
            for rank in ranks:
                if rank == "species": name_string+="{}".format(names[rank])
                else: name_string+="{}\t".format(names[rank])
            
            fhout.write("{}\t{}\t{}\t{}\n".format(seqid, taxid, seq_name, name_string))
            
    fhout.close()
    
    return summary1,summary2

def main():
    parser = ArgumentParser()
    parser.add_argument("--seqid2taxidmap", required=True,
                        help="Sequence id to taxid mapping file.")
    parser.add_argument("--rank", default="species",
                        help="Summarize sequences and sizes at this rank (default = species)")
    parser.add_argument("--mapTaxafile", default="seqid2taxa.map",
                        help="Add taxa info into seqid2taxidmap file (default = seqid2taxa.map)")
    parser.add_argument("--taxdbfile", required=True,
                        help="Ete3 sqlite database file.")
                        
    args = parser.parse_args()
    summary_tax, summary_rank = summarize_to_ranks(args.seqid2taxidmap, args.rank, args.mapTaxafile, args.taxdbfile)

    print("{}\t{}\t{}".format("Superkingdom","Taxa number", "Sequence numer"))
    for kingdom in sorted(summary_tax.keys()):
        print("{}\t{}\t{}".format(kingdom, len(summary_tax[kingdom]), sum(summary_tax[kingdom].values()) ))
    
    print("{}\t{}\t{}".format("Superkingdom",args.rank + " number", "Sequence numer"))
    for kingdom in sorted(summary_rank.keys()):
        print("{}\t{}\t{}".format(kingdom, len(summary_rank[kingdom]), sum(summary_rank[kingdom].values()) ))
    

if __name__ == '__main__':
    main()