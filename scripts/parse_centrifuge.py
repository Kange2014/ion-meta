#!/usr/bin/env python

from argparse import ArgumentParser
from ete3 import NCBITaxa
import pandas as pd
import sys
import re

def read_centrifuge_report(reportfile):
    df = pd.read_csv(reportfile, header=0, sep="\t")
    tmp = df.set_index(['taxID'])
    genomesizes = tmp.to_dict()["genomeSize"]
    return df, genomesizes


def read_centrifuge_result(args):
    taxa_counts = {}
    df = pd.read_csv(args.infile, header=0, sep="\t")
    if args.host_taxid: df = df.loc[df.taxID != int(args.host_taxid)]
    if args.unique:     df = df.loc[df.numMatches == 1]
    if args.min_length: df = df.loc[df.hitLength >= int(args.min_length)]
    if args.min_score:  df = df.loc[df.score >= int(args.min_score)]
    
    df['count'] = 1/df.numMatches
    df['aligned'] = df.hitLength/df.numMatches
    
    taxids = list(set([ taxid for taxid in df.taxID]))
    
    for taxID in taxids:
        if not taxID in taxa_counts.keys():
                taxa_counts[taxID] = {"count": 0, "aligned": 0, "score":0, "readID": "" }
        
        taxa_counts[taxID]["count"] = int(sum(df.loc[df.taxID == taxID,'count']))
        taxa_counts[taxID]["aligned"] = sum(df.loc[df.taxID == taxID,'aligned'])
        taxa_counts[taxID]["score"] = int(sum(df.loc[df.taxID == taxID,'score']))
        taxa_counts[taxID]["readID"] = ",".join(df.loc[df.taxID == taxID,'readID'])
        
    return taxa_counts

def generate_ete3db(taxdbfile):
    taxdb = NCBITaxa(taxdbfile)
    return taxdb

def get_lineage_names(lineage, taxdb):
    ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    names = {}
    for rank in ranks:
        names[rank] = "unclassified"
    for id, rank in lineage.items():
        if rank not in ranks:
            continue
        name = taxdb.translate_to_names([id])[0]
        if rank == "superkingdom": names[rank] = "k__" + name
        if rank == "phylum": names[rank] = "p__" + name
        if rank == "class": names[rank] = "c__" + name
        if rank == "order": names[rank] = "o__" + name
        if rank == "family": names[rank] = "f__" + name
        if rank == "genus": names[rank] = "g__" + name
        if rank == "species": names[rank] = "s__" + name
    name_string = ""
    for rank in ranks:
        if rank == "species": name_string+="{}".format(names[rank])
        else: name_string+="{}|".format(names[rank])
    return name_string

def generate_taxlineage_table(taxa_counts, reportdf, args, genomesizes):
    taxdb = generate_ete3db(args.taxdb)
    lineage_dict = {}
    with open(args.reportTaxalineage, 'w') as fhout:
        fhout.write("name\ttaxID\ttaxRank\tgenomeSize\tnumReads\tnumUniqueReads\tabundance\ttaxaLineage\ttotalScore\talignedNorm\n");
        for taxid in taxa_counts.keys():
            try:
                name_string = lineage_dict[taxid]
            except KeyError:
                try:
                    lineage = taxdb.get_rank(taxdb.get_lineage(taxid))
                    name_string = get_lineage_names(lineage, taxdb)
                except ValueError:
                    name_string = "|".join(["Unknown"]*7)
                lineage_dict[taxid] = name_string
            count = taxa_counts[taxid]["count"]
            aligned = taxa_counts[taxid]["aligned"]
            score = taxa_counts[taxid]['score']
            report = reportdf.loc[reportdf['taxID'] == int(taxid) ]
            report = report.values.tolist()
            report = report[0]
            if args.normalize and genomesizes:
                try:
                    genomesize = genomesizes[int(taxid)]
                except KeyError:
                    continue
                # Normalized is aligned bases per kb of genome
                try:
                    norm = aligned / float(genomesize) * 1000
                except ZeroDivisionError:
                    norm = 0
                fhout.write("{}\t{}\t{}\t{}\n".format("\t".join(map(str,report)), name_string,score,norm))
            else:
                fhout.write("{}\t{}\t{}\t{}\n".format("\t".join(map(str,report)), name_string,score,count))
            
            # for virus hits, extract mapped read ids
            if re.search(r"k__Viruses",name_string):
                #if taxa_counts[taxid]['count'] >= 3:
                    with open(args.taxidHitFolder + "/" + str(taxid) + ".txt","w") as readIDFile:
                        readIDFile.write(taxa_counts[taxid]['readID'].replace(",","\n") )
                        readIDFile.write("\n")
            
def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str, required=True,
                        help="Centrifuge output.")
    parser.add_argument("-r", "--reportfile", type=str,
                        help="Centrifuge report file.")
    parser.add_argument("--reportTaxalineage", type=str,
                        help="Produce centrifuge report file with taxa lineage")
    parser.add_argument("--normalize", action="store_true",
                        help="Normalize read counts by genome size")
    parser.add_argument("--unique", action="store_true",
                        help="Only count reads mapping uniquely")
    parser.add_argument("--min_score", type=int,
                        help="Require a minimum score for reads to be counted")
    parser.add_argument("--min_length", type=int,
                        help="Require a minimum alignment length to the read")
    parser.add_argument("--taxdb", type=str,
                        help="Ete3 sqlite database file or path to one that will be created.")
    parser.add_argument("--host_taxid", type=str,
                        help="Host taxonomy id that will be excluded.")
    parser.add_argument("--taxidHitFolder", type=str,
                        help="A filefoler to put readIDs for each taxid hit.")
    #parser.add_argument("--ranksplit", action="store_true",
    #                    help="")

    args = parser.parse_args()

    sys.stderr.write("Reading centrifuge results\n")
    taxa_counts = read_centrifuge_result(args)
    sys.stderr.write("{} taxa parsed\n".format(len(taxa_counts)))

    if args.reportfile:
        reportdf, genomesizes = read_centrifuge_report(args.reportfile)
    else:
        genomesizes = False

    if args.reportTaxalineage:
        generate_taxlineage_table(taxa_counts, reportdf, args, genomesizes)


if __name__ == '__main__':
    main()