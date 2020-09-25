#!/usr/bin/env python

from argparse import ArgumentParser
import subprocess
from subprocess import check_output
import os
import sys

def get_hit_fq(fastqInput,readIDFile,fastqOutput):
    cmd = "seqtk subseq " + fastqInput + " " + readIDFile + ">" + fastqOutput
    try:
        seqtk_out = check_output(cmd,shell=True)
    except subprocess.CalledProcessError as e:
        print(e.output)
    
    return fastqOutput;

def do_assembly(args):
    fastq = get_hit_fq(args.fastqInput, args.readIDFile, args.fastqOutput)
    if args.spades:
        # In metagenemomics dataset, generally assigned reads are with high GC, low or uneven coverage
        # spades assembly contigs are in "scaffolds.fasta" 
        cmd = "spades.py -t " + str(args.threads) + " --iontorrent --sc -s " + fastq + " --careful -k 21,33,55 -o " + args.outdir
        try:
            spades_out = check_output(cmd,shell=True)
        except subprocess.CalledProcessError as e:
            print(e.returncode)
            open(args.outdir + "/scaffolds.fasta", 'a').close()
        os.rename(args.outdir + "/scaffolds.fasta", args.outdir + "/contigs.fa")
    
    else:
        # use megahit, and assembly contigs are in "final.contigs.fa"
        cmd = "megahit -t " + str(args.threads) + " -r " + fastq + " -o " + args.outdir
        try:
            megahit_out = check_output(cmd,shell=True);
        except subprocess.CalledProcessError as e:
            print(e.returncode)
            open(args.outdir + "/final.contigs.fa", 'a').close()
        os.rename(args.outdir + "/final.contigs.fa", args.outdir + "/contigs.fa")

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--fastqInput", type=str, required=True,
                        help="FastqInput for reads extracting.")
    parser.add_argument("-o", "--fastqOutput", type=str,
                        help="Extracted hit reads.")
    parser.add_argument("--readIDFile", type=str,
                        help="Read ids hitting to one specific taxa.")
    parser.add_argument("--spades", action="store_true",
                        help="Use SPAdes to assembly hit reads.")
    parser.add_argument("--outdir", type=str,
                        help="Assembly output directory.")
    parser.add_argument("--threads", type=int,
                        help="Number of CPU threads.")
    #parser.add_argument("--taxidHitFolder", type=str,
    #                    help="A filefoler to put readIDs for each taxid hit.")

    args = parser.parse_args()

    sys.stderr.write("Doing assembly for each target taxa.\n")
    
    do_assembly(args)


if __name__ == '__main__':
    main()