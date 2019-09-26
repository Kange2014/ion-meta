from os.path import join as opj
import pandas as pd
import csv
from pathlib import Path
from snakemake.utils import min_version, validate

##### set minimum snakemake version #####
min_version("4.4.0")

# Disallow slashes in our sample names during Snakemake's wildcard evaluation.
# Slashes should always be interpreted as directory separators.

wildcard_constraints:
    sample="[^/]+"
  

shell.prefix("")

# Load config file
if not config:
    raise SystemExit(
        "No config file specified. Run `python scripts/init.py` to generate a "
        "config file, and specify with --configfile")

## Change your workdir
workdir: str(config['workdir'])

# Setting up config files and samples

def load_sample_list(samplelist_fp, paired_end):
    """
    Build a list of samples from a sample list file.
    :param samplelist_fp: a Path to a comma-delimited samplelist file,
       where the first entry is the sample name.
    :returns: A dictionary of samples with sample name and associated file(s)
    """
    Samples = {}
    with open(str(samplelist_fp)) as f:
        reader = csv.DictReader(f, fieldnames=['sample','1','2'])
        for row in reader:
            sample = row['sample']
            try:
                r1 = _verify_path(row['1'])
            except ValueError:
                raise ValueError("Associated file for {} not found.".format(
                    sample))
            r2 = ''
            if paired_end:
                try:
                    r2 = _verify_path(row['2'])
                except ValueError:
                    raise ValueError(
                        "Paired-end files specified, but mate pair for '{}' "
                        "is missing or does not exist.".format(sample))
            Samples[sample] = {'1': r1, '2': r2}
    return Samples

def _verify_path(fp):
    if not fp:
        raise ValueError("Missing filename")
    path = Path(fp)
    if not path.is_file():
        raise ValueError("File not found")
    return str(path.resolve())

Samples = load_sample_list(config['samplelist_fp'], config['paired_end'])

#Pairs = ['1', '2'] if config['paired_end'] else ['1']

# ---- Targets rules

rule all:
    input:
        expand(opj(config["report_path"],"bowtie2","{sample}.report.html"),sample=Samples.keys() )

# ---- Preprocess rules
include: "rules/Preprocess/preprocess.rules"

# ---- Classify rules
include: "rules/Classify/centrifuge.rules"

# ---- Filter rules
include: "rules/Filter/sourmash.rules"

# ---- Map rules
include: "rules/Map/bowtie2.rules"

# ---- Report rules
include: "rules/Report/report.rules"
