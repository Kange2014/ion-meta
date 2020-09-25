# ion-meta
A workflow using snakemake to analyze ion torrent based metagenomics data

## Usage Simple
#### Step 1: Install conda
First, you have to install the Miniconda Python3 distribution. See [here](https://docs.conda.io/en/latest/miniconda.html#installing) for installation instructions. 

#### Step 2: Clone the repository of ion-meta
    git clone https://github.com/Kange2014/ion-meta
    cd ion-meta

#### Step 3: Install the required software
    mkdir envs/ion-meta
    conda env create -f envs/environment.yaml -p envs/ion-meta
    conda config --add envs_dirs $(pwd)/envs/
    conda activate envs/ion-meta

Open R, and install R package pavian and Rsamtools for reporting

    > if (!require(remotes)) { install.packages("remotes") }
    > remotes::install_github("fbreitwieser/pavian")
    > source("https://bioconductor.org/biocLite.R")
    > biocLite("Rsamtools")

If finding errors in remotes::install_github("fbreitwieser/pavian"), e.g., about stringi, pls. re-install stringi package in R cmd window and then install pavian:

    > install.packages ("stringi")
    > remotes::install_github("fbreitwieser/pavian")

If install Rsamtools in R 3.5 or greater, install Bioconductor and packages first by typing the following in an R command window:

    > if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    > BiocManager::install("Rsamtools")

If finding errors in the above process about Rhtslib, pls. install Rhtslib package manually and set CPPFLAGS and LDFLAGS in the makefiles:

    wget https://bioconductor.org/packages/release/bioc/src/contrib/Rhtslib_1.16.1.tar.gz
    tar xvzf Rhtslib_1.16.1.tar.gz
    cd Rhtslib/src/htslib-1.7/

In this directory, change the flags in the Makefile and Makefile.Rhtslib:
    *Comment the lines CPPFLAGS = and LDFLAGS =*
    *Change the CFLAGS = to CFLAGS +=*

Then, re-tar it:

    tar -czvf Rhtslib_1.16.1.tar.gz ./Rhtslib/

and R CMD INSTALL as described: 
    
    R CMD INSTALL Rhtslib_1.16.1.tar.gz

Then, in an R command window: 

    > BiocManager::install("Rsamtools")

#### Step 4: Initiate and Configure workflow
Initiate and Configure the workflow according to your needs via editing the file config.yaml.

    python scripts/init.py --data_fp /path/to/bam/files /path/to/my_project --single_end --format {sample}.bam

In the generated directory, a new config file and a new sample list were created (by default named config.yaml and samplelist.csv, respectively). Edit the config file in your favorite text editor, in particular for below to ensure they match your case:

    taxdb: /results/luze/ion-meta/resources
    centrifuge_dir: /results/luze/ion-meta/resources/classifier_db
    centrifuge_base: bacteria.archaea.viral.fungi.protozoa.human

#### Step 5: Execute workflow
Test your configuration by performing a run via
    
    snakemake --configfile /path/to/my_project/config.yaml -j 10 --use-conda

## Appendix:
#### 1. update referenge genome database, taxonomy database and build index

ion-meta uses [centrifuge](https://ccb.jhu.edu/software/centrifuge/) as its primary reads classification. So, the database indexes can be built with arbritary sequences like centrifuge. Standard choices are all of the complete bacteria, archaea, viral, fungi, and human genomes, or using the sequences that are part of the BLAST nt database.

Building index on all complete bacterial, archaea, viral, fungi, and human genomes
    
    # first use centrifuge-download to download genomes and taxonomy information from NCBI. 
    $ centrifuge-download -o taxonomy taxonomy 
    $ centrifuge-download -o library -P 8 -m -d "archaea,bacteria,viral,fungi,protozoa" refseq > seqid2taxid.map
    $ centrifuge-download -o library -d "vertebrate_mammalian" -a "Chromosome" -t 9606 -c 'reference genome'  refseq >> seqid2taxid.map

    # to build the index, first concatenate all downloaded sequences into a single file, and then run centrifuge-build:
    $ cat library/*/*.fna > input-sequences.fna

    # build centrifuge index with 8 threads, which results in three index files named bacteria.archaea.viral.fungi.protozoa.human.*.[1234].cf index files
    $ centrifuge-build -p 8 --bmax 1342177280 --conversion-table seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp input-sequences.fna bacteria.archaea.viral.fungi.protozoa.human

If you want to get summary statistics info about the downloaded database, you could run:
    
    $ python resources/summarize_centrifuge_db.py -h

#### 2. update virus host annotation and KEGG pathogen database information:

To update virus host annotation:
    
    $ wget ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv
    $ mv virushostdb.tsv resources/

To update KEGG pathogen database:

    $ python resources/download_kegg_genomeinfo.py
    $ mv kegg_genomeinfo.tsv resources/

#### 3. install deeparg and database

For the first time to use ion-meta, pls. run below command, which will only install the required deeparg and the database without running the full workflow. Subsequent runs with --use-conda will make use of the local environments without requiring internet access. This is because deeparg has a different running environment from ion-meta.

    $ snakemake --use-conda --conda-create-envs-only
    
    
    
