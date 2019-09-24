# ion-meta
A workflow using snakemake to analyze ion torrent based metagenomics data

## Usage
#### Step 1: Install conda
First, you have to install the Miniconda Python3 distribution. See here () for installation instructions. 

#### Step 2: Clone the repository of ion-meta
git clone https://github.com/Kange2014/ion-meta
cd ion-meta

#### Step 3: Install the required software
mkdir envs/ion-meta
conda env create -f envs/environment.yaml -p envs/ion-meta
conda config --add envs_dirs $(pwd)/envs/
conda activate envs/ion-meta

Open R, and install R package pavian and Rsamtolls for reporting

> if (!require(remotes)) { install.packages("remotes") }
> remotes::install_github("fbreitwieser/pavian")
> source("https://bioconductor.org/biocLite.R")
> biocLite("Rsamtools")

#### Step 4: Initiate and Configure workflow
Initiate and Configure the workflow according to your needs via editing the file config.yaml.

python scripts/init.py --data_fp /path/to/bam/files /path/to/my_project --single_end --format {sample}.bam

In the generated directory, a new config file and a new sample list were created (by default named config.yaml and samplelist.csv, respectively). Edit the config file in your favorite text editor, in particular for below to ensure they match your case:
taxdb: /results/luze/ion-meta/resources
centrifuge_dir: /results/luze/ion-meta/resources/classifier_db
centrifuge_base: bacteria.archaea.viral.fungi.protozoa.human

#### Step 5: Execute workflow
Test your configuration by performing a run via
snakemake --configfile /path/to/my_project/config.yaml -j 10
