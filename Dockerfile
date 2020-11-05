FROM continuumio/miniconda3

LABEL maintainer="Lucius Zheng" email=lulu.zheng@thermofisher.com
LABEL description="Docker image for ion-meta workflow for metagenomics"

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set workdir
WORKDIR /ion-meta

# Add project files
COPY envs ./envs/
COPY config ./config/
COPY rules ./rules/
COPY scripts ./scripts/
COPY resources ./resources/
COPY Snakefile .

# Install environment into base
RUN apt-get install zip && apt-get clean
RUN conda env update -n base -f envs/environment.yaml \
&& conda clean -a

# install other tools not available on conda
RUN R -e "if (!require(remotes)) { install.packages('remotes', repos='https://cloud.r-project.org') } " \
 && R -e "remotes::install_github('fbreitwieser/pavian') "\
 && R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org') "\
 && R -e "BiocManager::install('Rsamtools') "

# Run workflow
#CMD python scripts/init.py --data_fp example example --single_end --format {sample}.bam
#CMD snakemake --configfile example/config.yaml --use-conda -j 10 -rp
CMD snakemake -h
