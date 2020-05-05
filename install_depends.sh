#!/bin/bash

wget -O - https://www.anaconda.com/distribution/ 2>/dev/null | sed -ne 's@.*\(https:\/\/repo\.anaconda\.com\/archive\/Anaconda3-.*-Linux-x86_64\.sh\)\">64-Bit (x86) Installer.*@\1@p' | xargs wget

find . -name "Anacond*" -exec bash {} \;
find . -name "Anacond*" | xargs rm
conda update -y -n root conda
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible
conda install -y -c bioconda sra-tools
conda install -y -c bioconda entrez-direct
conda install -y -c bioconda fastqc
conda install -y -c bioconda cutadapt
conda install -y -c bioconda samtools
conda install -y -c bioconda hisat2
conda install -y -c bioconda stringtie
conda install -y -c bioconda salmon
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'
chmod +x datasets

# conda config --set channel_priority strict
conda install -y -c conda-forge r-base
conda install -y -c conda-forge r-tidyr
conda install -y -c conda-forge r-stringr
conda install -y -c conda-forge r-dplyr
conda install -y -c conda-forge r-fastqcr
conda install -y -c conda-forge r-devtools

conda update -y --all
# conda config --show-sources
