#!/bin/bash

wget -O - https://www.anaconda.com/distribution/ 2>/dev/null | sed -ne 's@.*\(https:\/\/repo\.anaconda\.com\/archive\/Anaconda3-.*-Linux-x86_64\.sh\)\">64-Bit (x86) Installer.*@\1@p' | xargs wget

find . -name "Anacond*" -exec bash {} \;
conda config --add channels conda-forge bioconda
conda config --set channel_priority strict
conda install -c bioconda sra-tools entrez-direct fastqc cutadapt samtools hisat2 stringtie salmon
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets
chmod +x datasets

conda install -c conda-forge r-base
conda install -c conda-forge r-fastqcr r-tidyr r-stringr r-dplyr
