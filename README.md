
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BP4RNAseq

The assessment of gene expression is central to uncovering the functions
of the genome, understanding the regulation of development, and
investigating the molecular mechanisms that underlie cancer and other
diseases. RNA-sequencing (RNA-seq) now is the routine to assess the
genome wide gene expression due to its high speed, accuracy and
reproducibility, and low cost. An enormous volume of RNA-seq data that
have been accumulating and deposited in public data repositories, such
as the Gene Expression Omnibus (GEO) and the Sequence Read Archive
(SRA). Retrospectively analyzing these data or conducting a brand new
RNA-seq study is fundamentally important for researchers. However,
processing raw reads of RNA-seq data, no matter public or newly
sequenced data, involves a lot of specialized tools and technical
configurations that are often unfamiliar and time-consuming to learn for
non-bioinformatics researchers. The goal of BP4RNAseq is to make the
RNA-seq analysis smooth and easy. The package integrates the
state-of-art tools from both alignment-based and alignment-free
quantification workflows. The BP4RNAseq package uses an optimized
pipeline, applies to both retrospective and newly generated RNA-seq data
analyses and can take only two nontechnical parameters and output two
formatted gene expression quantification at gene and transcript levels
when working with local FASTQ files, therefore, facilitating the
application of RNA-seq.

### Operating System Requirements

BP4RNAseq runs in Windows (Subsystem for Linux), Linux and macOS.

### Dependencies

The BP4RNAseq requires the following utilities:

  - conda=4.83
  - SRA Toolkit=2.10.3
  - Entrez Direct=13.3
  - FastQC=v0.11.9
  - Cutadapt=2.10
  - datasets(Beta)
  - SAMtools=1.9
  - HISAT2=2.2.0
  - StringTie=2.1.1
  - Salmon=1.2.1  
  - R (\>= 3.5.0)

Users can install all the dependencies in the working directory with the
following commands:

``` r
cd WorkingDirectory ### change 'WorkingDirectory' to the the actual folder that you want to work in
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/sunshanwen/BP4RNAseq/master/install_depends.sh
./install_depends.sh
```

### Installation

<!-- You can install the released version of BP4RNAseq from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- #install.packages("BP4RNAseq") # remove comments later -->

<!-- ``` -->

<!-- And the development version from [GitHub](https://github.com/) with: -->

<!-- ``` r -->

<!-- # install.packages("devtools") -->

<!-- devtools::install_github("sunshanwen/BP4RNAseq") -->

<!-- ``` -->

You can install BP4RNAseq from [GitHub](https://github.com/) with:

``` r
devtools::install_github("sunshanwen/BP4RNAseq")
```

### Usage

The functions in BP4RNAseq are integrated into two main functions:
down2quan for public RNA-seq data, fastq2quan for newly generated
RNA-seq data.

down2quan can receive only two nontechnical parameters. The parameter
‘accession’ specifies the accession id of the target public RNA-seq
data in the Gene Expression Omnibus (GEO) or the Sequence Read Archive
(SRA). The accession id can be of a whole ‘BioProject’ or multiple
‘BioSample’. The parameter ‘taxa’ offers the scientific or common name
of the organism investigated. A simple example

``` r
library(BP4RNAseq)
down2quan(accession=c("SRR11486115","SRR11486114"), taxa="Drosophila melanogaster")
```

will download the public RNA-seq data from two ‘BioSample’ with
accession id “SRR11486115” and “SRR11486114”, respectively, and the
latest reference genome, transcript and annotation data of Drosophila
melanogaster, do the quality control (filter out the poor-quality reads
and contaminations), reads alignments and gene expression quantification
based on both alignment-free and alignment-based workflows in the
working directory. During the quality control procedure, if the
contamination of the adapter exists the program will automatically
detect the adapter sequence. However, an option is given to the users to
provide the adapter sequence to trim before the trimming process.

fastq2quan works with local RNA-seq data in fastq formats. It needs two
nontechnical parameters at a minimum, i.e., ‘taxa’ as explained above
and ‘pair’ which specifies the sequencing type with ‘single’ for
single-end (SE) reads or ‘paired’ for paired-end (PE) reads. Users
should place all the fastq files in the working directory. A simple
example

``` r
library(BP4RNAseq)
fastq2quan(taxa="Drosophila melanogaster", pair = "single")
```

will download the latest reference genome, transcript and annotation
data of Drosophila melanogaster, do the quality control, reads
alignments and gene expression quantification based on both
alignment-free and alignment-based workflows in the working directory as
the program down2quan do.

Both programs can do the parallel computing, which is specified by the
‘threads’ parameter.
