
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BP4RNAseq

The assessment of gene expression is central to uncovering the functions
of the genome, understanding the regulation of development and
investigating the molecular mechanisms that underlie cancer and other
diseases. RNA-sequencing (RNA-seq) now is the routine to assess the
genome wide gene expression due to its high speed, accuracy and
reproducibility, and low cost. An enormous volume of RNA-seq data have
been accumulating and deposited in public data repositories, such as the
Gene Expression Omnibus (GEO) and the Sequence Read Archive (SRA).
Retrospectively analyzing these data or conducting a brand new RNA-seq
study is fundamentally important for researchers. However, processing
raw reads of RNA-seq data, no matter public or newly sequenced data,
involves a lot of specialized tools and technical configurations that
are often unfamiliar and time-consuming to learn for non-bioinformatics
researchers. For example, when working with public RNA-seq data,
researchers need to download the RNA-seq data, convert data to FASTQ
format, check the sequencing type (i.e., single-end or pair-end), do the
quality control (when needed, trim adapters and poor quality reads),
download the reference genome, transcript and annotation file, align
reads to the reference genome or transcript and quantify gene
expression, etc. These steps and the details that they involve are even
tedious for bioinformatic scientist. The goal of BP4RNAseq is to make
the RNA-seq analysis smooth and easy and to minimize efforts from
researchers. The package offers several benefits to researchers. First,
the package is a highly automated tool. It can take only two
nontechnical parameters and output six formatted gene expression
quantification at gene and transcript levels. Second, it improves the
accuracy and sensitivity of RNA-seq analyses by using an optimized
pipeline. Third, it offers individual tools to provide users full
control to fine tune precisely how individual steps are optimized. This
can allow experienced users to further improve the accuracy and
sensitivity of RNA-seq analyses. Users can also use the package as a
toolbox to run the exact tools that suit their needs. Last but not
least, the package applies to both retrospective and newly generated
bulk RNA-seq data analyses and is also applicable for single-cell
RNA-seq analyses based on the Alevin algorithm \[1\].

### Operating System Requirements

BP4RNAseq runs in Windows (Subsystem for Linux), Linux and macOS.

### Dependencies

The BP4RNAseq requires the following utilities:

  - [SRA Toolkit
    \>= 2.10.3](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc)
  - [Entrez Direct
    \>= 13.3](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
  - [FastQC \>=
    v0.11.9](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  - [Cutadapt \>= 2.10](https://cutadapt.readthedocs.io/en/stable/)
  - [SAMtools \>= 1.9](http://www.htslib.org/)
  - [HISAT2 \>= 2.2.0](http://daehwankimlab.github.io/hisat2/)
  - [StringTie \>= 2.1.1](https://ccb.jhu.edu/software/stringtie/)
  - [Salmon \>= 1.2.1](https://combine-lab.github.io/salmon/)  
  - [jq \>= 1.6](https://stedolan.github.io/jq/)
  - [R \>= 4.0.0](https://www.r-project.org/)

Users can install these dependencies manually.

Alternatively, we provide a bash script to aid users to install all the
dependencies based on [conda](https://docs.conda.io/en/latest/). The
script uses Wget, which is pre-installed on most Linux distributions
such as Windows Subsystem for Linux, to download conda. If wget is not
installed, users can easily install it with the following commands.

  - Installing Wget on Ubuntu and Debian or Windows Subsystem
    equivalents

<!-- end list -->

``` r
sudo apt-get update 
sudo apt-get install -y wget
```

  - Installing Wget on CentOS and Fedora or Windows Subsystem
    equivalents

<!-- end list -->

``` r
sudo yum install wget
```

  - Installing Wget on macOS with Homebrew

<!-- end list -->

``` r
brew install wget
```

With Wget installed, users can install all the dependencies with the
following commands:

``` r
wget https://raw.githubusercontent.com/sunshanwen/BP4RNAseq/master/install_depends.sh
chmod +x install_depends.sh
./install_depends.sh
./use_conda.sh
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

#### Bulk RNA-seq analyses

The functions in BP4RNAseq are integrated into two main functions:
down2quan for public RNA-seq data, fastq2quan for newly generated
RNA-seq data.

down2quan requires no input data and can receive only two nontechnical
parameters. The parameter `accession` specifies the accession id of the
target public RNA-seq data in the Gene Expression Omnibus (GEO) or the
Sequence Read Archive (SRA). The accession id can be of a whole
‘BioProject’ or multiple ‘BioSample’. The parameter `taxa` offers the
scientific or common name of the organism investigated. A simple example

``` r
library(BP4RNAseq)
down2quan(accession=c("SRR11486115","SRR11486114"), taxa="Drosophila melanogaster")
```

will download the public RNA-seq data of two ‘BioSample’ with accession
id “SRR11486115” and “SRR11486114”, respectively, and the latest
reference genome, transcript and annotation data of Drosophila
melanogaster, do the quality control (filter out the poor-quality reads
and contaminations), reads alignments and gene expression quantification
based on both alignment-free and alignment-based workflows in the work
directory. During the quality control procedure, if the contamination of
the adapter exists, the program will automatically detect the adapter
sequence to trim. However, an option is given to the users to provide
the adapter sequence if they know it.

fastq2quan works with local RNA-seq data in fastq formats. It needs two
nontechnical parameters at a minimum, i.e., `taxa` as explained above
and `pair`which specifies the sequencing type with `single` for
single-end (SE) reads or `paired` for paired-end (PE) reads. Users
should place all the fastq files in the work directory. A simple example

``` r
library(BP4RNAseq)
fastq2quan(taxa="Drosophila melanogaster", pair = "single")
```

will download the latest reference genome, transcript and annotation
data of Drosophila melanogaster, and do the quality control, reads
alignments and gene expression quantification using the local RNA-seq
data based on both alignment-free and alignment-based workflows as the
program down2quan do.

Both programs support the parallel computing, which is specified by the
`threads` parameter.

Outputs from both functions are two gene count matrixes and two
transcript count matrixes based on the alignment-based workflow and the
alignment-free workflow, and corresponding average matrixes over two
workflows. These outputs can be directly processed with
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html),
[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
or
[limma](https://bioconductor.org/packages/release/bioc/html/limma.html).
Researchers may use the averages for downstream analyses \[2\].
Alternatively, we recommend to decide the type of data to use based on
their consistencies with qPCR results if available or/and the results
from the downstream analyses.

#### Single-cell RNA-seq analyses

down2quan and fastq2quan can also be extended to process single-cell
RNA-seq data by setting the `scRNA` parameter to be ‘TRUE’ and
specifying the protocols. Currently, dropseq, chromium and chromiumV3
are supported protocols. A simple example

``` r
library(BP4RNAseq)
down2quan(accession=c("SRR11402955","SRR11402974"), taxa="Homo sapiens", scRNA = TRUE, protocol = "dropseq")
```

will download the public single-cell RNA-seq data from two ‘BioSample’
with accession id “SRR11402955” and “SRR11402974”, respectively, and the
latest reference genome, transcript and annotation data of Homo sapiens,
do the quality control, reads alignments and gene expression
quantification based on the Alevin workflow.

Alternatively,

``` r
library(BP4RNAseq)
fastq2quan(taxa="Homo sapiens", scRNA = TRUE, protocol = "dropseq")
```

can preprocess local single-cell RNA-seq data in fastq formats. The data
are paired-end reads with one read containing cellular barcode and
unique molecule identifier (UMI) and the other read being the RNA
sequence.

The outputs of down2quan and fastq2quan are gene count matrix compressed
in binary format, and gene ids, barcode + UMI and tier categorization in
three separate files. These outputs can be further processed with
[tximport](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html)
and [Seurat](https://satijalab.org/seurat/).

#### Using BP4RNAseq as a toolbox for RNA-seq analyses and customizing gene expression quantification to improve the sensitivity and accuracy of RNA-seq analyses.

BP4RNAseq offers individual tools to users. These tools can allow users
to run the exact tools that suit their needs. Specifically, these tools
are:

  - sra\_download(), which can download RNA-seq data from public
    repositories.
  - sra2fastq(), which converts SRA files into FASTQ format.
  - qc\_test(), which assesses the quality of reads.
  - adapter\_trim(), which trims adapter.
  - quality\_trim(), which trims poor quality reads.
  - down\_Ref(), which downloads the latest reference genome, transcript
    and annotation files.
  - align\_based\_quan(), which performs the alignment-based gene
    expression quantification, including reads alignment and assembly.
  - align\_free\_quan(), which performs the alignment-free gene
    expression quantification.
  - scRNA\_quan(), which performs the single-cell RNA-seq gene
    expression quantification.

Additionally, these individual tools provide users full control to fine
tune precisely how individual steps are optimized. This can allow
experienced users to further improve the accuracy and sensitivity of
RNA-seq analyses. For example, setting the optional parameter
`salmon_quan_add` of align\_free\_quan() as `salmon_quan_add = "--useEM
--gcBias"` will allow users to apply the standard EM algorithm to
optimize abundance estimates and in the mean time to correct for
fragment-level GC biases in the input data when performing the
alignment-free workflow. Details about the optional customizing setting
in each tool can be found in package help page.

### References

1.  Srivastava, A., et al. Alevin efficiently estimates accurate gene
    abundances from dscRNA-seq data. Genome Biol 2019;20:16.

2.  Lachmann, A., et al. Interoperable RNA-Seq analysis in the cloud.
    Biochim. Biophys. Acta-Gene Regul. Mech. 2020;1863(6):1-11.
