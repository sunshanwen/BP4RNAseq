
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BP4RNAseq

<!-- badges: start -->

<!-- badges: end -->

The goal of BP4RNAseq is to make the RNA-seq analysis smooth and easy.
The package integrates the state-of-art tools from both alignment-based
and alignment-free quantification workflows. The BP4RNAseq package uses
an optimized pipeline, applies to both retrospective and newly generated
RNA-seq data analyses and can take only two nontechnical parameters and
output two formatted gene expression quantification at gene and
transcript levels when working with local FASTQ files, therefore,
facilitating the application of RNA-seq.

Many of its dependencies only works under Unix environment, therefore,
so does the package.

## Installation

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
# install.packages("devtools")
devtools::install_github("sunshanwen/BP4RNAseq")
```

and its dependencies in the working directory with:

``` r
./install_depends.sh
```

## Example

This is a basic example which shows you how to use the package:

``` r
# library(BP4RNAseq)
## basic example code to work with public RNA-seq data. This will download about 86 GB RNA-seq data and quantify the transcriptional consequences of the deletions of sulfur metabolism genes in Drosophila melanogaster
# down2quan(accession="PRJNA623275", taxa="Drosophila melanogaster")

# basic example code to work with local RNA-seq data
# fastq2quan(taxa="Drosophila melanogaster", pair = "paired")
```

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- ```{r cars} -->

<!-- summary(cars) -->

<!-- ``` -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. -->

<!-- You can also embed plots, for example: -->

<!-- ```{r pressure, echo = FALSE} -->

<!-- plot(pressure) -->

<!-- ``` -->

<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub! -->
