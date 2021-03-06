% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/align_free_quan.R
\name{align_free_quan}
\alias{align_free_quan}
\title{Alignment-free expression quantification with Salmon at gene and transcript levels.}
\usage{
align_free_quan(
  pair,
  genome,
  transcript,
  annotation,
  threads = 4,
  salmon_index_add = NULL,
  salmon_quan_add = NULL
)
}
\arguments{
\item{pair}{'single' for single-end (SE) or 'paired' for paired-end (PE) reads.}

\item{genome}{the reference genome.}

\item{transcript}{the reference transcript}

\item{annotation}{the annotation file.}

\item{threads}{the number of threads to be used. Default is 4.}

\item{salmon_index_add}{additional parameters to customize salmon index command. Default is NULL.}

\item{salmon_quan_add}{additional parameters to customize salmon quan command. Default is NULL.}
}
\value{
None
}
\description{
Alignment-free expression quantification with Salmon at gene and transcript levels.
}
\examples{

align_free_quan(
pair = 'paired', 'test_Homo_sapiens.fna', 
'test_transcript_Homo_sapiens.fna','test_Homo_sapiens.gff'
)



}
