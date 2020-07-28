#' Produce the quantifications of gene expression based on the fastq files.
#'
#' Produce the quantifications of gene expression based on the fastq files with alignment-based and alignment-free workflows.
#'
#'@param threads the number of threads to be used. Default is 4.
#'@param dir the working directory. Default is the current working directory.
#'@param pair "single" for single-end (SE) reads or "paired" for paired-end (PE) reads.
#'@param taxa the scientific or common name of the organism.
#'@param novel_transcript logic, whether identifying novel transcripts is expected or not. Default is FALSE.
#'@param scRNA logic, whether single-cell RNA-seq is quantified or not. Default is FALSE.
#'@param protocol the single-cell RNA sequencing protocol: dropseq, chromium, or chromiumV3.
#'@export fastq2quan
#'@return None
#'@examples
#'\dontrun{
#'
#'down2quan(threads = 4, dir = getwd(), "sesame", novel_transcript = FALSE, scRNA = FALSE)
#'}

fastq2quan <- function(threads = 4, dir = getwd(), pair = "paired", taxa, novel_transcript = FALSE, scRNA = FALSE, protocol) {
  setwd(dir)

  quality <- qc_test(threads, scRNA)

  if(length(quality[[1]]$sample)) cat("Adapter exists!\n")
  if(length(quality[[2]]$sample)) cat("Per base sequence quality failed!\n")
  if(length(quality[[3]]$sample)) cat("Per sequence quality scores failed!\n")

  ### quality trimming
  if(length(quality[[2]]$sample)|length(quality[[3]]$sample)) {
    print("Trim low quality bases.")
    quality_trim(quality[[2]]$sample, quality[[3]]$sample, pair, scRNA)### consider if add directory parameter.
    if(length(quality[[1]]$sample)) {
      print("Trim the adapter.")
      adapter_trim(quality[[1]]$sample, pair, scRNA)
    }
  }

  if(length(quality[[1]]$sample)) {
    print("Trim the adapter.")
    quality_trim(quality[[1]]$sample, quality[[1]]$sample, pair, scRNA )### consider if add directory parameter.
    adapter_trim(quality[[1]]$sample, pair, scRNA)
  }
  files <- list.files(dir, pattern = "fastqc", full.names = F)
  unlink(files)

  status <- down_Ref(taxa)
  if(status == 1)
  {
    stop("The download of reference genome and annotation files failed! Terminate the program.")
    break
  }
  # reference <- extract_genome(taxa)
  taxa_tmp <- gsub("\\s", "_", taxa)
  genome <- paste0(taxa_tmp, ".fna")
  transcript <- paste0("transcript_", taxa_tmp, ".fna")
  annotation <- paste0(taxa_tmp, ".gff")
  if(scRNA == FALSE){
    align_ge(pair, taxa, genome, annotation)
    trans_ass(novel_transcript)
    trans_quan()
    align_free_quan(pair, genome, transcript, annotation)
    average()
  } else if (scRNA == TRUE){
    scRNA(transcript, protocol)
  }
}
