
index_build <- function(taxa, genome, annotation)
{
  cmd1 <- paste("awk \'{if ($3==\"exon\") {print $1\"\\t\"$4-1\"\\t\"$5-1}}\'", annotation, "> exonsFile.table")
  system(cmd1)
  cmd2 <- paste("awk \'{if ($3==\"intron\") {print $1\"\\t\"$4-1\"\\t\"$5-1\"\\t\"$7}}\'", annotation, "> ssFile.table")
  system(cmd2)
  cmd3 <- paste("hisat2-build -f", genome, taxa, "--ss ssFile.table", "--exon exonsFile.table")
  system(cmd3)
}

#' Aligning the RNA-seq data to the reference genome with HISAT2.
#'
#' Aligning the RNA-seq data to the reference genome with HISAT2.
#' @param pair "single" for single-end (SE) or "paired" for paired-end (PE) reads.
#' @param taxa the scientific or common name of the organism.
#' @param genome the reference genome.
#' @param annotation the annotation file.
#' @export align_ge
#' @return None
#' @examples
#' \dontrun{
#'
#' align_ge("paired", "sesame", "sesame.fna", "sesame.gff")
#'}
#'
align_ge <- function(pair, taxa, genome, annotation)
{
  index_build(taxa, genome, annotation)
  index <- list.files(pattern = "ht2$", recursive = T, full.names = T)
  if(pair == "paired")
  {
    read <- list.files(pattern = "^Trimmed.*1\\.fastq$", full.names = F)
    if(length(read) == 0){
      read <- list.files(pattern = ".*1\\.fastq$", full.names = F)
    }
    # read2 <- list.files(pattern = "^Trimmed.*2\\.fastq$", full.names = F)
    for(f in read)
    {
      # read1 <- paste(read1, sep = ",", collapse = ',')
      # read2 <- paste(read2, sep = ",", collapse = ',')
      name <- gsub("_1.fastq", "", f)
      out_bam <- paste0(name, ".bam")
      read1 <- paste0(name, "_1.fastq")
      read2 <- paste0(name, "_2.fastq")
      cmd4 <- paste("hisat2 --dta -x", taxa, "-1", read1, "-2", read2,
                    "| samtools view -bh - | samtools sort - >", out_bam)
      # cat(cmd4, "\n")
      system(cmd4, intern = T)
    }

  } else if(pair == "single"){
    read <- list.files(pattern = "^Trimmed.*\\.fastq$", full.names = F)
    if(length(read) == 0){
      read <- list.files(pattern = ".*\\.fastq$", full.names = F)
    }
    # read <- paste(read, sep = ",", collapse = ',')
    for(f in read)
    {
      name <- gsub(".fastq", "", f)
      out_bam <- paste0(name, ".bam")
      cmd4 <- paste("hisat2 --dta -x", taxa,"-U", f,
                    "| samtools view -bh - | samtools sort - >", out_bam)
      # cat(cmd4, "\n")
      system(cmd4, intern = T)
    }

  } else stop("Paired-end and single-end mix. Please check the data source!")
}
