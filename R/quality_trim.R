#' Trim poor quality samples using cutadapt.
#'
#' Trim poor quality samples using cutadapt based on the results from quality_test().
#' @param per_base samples with poor per base sequence quality.
#' @param pair "single" for single-end (SE) or "paired" for paired-end (PE) reads.
#' @param per_seq samples with poor per sequence quality scores.
#' @param threads the number of threads to be used. Default is 4.
#' @export quality_trim
#' @return None
#' @examples
#' \dontrun{
#'
#' quality_trim(per_base, per_seq, pair, threads = 4)
#' }

quality_trim <- function(per_base, per_seq, pair, threads = 4)
{
  files <- unique(gsub("\\.fastq","",unique(cbind(per_base, per_seq))))
  for(f in files)
  {
    if(pair == "paired"){
      file <- c(paste0(f, "_1_fastqc"), paste0(f, "_2_fastqc"))
      qc_reports <- list.files(getwd(), pattern = "fastqc.zip$", full.names = F)
      qc_collection <- fastqcr::qc_read_collection(qc_reports, sample_names = gsub(".zip","",qc_reports), modules = "all", verbose = TRUE)
      Encoding <- unique(qc_collection$basic_statistics[qc_collection$basic_statistics$Measure == "Encoding"&qc_collection$basic_statistics$sample == file, 3])

      infile1 <- paste0(f, "_1.fastq")
      infile2 <- paste0(f, "_2.fastq")

      outfile1 <- paste0("Trimmed_",f, "_1.fastq")
      outfile2 <- paste0("Trimmed_",f, "_2.fastq")

      if(Encoding == "Sanger / Illumina 1.9"){
        cmd = paste("cutadapt -q 10", "-o", outfile1, "-p", outfile2,"-m 20","-j ", threads, infile1, infile2)
        # cat(cmd,"\n")#print the current command
        system(cmd)
      }else if(Encoding == "Illumina 1.3" | Encoding == "Illumina 1.5"){
        cmd = paste("cutadapt -q 10", "--quality-base=64", "-o", outfile1, "-p", outfile2,"-m 20", "-j ", threads, infile1, infile2)
        # cat(cmd,"\n")#print the current command
        system(cmd)
      }else stop("Wrong Encoding. Exit.")
    }else if(pair == "single"){
      file <- paste0(f, "_fastqc")
      qc_reports <- list.files(getwd(), pattern = "fastqc.zip$", full.names = F)
      qc_collection <- fastqcr::qc_read_collection(qc_reports, sample_names = gsub(".zip","",qc_reports), modules = "all", verbose = TRUE)
      Encoding <- unique(qc_collection$basic_statistics[qc_collection$basic_statistics$Measure == "Encoding"&qc_collection$basic_statistics$sample == file, 3])

      infile <- paste0(f, ".fastq")
      outfile <- paste0("Trimmed_",f, ".fastq")

      if(Encoding == "Sanger / Illumina 1.9"){
        cmd = paste("cutadapt -q 10", "-o", outfile, "-m 20","-j ", threads, infile)
        # cat(cmd,"\n")#print the current command
        system(cmd)
      }else if(Encoding == "Illumina 1.3" | Encoding == "Illumina 1.5"){
        cmd = paste("cutadapt -q 10", "--quality-base=64", "-o", outfile,"-m 20", "-j ", threads, infile)
        # cat(cmd,"\n")#print the current command
        system(cmd)
      }else stop("Wrong Encoding. Exit.")
    }
  }
}
