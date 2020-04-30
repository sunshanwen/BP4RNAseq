#' @order 2
### produce fastqc reports
fastqc_r <- function(threads, fq.dir, qc.dir)
{
  dir.create(file.path(qc.dir), showWarnings = FALSE) ### create the folder
  files_fastq<-list.files(fq.dir, pattern = ".fastq$", recursive = T, full.names = T)
  for(f in files_fastq)
  {
    cmd = paste("fastqc ", f, " --outdir=", fq.dir, " --threads ", threads)
    # cat(cmd,"\n")#print the current command
    system(cmd) # invoke command
  }
}


#' Quality assessments of raw sequence data.
#'
#' Per base sequence quality, per sequence quality scores and adapter content of raw sequence data are assessed and tested based on FastQC.
#' @param threads the number of threads to be used. Default is 4.
#' @export qc_test
#' @order 1
#' @return Problematic samples' names
#' @examples
#' \dontrun{
#'
#' qc_test()
#' }

qc_test <- function(threads = 4)
{
  fq.dir = getwd()
  qc.dir = getwd()
  fastqc_r(threads, fq.dir, qc.dir)
  qc <- fastqcr::qc_aggregate(qc.dir, progressbar = FALSE)
  index_ad <- subset(qc, module == "Adapter Content" & status == "Fail", select = sample)
  index_sq <- subset(qc, module == "Per base sequence quality" & status == "Fail", select = sample)
  index_sqs <- subset(qc, module == "Per sequence quality scores" & status == "Fail", select = sample)
  index <- list(index_ad, index_sq, index_sqs)
  return (index)
}
