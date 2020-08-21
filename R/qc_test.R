#' @order 2
### produce fastqc reports
fastqc_r <- function(threads, fq.dir, scRNA, fastqc_add)
{
  if (scRNA == FALSE) {
    files_fastq <-
      list.files(
        fq.dir,
        pattern = ".fastq$",
        recursive = FALSE,
        full.names = TRUE
      )
    if(length(list.files)){
      for (f in files_fastq)
        {
          cmd = paste("fastqc ", f, " --threads ", threads, fastqc_add)
          # cat(cmd,"\n")#print the current command
          system(cmd) # invoke command
        }
    } else print("No fastq files are found in the work directory.")
  } else if (scRNA == TRUE) {
    read <- list.files(pattern = ".*1\\.fastq$", full.names = FALSE)
    read_seq <- c()
    barcode_seq <- c()
    if(length(read)){
      for (f in read)
      {
        name <- gsub("_1.fastq", "", f)
        read1 <- paste0(name, "_1.fastq")
        read2 <- paste0(name, "_2.fastq")
        read_seq <- c(read_seq, read2)
        barcode_seq <- c(barcode_seq, read1)
        cmd1 <-
          paste0("head -n 1 ", read1, " | grep -o length=[0-9]* | cut -d '=' -f 2")
        leg_1 <- as.numeric(system(cmd1, intern = TRUE))
        cmd2 <-
          paste0("head -n 1 ", read2, " | grep -o length=[0-9]* | cut -d '=' -f 2")
        leg_2 <- as.numeric(system(cmd2, intern = TRUE))
        if (leg_1 > leg_2) {
          read_seq <- c(read_seq, read1)
          barcode_seq <- c(barcode_seq, read2)
        }
      }
      for (d in read_seq)
      {
        cmd = paste("fastqc ", d, " --threads ", threads, fastqc_add)
        # cat(cmd,"\n")#print the current command
        system(cmd) # invoke command
      }
    } else print("No fastq files are found in the work directory.")
  }
}

utils::globalVariables(c("module", "status"))


#' Quality assessments of raw sequence data.
#'
#' Per base sequence quality, per sequence quality scores and adapter content of raw sequence data are assessed and tested based on FastQC.
#' @param threads the number of threads to be used. Default is 4.
#' @param scRNA logic, whether single-cell RNA-seq is quantified or not. Default is FALSE.
#' @param fastqc_add additional parameters to customize FASTQC for quality control. Default is NULL.
#' @export qc_test
#' @order 1
#' @return Problematic samples' names
#' @examples
#'
#' qc_test(scRNA = TRUE)
#'
#'

qc_test <- function(threads = 4, scRNA = FALSE, fastqc_add = NULL)
{
  status <- tryCatch(
    system2(command = "which", args = "fastqc", stdout = FALSE, stderr = FALSE),
    error = function(err){
      1
    },
    warning = function(war){
      2
    }
  )
  if(status == 0){
    fq.dir = getwd()
    qc.dir = getwd()
    fastqc_r(threads, fq.dir, scRNA, fastqc_add)
    qc <- fastqcr::qc_aggregate(qc.dir, progressbar = FALSE)
    if(length(qc)){
      index_ad <-
        subset(qc, module == "Adapter Content" &
                 status == "FAIL", select = sample)
      index_sq <-
        subset(qc,
               module == "Per base sequence quality" &
                 status == "FAIL",
               select = sample)
      index_sqs <-
        subset(qc,
               module == "Per sequence quality scores" &
                 status == "FAIL",
               select = sample)
      index <- list(index_ad, index_sq, index_sqs)
      return (index)
    }
  } else print("FASTQC is not found. Please install it.")
}
