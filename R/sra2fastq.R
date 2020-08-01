

#' @order 2
check_paired <- function(files)
{
  temp <- data.frame(samples = NA, paired_or_single = NA)
  for (f in files) {
    cmd = paste("fastq-dump", "-X", "1", "-Z", "--split-spot", f, "| wc -l")
    # cat(cmd,"\n")#print the current command
    dd <- system(cmd, intern = TRUE) # invoke command\
    # print(dd)
    if (dd == "8") {
      temp <- rbind(temp, c(f, "paired"))
    } else if (dd == "4") {
      temp <- rbind(temp, c(f, "single"))
    } else
      temp <- rbind(temp, c(f, "single"))
  }
  temp <- temp[-1, ]
  temp$paired_or_single <- as.factor(temp$paired_or_single)
  # print(temp)
  if (levels(temp$paired_or_single) == "paired") {
    return("paired")
  } else if (levels(temp$paired_or_single) == "single") {
    return("single")
  } else
    stop("Paired-end and single-end mix. Please check the data source!")
}

#' Convert SRA format to fastq format and return the sequencing type, i.e., single-end (SE) or paired-end (PE) reads.
#'
#' Convert all SRA files to fastq files and return the sequencing type, i.e., single-end (SE) or paired-end (PE) reads.
#' @param threads the number of threads to be used. Default is 4.
#' @import stringr
#' @export sra2fastq
#' @return A string indicates single-end (SE) or paired-end (PE) reads.
#' @order 1
#' @examples
#'
#' sra_download(accession = "SRR11427582")
#' sra2fastq()
#'
sra2fastq <- function(threads = 4)
{
  sra.dir = getwd() ### get the directory to search SRA files
  fq.dir = getwd()
  files <-
    list.files(
      sra.dir,
      pattern = ".sra$",
      recursive = FALSE,
      full.names = FALSE
    )
  if(length(files)){
    pair <- check_paired(files)
    # files_sra <- stringr::str_subset(files,".sra")
    for (f in files) {
      cmd = paste("fasterq-dump", f, "-O", fq.dir, "-e", threads)
      # cat(cmd,"\n")#print the current command
      system(cmd) # invoke command
      unlink(f)
    }

    fastq_files <-
      list.files(
        sra.dir,
        pattern = "fastq$",
        recursive = FALSE,
        full.names = FALSE
      )
    file.rename(fastq_files, gsub(".sra", "", fastq_files))
    unlink("sra", recursive = TRUE)
    return(pair)
  } else print("No SRA files. Please download them.")
}
