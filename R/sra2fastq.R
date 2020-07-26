
#' @order 2
check_paired <- function(sra.dir)
{
  files<-list.files(sra.dir, pattern = ".sra$", recursive = F, full.names = F)
  temp <- data.frame(samples = NA, paired_or_single = NA)
  for(f in files) {
    cmd = paste("fastq-dump","-X","1","-Z","--split-spot", f, "| wc -l")
    # cat(cmd,"\n")#print the current command
    dd <- system(cmd, intern = TRUE) # invoke command\
    # print(dd)
    if(dd == "8") {
      temp <- rbind(temp, c(f, "paired"))
    } else if(dd == "4") {
      temp <- rbind(temp, c(f, "single"))
    } else temp <- rbind(temp, c(f, "single"))
  }
  temp <- temp[-1,]
  temp$paired_or_single <- as.factor(temp$paired_or_single)
  # print(temp)
  if(levels(temp$paired_or_single) == "paired"){
    return("paired")
  } else if(levels(temp$paired_or_single) == "single"){
    return("single")
  } else stop("Paired-end and single-end mix. Please check the data source!")
}

#' Convert SRA format to fastq format and return the sequencing type, i.e., single-end (SE) or paired-end (PE) reads.
#'
#' Convert all SRA files to fastq files and return the sequencing type, i.e., single-end (SE) or paired-end (PE) reads.
#' @param threads the number of threads to be used. Default is 4.
#'
#' @export sra2fastq
#' @return A string indicates single-end (SE) or paired-end (PE) reads.
#' @order 1
#' @examples
#' \dontrun{
#'
#' sra2fastq()
#' }
sra2fastq <- function(threads = 4)
{
  sra.dir=getwd() ### get the directory to search SRA files
  fq.dir = getwd()

  pair <- check_paired(sra.dir)
  #dir.create(file.path(fq.dir), showWarnings = FALSE) ### create the folder no matter it
  files<-list.files(sra.dir, pattern = ".sra$", recursive = F, full.names = F)
  # files_sra <- stringr::str_subset(files,".sra")
  for(f in files) {
    cmd = paste("fasterq-dump", f, "-O", fq.dir, "-e", threads)
    # cat(cmd,"\n")#print the current command
    system(cmd) # invoke command
    unlink(f)
  }

  fastq_files<-list.files(sra.dir, pattern = "fastq$", recursive = F, full.names = F)
  file.rename(fastq_files, gsub(".sra","",fastq_files))
  unlink("sra", recursive = TRUE)
  return(pair)
}
