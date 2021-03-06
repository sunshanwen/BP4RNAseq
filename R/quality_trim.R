#' Trim poor quality samples using cutadapt.
#'
#' Trim poor quality samples using cutadapt based on the results from quality_test().
#' @param per_base samples with poor per base sequence quality.
#' @param pair 'single' for single-end (SE) or 'paired' for paired-end (PE) reads.
#' @param per_seq samples with poor per sequence quality scores.
#' @param threads the number of threads to be used. Default is 4.
#' @param scRNA logic, whether single-cell RNA-seq is quantified or not. Default is FALSE.
#' @param cutadapt_add additional parameters to customize cutadapt to trim poor quality reads. Default is NULL
#' @export quality_trim
#' @return None
#' @examples
#'
#' quality_trim(per_base = 'test_read_2.fastq', per_seq = 'test_read_2.fastq',
#'              pair = 'paired', threads = 4, scRNA = TRUE)
#'

quality_trim <- function(per_base, per_seq, pair, threads = 4, scRNA = FALSE, cutadapt_add = NULL) {
    if (scRNA == FALSE) {
        files <- unique(gsub("\\.fastq", "", unique(cbind(per_base, per_seq))))
        files <- unique(gsub("_1$", "", files))
        files <- unique(gsub("_2$", "", files))
        files <- unique(files)
        for (f in files) {
            if (pair == "paired") {
                file <- c(paste0(f, "_1_fastqc"), paste0(f, "_2_fastqc"))
                pattern <- paste0(f, "_1_fastqc.zip")
                qc_reports <- list.files(pattern = pattern, full.names = FALSE, recursive = FALSE)
                qc_collection <- fastqcr::qc_read_collection(qc_reports, sample_names = gsub(".zip", "", qc_reports), modules = "all", verbose = TRUE)
                Encoding <- unique(qc_collection$basic_statistics[qc_collection$basic_statistics$Measure == "Encoding" & qc_collection$basic_statistics$sample == 
                  file[1], 3])
                
                infile1 <- paste0(f, "_1.fastq")
                infile2 <- paste0(f, "_2.fastq")
                
                outfile1 <- paste0("Trimmed_", f, "_1.fastq")
                outfile2 <- paste0("Trimmed_", f, "_2.fastq")
                
                if (length(Encoding)) {
                  if (Encoding == "Sanger / Illumina 1.9") {
                    cmd = paste("-q 10", "-o", outfile1, "-p", outfile2, "-m 20", "-j ", threads, infile1, infile2)
                    system2(command = "cutadapt", args = cmd)
                  } else if (Encoding == "Illumina 1.3" | Encoding == "Illumina 1.5") {
                    cmd = paste("-q 10", "--quality-base=64", "-o", outfile1, "-p", outfile2, "-m 20", "-j ", threads, infile1, infile2)
                    system2(command = "cutadapt", args = cmd)
                  } else {
                    stop("Wrong Encoding. Exit.")
                  }
                }
            } else if (pair == "single") {
                file <- paste0(f, "_fastqc")
                pattern <- paste0(f, "_fastqc.zip")
                qc_reports <- list.files(pattern = pattern, full.names = FALSE, recursive = FALSE)
                qc_collection <- fastqcr::qc_read_collection(qc_reports, sample_names = gsub(".zip", "", qc_reports), modules = "all", verbose = TRUE)
                Encoding <- unique(qc_collection$basic_statistics[qc_collection$basic_statistics$Measure == "Encoding" & qc_collection$basic_statistics$sample == 
                  file, 3])
                
                infile <- paste0(f, ".fastq")
                outfile <- paste0("Trimmed_", f, ".fastq")
                
                if (length(Encoding)) {
                  if (Encoding == "Sanger / Illumina 1.9") {
                    cmd = paste("-q 10", "-o", outfile, "-m 20", "-j ", threads, infile)
                    system2(command = "cutadapt", args = cmd)
                  } else if (Encoding == "Illumina 1.3" | Encoding == "Illumina 1.5") {
                    cmd = paste("-q 10", "--quality-base=64", "-o", outfile, "-m 20", "-j ", threads, infile)
                    system2(command = "cutadapt", args = cmd)
                  } else {
                    stop("Wrong Encoding. Exit.")
                  }
                } else {
                  print("Could not get the Encoding.")
                }
            }
        }
        unlink(files)
    } else if (scRNA == TRUE) {
        files <- unique(gsub("\\.fastq", "", unique(c(per_base, per_seq))))
        for (f in files) {
            pattern <- paste0(f, "_fastqc.zip")
            qc_reports <- list.files(pattern = pattern, full.names = FALSE, recursive = FALSE)
            qc_collection <- fastqcr::qc_read_collection(qc_reports, sample_names = gsub(".zip", "", qc_reports), modules = "all", verbose = TRUE)
            file <- paste0(f, "_fastqc")
            Encoding <- unique(qc_collection$basic_statistics[qc_collection$basic_statistics$Measure == "Encoding" & qc_collection$basic_statistics$sample == 
                file, 3])
            infile <- paste0(f, ".fastq")
            outfile <- paste0("Trimmed_", f, ".fastq")
            
            if (length(Encoding)) {
                if (Encoding == "Sanger / Illumina 1.9") {
                  cmd = paste("-q 10", "-o", outfile, "-m 20", "-j ", threads, infile)
                  system2(command = "cutadapt", args = cmd)
                } else if (Encoding == "Illumina 1.3" | Encoding == "Illumina 1.5") {
                  cmd = paste("-q 10", "--quality-base=64", "-o", outfile, "-m 20", "-j ", threads, infile)
                  system2(command = "cutadapt", args = cmd)
                }
            } else {
                print("Could not get the Encoding.")
            }
        }
    }
}
