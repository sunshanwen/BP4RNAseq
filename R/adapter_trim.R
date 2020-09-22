## Get the adapter based on FASTQC report.

get_adapter <- function(sample, pair) {
    adapter <- NULL
    adap_3or5 <- NULL
    answer <- readline(prompt = "Do you know the adapter (y/n)? ")
    if(tolower(answer) == "y") {
        if (pair == "paired") {
            sequence1 <- readline(prompt = "Please enter the adapter sequence for read1: ")
            sequence2 <- readline(prompt = "Please enter the adapter sequence for read2: ")
            adapter <- c(sequence1, sequence2)
        } else if (pair == "single") {
            sequence <- readline(prompt = "Please enter the adapter sequence: ")
            adap_3or5 <- readline(prompt = "Please enter the adapter type (3 for 3' adapter types, 5 for 5' adapter types, 0 for Don't know): ")
            adapter <- sequence
        }
    } else {
        adapter_list = list(
            "Illumina Universal Adapter" = "AGATCGGAAGAG",
            "Illumina Small RNA 3' Adapter" = "TGGAATTCTCGG",
            "Nextera Transposase Sequence" = "CTGTCTCTTATA",
            "SOLID Small RNA Adapter" = "CGCCTTGGCCGT",
            "Illumina Small RNA 5' Adapter" = "GATCGTCGGACT")
        pattern <- paste0(gsub("\\.fastq", "", sample), "_fastqc.zip")
        qc_reports <- list.files(pattern = pattern, full.names = FALSE)
        if(length(qc_reports)) {
            qc_collection <- fastqcr::qc_read_collection(
                qc_reports, sample_names = gsub(".zip", "", qc_reports),
                modules = "all", verbose = TRUE
            )
            if (length(sample)) {
                dat <- qc_collection$adapter_content
                tt <- apply(dat[, !(names(dat) %in% c("sample", "Position"))] > 10,
                            2, sum) > 0
                adapter_name <- names(tt)[tt]
                adapter <- adapter_list[[adapter_name]]
                if (adapter_name == "Illumina Small RNA 5' Adapter") {
                    adap_3or5 <- "5"
                } else {
                    adap_3or5 <- "3"
                }
            }
        }
    }
    return(list(adapter = adapter, adap_3or5 = adap_3or5))
}

#' Trim samples with adapter using cutadapt.
#'
#' Trim samples with adapter using cutadapt based on the results from quality_test().
#' @param sample samples with adapter.
#' @param pair 'single' for single-end (SE) or 'paired' for paired-end (PE) reads.
#' @param threads the number of threads to be used. Default is 4.
#' @param scRNA logic, whether single-cell RNA-seq is quantified or not. Default is FALSE.
#' @param cutadapt_add additional parameters to customize cutadapt to trim adapter. Default is NULL.
#' @export adapter_trim
#' @return None
#' @examples
#'
#' adapter_trim(sample = 'test_read_2.fastq', pair = 'paired', scRNA = TRUE)
#'
adapter_trim <- function(sample, pair, threads = 4, scRNA = FALSE, cutadapt_add = NULL) {
    existence <- check_dep(dependancy = "cutadapt")
    if (existence == TRUE) {
        if (scRNA == FALSE) {
            adapter_got = get_adapter(sample, pair)
            samples <- unique(
                gsub("_2$", "", gsub("_1$", "", gsub("\\.fastq", "", sample)))
            )
            if (pair == "paired") {
                if (length(adapter_got$adapter)) {
                    for (i in samples) {
                        outfile1 <- paste0("Adapter_trimmed_", i, "_1.fastq")
                        outfile2 <- paste0("Adapter_trimmed_", i, "_2.fastq")
                        infile1 <- paste0("Trimmed_", i, "_1.fastq")
                        infile2 <- paste0("Trimmed_", i, "_2.fastq")
                        if (length(adapter_got$adapter) == 1) {
                            cmd = paste(
                                "-a", adapter_got$adapter, "-A", adapter_got$adapter, 
                                "-o", outfile1, "-p", outfile2, 
                                "-j ", threads, infile1, infile2, cutadapt_add)
                            system2(command = "cutadapt", args = cmd)
                            unlink(c(infile1, infile2))
                        } else if (length(adapter_got$adapter) == 2) {
                            arguments = paste(
                                "-a", adapter_got$adapter[1], 
                                "-A", adapter_got$adapter[2], 
                                "-o", outfile1, 
                                "-p", outfile2, 
                                "-j ", threads, 
                                infile1, infile2, 
                                cutadapt_add
                            )
                            system2(command = "cutadapt", args = arguments)
                            unlink(c(infile1, infile2))
                        } else stop("Wrong adapter. Exit.")
                    }
                }
            } else if (pair == "single") {
                infile <- paste0("Trimmed_", i, ".fastq") 
                outfile <- paste0("Adapter_trimmed_", i, ".fastq")
                if (length(adapter_got$adap_3or5)) {
                    if (adapter_got$adap_3or5 == "3") {
                        cmd = paste(
                            "-a", adapter_got$adapter, 
                            "-o", outfile, 
                            "-j ", threads, infile, cutadapt_add
                        )
                        system2(command = "cutadapt", args = cmd)
                        unlink(infile)
                    } else if (adapter_got$adap_3or5 == "5") {
                        cmd = paste(
                            "-g", adapter_got$adapter, 
                            "-o", outfile, 
                            "-j ", threads, infile, cutadapt_add
                        )
                        system2(command = "cutadapt", args = cmd)
                        unlink(infile)
                    } else if (adapter_got$adap_3or5 == "0") {
                        cmd = paste(
                            "-b", adapter_got$adapter, 
                            "-o", outfile, 
                            "-j ", threads, infile, cutadapt_add
                        )
                        system2(command = "cutadapt", args = cmd)
                        unlink(infile)
                    } else stop("Wrong selection. Exit.")    
                }
            }
        } else if (scRNA == TRUE) {
            samples <- unique(gsub("\\.fastq", "", sample))
            adapter_got = get_adapter(sample, pair = 'single')
            for (i in samples) {
                infile <- paste0("Trimmed_", i, ".fastq")  
                outfile <- paste0("Adapter_trimmed_", i, ".fastq")
                if (length(adapter_got$adap_3or5)) {
                    if (adapter_got$adap_3or5 == "3") {
                        cmd = paste(
                            "-a", adapter_got$adapter, 
                            "-o", outfile, 
                            "-j ", threads, infile, cutadapt_add
                        )
                        system2(command = "cutadapt", args = cmd)
                        unlink(infile)
                    } else if (adapter_got$adap_3or5 == "5") {
                        cmd = paste(
                            "-g", adapter_got$adapter, 
                            "-o", outfile, 
                            "-j ", threads, infile, cutadapt_add
                        )
                        system2(command = "cutadapt", args = cmd)
                        unlink(infile)
                    } else if (adapter_got$adap_3or5 == "0") {
                        cmd = paste(
                            "-b", adapter_got$adapter, 
                            "-o", outfile, 
                            "-j ", threads, infile, cutadapt_add
                        )
                        system2(command = "cutadapt", args = cmd)
                        unlink(infile)
                    } else stop("Wrong selection. Exit.")                    
                }
            }
        }
    } else print("Cutadapt is not found. Please install it.")
}
