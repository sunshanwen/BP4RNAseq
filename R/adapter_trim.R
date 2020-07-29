#' Get the adapter based on FASTQC report.
#'
#' Get the adapter based on FASTQC report produced by quality_test().
#' @inheritParams adapter_trim
#' @return adapter name
#' @order 2
get_adapter <- function(sample)
{
  qc_reports <-
    list.files(getwd(), pattern = "fastqc.zip", full.names = FALSE)
  qc_collection <- fastqcr::qc_read_collection(
    qc_reports,
    sample_names = gsub(".zip", "", qc_reports),
    modules = "all",
    verbose = TRUE
  )
  if (length(sample))
  {
    dat <- qc_collection$adapter_content
    tt <- apply(dat[, -c(1, 2)] > 10, 2, sum) > 0
    adapter_found <- names(tt)[tt]
    return(adapter_found)
  }
}

#' Trim samples with adapter using cutadapt.
#'
#' Trim samples with adapter using cutadapt based on the results from quality_test().
#' @param sample samples with adapter.
#' @param pair "single" for single-end (SE) or "paired" for paired-end (PE) reads.
#' @param threads the number of threads to be used. Default is 4.
#' @export adapter_trim
#' @param scRNA logic, whether single-cell RNA-seq is quantified or not. Default is FALSE.
#' @return None
#' @order 1
#' @examples
#' \dontrun{
#' adapter_trim(sample, pair)
#'}
adapter_trim <- function(sample,
                         pair,
                         threads = 4,
                         scRNA = FALSE)
{
  if (scRNA == FALSE) {
    samples <- unique(gsub("\\.fastq", "", sample))
    samples <- gsub("_1$", "", samples)
    samples <- gsub("_2$", "", samples)
    samples <- unique(samples)
    answer3 <- readline(prompt = "Do you know the adapter (y/n)?")
    if (answer3 == "y")
    {
      for (i in samples)
      {
        if (pair[pair$samples == i, 2] == "pair")
        {
          outfile1 <- paste0("Adapter_trimmed_", i, "_1.fastq")
          outfile2 <- paste0("Adapter_trimmed_", i, "_2.fastq")
          infile1 <- paste0("Trimmed_", i, "_1.fastq")
          infile2 <- paste0("Trimmed_", i, "_2.fastq")
          #### to continue
          sequence1 <-
            readline(prompt = "Please enter the adapter sequence for read1: ")
          sequence2 <-
            readline(prompt = "Please enter the adapter sequence for read2: ")
          cmd = paste(
            "cutadapt -a",
            sequence1,
            "-A",
            sequence2,
            "-o",
            outfile1,
            "-p",
            outfile2,
            "-j ",
            threads,
            infile1,
            infile2
          )
          #cat(cmd,"\n")#print the current command
          system(cmd)
          unlink(c(infile1, infile2))
        } else if (pair[pair$samples == i, 2] == "single") {
          infile <- paste0("Trimmed_", i, ".fastq")### maybe need to add "_"
          outfile <-
            paste0("Adapter_trimmed_", i, ".fastq")### maybe need to add "_"
          sequence <-
            readline(prompt = "Please enter the adapter sequence: ")
          adap_3or5 <-
            readline(prompt = "Please enter the adapter type (3 for 3' adapter types, 5 for 5' adapter types, 0 for Don't know): ")
          if (adap_3or5 == "3") {
            cmd = paste("cutadapt -a",
                        sequence,
                        "-o",
                        outfile,
                        "-j ",
                        threads,
                        infile)
            system(cmd)
            unlink(c(infile))
          } else if (adap_3or5 == "5") {
            cmd = paste("cutadapt -g",
                        sequence,
                        "-o",
                        outfile,
                        "-j ",
                        threads,
                        infile)
            system(cmd)
            unlink(c(infile))
          } else if (adap_3or5 == "b") {
            cmd = paste("cutadapt -b",
                        sequence,
                        "-o",
                        outfile,
                        "-j ",
                        threads,
                        infile)
            system(cmd)
            unlink(c(infile))
          } else
            stop("Wrong selection. Exit.")
        }
      }
      #}
    } else if (answer3 == "n") {
      adapter = get_adapter(sample)
      for (i in samples)
      {
        if (pair[pair$samples == i, 2] == "pair")
        {
          outfile1 <- paste0("Adapter_trimmed_", i, "_1.fastq")
          outfile2 <- paste0("Adapter_trimmed_", i, "_2.fastq")
          infile1 <- paste0("Trimmed_", i, "_1.fastq")
          infile2 <- paste0("Trimmed_", i, "_2.fastq")
          #### to continue
          if (adapter == "Illumina Universal Adapter") {
            cmd = paste(
              "cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG",
              "-o",
              outfile1,
              "-p",
              outfile2,
              "-j ",
              threads,
              infile1,
              infile2
            )
            #cat(cmd,"\n")#print the current command
            system(cmd)
            unlink(c(infile1, infile2))
          } else if (adapter == "Illumina Small RNA 3' Adapter") {
            cmd = paste(
              "cutadapt -a TGGAATTCTCGG -A TGGAATTCTCGG",
              "-o",
              outfile1,
              "-p",
              outfile2,
              "-j ",
              threads,
              infile1,
              infile2
            )
            #cat(cmd,"\n")#print the current command
            system(cmd)
            unlink(c(infile1, infile2))
          } else if (adapter == "Nextera Transposase Sequence") {
            cmd = paste(
              "cutadapt -a CTGTCTCTTATA -A CTGTCTCTTATA",
              "-o",
              outfile1,
              "-p",
              outfile2,
              "-j ",
              threads,
              infile1,
              infile2
            )
            #cat(cmd,"\n")#print the current command
            system(cmd)
            unlink(c(infile1, infile2))
          } else if (adapter == "SOLID Small RNA Adapter") {
            cmd = paste(
              "cutadapt -a CGCCTTGGCCGT -A CGCCTTGGCCGT",
              "-o",
              outfile1,
              "-p",
              outfile2,
              "-j ",
              threads,
              infile1,
              infile2
            )
            #cat(cmd,"\n")#print the current command
            system(cmd)
            unlink(c(infile1, infile2))
          } else if (adapter == "Illumina Small RNA 5' Adapter") {
            cmd = paste(
              "cutadapt -a GATCGTCGGACT -A GATCGTCGGACT",
              "-o",
              outfile1,
              "-p",
              outfile2,
              "-j ",
              threads,
              infile1,
              infile2
            )
            #cat(cmd,"\n")#print the current command
            system(cmd)
            unlink(c(infile1, infile2))
          } else
            stop("Wrong adapter. Exit.")
        } else if (pair[pair$samples == i, 2] == "single") {
          infile <- paste0("Trimmed_", i, ".fastq")### maybe need to add "_"
          outfile <-
            paste0("Adapter_trimmed_", i, ".fastq")### maybe need to add "_"
          if (adapter == "Illumina Universal Adapter") {
            cmd = paste("cutadapt -a AGATCGGAAGAG",
                        "-o",
                        outfile,
                        "-j ",
                        threads,
                        infile)
            #cat(cmd,"\n")#print the current command
            system(cmd)
            unlink(c(infile))
          } else if (adapter == "Illumina Small RNA 3' Adapter") {
            cmd = paste("cutadapt -a TGGAATTCTCGG",
                        "-o",
                        outfile,
                        "-j ",
                        threads,
                        infile)
            #cat(cmd,"\n")#print the current command
            system(cmd)
            unlink(c(infile))
          } else if (adapter == "Nextera Transposase Sequence") {
            cmd = paste("cutadapt -a CTGTCTCTTATA",
                        "-o",
                        outfile,
                        "-j ",
                        threads,
                        infile)
            #cat(cmd,"\n")#print the current command
            system(cmd)
            unlink(c(infile))
          } else if (adapter == "SOLID Small RNA Adapter") {
            cmd = paste("cutadapt -a CGCCTTGGCCGT",
                        "-o",
                        outfile,
                        "-j ",
                        threads,
                        infile)
            #cat(cmd,"\n")#print the current command
            system(cmd)
            unlink(c(infile))
          } else if (adapter == "Illumina Small RNA 5' Adapter") {
            cmd = paste("cutadapt -g GATCGTCGGACT",
                        "-o",
                        outfile,
                        "-j ",
                        threads,
                        infile)
            #cat(cmd,"\n")#print the current command
            system(cmd)
            unlink(c(infile))
          } else
            stop("Wrong adapter. Exit.")
        }
      }
    }
  } else if (scRNA == TRUE) {
    samples <- unique(gsub("\\.fastq", "", sample))
    answer3 <- readline(prompt = "Do you know the adapter (y/n)?")
    if (answer3 == "y")
    {
      for (i in samples)
      {
        infile <- paste0("Trimmed_", i, ".fastq")### maybe need to add "_"
        outfile <-
          paste0("Adapter_trimmed_", i, ".fastq")### maybe need to add "_"
        sequence <-
          readline(prompt = "Please enter the adapter sequence: ")
        adap_3or5 <-
          readline(prompt = "Please enter the adapter type (3 for 3' adapter types, 5 for 5' adapter types, 0 for Don't know): ")
        if (adap_3or5 == "3") {
          cmd = paste("cutadapt -a",
                      sequence,
                      "-o",
                      outfile,
                      "-j ",
                      threads,
                      infile)
          system(cmd)
          unlink(c(infile))
        } else if (adap_3or5 == "5") {
          cmd = paste("cutadapt -g",
                      sequence,
                      "-o",
                      outfile,
                      "-j ",
                      threads,
                      infile)
          system(cmd)
          unlink(c(infile))
        } else if (adap_3or5 == "b") {
          cmd = paste("cutadapt -b",
                      sequence,
                      "-o",
                      outfile,
                      "-j ",
                      threads,
                      infile)
          system(cmd)
          unlink(c(infile))
        } else
          stop("Wrong selection. Exit.")

      }
      #}
    } else if (answer3 == "n") {
      adapter = get_adapter(sample)
      for (i in samples)
      {
        infile <- paste0("Trimmed_", i, ".fastq")### maybe need to add "_"
        outfile <-
          paste0("Adapter_trimmed_", i, ".fastq")### maybe need to add "_"
        if (adapter == "Illumina Universal Adapter") {
          cmd = paste("cutadapt -a AGATCGGAAGAG",
                      "-o",
                      outfile,
                      "-j ",
                      threads,
                      infile)
          #cat(cmd,"\n")#print the current command
          system(cmd)
          unlink(c(infile))
        } else if (adapter == "Illumina Small RNA 3' Adapter") {
          cmd = paste("cutadapt -a TGGAATTCTCGG",
                      "-o",
                      outfile,
                      "-j ",
                      threads,
                      infile)
          #cat(cmd,"\n")#print the current command
          system(cmd)
          unlink(c(infile))
        } else if (adapter == "Nextera Transposase Sequence") {
          cmd = paste("cutadapt -a CTGTCTCTTATA",
                      "-o",
                      outfile,
                      "-j ",
                      threads,
                      infile)
          #cat(cmd,"\n")#print the current command
          system(cmd)
          unlink(c(infile))
        } else if (adapter == "SOLID Small RNA Adapter") {
          cmd = paste("cutadapt -a CGCCTTGGCCGT",
                      "-o",
                      outfile,
                      "-j ",
                      threads,
                      infile)
          #cat(cmd,"\n")#print the current command
          system(cmd)
          unlink(c(infile))
        } else if (adapter == "Illumina Small RNA 5' Adapter") {
          cmd = paste("cutadapt -g GATCGTCGGACT",
                      "-o",
                      outfile,
                      "-j ",
                      threads,
                      infile)
          #cat(cmd,"\n")#print the current command
          system(cmd)
          unlink(c(infile))
        } else
          stop("Wrong adapter. Exit.")
      }
    }
  }
}
