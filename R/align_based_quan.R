
#### Reads alignment with Hisat2
.index_build <- function(taxa, genome, annotation, threads, hisat2_build_add) {
    if (file.exists(genome) && file.exists(annotation)) {
        cmd1 <- paste(
            "'{if ($3==\"exon\") {print $1\"\\t\"$4-1\"\\t\"$5-1}}'", 
                      annotation, "> exonsFile.table"
            )
        system2(command = "awk", args = cmd1)
        cmd2 <- paste(
            "'{if ($3==\"intron\") {print $1\"\\t\"$4-1\"\\t\"$5-1\"\\t\"$7}}'",
            annotation, "> ssFile.table"
            )
        system2(command = "awk", args = cmd2)
        if (file.size("exonsFile.table") * file.size("ssFile.table")) {
            cmd3 <- paste(
                "-f", genome, taxa, "--ss ssFile.table", 
                "--exon exonsFile.table", "-p", threads, hisat2_build_add
                )
            system2(command = "hisat2-build", args = cmd3)
        } else {
            cmd3 <- paste("-f", genome, taxa, "-p", threads, hisat2_build_add)
            system2(command = "hisat2-build", args = cmd3)
        }
    } else print("The reference genome and annotation are missing. Please download them first.")
}

# Aligning the RNA-seq data to the reference genome with HISAT2.


.align_ge <- function(pair, taxa, genome, annotation, threads, hisat2_build_add, hisat2_add) {
    
    taxa <- gsub("\\s", "_", taxa)
    
    .index_build(taxa, genome, annotation, threads, hisat2_build_add)
    index <- list.files(pattern = "ht2$", recursive = TRUE, full.names = TRUE)
    if (length(index)) {
        if (pair == "paired") {
            read <- list.files(pattern = "^Trimmed.*1\\.fastq$", full.names = FALSE)
            if (length(read) == 0) {
                read <- list.files(pattern = ".*1\\.fastq$", full.names = FALSE)
            }
            if (length(read)) {
                for (f in read) {
                  name <- gsub("_1.fastq", "", f)
                  out_bam <- paste0(name, ".bam")
                  read1 <- paste0(name, "_1.fastq")
                  read2 <- paste0(name, "_2.fastq")
                  cmd4 <- paste(
                      "-p", threads, hisat2_add, "--dta -x", taxa, 
                      "-1", read1, "-2", read2, 
                      "| samtools view -bh - | samtools sort - >", out_bam
                      )
                  system2(command = "hisat2", args = cmd4)
                }
            }
        } else if (pair == "single") {
            read <- list.files(pattern = "^Trimmed.*\\.fastq$", full.names = FALSE)
            if (length(read) == 0) {
                read <- list.files(pattern = ".*\\.fastq$", full.names = FALSE)
            }
            if (length(read)) {
                for (f in read) {
                  name <- gsub(".fastq", "", f)
                  out_bam <- paste0(name, ".bam")
                  cmd4 <- paste(
                      "-p", threads, hisat2_add, "--dta -x", taxa, "-U", f, 
                      "| samtools view -bh - | samtools sort - >", out_bam
                      )
                  system2(command = "hisat2", args = cmd4)
                }
            }
        }
    }
}


# Transcript assembly with StringTie.

.trans_ass <- function(novel_transcript = FALSE, threads, stringtie_add) {
    aligned_bam <- list.files(pattern = "*\\.bam$")
    gff <- list.files(pattern = "gff$", recursive = TRUE, full.names = TRUE)
    if (length(aligned_bam) && length(gff)) {
        for (f in aligned_bam) {
            taxa <- gsub("\\.bam", "", f)
            output <- paste0("ballgown/", taxa, "/", taxa, ".gtf")
            if (novel_transcript == TRUE) {
                cmd1 <- paste(
                    f, "-b ballgown -G", gff, "-o", output, "-p", 
                    threads, stringtie_add
                    )
                system2(command = "stringtie", args = cmd1)
            } else {
                cmd1 <- paste(
                    f, "-G", gff, "-e -b ballgown", "-o", output,
                    "-p", threads, stringtie_add
                    )
                system2(command = "stringtie", args = cmd1)
            }
        }
    }
}


#### Gene expression quantification based on the results from stringtie.

## extract the quantification data in gtf files produced by prepDE.py and save them to csv files.
.extract <- function() {
    files <- list.files(pattern = ".*\\.gtf$", recursive = TRUE, full.names = TRUE)
    outputs <- c()
    
    if (length(files)) {
        for (file in files) {
            output <- paste0("Tmp", gsub("\\.gtf$", "", gsub(".*/", "", file)), ".csv")
            outputs <- c(outputs, output)
            cmd = paste(
                "-F '\\t' '$3 == \"transcript\" {print $9}'", file, 
                "| awk -F ';' '{print $2, $5, $6}' | awk 'BEGIN{OFS = \",\"; print \"transcript_id\", \"FPKM\", \"TPM\"} {print $2, $4, $6}' >", 
                output
                )
            system2(command = "awk", args = cmd)
        }
    }
    return(outputs)
}

utils::globalVariables(c("transcript_id", "count", "gene_id"))

## comple all csv files into one file.
#' @importFrom magrittr '%>%'
#' @param files intermediate csv files

.compile <- function(files) {
    others <- data.frame(transcript_id = character(), FPKM = numeric(), TPM = numeric(), sample = character())
    for (f in files) {
        other <- utils::read.csv(f)
        other$sample <- gsub("\\.csv$", "", gsub("^Tmp", "", f))
        others <- rbind(others, other)
    }
    others <- stats::na.exclude(others)
    others <- others[!grepl("gene-", others$transcript_id), ]
    others$transcript_id <- sub("^.*?-", "", others$transcript_id)
    count <- utils::read.csv("transcript_count_matrix.csv")
    samples <- gsub("sorted_", "", names(count)[-1])
    names(count)[-1] <- samples
    count <- count[!grepl("gene-", count$transcript_id), ]
    count$transcript_id <- sub("^.*?-", "", count$transcript_id)
    tmp <- count %>% tidyr::gather(sample, count, -transcript_id)
    all <- merge(tmp, others, by = c("sample", "transcript_id"), all = TRUE)
    utils::write.csv(all, "transcript_alignment_based_quantification.csv")
}

# Expression quantification at gene and transcript levels.

#' @importFrom magrittr '%>%'
#'

.trans_quan <- function() {
    reticulate::source_python(system.file("prepDE_R.py", package = "BP4RNAseq"))
    outputs <- .extract()
    if (length(outputs)) {
        .compile(outputs)
        files <- list.files(pattern = "^Tmp.*csv$", recursive = TRUE, full.names = TRUE)
        unlink(files)
        files <- list.files(pattern = ".*bam$", recursive = TRUE, full.names = TRUE)
        unlink(files)
        files <- list.files(pattern = ".*ht2$", recursive = TRUE, full.names = TRUE)
        unlink(files)
        files <- list.files(pattern = ".*table$", recursive = TRUE, full.names = TRUE)
        unlink(files)
        gene_tmp <- utils::read.csv("gene_count_matrix.csv")
        gene_tmp <- gene_tmp %>% tidyr::gather(sample, count, -gene_id)
        gene_tmp$gene_id <- gsub("gene-", "", gene_tmp$gene_id)
        gene_tmp <- gene_tmp[, c("sample", "gene_id", "count")]
        utils::write.csv(gene_tmp, "gene_alignment_based_quantification.csv")
        unlink("gene_count_matrix.csv")
        unlink("transcript_count_matrix.csv")
        unlink("ballgown", recursive = TRUE)
    }
}

######## alignment-based workflow
#' alignment-based workflow
#'
#' @param pair 'single' for single-end (SE) or 'paired' for paired-end (PE) reads.
#' @param taxa the scientific or common name of the organism.
#' @param genome the reference genome.
#' @param annotation the annotation file.
#' @param novel_transcript logic, whether identifying novel transcripts is expected or not. Default is FALSE.
#' @param threads the number of threads to be used. Default is 4.
#' @param hisat2_build_add additional parameters to customize hisat2 build command. Default is NULL.
#' @param hisat2_add additional parameters to customize hisat2 command. Default is NULL.
#' @param stringtie_add additional parameters to customize stringtie command. Default is NULL.
#' @export align_based_quan
#' @return None
#' @examples
#' \dontrun{
#' test_Homo_sapiens.fna <- system.file('extdata', 'test_Homo_sapiens.fna', package = 'BP4RNAseq')
#' test_Homo_sapiens.gff <- system.file('extdata', 'test_Homo_sapiens.gff', package = 'BP4RNAseq')
#' align_based_quan(pair = 'paired', taxa = 'Homo sapiens', 
#' genome = test_Homo_sapiens.fna, 
#' annotation = test_Homo_sapiens.gff, novel_transcript = FALSE)
#'}
#'

align_based_quan <- function(pair, taxa, genome, annotation, novel_transcript = FALSE, threads = 4, hisat2_build_add = NULL, hisat2_add = NULL, stringtie_add = NULL) {
    existence <- check_dep(dependancy = "hisat2")
    if (existence == TRUE) {
        existence <- check_dep(dependancy = "stringtie")
        if (existence == TRUE) {
            .align_ge(pair, taxa, genome, annotation, threads, hisat2_build_add, hisat2_add)
            .trans_ass(novel_transcript, threads, stringtie_add)
            .trans_quan()
        }
    } else print("Hisat2 is not found. Please install it.")
}
