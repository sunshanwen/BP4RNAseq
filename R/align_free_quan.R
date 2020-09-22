convert_data <- function() {
    files <- list.files(pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
    if (length(files)) {
        others <- data.frame(transcript_id = character(), length = numeric(), TPM = numeric(), count = numeric(), sample = character())
        for (f in files) {
            other <- utils::read.table(f, header = TRUE, sep = "\t")
            other$sample <- gsub("^\\./", "", gsub("_trans.*", "", f))
            others <- rbind(others, other[, -2])
        }
        names(others) <- c("transcript_id", "length", "TPM", "count", "sample")
        others <- others[, c("sample", "transcript_id", "count", "TPM", "length")]
        utils::write.csv(others, "transcript_alignment_free_quantification.csv", row.names = FALSE)
    }
}

tx2gene <- function() {
    annotation <- list.files(pattern = "gff$", recursive = TRUE, full.names = TRUE)
    if (length(annotation)) {
        cmd1 <- paste("-v '^#|^$'", annotation, "| cut -f 9 | grep ID=rna | awk -F ';' 'BEGIN{OFS = \"=\";} {print $1, $2;}' | awk -F '=' 'BEGIN{OFS = \",\"} {print $NF, $2}' > raw_tx2gene.csv")
        system2(command = "egrep", args = cmd1)
        ### add support for Ensembl annotation file
        if (file.info("raw_tx2gene.csv")$size == 0) {
            cmd1 <- paste("-v '^#|^$'", annotation, "| cut -f 9 | grep ID=transcript | awk -F ';' 'BEGIN{OFS = \":\";} {print $1, $2;}' | awk -F ':' 'BEGIN{OFS = \",\"} {print $NF, $2}' > raw_tx2gene.csv")
            system2(command = "egrep", args = cmd1)
        }
        tx2gene <- utils::read.csv("raw_tx2gene.csv", header = FALSE)
        tx2gene[] <- lapply(tx2gene, as.character)
        index_to_be_changed <- which(tx2gene[, 1] %in% tx2gene[, 2])
        b <- length(index_to_be_changed)
        if (b > 0) {
            for (i in seq_len(b)) {
                tx2gene[index_to_be_changed[i], 1] <- tx2gene[tx2gene[, 2] == tx2gene[index_to_be_changed[i], 1], 1]
            }
        }
        tx2gene[, 1] <- gsub("gene-", "", tx2gene[, 1])
        tx2gene[, 2] <- gsub("rna-", "", tx2gene[, 2])
        tx2gene <- tx2gene[, c(2, 1)]
        colnames(tx2gene) <- c("transcript_id", "gene_id")
        utils::write.csv(tx2gene, row.names = FALSE, "tx2gene.csv")
    }
}

utils::globalVariables("TPM")

gene_quan <- function() {
    if (file.exists("tx2gene.csv") && file.exists("transcript_alignment_free_quantification.csv")) {
        tx2gene_file <- utils::read.csv("tx2gene.csv")
        transcript <- utils::read.csv("transcript_alignment_free_quantification.csv")
        all <- merge(tx2gene_file, transcript, by = "transcript_id", all = TRUE)
        all <- stats::na.omit(all)
        tmp1 <- all %>% dplyr::group_by(sample, gene_id) %>% dplyr::summarise(count = sum(count))
        gene_quantification <- tmp1[, c("sample", "gene_id", "count")]
        utils::write.csv(gene_quantification, "gene_alignment_free_quantification.csv", row.names = FALSE)
    }
}

#' Alignment-free expression quantification with Salmon at gene and transcript levels.
#'
#' Alignment-free expression quantification with Salmon at gene and transcript levels.
#' @param pair 'single' for single-end (SE) or 'paired' for paired-end (PE) reads.
#' @param genome the reference genome.
#' @param transcript the reference transcript
#' @param annotation the annotation file.
#' @param threads the number of threads to be used. Default is 4.
#' @param salmon_index_add additional parameters to customize salmon index command. Default is NULL.
#' @param salmon_quan_add additional parameters to customize salmon quan command. Default is NULL.
#' @export align_free_quan
#' @return None
#' @examples
#'
#' align_free_quan(
#' pair = 'paired', 'test_Homo_sapiens.fna', 
#' 'test_transcript_Homo_sapiens.fna','test_Homo_sapiens.gff'
#' )
#'
#'
#'
align_free_quan <- function(pair, genome, transcript, annotation, threads = 4, salmon_index_add = NULL, salmon_quan_add = NULL) {
    existence <- check_dep(dependancy = "salmon")
    if (existence == TRUE) {
        if (file.exists(genome) && file.exists(transcript) && file.exists(annotation)) {
            cmd3 <- paste("index -t", transcript, "-i salmon_index", "-p", threads, salmon_index_add)
            system2(command = "salmon", args = cmd3)
            if (pair == "paired") {
                read <- list.files(pattern = "^Trimmed.*1\\.fastq$", full.names = FALSE)
                if (length(read) == 0) {
                  read <- list.files(pattern = ".*1\\.fastq$", full.names = FALSE)
                }
                if (length(read)) {
                  for (f in read) {
                    name <- gsub("_1.fastq", "", f)
                    out <- paste0(name, "_transcripts_quant")
                    read1 <- paste0(name, "_1.fastq")
                    read2 <- paste0(name, "_2.fastq")
                    cmd4 <- paste("quant -p", threads, "-i salmon_index -l A", "-1", read1, "-2", read2, "--validateMappings -o", out, salmon_quan_add)
                    system2(command = "salmon", args = cmd4)
                  }
                } else print("No fastq files are found in the work directory.")
            } else if (pair == "single") {
                read <- list.files(pattern = "^Trimmed.*\\.fastq$", full.names = FALSE)
                if (length(read) == 0) {
                  read <- list.files(pattern = ".*\\.fastq$", full.names = FALSE)
                }
                if (length(read)) {
                  for (f in read) {
                    name <- gsub(".fastq", "", f)
                    out <- paste0(name, "_transcripts_quant")
                    cmd4 <- paste("quant -p", threads, "-i salmon_index -l A -r", f, "--validateMappings -o", out, salmon_quan_add)
                    system2(command = "salmon", args = cmd4)
                  }
                }
            } else stop("Paired-end and single-end mix. Please check the data source!")
            convert_data()  #### rna quantification to use for quantifying genes.
            tx2gene()
            gene_quan()
            unlink("raw_tx2gene.csv")
            unlink("tx2gene.csv")
            folders <- dir(pattern = "transcripts_quant$")
            unlink(folders, recursive = TRUE)
            unlink("salmon_index", recursive = TRUE)
        } else print("No reference genome and annotation files are found. Please download them first")
    } else print("Salmon is not found. Please install it.")
}
