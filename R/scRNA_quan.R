tx2gene_scRNA <- function() {
    annotation <- list.files(pattern = "gff$", recursive = TRUE, full.names = TRUE)
    existence <- check_dep(dependancy = "egrep")
    if(existence == TRUE){
        cmd1 <- paste(
            "-v '^#|^$'", annotation, 
            "| cut -f 9 | grep ID=rna | awk -F ';' 'BEGIN{OFS = \"=\";} {print $1, $2;}' | awk -F '=' 'BEGIN{OFS = \",\"} {print $NF, $2}' > raw_tx2gene.csv"
        )
        system2(command = "egrep", args = cmd1)
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
        utils::write.table(tx2gene, sep = "\t", col.names = FALSE, row.names = FALSE, "tx2gene.tsv", quote = FALSE)
        unlink("raw_tx2gene.csv")        
    }
}


#' scRNA expression quantification with Salmon.
#' @param transcript the reference transcript
#' @param protocol the single-cell RNA sequencing protocol: dropseq, chromium, or chromiumV3.
#' @param threads the number of threads to be used. Default is 4.
#' @param salmon_index_add additional parameters to customize salmon index command. Default is NULL.
#' @param salmon_alevin_add additional parameters to customize salmon alevin command. Default is NULL.
#' @export scRNA_quan
#' @return None
#' @examples
#'
#' scRNA_quan(transcript = 'test_transcript_Homo_sapiens.fna', protocol = 'dropseq')
#'
scRNA_quan <- function(transcript, protocol, threads = 4, salmon_index_add = NULL, salmon_alevin_add = NULL) {
    existence <- check_dep(dependancy = "salmon")
    if (existence == TRUE) {
        tx2gene_scRNA()
        cmd3 <- paste(
            "index -t", transcript, 
            "-i salmon_index", "-p", threads, salmon_index_add
            )
        system2(command = "salmon", args = cmd3)
        read <- list.files(pattern = "^Trimmed.*1\\.fastq$", full.names = FALSE)
        if (length(read) == 0) {
            read <- list.files(pattern = ".*1\\.fastq$", full.names = FALSE)
        }
        for (f in read) {
            name <- gsub("_1.fastq", "", f)
            out <- paste0("scRNA_", name)
            read1 <- paste0(name, "_1.fastq")
            read2 <- paste0(name, "_2.fastq")
            read1 <- paste0(name, "_1.fastq")
            read2 <- paste0(name, "_2.fastq")
            read_seq <- read2
            barcode_seq <- read1
            cmd1 <- paste0("-n 1 ", read1, " | grep -o length=[0-9]* | cut -d '=' -f 2")
            leg_1 <- as.numeric(system2(command = "head", args = cmd1))
            cmd2 <- paste0("-n 1 ", read2, " | grep -o length=[0-9]* | cut -d '=' -f 2")
            leg_2 <- as.numeric(system2(command = "head", args = cmd2))
            if (leg_1 > leg_2) {
                read_seq <- read1
                barcode_seq <- read2
            }
            cmd4 <- paste(
                "alevin -p", threads, "-i salmon_index -l ISR --", 
                protocol, " -1 ", barcode_seq, " -2 ", read_seq, 
                " --tgMap tx2gene.tsv -o ", 
                out, salmon_alevin_add
                )
            system2(command = "salmon", args = cmd4)
        }
        unlink("tx2gene.tsv")
        unlink("salmon_index", recursive = TRUE)
    } else print("Salmon is not found. Please install it.")
}
