#' Produce the quantifications of gene expression based on the fastq files.
#'
#' Produce the quantifications of gene expression based on the fastq files with alignment-based and alignment-free workflows.
#'
#'@param threads the number of threads to be used. Default is 4.
#'@param dir the working directory. Default is the current working directory.
#'@param pair 'single' for single-end (SE) reads or 'paired' for paired-end (PE) reads.
#'@param taxa the scientific or common name of the organism.
#'@param novel_transcript logic, whether identifying novel transcripts is expected or not. Default is FALSE.
#'@param scRNA logic, whether single-cell RNA-seq is quantified or not. Default is FALSE.
#'@param protocol the single-cell RNA sequencing protocol: dropseq, chromium, or chromiumV3.
#'@export fastq2quan
#'@return None
#'@examples
#'
#'sra_download(accession = 'SRR11427582')
#'sra2fastq()
#'fastq2quan(pair = 'single', taxa = 'Drosophila melanogaster')


fastq2quan <- function(threads = 4, dir = getwd(), pair, taxa, novel_transcript = FALSE, scRNA = FALSE, protocol = NULL) {
    setwd(dir)
    fastq.files <- list.files(pattern = "fastq$", full.names = FALSE)
    
    if (length(fastq.files)) {
        quality <- qc_test(threads, scRNA)
        
        if (length(quality$adapter$sample)) 
            cat("Adapter exists!\n")
        if (length(quality$base_quality$sample)) 
            cat("Per base sequence quality failed!\n")
        if (length(quality$sequence_score$sample)) 
            cat("Per sequence quality scores failed!\n")
        
        ### quality trimming
        if (length(quality$base_quality$sample) | length(quality$sequence_score$sample)) {
            print("Trim low quality bases.")
            quality_trim(quality$base_quality$sample, quality$sequence_score$sample, pair, scRNA)  ### consider if add directory parameter.
            if (length(quality$adapter$sample)) {
                print("Trim the adapter.")
                adapter_trim(quality$adapter$sample, pair, scRNA)
            }
        }
        
        if (length(quality$adapter$sample)) {
            print("Trim the adapter.")
            quality_trim(quality$adapter$sample, quality$adapter$sample, pair, scRNA)  ### consider if add directory parameter.
            adapter_trim(quality$adapter$sample, pair, scRNA)
        }
        files <- list.files(dir, pattern = "fastqc", full.names = FALSE)
        unlink(files)
        
        status <- down_Ref(taxa)
        if (status == 0) {
            # reference <- extract_genome(taxa)
            taxa_tmp <- gsub("\\s", "_", taxa)
            genome <- paste0(taxa_tmp, ".fna")
            transcript <- paste0("transcript_", taxa_tmp, ".fna")
            annotation <- paste0(taxa_tmp, ".gff")
            if (scRNA == FALSE) {
                align_based_quan(pair, taxa, genome, annotation, novel_transcript, threads)
                
                align_free_quan(pair, genome, transcript, annotation, threads)
                
                ### renaming results from trans_quan and align_free_quan first gene count
                align_based_gene <- utils::read.csv("gene_alignment_based_quantification.csv")
                align_based_gene <- align_based_gene[, c("sample", "gene_id", "count")]
                align_free_gene <- utils::read.csv("gene_alignment_free_quantification.csv")
                align_free_gene <- align_free_gene[, c("sample", "gene_id", "count")]
                utils::write.csv(align_based_gene, "gene_alignment_based_quantification.csv", row.names = FALSE)
                utils::write.csv(align_free_gene, "gene_alignment_free_quantification.csv", row.names = FALSE)
                ### then transcript count
                align_based_transcript <- utils::read.csv("transcript_alignment_based_quantification.csv")
                align_based_transcript <- align_based_transcript[, c("sample", "transcript_id", "count")]
                align_free_transcript <- utils::read.csv("transcript_alignment_free_quantification.csv")
                align_free_transcript <- align_free_transcript[, c("sample", "transcript_id", "count")]
                utils::write.csv(align_based_transcript, "transcript_alignment_based_quantification.csv", row.names = FALSE)
                utils::write.csv(align_free_transcript, "transcript_alignment_free_quantification.csv", row.names = FALSE)
                average(align_based_gene, align_free_gene, align_based_transcript, align_free_transcript)
            } else if (scRNA == TRUE) {
                scRNA_quan(transcript, protocol, threads)
            }
        }
    } else print("No fastq files are found in the work directory.")
}
