#' Download the RNA-seq samples from the NCBI and produce the quantifications of gene expression.
#'
#'Download the RNA-seq samples from the NCBI and produce the quantifications of gene expression with alignment-based and alignment-free workflows.
#'
#'@param accession the bioproject accession code for the RNA-seq samples deposited in the NCBI.
#'@param dir the working directory. Default is the current working directory.
#'@param taxa the scientific or common name of the organism.
#'@param novel_transcript logic, whether identifying novel transcripts is expected or not. Default is FALSE.
#'@param threads the number of threads to be used. Default is 4.
#'@param scRNA logic, whether single-cell RNA-seq is quantified or not. Default is FALSE.
#'@param protocol the single-cell RNA sequencing protocol: dropseq, chromium, or chromiumV3.
#'@export down2quan
#'@return None
#'@examples
#'
#' down2quan(accession = "SRR11427582", dir = getwd(), taxa = "Drosophila melanogaster", novel_transcript = FALSE, scRNA = FALSE)
#'

down2quan <-
  function(accession,
           dir = getwd(),
           taxa,
           novel_transcript = FALSE,
           threads = 4,
           scRNA = FALSE,
           protocol) {
    setwd(dir)
    status <- sra_download(accession, dir)
    if(status == 0){
      pair <- sra2fastq(threads)
      quality <- qc_test(threads, scRNA)

      if (length(quality[[1]]$sample))
        cat("Adapter exists!\n")
      if (length(quality[[2]]$sample))
        cat("Per base sequence quality failed!\n")
      if (length(quality[[3]]$sample))
        cat("Per sequence quality scores failed!\n")

      ### quality trimming
      if (length(quality[[2]]$sample) | length(quality[[3]]$sample)) {
        print("Trim low quality bases.")
        quality_trim(quality[[2]]$sample, quality[[3]]$sample, pair, scRNA)### consider if add directory parameter.
        if (length(quality[[1]]$sample)) {
          print("Trim the adapter.")
          adapter_trim(quality[[1]]$sample, pair, scRNA)
        }
      }

      if (length(quality[[1]]$sample)) {
        print("Trim the adapter.")
        quality_trim(quality[[1]]$sample, quality[[1]]$sample, pair, scRNA)### consider if add directory parameter.
        adapter_trim(quality[[1]]$sample, pair, scRNA)
      }
      files <- list.files(dir, pattern = "fastqc", full.names = FALSE)
      unlink(files)

      status <- down_Ref(taxa)
      if (status > 0)
      {
        print(
          "The download of reference genome and annotation files failed! Please try it later."
        )
      } else if(status == 0){
        # reference <- extract_genome(taxa)
        taxa_tmp <- gsub("\\s", "_", taxa)
        genome <- paste0(taxa_tmp, ".fna")
        transcript <- paste0("transcript_", taxa_tmp, ".fna")
        annotation <- paste0(taxa_tmp, ".gff")
        if (scRNA == FALSE) {
          align_based_quan(pair, taxa, genome, annotation, novel_transcript)

          align_free_quan(pair, genome, transcript, annotation)

          ### renaming results from trans_quan and align_free_quan
          # first gene count
          align_based_gene <- utils::read.csv("gene_alignment_based_quantification.csv")
          align_based_gene <-
            align_based_gene[, c("sample", "gene_id", "count")]
          align_free_gene <-
            utils::read.csv("gene_alignment_free_quantification.csv")
          align_free_gene <-
            align_free_gene[, c("sample", "gene_id", "count")]
          utils::write.csv(align_based_gene,
                          "gene_alignment_based_quantification.csv",
                          row.names = FALSE)
          utils::write.csv(align_free_gene,
                          "gene_alignment_free_quantification.csv",
                          row.names = FALSE)
          # unlink("gene_quantification.csv")
          # unlink("salmon_gene_quantification.csv")

          ### then transcript count
          align_based_transcript <-
            utils::read.csv("transcript_alignment_based_quantification.csv")
          align_based_transcript <-
            align_based_transcript[, c("sample", "transcript_id", "count")]
          align_free_transcript <-
            utils::read.csv("transcript_alignment_free_quantification.csv")
          align_free_transcript <-
            align_free_transcript[, c("sample", "transcript_id", "count")]
          utils::write.csv(
            align_based_transcript,
            "transcript_alignment_based_quantification.csv",
            row.names = FALSE
          )
          utils::write.csv(
            align_free_transcript,
            "transcript_alignment_free_quantification.csv",
            row.names = FALSE
          )
          # unlink("transcript_quantifications.csv")
          # unlink("salmon_transcript_quantifications.csv")

          average(align_based_gene, align_free_gene, align_based_transcript, align_free_transcript)
        } else if (scRNA == TRUE) {
          scRNA_quan(transcript, protocol)
        }
    }
  }
}
