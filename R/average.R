#' Take averages with alignment-free and alignment-based quantification
#' @param based_gene gene count matrix assessed based on alignment-based workflow
#' @param free_gene gene count matrix assessed based on alignment-free workflow
#' @param based_transcript transcript count matrix assessed based on alignment-based workflow
#' @param free_transcript transcript count matrix assessed based on alignment-free workflow
#' @importFrom magrittr "%>%"
#' @export average
#' @return the average of gene expression quantification with alignment-free and alignment-based quantification workflow
#' @examples
#'
#' align_based_gene <- system.file("extdata", "test_gene_alignment_based_quantification.csv", package = "BP4RNAseq")
#' align_free_gene <- system.file("extdata", "test_gene_alignment_free_quantification.csv", package = "BP4RNAseq")
#' align_based_transcript <- system.file("extdata", "test_transcript_alignment_based_quantification.csv", package = "BP4RNAseq")
#' align_free_transcript <- system.file("extdata", "test_transcript_alignment_free_quantification.csv", package = "BP4RNAseq")
#
#' align_based_gene <- read.csv(align_based_gene, header = TRUE)
#' align_free_gene <- read.csv(align_free_gene, header = TRUE)
#' align_based_transcript <- read.csv(align_based_transcript, header = TRUE)
#' align_free_transcript <- read.csv(align_free_transcript, header = TRUE)

#' average(align_based_gene, align_free_gene, align_based_transcript, align_free_transcript)

average <-
  function(based_gene,
           free_gene,
           based_transcript,
           free_transcript)
  {
    ### gene count
    align_based_gene <- as.data.frame(based_gene)
    align_free_gene <- as.data.frame(free_gene)
    align_based_transcript <- as.data.frame(based_transcript)
    align_free_transcript <- as.data.frame(free_transcript)

    all_gene <- rbind(align_based_gene, align_free_gene)
    all_gene$sample <- as.factor(all_gene$sample)
    all_gene$gene_id <- as.factor(all_gene$gene_id)
    ave_gene <-
      all_gene %>% dplyr::group_by(sample, gene_id) %>% dplyr::summarise(count = mean(count, na.rm = TRUE))
    utils::write.csv(ave_gene, "average_gene_quantification.csv", row.names =
                       FALSE)

    ### transcript count

    all_transcript <-
      rbind(align_based_transcript, align_free_transcript)
    all_transcript$sample <- as.factor(all_transcript$sample)
    all_transcript$transcript_id <-
      as.factor(all_transcript$transcript_id)
    ave_transcript <-
      all_transcript %>% dplyr::group_by(sample, transcript_id) %>% dplyr::summarise(count = mean(count, na.rm = TRUE))
    utils::write.csv(ave_transcript,
                     "average_transcript_quantification.csv",
                     row.names = FALSE)
  }
