#' Take averages with alignment-free and alignment-based quantification
#' #'
#' @export average
#' @return the average of gene expression quantification with alignment-free and alignment-based quantification workflow
#' @examples
#' \dontrun{
#'
#' average()
#' }

average <- function()
{
  ### gene count
  align_based_gene <- utils::read.csv("gene_quantification.csv")
  align_based_gene <- align_based_gene[, c("sample", "gene_id", "count")]
  align_free_gene <- utils::read.csv("salmon_gene_quantification.csv")
  align_free_gene <- align_free_gene[, c("sample", "gene_id", "count")]
  utils::write.csv(align_based_gene, "gene_alignment_based_quantification.csv", row.names=FALSE)
  utils::write.csv(align_free_gene, "gene_alignment_free_quantification.csv", row.names=FALSE)
  unlink("gene_quantification.csv")
  unlink("salmon_gene_quantification.csv")
  all_gene <- rbind(align_based_gene, align_free_gene)
  ave_gene <- all_gene %>% dplyr::group_by(sample, gene_id) %>% dplyr::summarise(count = mean(count, na.rm = TRUE))
  utils::write.csv(ave_gene, "average_gene_quantification.csv", row.names=FALSE)

  ### transcript count
  align_based_transcript <- utils::read.csv("transcript_quantifications.csv")
  align_based_transcript <- align_based_transcript[, c("sample", "transcript_id", "count")]
  align_free_transcript <- utils::read.csv("salmon_transcript_quantifications.csv")
  align_free_transcript <- align_free_transcript[, c("sample", "transcript_id", "count")]
  utils::write.csv(align_based_transcript, "transcript_alignment_based_quantification.csv", row.names=FALSE)
  utils::write.csv(align_free_transcript, "transcript_alignment_free_quantification.csv", row.names=FALSE)
  unlink("transcript_quantifications.csv")
  unlink("salmon_transcript_quantifications.csv")
  all_transcript <- rbind(align_based_transcript, align_free_transcript)
  ave_transcript <- all_transcript %>% dplyr::group_by(sample, transcript_id) %>% dplyr::summarise(count = mean(count, na.rm = TRUE))
  utils::write.csv(ave_transcript, "average_transcript_quantification.csv", row.names=FALSE)
}
