#' Extract the reference genome, transcript and the corresponding annotation files.
#'
#' @param taxa the scientific or common name of the organism.
#' @export extract_genome
#' @return the extracted reference genome, transcript and annotation file names.
#' @examples
#' \dontrun{
#'
#' extract_genome("sesame")
#' }

extract_genome <- function(taxa)
{
  ref_seq <-
    list.files(pattern = "fna$",
               recursive = TRUE,
               full.names = TRUE)
  transcript <-
    list.files(pattern = "rna\\.fna",
               recursive = TRUE,
               full.names = TRUE)
  # ref_seq1 <- ref_seq[-which(ref_seq == transcript)]
  taxa <- gsub("\\s", "_", taxa)
  name1 <- paste0(taxa, ".fna")
  # tmp <- paste0("./", name1)
  # ref_seq1 <- paste(ref_seq[-which(ref_seq == tmp)], collapse = ' ')
  ref_seq1 <- paste(ref_seq, collapse = ' ')
  cmd2 <- paste("cat", ref_seq1, ">", name1)
  system(cmd2)
  gff <-
    list.files(pattern = "gff$",
               recursive = TRUE,
               full.names = TRUE)
  name2 <- paste0(taxa, ".gff")
  cmd3 <- paste("mv", gff, name2)
  system(cmd3)

  name3 <- paste0("transcript_", taxa, ".fna")
  cmd4 <- paste("mv", transcript, name3)
  system(cmd4)
  # results <- c(name1, name2, name3)
  unlink("./ncbi_dataset.zip")
  unlink("README.md")
  unlink("./dehydrated", recursive = TRUE)
  # return(results)
}
