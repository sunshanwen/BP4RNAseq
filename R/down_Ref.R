#' Download the latest reference genome, transcript and the corresponding annotation files from NCBI.
#'
#' Download the latest reference genome and the corresponding annotation files from NCBI using the NCBI datasets command-line tool.
#' @param taxa the scientific or common name of the organism.
#' @export down_Ref
#' @return None
#' @examples
#' \dontrun{
#'
#' down_Ref("sesame")
#'}
down_Ref <- function(taxa) {
  #cmd1 <- paste("./datasets assembly_descriptors tax_name", taxa, "-r | jq .datasets[].assembly_accession -r") ### change "_" to "-" according to the official documentation of datasets
  taxa <- paste0("\"",taxa,"\"")
  cmd1 <- paste("./datasets assembly-descriptors tax-name", taxa, "--refseq --assmaccs | jq .datasets[].assembly_accession -r") ### may change later
  # cat(cmd1, "\n")
  accession_id <- system(cmd1,  intern = TRUE)
  # print(accession_id)
  #cmd2 <- paste("./datasets download assembly", accession_id, "--include_gff3 --include_rna") ### change "_" to "-" according to the official documentation of datasets
  cmd2 <- paste("./datasets download assembly", accession_id, "--include-gff3 --include-rna")
  # cat(cmd2, "\n")
  system(cmd2)
}
