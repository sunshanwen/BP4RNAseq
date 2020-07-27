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
  datasets <- system.file("datasets", package = "BP4RNAseq")
  cmd0 <- system(paste("chmod +x", datasets))
  # cmd1 <- paste("./datasets assembly-descriptors tax-name", taxa, "--refseq --assmaccs | jq .datasets[].assembly_accession -r")
  cmd1 <- paste(datasets, "assembly-descriptors tax-name", taxa, "--refseq --assmaccs | jq .datasets[].assembly_accession -r")
  # cat(cmd1, "\n")
  accession_id <- system(cmd1,  intern = TRUE)

  continue <- TRUE
  while(continue){
    # print(accession_id)
    # cmd2 <- paste("./datasets download assembly", accession_id, "--include_gff3 --include_rna") ### change "_" to "-" according to the official documentation of datasets
    cmd2 <- paste(datasets, "download assembly", accession_id, "-g -r")
    # cat(cmd2, "\n")
    # cmd2 <- paste("./datasets download assembly", accession_id, "--unresolved")

    system(cmd2)
    file <- list.files(pattern = "^ncbi.*zip$")
    #file <- paste0("./", file)
    tmp <- try(utils::unzip(file, list = TRUE)$Name, silent = T)
    if(class(tmp) != "try-error"){
      files <- paste(tmp, collapse = " ")
      if(grepl("rna.fna", files) && grepl("genomic.gff", files)) {continue = FALSE}
      # files <- gsub("[0-9a-zA-Z.]*$","", unzip(file, list = TRUE)$Name)
      # print(files)
      # if("rna.fna" %in% files && "genomic.gff" %in% files) {continue = FALSE}
    }
  }
}
