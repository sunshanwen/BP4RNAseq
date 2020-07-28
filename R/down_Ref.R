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
  datasets <- list.files(pattern = "^datasets$", full.names = TRUE)
  if(length(datasets) == 0){
    ### switch datasets according to the platform
    if(Sys.info()['sysname'] == "Linux")
    {
      utils::download.file("https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets", destfile = "datasets", quit = TRUE)
      datasets <- list.files(pattern = "^datasets$", full.names = TRUE)
      # datasets <- system.file("datasets_L", package = "BP4RNAseq")
    } else if(Sys.info()['sysname'] == "Darwin"){
      utils::download.file("https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets", destfile = "datasets", quit = TRUE)
      datasets <- list.files(pattern = "^datasets$", full.names = TRUE)
      # datasets <- system.file("datasets_D", package = "BP4RNAseq")
    }
  }

  cmd0 <- system(paste("chmod +x", datasets))
  # cmd1 <- paste("./datasets assembly-descriptors tax-name", taxa, "--refseq --assmaccs | jq .datasets[].assembly_accession -r")
  cmd1 <- paste(datasets, "assembly-descriptors tax-name", taxa, "--refseq --assmaccs | jq .datasets[].assembly_accession -r")
  # cat(cmd1, "\n")
  accession_id <- system(cmd1,  intern = TRUE)

  cmd2 <- paste(datasets, "download assembly", accession_id, "-g -r --dehydrated --filename dehydrated.zip")
  system(cmd2)
  # cat(cmd2, "\n")
  # file <- list.files(pattern = "^dehydrated.zip$")
  utils::unzip("dehydrated.zip", list = FALSE, exdir = "dehydrated")

  cmd3 <- paste(datasets, "rehydrate --filename dehydrated")
  system(cmd3)
  unlink(datasets)
  unlink("dehydrated.zip")
}
