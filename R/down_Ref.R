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
  taxa <- taxa_raw
  taxa_tmp <- gsub("\\s", "_", taxa)
  genome <- paste0(taxa_tmp, ".fna")
  transcript <- paste0("transcript_", taxa_tmp, ".fna")
  annotation <- paste0(taxa_tmp, ".gff")
  if(file.exists(genome) && file.exists(transcript) && file.exists(annotation)) {
    stop("The reference genome and annotation files already exist")
  }

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

  system(paste("chmod +x", datasets))

  status = 0
  if(dir.exists("dehydrated")) {
    print("Downloading the reference genome and annotation files.")
    cmd3 <- paste(datasets, "rehydrate --filename dehydrated")
    status <- system(cmd3)
  } else {
    # cmd1 <- paste("./datasets assembly-descriptors tax-name", taxa, "--refseq --assmaccs | jq .datasets[].assembly_accession -r")
    cmd1 <- paste(datasets, "assembly-descriptors tax-name", taxa, "--refseq --assmaccs | jq .datasets[].assembly_accession -r")
    # cat(cmd1, "\n")
    accession_id <- system(cmd1,  intern = TRUE)

    cmd2 <- paste(datasets, "download assembly", accession_id, "-g -r --dehydrated --filename dehydrated.zip")
    system(cmd2, intern = T)
    # cat(cmd2, "\n")
    # file <- list.files(pattern = "^dehydrated.zip$")
    utils::unzip("dehydrated.zip", list = FALSE, exdir = "dehydrated")
    print("Downloading the reference genome and annotation files.")
    cmd3 <- paste(datasets, "rehydrate --filename dehydrated")
    status <- system(cmd3)
  }
  if(status == 1)
  {
    print("The internet connection is poor. Please retry later!")
  } else {unlink(datasets)}
  extract_genome(taxa_raw)
  unlink("dehydrated.zip")
  return(status)
}

dd <- down_Ref(taxa = "Homo sapiens")
