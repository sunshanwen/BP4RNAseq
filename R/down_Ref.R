#' @order 2
### extract the genome and annotation files
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



#' Download the latest reference genome, transcript and the corresponding annotation files from NCBI.
#'
#' Download the latest reference genome and the corresponding annotation files from NCBI using the NCBI datasets command-line tool.
#' @param taxa the scientific or common name of the organism.
#' @export down_Ref
#' @order 1
#' @return None
#' @examples
#' down_Ref("Drosophila melanogaster")
#'
down_Ref <- function(taxa) {
  #cmd1 <- paste("./datasets assembly_descriptors tax_name", taxa, "-r | jq .datasets[].assembly_accession -r") ### change "_" to "-" according to the official documentation of datasets
  taxa_raw <- taxa
  taxa_tmp <- gsub("\\s", "_", taxa)
  genome <- paste0(taxa_tmp, ".fna")
  transcript <- paste0("transcript_", taxa_tmp, ".fna")
  annotation <- paste0(taxa_tmp, ".gff")
  status = 0

  if (!(file.exists(genome) &&
      file.exists(transcript) && file.exists(annotation))) {
    taxa <- paste0("\"", taxa, "\"")
    datasets <- list.files(pattern = "^datasets$", full.names = TRUE)
    if (length(datasets) == 0) {
      ### switch datasets according to the platform
      if (Sys.info()['sysname'] == "Linux")
      {
        status <-
        tryCatch(
          utils::download.file(
          "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets",
          destfile = "datasets",
          quit = TRUE
        ),
          error = function(err) {1},
          warning = function(war){2}
        )
        datasets <-
          list.files(pattern = "^datasets$", full.names = TRUE)
        # datasets <- system.file("datasets_L", package = "BP4RNAseq")
      } else if (Sys.info()['sysname'] == "Darwin") {
        utils::download.file(
          "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets",
          destfile = "datasets",
          quit = TRUE
        )
        status <-
        tryCatch(
          utils::download.file(
          "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets",
          destfile = "datasets",
          quit = TRUE
        ),
          error = function(err) {1},
          warning = function(war){2}
        )
        datasets <-
          list.files(pattern = "^datasets$", full.names = TRUE)
        # datasets <- system.file("datasets_D", package = "BP4RNAseq")
      }
    }
    if(status == 0){
      system(paste("chmod +x", datasets))

      if (dir.exists("dehydrated")&&(length(dir(path = "dehydrated", all.files = FALSE)) > 0)) {
        print("Downloading the reference genome and annotation files.")
        cmd3 <- paste(datasets, "rehydrate --filename dehydrated")
        status <- system(cmd3)
      } else {
        # cmd1 <- paste("./datasets assembly-descriptors tax-name", taxa, "--refseq --assmaccs | jq .datasets[].assembly_accession -r")
        cmd1 <-
          paste(
            datasets,
            "assembly-descriptors tax-name",
            taxa,
            "--refseq --assmaccs | jq .datasets[].assembly_accession -r"
          )
        # cat(cmd1, "\n")
        accession_id <- system(cmd1,  intern = TRUE)
        cmd2 <-
          paste(
            datasets,
            "download assembly",
            accession_id,
            "-g -r --dehydrated --filename dehydrated.zip"
          )
        system(cmd2, intern = TRUE)
        # cat(cmd2, "\n")
        # file <- list.files(pattern = "^dehydrated.zip$")
        utils::unzip("dehydrated.zip", list = FALSE, exdir = "dehydrated")
        print("Downloading the reference genome and annotation files.")
        cmd3 <- paste(datasets, "rehydrate --filename dehydrated")
        status <- system(cmd3)
      }

      if (status == 1)
      {
        print("The internet connection is poor. Please retry later!")
      } else {
        unlink(datasets)
        extract_genome(taxa_raw)
      }

      unlink("dehydrated.zip")

    } else {
      print ("The reference genome and annotation files already exist")
    }
    }
  return(status)
}
