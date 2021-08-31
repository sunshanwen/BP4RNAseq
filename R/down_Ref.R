#' @order 2
### extract the genome and annotation files
extract_genome <- function(taxa) {
    transcript <- list.files(path = "./dehydrated", pattern = "rna\\.fna", recursive = TRUE, full.names = TRUE)
    ref_seq <- list.files(path = "./dehydrated", pattern = "fna$", recursive = TRUE, full.names = TRUE)
    ref_seq <- ref_seq[ref_seq != transcript]
    taxa <- gsub("\\s", "_", taxa)
    name1 <- paste0(taxa, ".fna")
    ref_seq1 <- paste(ref_seq, collapse = " ")
    cmd2 <- paste(ref_seq1, ">", name1)
    system2(command = "cat", args = cmd2)
    gff <- list.files(path = "./dehydrated", pattern = "gff$", recursive = TRUE, full.names = TRUE)
    name2 <- paste0(taxa, ".gff")
    cmd3 <- paste(gff, name2)
    system2(command = "mv", args = cmd3)
    name3 <- paste0("transcript_", taxa, ".fna")
    cmd4 <- paste(transcript, name3)
    system2(command = "mv", args = cmd4)
    unlink("./ncbi_dataset.zip")
    unlink("README.md")
    unlink("./dehydrated", recursive = TRUE)
}



#' Download the latest reference genome, transcript and the corresponding annotation files from NCBI.
#'
#' Download the latest reference genome and the corresponding annotation files from NCBI using the NCBI datasets command-line tool.
#' @param taxa the scientific or common name of the organism.
#' @export down_Ref
#' @order 1
#' @return None
#' @examples
#' down_Ref(taxa = 'Drosophila melanogaster')
#'
down_Ref <- function(taxa) {

    taxa_raw <- taxa
    taxa_tmp <- gsub("\\s", "_", taxa)
    genome <- paste0(taxa_tmp, ".fna")
    transcript <- paste0("transcript_", taxa_tmp, ".fna")
    annotation <- paste0(taxa_tmp, ".gff")
    
    if (!(file.exists(genome) && file.exists(transcript) && file.exists(annotation))) {
        taxa <- paste0("'", taxa, "'")
        datasets <- list.files(pattern = "^datasets$", full.names = TRUE)
        if (length(datasets) == 0) {
            print('Downloading ncbi datasets.')
            ### switch datasets according to the platform
            if (Sys.info()["sysname"] == "Linux") {
              status <- tryCatch(utils::download.file("https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets", destfile = "datasets", 
                quit = TRUE), error = function(err) {
                1
              }, warning = function(war) {
                2
              })
              datasets <- list.files(pattern = "^datasets$", full.names = TRUE)
            } else if (Sys.info()["sysname"] == "Darwin") {
              status <- tryCatch(utils::download.file("https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets", destfile = "datasets", 
                quit = TRUE), error = function(err) {
                1
              }, warning = function(war) {
                2
              })
              datasets <- list.files(pattern = "^datasets$", full.names = TRUE)
            }
            if (status){
                print("The internet connection is poor. The download of datasets is failed. Please retry later.")
                unlink(datasets)
                return(status)
            }
        }
        cmd0 <- paste("+x", datasets)
        system2(command = "chmod", args = cmd0)
        i = 0
        while (!file.exists("dehydrated/ncbi_dataset/fetch.txt")&& i <3) {
            unlink('dehydrated', recursive = TRUE)
            cmd2 <- paste("download genome taxon", taxa, "--reference --dehydrated --filename dehydrated.zip")
            print(cmd2)
            print(datasets)
            system2(command = datasets, args = cmd2, stdout = TRUE)
            if (Sys.info()["sysname"] == "Linux") {
                system2(command = 'unzip', args = 'dehydrated.zip -d dehydrated', stdout = TRUE)
            } else{
                utils::unzip("dehydrated.zip", list = FALSE, overwrite = TRUE, unzip = "internal", exdir = "dehydrated")
                }
            i = i + 1
        } 
        if(!file.exists("dehydrated/ncbi_dataset/fetch.txt")){
            status = 1
            print("The internet connection is poor and the download of reference genome and annotation files failed. Please retry later!")
            return(status)
        }
        print("Downloading the reference genome and annotation files.")
        cmd3 <- paste("rehydrate --directory ./dehydrated")
        # status <- system2(command = datasets, args = cmd3, stdout = FALSE, stderr = FALSE)
        status <- system2(command = datasets, args = cmd3)
        
        # print(status)
        if (status) {
            print("The internet connection is poor and the download of reference genome and annotation files failed. Please retry later!")
        } else {
          unlink(datasets)
          extract_genome(taxa_raw)
        }
        unlink("dehydrated.zip")
        
    } else {
        status = 0
        print("The reference genome and annotation files already exist")
    }
    return(status)
}
