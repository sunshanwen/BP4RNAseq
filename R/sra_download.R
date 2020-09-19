#' Download RNA-seq samples from the NCBI.
#'
#' Download RNA-seq samples from the NCBI using SRA Toolkit and Entrez Direct based on the accession code.
#' @param accession the bioproject accession code for the RNA-seq samples deposited in the NCBI.
#' @param dir the working directory. Default is the current working directory.
#' @export sra_download
#' @return the status of the download (0 for success, 1 for error).
#' @examples
#' sra_download(accession = 'SRR11427582')
#'
### Re-code the function to make examples runnable

sra_download <- function(accession, dir = getwd()) {
    if (Sys.info()["sysname"] == "Windows") {
        print("Please use Windows subsystem for linux.")
        status <- 1
    } else {
        ### check if sratools is available or not.
        existence <- check_dep(dependancy = "vdb-config")
        if (existence == TRUE) {status <- 0} else {status <- 1}
        if (status == 0) {
            system2(command = "vdb-config", args = "--interactive | pidof")
            ## Get the default download directory
            dir.download <- NULL
            
            root.dir <- system2(command = "vdb-config", args = "| grep '<root>' | cut -d '>' -f 2 | cut -d '<' -f 1", stdout = TRUE)
            root.exist <- NULL
            for (i in seq_len(length(root.dir))) {
                if (root.dir[i] != ".") {
                  dir.download = root.dir[i]
                  root.exist <- TRUE
                }
            }
            if (length(dir.download) == 0) {
                dir.download <- system2(command = "vdb-config", args = "| grep '<HOME>' | cut -d '>' -f 2 | cut -d '<' -f 1", stdout = TRUE)
                root.exist <- FALSE
            }
            
            ### check if Entrez Direct tool is avaliable
            existence <- check_dep(dependancy = "esearch")
            #### start downloading
            if (existence == TRUE) {
                file_path <- NULL
                status <- NULL
                if (root.exist == TRUE) {
                  for (f in accession) {
                    #### use Entrez Direct tool to get all samples accession code, then use prefetch to download all the data
                    cmd2 = paste("$(esearch -db sra -query", f, "| efetch --format runinfo | grep", f, "| cut -d \",\" -f 1)")
                    status <- tryCatch(system2("prefetch", args = cmd2, stdout = FALSE, stderr = FALSE), error = function(err) {
                      1
                    }, warning = function(war) {
                      2
                    })
                    if (status == 0) {
                      file_path <- paste0(dir.download, "/sra")
                      file <- list.files(file_path, pattern = ".sra$", recursive = FALSE, full.names = TRUE)
                      files <- paste(file, collapse = " ")
                      cmd4 <- paste(files, "-t", dir)
                      system2(command = "mv", args = cmd4)
                    }
                  }
                  unlink(file_path)
                } else {
                  for (f in accession) {
                    #### use Entrez Direct tool to get all samples accession code, then use prefetch to download all the data
                    cmd2 = paste("$(esearch -db sra -query", f, "| efetch --format runinfo | grep", f, "| cut -d \",\" -f 1)")
                    status <- tryCatch(system2("prefetch", args = cmd2, stdout = FALSE, stderr = FALSE), error = function(err) {
                      1
                    }, warning = function(war) {
                      2
                    })
                    if (status == 0) {
                      file_path <- paste0(dir.download, "/", f)
                      file <- list.files(file_path, pattern = ".sra$", recursive = FALSE, full.names = TRUE)
                      files <- paste(file, collapse = " ")
                      cmd4 <- paste(files, "-t", dir)
                      system2(command = "mv", args = cmd4)
                    }
                  }
                  unlink(file_path)
                }
                if (status != 0) 
                  print("Poor internet connection. The download of RNAseq data failed. Please try it later.")
            } else print("Entrez Direct is not found. Please install it.")
        } else print("SRA Toolkit is not found. Please install it.")
    }
    return(status)
}
