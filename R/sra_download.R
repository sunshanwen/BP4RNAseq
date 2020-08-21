#' Download RNA-seq samples from the NCBI.
#'
#' Download RNA-seq samples from the NCBI using SRA Toolkit and Entrez Direct based on the accession code.
#' @param accession the bioproject accession code for the RNA-seq samples deposited in the NCBI.
#' @param dir the working directory. Default is the current working directory.
#' @export sra_download
#' @return the status of the download (0 for success, 1 for error).
#' @examples
#' sra_download(accession = "SRR11427582")
#'
### Re-code the function to make examples runnable

sra_download <- function(accession, dir = getwd())
{
  if (Sys.info()['sysname'] == "Windows") {
    print("Please use Windows subsystem for linux.")
    status <- 1
  } else{
    ### check if sratools is available or not.
    status <-
      tryCatch(
        system2(
          "which",
          args = shQuote("vdb-config"),
          stdout = FALSE,
          stderr = FALSE
        ),
        error = function(err) {
          1
        },
        warning = function(war) {
          2
        }
      )
    if (status == 0) {
      system(
        paste(
          "/bin/bash -c",
          shQuote("vdb-config --interactive & read -t 3; kill $!")
        ),
        ignore.stdout = TRUE,
        ignore.stderr = TRUE
      )

      ## Get the default download directory
      dir.download <- NULL

      root.dir <-
        system("vdb-config | grep '<root>' | cut -d '>' -f 2 | cut -d '<' -f 1",
               intern = TRUE)
      root.exist <- NULL
      for (i in seq_len(length(root.dir))) {
        if (root.dir[i] != ".") {
          dir.download = root.dir[i]
          root.exist <- TRUE
        }
      }
      if (length(dir.download) == 0) {
        dir.download <-
          system("vdb-config | grep '<HOME>' | cut -d '>' -f 2 | cut -d '<' -f 1",
                 intern = TRUE)
        root.exist <- FALSE
      }

      ### check if Entrez Direct tool is avaliable
      status <-
        tryCatch(
          system2(
            "which",
            args = shQuote("esearch"),
            stdout = FALSE,
            stderr = FALSE
          ),
          error = function(err) {
            1
          },
          warning = function(war) {
            2
          }
        )
      #### start downloading
      if (status == 0) {
        file_path <- NULL
        status <- NULL
        if (root.exist == TRUE) {
          for (f in accession)
          {
            #### use Entrez Direct tool to get all samples accession code, then use prefetch to download all the data
            cmd2 = paste(
              "$(esearch -db sra -query",
              f,
              "| efetch --format runinfo | grep",
              f,
              "| cut -d \",\" -f 1)"
            )
            status <-
              tryCatch(
                system2(
                  "prefetch",
                  args = cmd2,
                  stdout = FALSE,
                  stderr = FALSE
                ),
                error = function(err) {
                  1
                },
                warning = function(war) {
                  2
                }
              )
            if (status == 0) {
              file_path <- paste0(dir.download, "/sra")
              file <-
                list.files(
                  file_path,
                  pattern = ".sra$",
                  recursive = FALSE,
                  full.names = TRUE
                )
              # file<-list.files("./sra", pattern = ".sra$", recursive = F, full.names = T)
              files <- paste(file, collapse = " ")
              cmd4 <- paste("mv", files, "-t", dir)
              system(cmd4)
            }
          }
          unlink(file_path)
        } else{
          for (f in accession)
          {
            #### use Entrez Direct tool to get all samples accession code, then use prefetch to download all the data
            cmd2 = paste(
              "$(esearch -db sra -query",
              f,
              "| efetch --format runinfo | grep",
              f,
              "| cut -d \",\" -f 1)"
            )
            status <-
              tryCatch(
                system2(
                  "prefetch",
                  args = cmd2,
                  stdout = FALSE,
                  stderr = FALSE
                ),
                error = function(err) {
                  1
                },
                warning = function(war) {
                  2
                }
              )
            if (status == 0) {
              file_path <- paste0(dir.download, "/", f)
              # print(file_path)
              file <-
                list.files(
                  file_path,
                  pattern = ".sra$",
                  recursive = FALSE,
                  full.names = TRUE
                )
              #print(file)
              files <- paste(file, collapse = " ")
              cmd4 <- paste("mv", files, "-t", dir)
              system(cmd4)
            }
          }
          unlink(file_path)
        }
        if (status != 0)
          print(
            "Poor internet connection. The download of RNAseq data failed. Please try it later."
          )
      } else
        print("Entrez Direct is not found. Please install it.")
    } else
      print("SRA Toolkit is not found. Please install it.")
  }
  return(status)
}




#
# sra_download <- function(accession, dir = getwd())
# {
#
#     ##check the version of sratoolkit
#     ver <-
#         system("vdb-config --version | cut -d '.' -f 3", intern = TRUE)
#     ver <- stats::na.exclude(as.numeric(as.character(ver)))
#     if (ver <= 3) {
#         system("mkdir -p $HOME/.ncbi/")
#         system("touch $HOME/.ncbi/user-settings.mkfg")
#         cmd1 <-
#             paste0(
#                 "echo \"/repository/user/main/public/root = \\\"",
#                 dir,
#                 "\\\"\" > $HOME/.ncbi/user-settings.mkfg"
#             )
#         system(cmd1)
#         system(
#             paste(
#                 "/bin/bash -c",
#                 shQuote("vdb-config -i & read -t 3; kill $!")
#             ),
#             ignore.stdout = TRUE,
#             ignore.stderr = TRUE
#         )
#
#         for (f in accession)
#         {
#             system(
#                 paste(
#                     "/bin/bash -c",
#                     shQuote("vdb-config -i & read -t 3; kill $!")
#                 ),
#                 ignore.stdout = TRUE,
#                 ignore.stderr = TRUE
#             )
#             #### use Entrez Direct tool to get all samples accession code, then use prefetch to download all the data
#             cmd2 = paste(
#                 "prefetch -O",
#                 dir,
#                 "-X 100000000 $(esearch -db sra -query",
#                 f,
#                 "| efetch --format runinfo | grep",
#                 f,
#                 "| cut -d \",\" -f 1)"
#             )
#             system(cmd2)
#         }
#     } else {
#         cmd1 <-
#             paste0(
#                 "echo \"/repository/user/main/public/root = \\\"",
#                 dir,
#                 "\\\"\" > $HOME/.ncbi/user-settings.mkfg"
#             )
#         # cat(cmd1)
#         system(cmd1)
#         # system("vdb-config --interactive", intern = T)
#         # cmd3 <- "vdb-config -i & read -t 3; kill $!"
#         # cat(cmd3)
#         # system2(cmd3, intern = T)
#         system(
#             paste(
#                 "/bin/bash -c",
#                 shQuote("vdb-config -i & read -t 3; kill $!")
#             ),
#             ignore.stdout = TRUE,
#             ignore.stderr = TRUE
#         )
#
# for (f in accession)
# {
#     #### use Entrez Direct tool to get all samples accession code, then use prefetch to download all the data
#     cmd2 = paste(
#         "prefetch $(esearch -db sra -query",
#         f,
#         "| efetch --format runinfo | grep",
#         f,
#         "| cut -d \",\" -f 1)"
#     )
#     system(cmd2)
#     path <- paste0(dir, "/sra")
#     file <-
#         list.files(
#             path,
#             pattern = ".sra$",
#             recursive = FALSE,
#             full.names = TRUE
#         )
#     # file<-list.files("./sra", pattern = ".sra$", recursive = F, full.names = T)
#     cmd4 <- paste("mv", file, dir)
#     system(cmd4)
# }
#     }
# }
