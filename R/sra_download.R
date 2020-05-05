#' Download RNA-seq samples from the NCBI.
#'
#' Download RNA-seq samples from the NCBI using SRA Toolkit and Entrez Direct based on the accession code.
#' @param accession the bioproject accession code for the RNA-seq samples deposited in the NCBI.
#' @param dir the working directory. Default is the current working directory.
#' @export sra_download
#' @return None
#' @examples
#' \dontrun{
#'
#' sra_download("PRJNA615381")
#' }




sra_download <- function(accession, dir)
{
    dir=getwd()
    cmd1 <- paste0("echo \"/repository/user/main/public/root = \\\"",dir,"\\\"\" > $HOME/.ncbi/user-settings.mkfg")
    # cat(cmd1)
    system(cmd1)
    for(f in accession)
    {
        #### use Entrez Direct tool to get all samples accession code, then use prefetch to download all the data
        cmd2 = paste("prefetch -O", dir, "-X 100000000 $(esearch -db sra -query", f, "| efetch --format runinfo |cut -d \",\" -f 1 | grep SRR)")
        #cmd2 = paste("prefetch -X 100000000 $(esearch -db sra -query", f, "| efetch --format runinfo |cut -d \",\" -f 1 | grep SRR)")
        # cat(cmd1)
        system(cmd2)

    }
    files <- list.files(pattern = ".sra$", recursive = TRUE, full.names = T)
    direc <- unique(gsub("[a-zA-Z0-9]*\\.sra$", replacement = "",files))
    if(direc !="./"){
        file <- paste(files, collapse=" ")
        cmd <- paste("mv", file, dir)
        system(cmd)
        #print(direc)
        unlink(direc, recursive = TRUE)
    }
}


