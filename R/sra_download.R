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
    ##check the version of sratoolkit
    ver <- system("vdb-config --version | cut -d '.' -f 3", intern = TRUE)
    ver <- na.exclude(as.numeric(as.character(ver)))
    if(ver <= 3){
        system("mkdir -p $HOME/.ncbi/")
        system("touch $HOME/.ncbi/user-settings.mkfg")
        cmd1 <- paste0("echo \"/repository/user/main/public/root = \\\"",dir,"\\\"\" > $HOME/.ncbi/user-settings.mkfg")
        system(cmd1)
        for(f in accession)
        {
            #### use Entrez Direct tool to get all samples accession code, then use prefetch to download all the data
            cmd2 = paste("prefetch -O", dir, "-X 100000000 $(esearch -db sra -query", f, "| efetch --format runinfo | grep",f, "| cut -d \",\" -f 1)")
            system(cmd2)
        }
    } else {
        cmd1 <- paste0("echo \"/repository/user/main/public/root = \\\"",dir,"\\\"\" > $HOME/.ncbi/user-settings.mkfg")
        # cat(cmd1)
        system(cmd1)
        # system("vdb-config --interactive", intern = T)
        # cmd3 <- "vdb-config -i & read -t 3; kill $!"
        # cat(cmd3)
        # system2(cmd3, intern = T)
        system(paste("/bin/bash -c", shQuote("vdb-config -i & read -t 3; kill $!")), ignore.stdout = TRUE, ignore.stderr = TRUE)

        for(f in accession)
        {
            #### use Entrez Direct tool to get all samples accession code, then use prefetch to download all the data
            cmd2 = paste("prefetch $(esearch -db sra -query", f, "| efetch --format runinfo | grep",f, "| cut -d \",\" -f 1)")
            system(cmd2)
            system("cd sra")
            file<-list.files(sra.dir, pattern = ".sra$", recursive = F, full.names = F)
            cmd4 <- paste("mv", file, dir)
            system(cmd4)
        }
    }


}


