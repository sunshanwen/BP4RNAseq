#' Download RNA-seq samples from the NCBI.
#'
#' Download RNA-seq samples from the NCBI using SRA Toolkit and Entrez Direct based on the accession code.
#' @param accession the bioproject accession code for the RNA-seq samples deposited in the NCBI.
#' @export sra_download
#' @return None
#' @examples
#' \dontrun{
#'
#' sra_download("PRJNA615381")
#' }




sra_download <- function(accession)
{
    dir=getwd()
    for(f in accession)
    {
        #### use Entrez Direct tool to get all samples accession code, then use prefetch to download all the data
        # cmd1 = paste("prefetch -O", dir, "--max-size 100000000 $(esearch -db sra -query", f, "| efetch --format runinfo |cut -d \",\" -f 1 | grep SRR)")
        cmd1 = paste("prefetch $(esearch -db sra -query", f, "| efetch --format runinfo |cut -d \",\" -f 1 | grep SRR)")
        # cat(cmd1)
        system(cmd1)
    }
}


