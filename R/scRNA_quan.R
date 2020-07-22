tx2gene_scRNA <- function()
{
  annotation <- list.files(pattern = "gff$", recursive = T, full.names = TRUE)
  #cmd1 <- paste("egrep -v '^#|^$'", annotation, "| awk -F '\t' '$3 ~ /RNA/ {print $9}' | awk -F ';' 'BEGIN{OFS = \"=\";} {print $1, $2;}' | awk -F '=' 'BEGIN{OFS = \"-\";}{print $NF, $2;}'| grep 'gene' | awk -F '-' 'BEGIN{OFS = \",\";print \"gene_id\", \"transcript_id\"}{print $2, $4}' > tx2gene.csv")
  cmd1 <- paste("egrep -v '^#|^$'", annotation, "| cut -f 9 | grep ID=rna | awk -F ';' 'BEGIN{OFS = \"=\";} {print $1, $2;}' | awk -F '=' 'BEGIN{OFS = \",\"} {print $NF, $2}' > raw_tx2gene.csv")
  #cat(cmd1)
  system(cmd1)
  tx2gene <- read.csv("raw_tx2gene.csv", header = FALSE)
  index_to_be_changed <- which(tx2gene[,1] %in% tx2gene[,2])
  b <- length(index_to_be_changed)
  for (i in 1:b)
  {
    tx2gene[index_to_be_changed[i],1] <- tx2gene[tx2gene[,2] == tx2gene[index_to_be_changed[i],1] ,1]
  }
  tx2gene[,1] <- gsub("gene-","", tx2gene[,1])
  tx2gene[,2] <- gsub("rna-","", tx2gene[,2])
  tx2gene[, c(2,1)]
  write.table(tx2gene, sep = "\t", col.names = FALSE,row.names=FALSE, "tx2gene.tsv", quote=FALSE)
}


#' scRNA expression quantification with Salmon.
#' @param transcript the reference transcript
#' @param protocol the single-cell RNA sequencing protocol: dropseq, chromium, or chromiumV3.
#' @export scRNA_quan
#' @return None
#' @examples
#' \dontrun{
#'
#' align_free_quan("paired", "sesame.fna", "transcript_sesame.fna","sesame.gff")
#'}
#'
scRNA_quan <- function(transcript, protocol)
{
  tx2gene_scRNA()
  cmd3 <- paste("salmon index --gencode -t", transcript, "-i salmon_index")
  # cat(cmd3, "\n")
  system(cmd3)
  read <- list.files(pattern = "^Trimmed.*1\\.fastq$", full.names = F)
  if(length(read) == 0){
    read <- list.files(pattern = ".*1\\.fastq$", full.names = F)
  }
  for(f in read)
  {
    name <- gsub("_1.fastq", "", f)
    out <- paste0("scRNA_",name)
    read1 <- paste0(name, "_1.fastq")
    read2 <- paste0(name, "_2.fastq")
    read1 <- paste0(name, "_1.fastq")
    read2 <- paste0(name, "_2.fastq")
    read_seq <- read2
    barcode_seq <- read1
    cmd1 <- paste0("head -n 1 ", read1, " | grep -o length=[0-9]* | cut -d '=' -f 2")
    leg_1 <- as.numeric(system(cmd1, intern = TRUE))
    cmd2 <- paste0("head -n 1 ", read2, " | grep -o length=[0-9]* | cut -d '=' -f 2")
    leg_2 <- as.numeric(system(cmd2, intern = TRUE))
    if(leg_1 > leg_2){
      read_seq <- read1
      barcode_seq <- read2
    }
    # cmd4 <- paste("salmon quant -i salmon_index -l A", gentrome.fna, "-1", read1, "-2", read2, "--validateMappings -o", out)
    cmd4 <- paste("salmon alevin -i salmon_index -l ISR", protocol, "-1", barcode_seq, "-2", read_seq, "--tgMap tx2gene.tsv -o", out)
    # cat(cmd4, "\n")
    system(cmd4, intern = T)
  }
  unlink("tx2gene.tsv")
  unlink("salmon_index", recursive = TRUE)
}
