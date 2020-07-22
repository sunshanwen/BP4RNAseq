
convert_data <- function()
{
  files <- list.files(pattern = "quant.sf", recursive = T, full.names = TRUE)

  others <- data.frame(transcript_id = character(), length = numeric(), TPM = numeric(), count = numeric(), sample = character())
  for (f in files)
  {
    other <- utils::read.table(f, header = TRUE, sep = "\t")
    other$sample <- gsub("^\\./", "", gsub("_trans.*", "", f))
    others <- rbind(others, other[,-2])
  }
  names(others) <- c("transcript_id", "length", "TPM", "count", "sample")
  others <- others[, c("sample", "transcript_id", "count", "TPM", "length")]
  utils::write.csv(others, "salmon_transcript_quantifications.csv")
}

tx2gene <- function()
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
  colnames(tx2gene) <- c("transcript_id", "gene_id")
  write.csv(tx2gene, row.names=FALSE, "tx2gene.csv")
}

utils::globalVariables("TPM")

gene_quan <- function()
{
  tx2gene_file <- utils::read.csv("tx2gene.csv")
  transcript <- utils::read.csv("salmon_transcript_quantifications.csv")
  transcript <- transcript[, -1]
  all <- merge(tx2gene_file, transcript, by = "transcript_id", all = TRUE)
  all <- stats::na.omit(all)

  # tmp1 <- all %>% dplyr::group_by(sample, gene_id) %>% dplyr::summarise(abundance = sum(TPM), count = sum(count), length = sum(length*TPM)/sum(TPM))
  # # pre-calculate a simple average transcript length
  # # for the case the abundances are all zero for a gene.
  # mean_TPM <- all %>% dplyr::group_by(gene_id) %>% dplyr::summarise(length = mean(length))
  #
  # missing <- unique(tmp1$gene_id[is.nan(tmp1$length)])
  # for(i in missing)
  # {
  #   tmp1[tmp1$gene_id == i, ]$length <- mean_TPM[mean_TPM$gene_id == i, ]$length
  # }
  # gene_quantification <- tmp1[, c("sample", "gene_id", "count", "abundance", "length")]
  tmp1 <- all %>% dplyr::group_by(sample, gene_id) %>% dplyr::summarise(count = sum(count))
  gene_quantification <- tmp1[, c("sample", "gene_id", "count")]
  utils::write.csv(gene_quantification, "salmon_gene_quantification.csv")
}

#' Alignment-free expression quantification with Salmon at gene and transcript levels.
#'
#' Alignment-free expression quantification with Salmon at gene and transcript levels.
#' @param pair "single" for single-end (SE) or "paired" for paired-end (PE) reads.
#' @param genome the reference genome.
#' @param transcript the reference transcript
#' @param annotation the annotation file.
#' @export align_free_quan
#' @return None
#' @examples
#' \dontrun{
#'
#' align_free_quan("paired", "sesame.fna", "transcript_sesame.fna","sesame.gff")
#'}
#'
align_free_quan <- function(pair, genome, transcript, annotation)
{
  cmd3 <- paste("salmon index --gencode -t", transcript, "-i salmon_index")
  # cat(cmd3, "\n")
  system(cmd3)
  if(pair == "paired")
  {
    read <- list.files(pattern = "^Trimmed.*1\\.fastq$", full.names = F)
    if(length(read) == 0){
      read <- list.files(pattern = ".*1\\.fastq$", full.names = F)
    }
    # read2 <- list.files(pattern = "^Trimmed.*2\\.fastq$", full.names = F)
    for(f in read)
    {
      # read1 <- paste(read1, sep = ",", collapse = ',')
      # read2 <- paste(read2, sep = ",", collapse = ',')
      name <- gsub("_1.fastq", "", f)
      out <- paste0(name, "_transcripts_quant")
      read1 <- paste0(name, "_1.fastq")
      read2 <- paste0(name, "_2.fastq")
      # cmd4 <- paste("salmon quant -i salmon_index -l A", gentrome.fna, "-1", read1, "-2", read2, "--validateMappings -o", out)
      cmd4 <- paste("salmon quant -i salmon_index -l A", "-1", read1, "-2", read2, "--validateMappings -o", out)
      # cat(cmd4, "\n")
      system(cmd4, intern = T)
    }
  } else if(pair == "single"){
    read <- list.files(pattern = "^Trimmed.*\\.fastq$", full.names = F)
    if(length(read) == 0){
      read <- list.files(pattern = ".*\\.fastq$", full.names = F)
    }
    # read <- paste(read, sep = ",", collapse = ',')
    for(f in read)
    {
      name <- gsub(".fastq", "", f)
      out <- paste0(name, "_transcripts_quant")
      # cmd4 <- paste("salmon quant -i salmon_index -l A -1", gentrome.fna, "-r", f, "--validateMappings -o", out)
      cmd4 <- paste("salmon quant -i salmon_index -l A -r", f, "--validateMappings -o", out)

      # cat(cmd4, "\n")
      system(cmd4, intern = T)
    }

  } else stop("Paired-end and single-end mix. Please check the data source!")
  convert_data() #### rna quantification to use for quantifying genes.
  tx2gene()
  gene_quan()
  unlink("tx2gene.csv")
  folders <- dir(pattern = "transcripts_quant$")
  unlink(folders, recursive = TRUE)
  unlink("salmon_index", recursive = TRUE)
}

