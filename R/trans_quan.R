


## extract the quantification data in gtf files produced by prepDE.py and save them to csv files.
extract <- function()
{
  files <-
    list.files(pattern = ".*\\.gtf$",
               recursive = TRUE,
               full.names = TRUE)
  outputs <- c()
  for (file in files)
  {
    output <-
      paste0("Tmp", gsub("\\.gtf$", "", gsub(".*/", "", file)), ".csv")
    # print(output)
    outputs <- c(outputs, output)
    cmd = paste(
      "awk -F \'\\t\' \'$3 == \"transcript\" {print $9}\'",
      file,
      "| awk -F \';\' \'{print $2, $5, $6}\' | awk \'BEGIN{OFS = \",\"; print \"transcript_id\", \"FPKM\", \"TPM\"} {print $2, $4, $6}\' >",
      output
    )
    # cat(cmd, "\n")
    system(cmd)
  }
  return(outputs)
}

utils::globalVariables(c("transcript_id", "count", "gene_id"))

## comple all csv files into one file.
#' @importFrom magrittr "%>%"

compile <- function(files)
{
  others <-
    data.frame(
      transcript_id = character(),
      FPKM = numeric(),
      TPM = numeric(),
      sample = character()
    )
  for (f in files)
  {
    other <- utils::read.csv(f)
    other$sample <- gsub("\\.csv$", "", gsub("^Tmp", "", f))
    others <- rbind(others, other)
  }
  others <- stats::na.exclude(others)
  others <- others[!grepl("gene-", others$transcript_id),]
  others$transcript_id <- sub("^.*?-", "", others$transcript_id)
  count <- utils::read.csv("transcript_count_matrix.csv")
  samples <- gsub("sorted_", "", names(count)[-1])
  names(count)[-1] <- samples
  count <- count[!grepl("gene-", count$transcript_id),]
  count$transcript_id <- sub("^.*?-", "", count$transcript_id)
  tmp <- count %>% tidyr::gather(sample, count,-transcript_id)
  all <-
    merge(tmp,
          others,
          by = c("sample", "transcript_id"),
          all = TRUE)
  # all$transcript_id <- gsub("gene-", "", all$transcript_id)
  utils::write.csv(all, "transcript_quantifications.csv")
}

#' Expression quantification at gene and transcript levels.
#'
#' Expression quantification at gene and transcript levels.
#' @return None
#' @export trans_quan
#' @examples
#' \dontrun{
#'
#' trans_quan()
#' }
#' @importFrom magrittr "%>%"

trans_quan <- function()
{
  # cmd <- paste("prepDE.py")
  # system(cmd)
  reticulate::source_python(system.file("prepDE_R.py", package = "BP4RNAseq"))
  outputs <- extract()
  compile(outputs)
  files <-
    list.files(pattern = "^Tmp.*csv$",
               recursive = TRUE,
               full.names = TRUE)
  unlink(files)
  files <-
    list.files(pattern = ".*bam$",
               recursive = TRUE,
               full.names = TRUE)
  unlink(files)
  files <-
    list.files(pattern = ".*ht2$",
               recursive = TRUE,
               full.names = TRUE)
  unlink(files)
  files <-
    list.files(pattern = ".*table$",
               recursive = TRUE,
               full.names = TRUE)
  unlink(files)
  gene_tmp <- utils::read.csv("gene_count_matrix.csv")
  gene_tmp <- gene_tmp %>% tidyr::gather(sample, count,-gene_id)
  gene_tmp$gene_id <- gsub("gene-", "", gene_tmp$gene_id)
  gene_tmp <- gene_tmp[, c("sample", "gene_id", "count")]
  utils::write.csv(gene_tmp, "gene_quantification.csv")
  unlink("gene_count_matrix.csv")
  unlink("transcript_count_matrix.csv")
  unlink("ballgown", recursive = TRUE)
}
