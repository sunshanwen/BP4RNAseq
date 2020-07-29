#' Transcript assembly with StringTie.
#'
#' Transcript assembly with StringTie.
#' @param novel_transcript logic, whether identifying novel transcripts is expected or not. Default is FALSE.
#' @export trans_ass
#' @return None
#' @examples
#' \dontrun{
#'
#' trans_ass(novel_transcript = FALSE)
#' }


trans_ass <- function(novel_transcript = FALSE)
{
  aligned_bam <- list.files(pattern = "*\\.bam$")
  gff <-
    list.files(pattern = "gff$",
               recursive = TRUE,
               full.names = TRUE)
  for (f in aligned_bam)
  {
    # taxa <- gsub("^sorted", "", gsub("\\.bam", "", f))
    taxa <- gsub("\\.bam", "", f)
    output <- paste0("ballgown/", taxa, "/", taxa, ".gtf")
    if (novel_transcript == TRUE) {
      cmd1 <- paste("stringtie", f, "-b ballgown -G", gff, "-o", output)
      # cat(cmd1, "\n")
      system(cmd1)
      # cmd2 <- paste("stringtie", aligned_bam, "-G", taxa, "-eB -o", taxa)
      # # cat(cmd1, "\n")
      # system(cmd2)
      #### consider to add merge
    } else {
      cmd1 <-
        paste("stringtie", f, "-G", gff, "-e -b ballgown", "-o", output)
      # cat(cmd1, "\n")
      system(cmd1)
      #### consider to add merge
    }
  }
}
