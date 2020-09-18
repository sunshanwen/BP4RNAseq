#' Check the existence of dependancies.
#'
#' Trim samples with adapter using cutadapt based on the results from quality_test().
#' @param dependancy the dependancy to check.
#' @export check_dep
#' @return TRUE if the dependancy exists, otherwise FALSE. 
#' @examples
#'
#' check_dep(dependancy = 'cutadapt')
#'

check_dep <- function(dependancy) {
    status <- tryCatch(system2(command = "which", 
                               args = dependancy, 
                               stdout = FALSE, stderr = FALSE), 
                       error = function(err) {1},
                       warning = function(war) {2}
                       )
    if (status == 0) 
        return(TRUE) else return(FALSE)
}
