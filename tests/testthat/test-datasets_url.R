test_that("ncbi datasets url works", {
    status <- 0
    if (Sys.info()["sysname"] == "Linux") {
        status <- tryCatch(utils::download.file("https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets", destfile = "datasets", 
                                                quit = TRUE), error = function(err) {
                                                    1
                                                }, warning = function(war) {
                                                    2
                                                })
        datasets <- list.files(pattern = "^datasets$", full.names = TRUE)
    } else if (Sys.info()["sysname"] == "Darwin") {
        status <- tryCatch(utils::download.file("https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets", destfile = "datasets", 
                                                quit = TRUE), error = function(err) {
                                                    1
                                                }, warning = function(war) {
                                                    2
                                                })
        datasets <- list.files(pattern = "^datasets$", full.names = TRUE)
    } else if (Sys.info()["sysname"] == "Windows") print("Please use Windows Subsystem for Linux.")
    expect_equal(status, 0)
})
