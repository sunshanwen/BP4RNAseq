test_that("get_adapter works", {
    a <- get_adapter("test_read_2_fastqc", "pair")
    expect_output(str(a), "List of 2")
})
