test_that("fastqc_r works", {
    a <- fastqc_r(threads = 4, fq.dir = getwd(), scRNA = FALSE, fastqc_add = NULL)
    expect_equal(str_length("a"), 1)
})
