test_that("evaluateNMF returns loss value", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)
    sce <- runNMFscape(sce, k = 5, verbose = FALSE)

    loss <- evaluateNMF(sce)

    expect_true(is.numeric(loss))
    expect_length(loss, 1)
    expect_true(loss > 0)
})

test_that("evaluateNMF errors on missing model", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    expect_error(evaluateNMF(sce), "not found in metadata")
})
