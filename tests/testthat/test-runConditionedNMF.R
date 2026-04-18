test_that("runConditionedNMF works with categorical condition", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)
    sce$batch <- sample(c("A", "B"), ncol(sce), replace = TRUE)

    sce <- runConditionedNMF(sce, k = 5, condition_col = "batch",
                             verbose = FALSE)

    expect_true("CondNMF" %in% reducedDimNames(sce))
    expect_equal(nrow(reducedDim(sce, "CondNMF")), ncol(sce))
    expect_equal(ncol(reducedDim(sce, "CondNMF")), 5)

    expect_true("CondNMF_basis" %in% names(metadata(sce)))
    expect_equal(nrow(metadata(sce)$CondNMF_basis), nrow(sce))

    expect_true("CondNMF_model" %in% names(metadata(sce)))
    expect_true("CondNMF_condition" %in% names(metadata(sce)))
    expect_equal(metadata(sce)$CondNMF_condition$col, "batch")
})

test_that("runConditionedNMF errors on missing condition_col", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    expect_error(runConditionedNMF(sce, k = 3, condition_col = "nonexistent"),
                 "not found in colData")
})

test_that("runConditionedNMF errors on NA in condition", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)
    sce$batch <- c(rep("A", 25), rep(NA, 25))

    expect_error(runConditionedNMF(sce, k = 3, condition_col = "batch"),
                 "contains NA")
})
