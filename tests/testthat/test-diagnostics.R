test_that("selectDistribution works", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    result <- selectDistribution(sce, k = 3, verbose = FALSE)

    expect_type(result, "list")
    expect_true("best" %in% names(result))
    expect_true("comparison" %in% names(result))
    expect_true(result$best %in% c("mse", "gp", "nb"))
    expect_s3_class(result$comparison, "data.frame")
    expect_true(all(c("distribution", "bic", "selected") %in%
                    colnames(result$comparison)))
})

test_that("selectRank works with single rep", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    cv <- selectRank(sce, k = 2:5, verbose = FALSE)

    expect_s3_class(cv, "data.frame")
    expect_equal(nrow(cv), 4)
    expect_true(all(c("k", "train_error", "test_error") %in% colnames(cv)))
    expect_equal(cv$k, 2:5)
    expect_true(all(cv$test_error > 0))
})

test_that("selectRank works with multiple reps", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    cv <- selectRank(sce, k = 2:4, n_reps = 2, seed = 42, verbose = FALSE)

    expect_s3_class(cv, "data.frame")
    expect_equal(nrow(cv), 3)
    expect_true(all(c("k", "train_error", "test_error",
                      "sd_train_error", "sd_test_error") %in% colnames(cv)))
})

test_that("plotRankSelection returns ggplot", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    cv <- selectRank(sce, k = 2:5, verbose = FALSE)
    p <- plotRankSelection(cv)

    expect_s3_class(p, "ggplot")
})
