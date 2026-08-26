test_that("character subset_row names the basis rows correctly", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)
    hvg <- rownames(sce)[seq_len(30)]

    sce <- runNMFscape(sce, k = 3, subset_row = hvg, verbose = FALSE)

    basis <- getBasis(sce)
    expect_equal(nrow(basis), 30)
    expect_equal(rownames(basis), hvg)
    expect_false(any(is.na(rownames(basis))))

    # top features must be real gene names, not NA
    top <- getTopFeatures(sce, n = 5)
    expect_false(any(is.na(unlist(top))))
    expect_true(all(unlist(top) %in% hvg))
})

test_that("integer and logical subset_row still work", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    int_res <- runNMFscape(sce, k = 3, subset_row = seq_len(30),
                           verbose = FALSE)
    expect_equal(rownames(getBasis(int_res)), rownames(sce)[seq_len(30)])

    keep <- rep(c(TRUE, FALSE), length.out = nrow(sce))
    log_res <- runNMFscape(sce, k = 3, subset_row = keep, verbose = FALSE)
    expect_equal(rownames(getBasis(log_res)), rownames(sce)[keep])
})

test_that("unknown feature names in subset_row error informatively", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    expect_error(
        runNMFscape(sce, k = 3, subset_row = c("not_a_gene"), verbose = FALSE),
        "not found in rownames"
    )
})

test_that("character subset_row works across the other wrappers", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)
    sce$batch <- rep(c("A", "B"), length.out = ncol(sce))
    hvg <- rownames(sce)[seq_len(30)]

    cons <- consensusNMF(sce, k = 3, subset_row = hvg, n_runs = 5,
                         verbose = FALSE)
    expect_equal(rownames(getBasis(cons, "cNMF")), hvg)

    cond <- runConditionedNMF(sce, k = 3, condition_col = "batch",
                              subset_row = hvg, verbose = FALSE)
    expect_equal(rownames(getBasis(cond, "CondNMF")), hvg)

    deep <- runDeepNMF(sce, k = c(10, 3), subset_row = hvg, verbose = FALSE)
    expect_equal(rownames(getBasis(deep, "DeepNMF")), hvg)
})
