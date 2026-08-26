setup_diag_sce <- function(ngenes = 200, ncells = 100) {
    sce <- scuttle::mockSCE(ngenes = ngenes, ncells = ncells)
    scuttle::logNormCounts(sce)
}

test_that("diagnoseNMF returns the four diagnostics", {
    skip_if_not_installed("scuttle")
    sce <- setup_diag_sce()
    fit <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)
    res <- diagnoseNMF(fit)

    expect_equal(res$k, 4)
    expect_s3_class(res$sparsity, "data.frame")
    expect_true(all(c("dispersion", "zero_inflation", "distribution") %in%
                    names(res)))
})

test_that("diagnoseNMF explains itself on FactorNet results", {
    skip_if_not_installed("scuttle")
    sce <- setup_diag_sce()
    sce$batch <- rep(c("A", "B"), length.out = ncol(sce))
    cond <- runConditionedNMF(sce, k = 4, condition_col = "batch",
                              verbose = FALSE)
    expect_error(diagnoseNMF(cond, name = "CondNMF"), "FactorNet recipe")
})
