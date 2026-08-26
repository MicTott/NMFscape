setup_diag_sce <- function(ngenes = 200, ncells = 100) {
    sce <- scuttle::mockSCE(ngenes = ngenes, ncells = ncells)
    scuttle::logNormCounts(sce)
}

test_that("diagnoseNMF returns all four diagnostics", {
    skip_if_not_installed("scuttle")
    sce <- setup_diag_sce()
    fit <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)

    res <- diagnoseNMF(fit)

    expect_type(res, "list")
    expect_equal(names(res), c("dispersion", "zero_inflation", "distribution",
                               "sparsity", "k"))
    expect_equal(res$k, 4)

    expect_true(res$dispersion$mode %in% c("global", "per_row", "per_col"))
    expect_type(res$dispersion$global_phi, "double")
    expect_true(all(c("row_cv", "col_cv") %in% names(res$dispersion)))

    expect_type(res$zero_inflation$has_zi, "logical")
    expect_true(res$zero_inflation$zi_mode %in% c("none", "row", "col"))
    expect_length(res$zero_inflation$row_excess, nrow(sce))
    expect_length(res$zero_inflation$col_excess, ncol(sce))

    expect_s3_class(res$distribution$scores, "data.frame")
    expect_equal(colnames(res$distribution$scores),
                 c("power", "T_stat", "abs_T", "distribution"))
    expect_true(res$distribution$best_power %in% res$distribution$scores$power)
    expect_type(res$distribution$best_distribution, "character")

    expect_s3_class(res$sparsity, "data.frame")
    expect_equal(colnames(res$sparsity), c("factor", "sparsity", "model"))
    expect_setequal(res$sparsity$model, c("w", "h"))
    expect_equal(nrow(res$sparsity), 8)
    expect_true(all(res$sparsity$sparsity >= 0 & res$sparsity$sparsity <= 1))
})

test_that("diagnoseNMF honors its tuning arguments", {
    skip_if_not_installed("scuttle")
    sce <- setup_diag_sce()
    fit <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)

    # a CV threshold above every observed CV forces the global mode
    lax <- diagnoseNMF(fit, cv_threshold = 1e6)
    expect_identical(lax$dispersion$mode, "global")
    strict <- diagnoseNMF(fit, cv_threshold = 0)
    expect_false(identical(strict$dispersion$mode, "global"))

    # an unreachable excess-zero threshold rules zero inflation out
    expect_false(diagnoseNMF(fit, zi_threshold = 1)$zero_inflation$has_zi)

    subset_powers <- diagnoseNMF(fit, powers = c(0, 1), test_nb = FALSE)
    expect_equal(subset_powers$distribution$scores$power, c(0, 1))
    expect_null(subset_powers$distribution$nb_diagnostic)
})

test_that("diagnoseNMF errors for FactorNet-backed results", {
    skip_if_not_installed("scuttle")
    sce <- setup_diag_sce(ngenes = 120, ncells = 60)

    deep <- runDeepNMF(sce, k = c(10, 4), verbose = FALSE)
    expect_error(diagnoseNMF(deep, name = "DeepNMF"), "FactorNet recipe")
    expect_error(diagnoseNMF(deep, name = "DeepNMF"), "total_loss")

    sce$batch <- rep(c("A", "B"), length.out = ncol(sce))
    cond <- runConditionedNMF(sce, k = 4, condition_col = "batch",
                              verbose = FALSE)
    expect_error(diagnoseNMF(cond, name = "CondNMF"), "FactorNet recipe")
})

test_that("diagnoseNMF validates its inputs", {
    skip_if_not_installed("scuttle")
    sce <- setup_diag_sce()
    fit <- runNMFscape(sce, k = 3, seed = 1, verbose = FALSE)

    expect_error(diagnoseNMF("not an sce"), "SingleCellExperiment")
    expect_error(diagnoseNMF(fit, assay = "nope"), "not found")
    expect_error(diagnoseNMF(fit, name = "missing"), "not found in metadata")

    # a fit on a gene subset cannot be diagnosed against the full assay
    hvgs <- rownames(sce)[seq_len(100)]
    subset_fit <- runNMFscape(sce, k = 3, subset_row = hvgs, seed = 1,
                              verbose = FALSE)
    expect_error(diagnoseNMF(subset_fit), "subset_row")
})

test_that("diagnoseNMF works on refineNMF results", {
    skip_if_not_installed("scuttle")
    sce <- setup_diag_sce(ngenes = 120, ncells = 60)
    fit <- runNMFscape(sce, k = 3, seed = 1, verbose = FALSE)
    fit$celltype <- rep(c("A", "B"), length.out = ncol(fit))
    refined <- refineNMF(fit, label_col = "celltype", verbose = FALSE)

    res <- diagnoseNMF(refined, name = "NMF_refined")
    expect_equal(res$k, 3)
    expect_s3_class(res$sparsity, "data.frame")
})

test_that("diagnoseNMF explains that consensusNMF stores no model", {
    skip_if_not_installed("scuttle")
    sce <- setup_diag_sce(ngenes = 120, ncells = 60)
    cons <- consensusNMF(sce, k = 3, n_runs = 5, seed = 1, verbose = FALSE)

    expect_error(diagnoseNMF(cons, name = "cNMF"), "not found in metadata")
})
