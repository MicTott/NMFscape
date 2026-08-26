makeTuneSCE <- function(ngenes = 60, ncells = 40, seed = 42) {
    set.seed(seed)
    sce <- scuttle::mockSCE(ngenes = ngenes, ncells = ncells)
    scuttle::logNormCounts(sce)
}

makeMultiSCE <- function() {
    sce <- makeTuneSCE()
    set.seed(7)
    atac <- SingleCellExperiment(list(
        counts = abs(matrix(rnorm(20 * ncol(sce)), 20, ncol(sce)))
    ))
    rownames(atac) <- paste0("peak_", seq_len(20))
    colnames(atac) <- colnames(sce)
    altExp(sce, "ATAC") <- atac
    sce
}

test_that("standard recipe cross-validates k and transfers to runNMFscape", {
    sce <- makeTuneSCE()
    fit <- tuneNMF(sce, params = list(k = c(2, 4)), reps = 2, verbose = FALSE)

    expect_s3_class(fit, "nmfscape_tuning")
    expect_true(fit$best_params$k %in% c(2, 4))
    expect_true(all(c("k", "mean_test_loss") %in% names(fit$summary)))

    out <- do.call(runNMFscape, c(list(x = sce, verbose = FALSE),
                                  fit$best_params))
    expect_equal(ncol(reducedDim(out, "NMF")), fit$best_params$k)
})

test_that("conditioned and multimodal recipes run and transfer", {
    sce <- makeTuneSCE()
    sce$batch <- rep(c("A", "B"), length.out = ncol(sce))
    cond <- tuneNMF(sce, params = list(k = c(2, 3)), recipe = "conditioned",
                    condition_col = "batch", reps = 1, verbose = FALSE)
    out <- do.call(runConditionedNMF,
                   c(list(x = sce, condition_col = "batch", verbose = FALSE),
                     cond$best_params))
    expect_equal(ncol(reducedDim(out, "CondNMF")), cond$best_params$k)

    mm <- makeMultiSCE()
    tuned <- tuneNMF(mm, params = list(k = c(2, 3)), recipe = "multimodal",
                     alt_exps = c("main", "ATAC"),
                     modality_names = c("RNA", "ATAC"), reps = 1,
                     verbose = FALSE)
    expect_true(tuned$best_params$k %in% c(2, 3))
})

test_that("layer_fn overrides the recipe topology", {
    sce <- makeTuneSCE()
    fit <- tuneNMF(sce, params = list(k = c(2, 4)), reps = 1, verbose = FALSE,
                   layer_fn = function(p, inputs) {
                       RcppML::nmf_layer(inputs[[1]], k = p$k, name = "custom")
                   })
    expect_true(fit$best_params$k %in% c(2, 4))
    expect_error(tuneNMF(sce, params = list(k = 2), layer_fn = "nope"),
                 "layer_fn must be a function")
})

test_that("plotTuning returns a ggplot and bad input errors", {
    sce <- makeTuneSCE()
    one <- tuneNMF(sce, params = list(k = c(2, 4)), reps = 1, verbose = FALSE)
    expect_s3_class(plotTuning(one), "ggplot")

    two <- tuneNMF(sce, params = list(k = c(2, 4), L1 = c(0, 0.01)),
                   reps = 1, verbose = FALSE)
    expect_s3_class(plotTuning(two), "ggplot")

    expect_error(tuneNMF(sce, params = list(), verbose = FALSE), "non-empty")
    expect_error(plotTuning(sce))
})
