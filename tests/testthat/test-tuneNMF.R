library(testthat)
library(NMFscape)
library(SingleCellExperiment)

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

test_that("standard recipe cross-validates k and returns a usable result", {
    sce <- makeTuneSCE()
    tuned <- tuneNMF(sce, params = list(k = c(2, 4, 6)), reps = 2,
                     verbose = FALSE)

    expect_s3_class(tuned, "nmfscape_tuning")
    expect_s3_class(tuned, "factor_net_cv")
    expect_true(is.data.frame(tuned$summary))
    expect_true(all(c("k", "mean_test_loss", "se_test_loss") %in%
                        colnames(tuned$summary)))
    expect_equal(nrow(tuned$summary), 3L)
    expect_equal(nrow(tuned$results), 6L)
    expect_true(is.list(tuned$best_params))
    expect_true(tuned$best_params$k %in% c(2, 4, 6))
    expect_identical(tuned$n_fits, 6L)
    # summary is ranked best-first
    expect_false(is.unsorted(tuned$summary$mean_test_loss))
    # held-out losses are real numbers, not the all-zero upstream failure mode
    expect_true(all(tuned$summary$mean_test_loss > 0))
})

test_that("standard best_params feed straight into runNMFscape", {
    sce <- makeTuneSCE()
    tuned <- tuneNMF(sce, params = list(k = c(3, 5)), reps = 1,
                     verbose = FALSE)
    out <- do.call(runNMFscape,
                   c(list(x = sce, verbose = FALSE), tuned$best_params))
    expect_equal(ncol(reducedDim(out, "NMF")), tuned$best_params$k)
})

test_that("L1 is tuned alongside k and widened for runNMFscape", {
    sce <- makeTuneSCE()
    tuned <- tuneNMF(sce, params = list(k = c(3, 5), L1 = c(0, 0.05)),
                     reps = 1, verbose = FALSE)

    expect_equal(nrow(tuned$summary), 4L)
    expect_true(all(c("k", "L1") %in% colnames(tuned$summary)))
    # runNMFscape() wants a length-2 c(w, h) penalty
    expect_length(tuned$best_params$L1, 2L)
    out <- do.call(runNMFscape,
                   c(list(x = sce, verbose = FALSE), tuned$best_params))
    expect_true("NMF" %in% reducedDimNames(out))
})

test_that("conditioned recipe cross-validates and transfers to runConditionedNMF", {
    sce <- makeTuneSCE()
    sce$batch <- rep(c("A", "B"), length.out = ncol(sce))

    tuned <- tuneNMF(sce, params = list(k = c(2, 4)), recipe = "conditioned",
                     condition_col = "batch", reps = 2, verbose = FALSE)

    expect_equal(tuned$recipe, "conditioned")
    expect_equal(nrow(tuned$summary), 2L)
    expect_true(all(tuned$summary$mean_test_loss > 0))

    out <- do.call(runConditionedNMF,
                   c(list(x = sce, condition_col = "batch", verbose = FALSE),
                     tuned$best_params))
    expect_equal(ncol(reducedDim(out, "CondNMF")), tuned$best_params$k)
})

test_that("conditioned recipe requires condition_col", {
    sce <- makeTuneSCE()
    expect_error(
        tuneNMF(sce, params = list(k = c(2, 4)), recipe = "conditioned",
                reps = 1, verbose = FALSE),
        "condition_col is required"
    )
    expect_error(
        tuneNMF(sce, params = list(k = 2), recipe = "conditioned",
                condition_col = "nope", reps = 1, verbose = FALSE),
        "not found in colData"
    )
})

test_that("multimodal recipe cross-validates altExps and transfers", {
    sce <- makeMultiSCE()

    tuned <- tuneNMF(sce, params = list(k = c(2, 4)), recipe = "multimodal",
                     alt_exps = c("main", "ATAC"),
                     modality_names = c("RNA", "ATAC"),
                     reps = 2, verbose = FALSE)

    expect_equal(tuned$recipe, "multimodal")
    expect_equal(nrow(tuned$summary), 2L)
    expect_true(all(tuned$summary$mean_test_loss > 0))

    out <- do.call(runMultiModalNMF,
                   c(list(x = sce, alt_exps = c("main", "ATAC"),
                          modality_names = c("RNA", "ATAC"), verbose = FALSE),
                     tuned$best_params))
    expect_equal(ncol(reducedDim(out, "MultiNMF")), tuned$best_params$k)
})

test_that("multimodal recipe accepts several assays and rejects one modality", {
    sce <- makeTuneSCE()
    tuned <- tuneNMF(sce, params = list(k = c(2, 3)), recipe = "multimodal",
                     assay = c("counts", "logcounts"), reps = 1,
                     verbose = FALSE)
    expect_equal(nrow(tuned$summary), 2L)

    expect_error(
        tuneNMF(sce, params = list(k = 2), recipe = "multimodal",
                reps = 1, verbose = FALSE),
        "two or more assay names"
    )
    expect_error(
        tuneNMF(sce, params = list(k = 2), recipe = "multimodal",
                alt_exps = "main", reps = 1, verbose = FALSE),
        "at least 2 modalities"
    )
})

test_that("deep recipe tunes layer by layer and transfers to runDeepNMF", {
    sce <- makeTuneSCE()

    tuned <- tuneNMF(sce, recipe = "deep",
                     params = list(k = list(c(6, 10), c(2, 3))),
                     reps = 2, verbose = FALSE)

    expect_equal(tuned$recipe, "deep")
    expect_true("layer" %in% colnames(tuned$summary))
    expect_setequal(unique(tuned$summary$layer), c(1, 2))
    expect_true(all(tuned$summary$mean_test_loss > 0))

    expect_length(tuned$best_params$k, 2L)
    expect_true(is.integer(tuned$best_params$k))
    expect_true(tuned$best_params$k[1] %in% c(6, 10))
    expect_true(tuned$best_params$k[2] %in% c(2, 3))

    out <- do.call(runDeepNMF,
                   c(list(x = sce, verbose = FALSE), tuned$best_params))
    expect_equal(ncol(reducedDim(out, "DeepNMF")), tuned$best_params$k[2])
})

test_that("deep recipe holds non-k parameters fixed after layer 1", {
    sce <- makeTuneSCE()
    tuned <- tuneNMF(sce, recipe = "deep",
                     params = list(k = list(c(6, 8), c(2, 3)),
                                   L2 = c(0, 0.05)),
                     reps = 1, verbose = FALSE)

    expect_length(tuned$best_params$L2, 1L)
    expect_true(tuned$best_params$L2 %in% c(0, 0.05))
    # layer 2 only searched k, at the layer-1 penalty
    layer2 <- tuned$summary[tuned$summary$layer == 2, , drop = FALSE]
    expect_equal(nrow(layer2), 2L)
    expect_equal(unique(layer2$L2), tuned$best_params$L2)
})

test_that("deep recipe rejects grids that are not per-layer lists", {
    sce <- makeTuneSCE()
    expect_error(
        tuneNMF(sce, recipe = "deep", params = list(k = c(10, 3)),
                reps = 1, verbose = FALSE),
        "list with one candidate grid per layer"
    )
    expect_error(
        tuneNMF(sce, recipe = "deep", params = list(k = list(c(10, 5))),
                reps = 1, verbose = FALSE),
        "at least 2 layers"
    )
    expect_error(
        tuneNMF(sce, recipe = "deep", params = list(L1 = 0),
                reps = 1, verbose = FALSE),
        "must include 'k'"
    )
})

test_that("layer_fn overrides the recipe topology", {
    sce <- makeTuneSCE()

    custom <- function(p, inputs) {
        RcppML::nmf_layer(inputs[[1]], k = p$k, L1 = p$L1, name = "custom")
    }
    tuned <- tuneNMF(sce, params = list(k = c(2, 4), L1 = c(0, 0.1)),
                     layer_fn = custom, reps = 1, verbose = FALSE)

    expect_equal(nrow(tuned$summary), 4L)
    expect_true(all(tuned$summary$mean_test_loss > 0))
    # escape hatch leaves penalties as the scalars nmf_layer() wants
    expect_length(tuned$best_params$L1, 1L)

    expect_error(tuneNMF(sce, params = list(k = 2), layer_fn = "nope"),
                 "layer_fn must be a function")
})

test_that("layer_fn may be supplied for a multi-input graph", {
    sce <- makeMultiSCE()
    custom <- function(p, inputs) {
        RcppML::nmf_layer(do.call(RcppML::factor_shared, inputs), k = p$k,
                          name = "shared_nmf")
    }
    tuned <- tuneNMF(sce, params = list(k = c(2, 3)), recipe = "multimodal",
                     alt_exps = c("main", "ATAC"), layer_fn = custom,
                     reps = 1, verbose = FALSE)
    expect_equal(nrow(tuned$summary), 2L)
})

test_that("invalid inputs error clearly", {
    sce <- makeTuneSCE()

    expect_error(tuneNMF(matrix(1:4, 2), params = list(k = 2)),
                 "SingleCellExperiment")
    expect_error(tuneNMF(sce, params = list(k = 2), recipe = "bogus"),
                 "should be one of")
    expect_error(tuneNMF(sce, params = list()), "non-empty named list")
    expect_error(tuneNMF(sce, params = list(2, 3)), "unique, non-empty names")
    expect_error(tuneNMF(sce, params = list(k = numeric(0))),
                 "non-empty atomic vector")
    expect_error(tuneNMF(sce, params = list(k = list(2, 3))),
                 "non-empty atomic vector")
    expect_error(tuneNMF(sce, params = list(k = 2), reps = 0),
                 "positive integer")
    expect_error(tuneNMF(sce, params = list(k = 2), test_fraction = 0),
                 "strictly between 0 and 1")
    expect_error(tuneNMF(sce, params = list(k = 2), assay = "missing"),
                 "not found in x")
})

test_that("random strategy caps the number of combinations", {
    sce <- makeTuneSCE()
    tuned <- tuneNMF(sce, params = list(k = c(2, 3, 4, 5), L1 = c(0, 0.1)),
                     strategy = "random", n_random = 3, reps = 1,
                     verbose = FALSE)
    expect_equal(nrow(tuned$summary), 3L)
    expect_identical(tuned$n_fits, 3L)
    expect_equal(tuned$strategy, "random")
})

test_that("verbose reports the fit budget and warns when it is large", {
    sce <- makeTuneSCE()
    expect_message(
        tuneNMF(sce, params = list(k = c(2, 3)), reps = 1, verbose = TRUE),
        "2 model fits"
    )
    expect_warning(
        tuneNMF(sce, params = list(k = seq(2, 12)), reps = 10,
                strategy = "grid", verbose = FALSE),
        "model fits requested"
    )
})

test_that("plotTuning returns a ggplot for one varying parameter", {
    sce <- makeTuneSCE()
    tuned <- tuneNMF(sce, params = list(k = c(2, 4, 6)), reps = 2,
                     verbose = FALSE)
    p <- plotTuning(tuned)
    expect_s3_class(p, "ggplot")
    expect_true(is.data.frame(p$data))
    expect_true(all(c("xval", "mean_test_loss") %in% colnames(p$data)))
    expect_true(nrow(ggplot2::ggplot_build(p)$data[[1]]) > 0)

    expect_s3_class(plotTuning(tuned, show_se = FALSE,
                               highlight_best = FALSE), "ggplot")
})

test_that("plotTuning returns a heatmap for two varying parameters", {
    sce <- makeTuneSCE()
    tuned <- tuneNMF(sce, params = list(k = c(2, 4), L1 = c(0, 0.05)),
                     reps = 1, verbose = FALSE)
    p <- plotTuning(tuned)
    expect_s3_class(p, "ggplot")
    expect_true(all(c("xval", "yval") %in% colnames(p$data)))
    expect_equal(nrow(p$data), 4L)
    expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 4L)
})

test_that("plotTuning falls back to a ranked dot plot for three parameters", {
    sce <- makeTuneSCE()
    tuned <- tuneNMF(sce, params = list(k = c(2, 4), L1 = c(0, 0.05),
                                        L2 = c(0, 0.05)),
                     reps = 1, verbose = FALSE)
    p <- plotTuning(tuned, top_n = 5)
    expect_s3_class(p, "ggplot")
    expect_true("combo_label" %in% colnames(p$data))
    expect_equal(nrow(p$data), 5L)
    expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 5L)
})

test_that("plotTuning facets deep results by layer", {
    sce <- makeTuneSCE()
    tuned <- tuneNMF(sce, recipe = "deep",
                     params = list(k = list(c(6, 10), c(2, 3))),
                     reps = 1, verbose = FALSE)
    p <- plotTuning(tuned)
    expect_s3_class(p, "ggplot")
    expect_true("layer" %in% colnames(p$data))
    expect_equal(nlevels(p$data$layer), 2L)
    expect_length(unique(ggplot2::ggplot_build(p)$data[[1]]$PANEL), 2L)
})

test_that("plotTuning rejects non-tuning input and constant grids", {
    sce <- makeTuneSCE()
    expect_error(plotTuning(data.frame(k = 1:3)), "nmfscape_tuning object")

    tuned <- tuneNMF(sce, params = list(k = 4), reps = 1, verbose = FALSE)
    expect_error(plotTuning(tuned), "No parameter varies")
})
