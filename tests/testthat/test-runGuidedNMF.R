# Adjusted Rand Index, computed directly so the tests add no dependency.
.testARI <- function(a, b) {
    tab <- table(a, b)
    n <- sum(tab)
    sum_a <- sum(choose(rowSums(tab), 2))
    sum_b <- sum(choose(colSums(tab), 2))
    sum_ij <- sum(choose(tab, 2))
    expected <- sum_a * sum_b / choose(n, 2)
    (sum_ij - expected) / ((sum_a + sum_b) / 2 - expected)
}

# A noisy k-program simulation where the unguided fit is good but imperfect,
# leaving room for guidance to help or hurt measurably.
.testGuidedSCE <- function(seed = 1, n_genes = 200, n_cells = 400,
                           n_prog = 5, noise = 5, overlap = 0.5) {
    set.seed(seed)
    genes_per <- n_genes %/% n_prog
    group <- rep(seq_len(n_prog), length.out = n_cells)

    w_mat <- matrix(runif(n_genes * n_prog, 0, overlap), n_genes, n_prog)
    for (i in seq_len(n_prog)) {
        idx <- seq((i - 1) * genes_per + 1, i * genes_per)
        w_mat[idx, i] <- w_mat[idx, i] + runif(genes_per, 1, 3)
    }
    h_mat <- matrix(runif(n_prog * n_cells, 0, overlap), n_prog, n_cells)
    for (j in seq_len(n_cells)) {
        h_mat[group[j], j] <- h_mat[group[j], j] + runif(1, 1, 2)
    }
    mat <- w_mat %*% h_mat +
        matrix(abs(rnorm(n_genes * n_cells, 0, noise)), n_genes, n_cells)

    rownames(mat) <- sprintf("Gene_%03d", seq_len(n_genes))
    colnames(mat) <- sprintf("Cell_%03d", seq_len(n_cells))

    sce <- SingleCellExperiment::SingleCellExperiment(list(logcounts = mat))
    sce$program <- paste0("P", group)
    sce
}

.argmaxLabels <- function(sce, name) {
    apply(reducedDim(sce, name), 1, which.max)
}


test_that("enrichment improves agreement with labels", {
    sce <- .testGuidedSCE()

    sce <- runGuidedNMF(sce, k = 5, label_col = "program", mode = "enrich",
                        strength = 1, name = "Enriched",
                        init_name = "Base", seed = 42, verbose = FALSE)

    ari_base <- .testARI(.argmaxLabels(sce, "Base"), sce$program)
    ari_enrich <- .testARI(.argmaxLabels(sce, "Enriched"), sce$program)

    expect_gt(ari_enrich, ari_base)
    expect_gt(ari_enrich, 0.95)
})

test_that("runGuidedNMF stores results in the standard slots", {
    sce <- .testGuidedSCE(n_genes = 60, n_cells = 80, noise = 2)

    sce <- runGuidedNMF(sce, k = 4, label_col = "program", verbose = FALSE)

    expect_true("GuidedNMF" %in% reducedDimNames(sce))
    expect_equal(dim(reducedDim(sce, "GuidedNMF")), c(ncol(sce), 4L))
    expect_equal(rownames(reducedDim(sce, "GuidedNMF")), colnames(sce))
    expect_equal(colnames(reducedDim(sce, "GuidedNMF")),
                 paste0("NMF_", seq_len(4)))

    basis <- metadata(sce)$GuidedNMF_basis
    expect_equal(dim(basis), c(nrow(sce), 4L))
    expect_equal(rownames(basis), rownames(sce))

    expect_s4_class(metadata(sce)$GuidedNMF_model, "nmf")

    guidance <- metadata(sce)$GuidedNMF_guidance
    expect_equal(guidance$label_col, "program")
    expect_equal(guidance$mode, "enrich")
    expect_equal(guidance$target_lambda, 0.5)
    expect_true(guidance$whiten)
    expect_equal(guidance$n_guided, ncol(sce))
})

test_that("mode = 'remove' is refused, and says what to use instead", {
    sce <- .testGuidedSCE(n_genes = 60, n_cells = 80, noise = 2)

    # RcppML's adversarial target removal collapses the embedding onto a
    # single program at every strength from 0.5 to 10, leaving most programs
    # identically zero, so the mode is disabled rather than shipped.
    expect_error(
        runGuidedNMF(sce, k = 4, label_col = "program", mode = "remove",
                     strength = 0.75, verbose = FALSE),
        "disabled"
    )
    expect_error(
        runGuidedNMF(sce, k = 4, label_col = "program", mode = "remove",
                     verbose = FALSE),
        "runConditionedNMF"
    )
})

test_that("partially labelled (NA) cells are retained and left unguided", {
    sce <- .testGuidedSCE()

    partial <- sce$program
    set.seed(7)
    partial[sample(ncol(sce), round(ncol(sce) * 0.6))] <- NA
    sce$partial <- partial

    sce <- runGuidedNMF(sce, k = 5, label_col = "partial", mode = "enrich",
                        strength = 2, name = "Partial", init_name = "Base",
                        seed = 42, verbose = FALSE)

    # no cells dropped
    expect_equal(nrow(reducedDim(sce, "Partial")), ncol(sce))
    expect_equal(metadata(sce)$Partial_guidance$n_guided, sum(!is.na(partial)))

    ari_base <- .testARI(.argmaxLabels(sce, "Base"), sce$program)
    ari_partial <- .testARI(.argmaxLabels(sce, "Partial"), sce$program)

    # guidance from 40% of the cells still improves the whole embedding
    expect_gt(ari_partial, ari_base)
})

test_that("init_name stores the unguided base fit", {
    sce <- .testGuidedSCE(n_genes = 60, n_cells = 80, noise = 2)

    sce <- runGuidedNMF(sce, k = 4, label_col = "program",
                        init_name = "BaseFit", seed = 3, tol = 1e-5,
                        maxit = 100, verbose = FALSE)

    expect_true(all(c("GuidedNMF", "BaseFit") %in% reducedDimNames(sce)))
    expect_true("BaseFit_basis" %in% names(metadata(sce)))
    expect_s4_class(metadata(sce)$BaseFit_model, "nmf")
    expect_false("BaseFit_guidance" %in% names(metadata(sce)))

    # the base fit must match a plain runNMFscape() with the same seed
    ref <- runNMFscape(sce, k = 4, name = "Ref", seed = 3, verbose = FALSE)
    expect_equal(reducedDim(sce, "BaseFit"), reducedDim(ref, "Ref"))

    # and the guided fit must differ from it
    expect_false(isTRUE(all.equal(reducedDim(sce, "BaseFit"),
                                  reducedDim(sce, "GuidedNMF"),
                                  check.attributes = FALSE)))
})

test_that("standard accessors work on guided results", {
    sce <- .testGuidedSCE(n_genes = 60, n_cells = 80, noise = 2)
    sce <- runGuidedNMF(sce, k = 4, label_col = "program",
                        name = "Guided", verbose = FALSE)

    expect_equal(dim(getBasis(sce, "Guided")), c(nrow(sce), 4L))
    expect_equal(dim(getCoefficients(sce, "Guided")), c(ncol(sce), 4L))
    expect_s4_class(getModel(sce, "Guided"), "nmf")
    expect_length(getDiagonal(sce, "Guided"), 4L)

    top <- getTopFeatures(sce, name = "Guided", n = 5)
    expect_length(top, 4L)
    expect_true(all(lengths(top) == 5L))

    recon <- reconstructNMF(sce, name = "Guided")
    expect_equal(dim(recon), dim(sce))
})

test_that("runGuidedNMF respects subset_row", {
    sce <- .testGuidedSCE(n_genes = 60, n_cells = 80, noise = 2)
    keep <- rownames(sce)[seq_len(30)]

    sce <- runGuidedNMF(sce, k = 3, label_col = "program",
                        subset_row = keep, verbose = FALSE)

    expect_equal(rownames(metadata(sce)$GuidedNMF_basis), keep)
})

test_that("runGuidedNMF validates its arguments", {
    sce <- .testGuidedSCE(n_genes = 60, n_cells = 60, noise = 2)

    expect_error(runGuidedNMF(sce, k = 3, label_col = "nonexistent"),
                 "not found in colData")
    expect_error(runGuidedNMF(sce, k = 3, label_col = "program",
                              mode = "adversarial"),
                 "should be one of")
    expect_error(runGuidedNMF(sce, k = 3, label_col = "program",
                              strength = -0.5),
                 "single positive number")
    expect_error(runGuidedNMF(sce, k = 3, label_col = "program",
                              strength = 0),
                 "single positive number")
    expect_error(runGuidedNMF(sce, k = 3, label_col = "program",
                              strength = c(1, 2)),
                 "single positive number")
    expect_error(runGuidedNMF(sce, k = 3, label_col = "program",
                              assay = "missing"),
                 "not found in x")
    expect_error(runGuidedNMF("not an sce", k = 3, label_col = "program"),
                 "SingleCellExperiment")

    sce$numeric_label <- seq_len(ncol(sce))
    expect_error(runGuidedNMF(sce, k = 3, label_col = "numeric_label"),
                 "character or factor")

    sce$one_level <- "only"
    expect_error(runGuidedNMF(sce, k = 3, label_col = "one_level"),
                 "at least 2 non-NA levels")
})
