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

test_that("enrichment improves agreement with the labels", {
    sce <- .testGuidedSCE()
    sce <- runGuidedNMF(sce, k = 5, label_col = "program", mode = "enrich",
                        strength = 1, name = "Enriched", init_name = "Base",
                        seed = 42, verbose = FALSE)

    ari_base <- .testARI(.argmaxLabels(sce, "Base"), sce$program)
    ari_enrich <- .testARI(.argmaxLabels(sce, "Enriched"), sce$program)
    expect_gt(ari_enrich, ari_base)
    expect_gt(ari_enrich, 0.95)
})

test_that("guided results use the standard slots", {
    sce <- .testGuidedSCE(n_genes = 60, n_cells = 80, noise = 2)
    sce <- runGuidedNMF(sce, k = 4, label_col = "program", seed = 1,
                        verbose = FALSE)

    expect_equal(dim(reducedDim(sce, "GuidedNMF")), c(80, 4))
    expect_equal(nrow(getBasis(sce, "GuidedNMF")), 60)
    expect_s4_class(getModel(sce, "GuidedNMF"), "nmf")
})

test_that("partially labelled cells are retained", {
    sce <- .testGuidedSCE(n_genes = 60, n_cells = 80, noise = 2)
    sce$partial <- sce$program
    sce$partial[seq(1, 80, by = 2)] <- NA
    sce <- runGuidedNMF(sce, k = 4, label_col = "partial", seed = 1,
                        verbose = FALSE)
    expect_equal(nrow(reducedDim(sce, "GuidedNMF")), 80)
})

test_that("mode = 'remove' is refused, and says what to use instead", {
    sce <- .testGuidedSCE(n_genes = 60, n_cells = 80, noise = 2)
    # RcppML's adversarial removal collapses the embedding onto one program
    # at every strength tested, so the mode is disabled rather than shipped.
    expect_error(runGuidedNMF(sce, k = 4, label_col = "program",
                              mode = "remove", verbose = FALSE),
                 "runConditionedNMF")
})
