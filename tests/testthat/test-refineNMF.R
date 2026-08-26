test_that("refineNMF works with cell type labels", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)
    sce <- runNMFscape(sce, k = 5, verbose = FALSE)

    # Add labels
    sce$celltype <- sample(c("TypeA", "TypeB", "TypeC"), ncol(sce),
                           replace = TRUE)

    sce <- refineNMF(sce, label_col = "celltype", verbose = FALSE)

    expect_true("NMF_refined" %in% reducedDimNames(sce))
    expect_equal(nrow(reducedDim(sce, "NMF_refined")), ncol(sce))
    expect_equal(ncol(reducedDim(sce, "NMF_refined")), 5)

    expect_true("NMF_refined_basis" %in% names(metadata(sce)))
    expect_true("NMF_refined_model" %in% names(metadata(sce)))
    expect_s4_class(metadata(sce)$NMF_refined_model, "nmf")
})

test_that("refineNMF errors on missing label_col", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)
    sce <- runNMFscape(sce, k = 3, verbose = FALSE)

    expect_error(refineNMF(sce, label_col = "nonexistent"),
                 "not found in colData")
})

test_that("cycles controls whether the gene programs move", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 150, ncells = 120)
    sce <- logNormCounts(sce)
    sce$celltype <- rep(c("A", "B", "C"), length.out = ncol(sce))
    fit <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)
    basis0 <- getBasis(fit)

    # the default is an embedding-only operation: W comes back untouched
    kept <- refineNMF(fit, label_col = "celltype", verbose = FALSE)
    expect_identical(unname(getBasis(kept, "NMF_refined")), unname(basis0))
    expect_false(identical(reducedDim(kept, "NMF_refined"),
                           reducedDim(fit, "NMF")))

    # a positive cycles count refits the basis
    moved <- refineNMF(fit, label_col = "celltype", cycles = 5,
                       refined_name = "NMF_c5", verbose = FALSE)
    expect_gt(max(abs(getBasis(moved, "NMF_c5") - basis0)), 0)
})

test_that("refineNMF tolerates a Matrix-backed assay", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 120, ncells = 90)
    sce <- logNormCounts(sce)
    sce$celltype <- rep(c("A", "B"), length.out = ncol(sce))

    dense_m <- Matrix::Matrix(as.matrix(assay(sce, "logcounts")))
    expect_s4_class(dense_m, "dgeMatrix")
    assay(sce, "logcounts") <- dense_m

    fit <- runNMFscape(sce, k = 3, seed = 1, verbose = FALSE)
    expect_no_error(refineNMF(fit, label_col = "celltype", cycles = 5,
                              verbose = FALSE))
})
