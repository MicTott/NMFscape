test_that("runDeepNMF works with 2 layers", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    sce <- runDeepNMF(sce, k = c(15, 5), verbose = FALSE)

    # Check reducedDim (cell embedding)
    expect_true("DeepNMF" %in% reducedDimNames(sce))
    expect_equal(nrow(reducedDim(sce, "DeepNMF")), ncol(sce))
    expect_equal(ncol(reducedDim(sce, "DeepNMF")), 5)

    # Check effective basis (genes x k_final)
    expect_true("DeepNMF_basis" %in% names(metadata(sce)))
    eff_basis <- getBasis(sce, "DeepNMF")
    expect_equal(nrow(eff_basis), nrow(sce))
    expect_equal(ncol(eff_basis), 5)

    # Check per-layer data
    layers <- metadata(sce)$DeepNMF_layers
    expect_length(layers, 2)
    expect_equal(dim(layers$layer_1$W), c(100, 15))
    expect_equal(dim(layers$layer_2$W), c(50, 5))

    # Check raw model
    expect_true("DeepNMF_model" %in% names(metadata(sce)))
})

test_that("runDeepNMF errors on single k", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    expect_error(runDeepNMF(sce, k = 5), "length >= 2")
})

test_that("runDeepNMF warns on non-decreasing k", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    expect_warning(
        runDeepNMF(sce, k = c(5, 10), verbose = FALSE),
        "not strictly decreasing"
    )
})

test_that("getTopFeatures works with deep NMF", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)
    sce <- runDeepNMF(sce, k = c(15, 5), verbose = FALSE)

    top <- getTopFeatures(sce, name = "DeepNMF", n = 5)
    expect_type(top, "list")
    expect_length(top, 5)
    expect_length(top[[1]], 5)
})

test_that("runDeepNMF handles three or more layers", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 150, ncells = 80)
    sce <- logNormCounts(sce)

    for (k in list(c(20, 10, 5), c(30, 20, 10, 5))) {
        res <- runDeepNMF(sce, k = k, seed = 1, verbose = FALSE)
        basis <- getBasis(res, "DeepNMF")
        coeff <- reducedDim(res, "DeepNMF")

        k_final <- k[length(k)]
        expect_equal(dim(basis), c(nrow(sce), k_final))
        expect_equal(dim(coeff), c(ncol(sce), k_final))
        expect_equal(rownames(basis), rownames(sce))
        expect_true(all(basis >= 0))
        expect_true(all(coeff >= 0))
        expect_false(any(is.na(basis)))
        expect_length(metadata(res)$DeepNMF_layers, length(k))
    }
})

test_that("deep NMF reconstruction tracks the input", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 150, ncells = 80)
    sce <- logNormCounts(sce)
    orig <- as.matrix(assay(sce, "logcounts"))

    for (k in list(c(20, 5), c(20, 10, 5))) {
        res <- runDeepNMF(sce, k = k, seed = 1, verbose = FALSE)
        recon <- getBasis(res, "DeepNMF") %*% t(reducedDim(res, "DeepNMF"))
        expect_equal(dim(recon), dim(orig))
        expect_gt(cor(as.vector(orig), as.vector(recon)), 0.3)
    }
})
