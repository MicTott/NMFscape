test_that("runNMFscape works with SingleCellExperiment", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    # Test basic functionality
    sce <- runNMFscape(sce, k = 5, verbose = FALSE)

    # Check that results are stored correctly
    expect_true("NMF" %in% reducedDimNames(sce))
    expect_true("NMF_basis" %in% names(metadata(sce)))
    expect_true("NMF_model" %in% names(metadata(sce)))

    # Check dimensions
    expect_equal(nrow(reducedDim(sce, "NMF")), ncol(sce))
    expect_equal(ncol(reducedDim(sce, "NMF")), 5)
    expect_equal(nrow(getBasis(sce)), nrow(sce))
    expect_equal(ncol(getBasis(sce)), 5)

    # Check model is S4 nmf class
    model <- getModel(sce)
    expect_s4_class(model, "nmf")
    expect_equal(length(model@d), 5)

    # Test accessor functions
    coeffs <- getCoefficients(sce)
    basis <- getBasis(sce)
    top_features <- getTopFeatures(sce, n = 10)

    expect_true(is.matrix(coeffs))
    expect_true(is.matrix(basis))
    expect_type(top_features, "list")
    expect_length(top_features, 5)

    # Test getDiagonal
    d <- getDiagonal(sce)
    expect_true(is.numeric(d))
    expect_length(d, 5)

    # Test reconstruction
    reconstructed <- reconstructNMF(sce)
    expect_true(is.matrix(reconstructed))
    expect_equal(dim(reconstructed), c(nrow(sce), ncol(sce)))
})

test_that("runNMFscape d absorption works correctly", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    # Run with d absorbed (default) and not absorbed
    sce_absorbed <- runNMFscape(sce, k = 3, seed = 42, verbose = FALSE)
    sce_raw <- runNMFscape(sce, k = 3, seed = 42, absorb_d = FALSE,
                           name = "NMF_raw", verbose = FALSE)

    # Both should reconstruct to approximately the same matrix
    recon_absorbed <- reconstructNMF(sce_absorbed)
    recon_raw <- reconstructNMF(sce_raw, name = "NMF_raw")
    expect_equal(recon_absorbed, recon_raw, tolerance = 1e-6)
})

test_that("runNMFscape distribution parameter works", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    # Test with explicit distribution
    sce <- runNMFscape(sce, k = 3, distribution = "mse", verbose = FALSE)
    expect_true("NMF" %in% reducedDimNames(sce))
})
