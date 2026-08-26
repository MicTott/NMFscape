test_that("consensusNMF works with SingleCellExperiment", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    sce <- consensusNMF(sce, k = 3, n_runs = 5, verbose = FALSE)

    # Check reducedDim stored
    expect_true("cNMF" %in% reducedDimNames(sce))
    expect_equal(nrow(reducedDim(sce, "cNMF")), ncol(sce))
    expect_equal(ncol(reducedDim(sce, "cNMF")), 3)

    # Check basis stored
    expect_true("cNMF_basis" %in% names(metadata(sce)))
    expect_equal(nrow(metadata(sce)$cNMF_basis), nrow(sce))
    expect_equal(ncol(metadata(sce)$cNMF_basis), 3)

    # Check consensus result stored
    cres <- metadata(sce)$cNMF_consensus
    expect_type(cres, "list")
    expect_true("cophenetic" %in% names(cres))
    expect_true("clusters" %in% names(cres))
    expect_true("consensus" %in% names(cres))
    expect_true(cres$cophenetic >= 0 && cres$cophenetic <= 1)

    # Consensus clusters are over features (rows), not cells
    expect_equal(length(cres$clusters), nrow(sce))
})

test_that("consensusNMF stores a usable S4 model", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 120, ncells = 60)
    sce <- logNormCounts(sce)
    sce <- consensusNMF(sce, k = 3, n_runs = 5, verbose = FALSE)

    expect_true("cNMF_model" %in% names(metadata(sce)))
    expect_s4_class(getModel(sce, "cNMF"), "nmf")

    # the accessors that require an S4 nmf must now work on consensus results
    expect_length(getDiagonal(sce, "cNMF"), 3)
    expect_type(evaluateNMF(sce, name = "cNMF"), "double")
    expect_equal(dim(reconstructNMF(sce, "cNMF")), c(nrow(sce), ncol(sce)))
})
