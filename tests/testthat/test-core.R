test_that("runNMFscape stores a usable factorization", {
    sce <- makeSCE()
    sce <- runNMFscape(sce, k = 5, seed = 1, verbose = FALSE)

    expect_equal(dim(reducedDim(sce, "NMF")), c(120L, 5L))
    expect_equal(rownames(getBasis(sce)), rownames(sce))
    expect_s4_class(getModel(sce), "nmf")
    expect_length(getDiagonal(sce), 5)
    expect_true(all(getBasis(sce) >= 0))
    expect_equal(dim(reconstructNMF(sce)), dim(sce))
    expect_length(getTopFeatures(sce, n = 3), 5)
})

test_that("absorbing the diagonal keeps A ~ W t(H)", {
    sce <- makeSCE()
    orig <- as.matrix(assay(sce, "logcounts"))
    fit <- runNMFscape(sce, k = 5, seed = 1, verbose = FALSE)

    recon <- getBasis(fit) %*% t(reducedDim(fit, "NMF"))
    expect_gt(cor(as.vector(orig), as.vector(recon)), 0.5)
})

test_that("consensusNMF reports stability and stores a model", {
    sce <- makeSCE()
    sce <- consensusNMF(sce, k = 3, n_runs = 5, seed = 1, verbose = FALSE)

    expect_equal(ncol(reducedDim(sce, "cNMF")), 3L)
    expect_s4_class(getModel(sce, "cNMF"), "nmf")
    expect_type(metadata(sce)$cNMF_consensus$cophenetic, "double")
})

test_that("predictNMF projects new data onto a reference", {
    ref <- runNMFscape(makeSCE(), k = 4, seed = 1, verbose = FALSE)
    query <- makeSCE(ncells = 40, seed = 2)
    rownames(query) <- rownames(ref)

    proj <- predictNMF(query, reference = ref, verbose = FALSE)
    expect_equal(dim(reducedDim(proj, "NMF_projected")), c(40L, 4L))
    # a raw model works too
    expect_no_error(predictNMF(query, reference = getModel(ref),
                               verbose = FALSE))
})

test_that("refineNMF shifts the embedding toward the labels", {
    fit <- runNMFscape(makeSCE(), k = 4, seed = 1, verbose = FALSE)
    ref <- refineNMF(fit, label_col = "celltype", verbose = FALSE)

    expect_equal(dim(reducedDim(ref, "NMF_refined")), c(120L, 4L))
    expect_false(identical(reducedDim(ref, "NMF_refined"),
                           reducedDim(fit, "NMF")))
})
