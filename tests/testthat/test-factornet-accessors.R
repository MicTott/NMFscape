# FactorNet recipes store a factor_net_result rather than an S4 nmf; the
# accessors used to fail with "no applicable method for @".

make_sce <- function() {
    sce <- scuttle::mockSCE(ngenes = 120, ncells = 60)
    sce <- scuttle::logNormCounts(sce)
    sce$batch <- rep(c("A", "B"), length.out = ncol(sce))
    sce
}

test_that("getDiagonal works for FactorNet results", {
    skip_if_not_installed("scuttle")
    sce <- make_sce()

    cond <- runConditionedNMF(sce, k = 4, condition_col = "batch",
                              verbose = FALSE)
    d <- getDiagonal(cond, "CondNMF")
    expect_length(d, 4)
    expect_true(all(d >= 0))

    mm <- runMultiModalNMF(sce, k = 3, assays = c("counts", "logcounts"),
                           modality_names = c("A", "B"), verbose = FALSE)
    expect_length(getDiagonal(mm, "MultiNMF"), 3)
})

test_that("reconstructNMF works for FactorNet results", {
    skip_if_not_installed("scuttle")
    sce <- make_sce()

    cond <- runConditionedNMF(sce, k = 4, condition_col = "batch",
                              verbose = FALSE)
    recon <- reconstructNMF(cond, "CondNMF")
    expect_equal(dim(recon), c(nrow(sce), ncol(sce)))
    expect_false(any(is.na(recon)))
})

test_that("predictNMF projects single-layer FactorNet results", {
    skip_if_not_installed("scuttle")
    sce <- make_sce()
    query <- scuttle::logNormCounts(scuttle::mockSCE(ngenes = 120, ncells = 30))
    rownames(query) <- rownames(sce)

    cond <- runConditionedNMF(sce, k = 4, condition_col = "batch",
                              verbose = FALSE)
    proj <- predictNMF(query, reference = cond, ref_name = "CondNMF",
                       verbose = FALSE)
    expect_equal(dim(reducedDim(proj, "NMF_projected")), c(ncol(query), 4))
})

test_that("unsupported FactorNet projections explain why", {
    skip_if_not_installed("scuttle")
    sce <- make_sce()
    query <- scuttle::logNormCounts(scuttle::mockSCE(ngenes = 120, ncells = 30))
    rownames(query) <- rownames(sce)

    # a multi-layer graph is still reachable through a custom FactorNet
    mat <- SummarizedExperiment::assay(sce, "logcounts")
    inp <- RcppML::factor_input(mat, name = "input")
    chain <- RcppML::nmf_layer(RcppML::nmf_layer(inp, k = 10, name = "l1"),
                               k = 4, name = "l2")
    multi <- RcppML::fit(RcppML::factor_net(
        inputs = list(inp), output = chain,
        config = RcppML::factor_config(maxit = 10, seed = 1)))
    expect_error(predictNMF(query, reference = multi, verbose = FALSE),
                 "multi-layer")

    cond <- runConditionedNMF(sce, k = 4, condition_col = "batch",
                              verbose = FALSE)
    expect_error(evaluateNMF(cond, name = "CondNMF"), "total_loss")
})

test_that("S4 nmf results are unaffected", {
    skip_if_not_installed("scuttle")
    sce <- make_sce()
    query <- scuttle::logNormCounts(scuttle::mockSCE(ngenes = 120, ncells = 30))
    rownames(query) <- rownames(sce)

    res <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)
    expect_length(getDiagonal(res), 4)
    expect_type(evaluateNMF(res), "double")
    expect_equal(dim(reconstructNMF(res)), c(nrow(sce), ncol(sce)))
    proj <- predictNMF(query, reference = res, verbose = FALSE)
    expect_equal(dim(reducedDim(proj, "NMF_projected")), c(ncol(query), 4))
})
