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

    deep <- runDeepNMF(sce, k = c(10, 4), verbose = FALSE)
    expect_length(getDiagonal(deep, "DeepNMF"), 4)
})

test_that("reconstructNMF works for FactorNet results", {
    skip_if_not_installed("scuttle")
    sce <- make_sce()

    deep <- runDeepNMF(sce, k = c(10, 4), verbose = FALSE)
    recon <- reconstructNMF(deep, "DeepNMF")
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

    deep <- runDeepNMF(sce, k = c(10, 4), verbose = FALSE)
    expect_error(
        predictNMF(query, reference = deep, ref_name = "DeepNMF",
                   verbose = FALSE),
        "multi-layer"
    )
    expect_error(evaluateNMF(deep, name = "DeepNMF"), "total_loss")
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
