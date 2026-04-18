test_that("predictNMF projects new data onto existing model", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    ref_sce <- mockSCE(ngenes = 100, ncells = 50)
    ref_sce <- logNormCounts(ref_sce)
    ref_sce <- runNMFscape(ref_sce, k = 5, verbose = FALSE)

    query_sce <- mockSCE(ngenes = 100, ncells = 30)
    query_sce <- logNormCounts(query_sce)

    query_sce <- predictNMF(query_sce, reference = ref_sce, verbose = FALSE)

    expect_true("NMF_projected" %in% reducedDimNames(query_sce))
    expect_equal(nrow(reducedDim(query_sce, "NMF_projected")), 30)
    expect_equal(ncol(reducedDim(query_sce, "NMF_projected")), 5)
})

test_that("predictNMF works with raw model object", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    ref_sce <- mockSCE(ngenes = 100, ncells = 50)
    ref_sce <- logNormCounts(ref_sce)
    ref_sce <- runNMFscape(ref_sce, k = 3, verbose = FALSE)

    model <- getModel(ref_sce)

    query_sce <- mockSCE(ngenes = 100, ncells = 20)
    query_sce <- logNormCounts(query_sce)

    query_sce <- predictNMF(query_sce, reference = model, verbose = FALSE)

    expect_true("NMF_projected" %in% reducedDimNames(query_sce))
    expect_equal(ncol(reducedDim(query_sce, "NMF_projected")), 3)
})
