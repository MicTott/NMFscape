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
