test_that("assessNMF returns assessment metrics", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 200)
    sce <- logNormCounts(sce)
    sce <- runNMFscape(sce, k = 5, verbose = FALSE)

    sce$celltype <- sample(c("TypeA", "TypeB", "TypeC"), ncol(sce),
                           replace = TRUE)

    result <- assessNMF(sce, label_col = "celltype", min_class_size = 5)

    expect_s3_class(result, "nmf_assessment")
    expect_true("metrics" %in% names(result))
    expect_true("classification" %in% names(result))
    expect_true("ari" %in% names(result$metrics))
    expect_true("nmi" %in% names(result$metrics))
    expect_true("silhouette" %in% names(result$metrics))
})

test_that("assessNMF errors on missing label_col", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)
    sce <- runNMFscape(sce, k = 3, verbose = FALSE)

    expect_error(assessNMF(sce, label_col = "nonexistent"),
                 "not found in colData")
})
