.makeMultiModalSCE <- function() {
    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    # Create altExp for ATAC modality
    atac <- SingleCellExperiment(list(
        counts = abs(matrix(rnorm(50 * ncol(sce)), 50, ncol(sce)))
    ))
    rownames(atac) <- paste0("peak_", seq_len(50))
    colnames(atac) <- colnames(sce)
    altExp(sce, "ATAC") <- atac
    sce
}

test_that("runMultiModalNMF works with altExps", {
    skip_if_not_installed("scuttle")

    sce <- .makeMultiModalSCE()

    sce <- runMultiModalNMF(sce, k = 5,
                            alt_exps = c("main", "ATAC"),
                            modality_names = c("RNA", "ATAC"),
                            verbose = FALSE)

    # Shared cell embedding
    expect_true("MultiNMF" %in% reducedDimNames(sce))
    expect_equal(nrow(reducedDim(sce, "MultiNMF")), ncol(sce))
    expect_equal(ncol(reducedDim(sce, "MultiNMF")), 5)

    # Per-modality basis
    rna_basis <- getBasis(sce, "MultiNMF", modality = "RNA")
    expect_equal(nrow(rna_basis), 100)
    expect_equal(ncol(rna_basis), 5)

    atac_basis <- getBasis(sce, "MultiNMF", modality = "ATAC")
    expect_equal(nrow(atac_basis), 50)
    expect_equal(ncol(atac_basis), 5)

    # Modalities metadata
    expect_equal(metadata(sce)$MultiNMF_modalities, c("RNA", "ATAC"))
    expect_true("MultiNMF_model" %in% names(metadata(sce)))
})

test_that("getBasis errors helpfully for multi-modal without modality", {
    skip_if_not_installed("scuttle")

    sce <- .makeMultiModalSCE()
    sce <- runMultiModalNMF(sce, k = 3,
                            alt_exps = c("main", "ATAC"),
                            modality_names = c("RNA", "ATAC"),
                            verbose = FALSE)

    expect_error(getBasis(sce, "MultiNMF"), "requires a modality")
})

test_that("runMultiModalNMF errors on single modality", {
    skip_if_not_installed("scuttle")

    sce <- .makeMultiModalSCE()

    expect_error(runMultiModalNMF(sce, k = 5, alt_exps = "main"),
                 "at least 2 modalities")
})

test_that("getTopFeatures works with multi-modal", {
    skip_if_not_installed("scuttle")

    sce <- .makeMultiModalSCE()
    sce <- runMultiModalNMF(sce, k = 3,
                            alt_exps = c("main", "ATAC"),
                            modality_names = c("RNA", "ATAC"),
                            verbose = FALSE)

    top_rna <- getTopFeatures(sce, "MultiNMF", n = 5, modality = "RNA")
    expect_type(top_rna, "list")
    expect_length(top_rna, 3)

    top_atac <- getTopFeatures(sce, "MultiNMF", n = 5, modality = "ATAC")
    expect_type(top_atac, "list")
    expect_length(top_atac, 3)
})

test_that("runMultiModalNMF errors when neither assays nor alt_exps given", {
    skip_if_not_installed("scuttle")

    sce <- .makeMultiModalSCE()
    expect_error(runMultiModalNMF(sce, k = 5), "Specify either")
})
