test_that("runConditionedNMF factors out a covariate", {
    sce <- runConditionedNMF(makeSCE(), k = 4, condition_col = "batch",
                             seed = 1, verbose = FALSE)
    expect_equal(dim(reducedDim(sce, "CondNMF")), c(120L, 4L))
    expect_equal(rownames(getBasis(sce, "CondNMF")), rownames(sce))
})

test_that("runMultiModalNMF shares one embedding across modalities", {
    sce <- makeSCE()
    atac <- SingleCellExperiment::SingleCellExperiment(
        list(counts = abs(matrix(rnorm(40 * ncol(sce)), 40, ncol(sce)))))
    rownames(atac) <- paste0("peak_", seq_len(40))
    colnames(atac) <- colnames(sce)
    altExp(sce, "ATAC") <- atac

    sce <- runMultiModalNMF(sce, k = 4, alt_exps = c("main", "ATAC"),
                            modality_names = c("RNA", "ATAC"),
                            seed = 1, verbose = FALSE)

    expect_equal(ncol(reducedDim(sce, "MultiNMF")), 4L)
    expect_equal(nrow(getBasis(sce, "MultiNMF", modality = "RNA")), 150L)
    expect_equal(nrow(getBasis(sce, "MultiNMF", modality = "ATAC")), 40L)
    # asking without a modality should say which are available
    expect_error(getBasis(sce, "MultiNMF"), "ATAC")
})

test_that("runGuidedNMF enriches toward the labels", {
    sce <- .testGuidedSCE()
    sce <- runGuidedNMF(sce, k = 5, label_col = "program", strength = 1,
                        name = "Enriched", init_name = "Base", seed = 42,
                        verbose = FALSE)

    base <- adjustedRandIndex(.argmaxLabels(sce, "Base"), sce$program)
    guided <- adjustedRandIndex(.argmaxLabels(sce, "Enriched"), sce$program)
    expect_gt(guided, base)
    expect_s4_class(getModel(sce, "Enriched"), "nmf")
})

test_that("runSpatialNMF uses the spatial graph, and lambda 0 does not", {
    spe <- makeSpatialGrid()
    plain <- runNMFscape(spe, k = 4, seed = 1, verbose = FALSE)
    spatial <- runSpatialNMF(spe, k = 4, seed = 1, verbose = FALSE)
    off <- runSpatialNMF(spe, k = 4, graph_lambda = 0, seed = 1,
                         verbose = FALSE)

    # switching the penalty off must recover plain NMF exactly
    expect_equal(unname(reducedDim(off, "SpatialNMF")),
                 unname(reducedDim(plain, "NMF")), tolerance = 1e-8)

    # and switching it on must recover the true domains better
    cl <- function(m) apply(m, 1, which.max)
    expect_gt(adjustedRandIndex(cl(reducedDim(spatial, "SpatialNMF")),
                                spe$domain),
              adjustedRandIndex(cl(reducedDim(plain, "NMF")), spe$domain))

    expect_equal(nrow(spatialAutocorrelation(spatial, "SpatialNMF")), 4L)
})

test_that("runFactorNet stores a custom graph result", {
    sce <- makeSCE()
    mat <- assay(sce, "logcounts")
    inp <- RcppML::factor_input(mat, name = "input")
    net <- RcppML::factor_net(
        inputs = list(inp),
        output = RcppML::nmf_layer(inp, k = 4, name = "my_nmf"),
        config = RcppML::factor_config(maxit = 20, seed = 1))

    sce <- runFactorNet(sce, RcppML::fit(net), layer_name = "my_nmf",
                        verbose = FALSE)
    expect_equal(ncol(reducedDim(sce, "FactorNet")), 4L)
})
