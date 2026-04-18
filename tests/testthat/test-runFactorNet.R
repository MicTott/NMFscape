test_that("runFactorNet stores single-input results", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    library(RcppML)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    mat <- as.matrix(assay(sce, "logcounts"))
    mat[mat < 0] <- 0
    inp <- factor_input(mat, name = "input")
    layer <- nmf_layer(inp, k = 5, name = "my_nmf")
    net <- factor_net(inputs = list(inp), output = layer,
                      config = factor_config(verbose = FALSE))
    result <- fit(net)

    sce <- runFactorNet(sce, result, layer_name = "my_nmf", verbose = FALSE)

    expect_true("FactorNet" %in% reducedDimNames(sce))
    expect_equal(ncol(reducedDim(sce, "FactorNet")), 5)
    expect_true("FactorNet_basis" %in% names(metadata(sce)))
    expect_true("FactorNet_model" %in% names(metadata(sce)))
})

test_that("runFactorNet errors on missing layer", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    library(RcppML)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    mat <- as.matrix(assay(sce, "logcounts"))
    mat[mat < 0] <- 0
    inp <- factor_input(mat)
    layer <- nmf_layer(inp, k = 3, name = "my_nmf")
    net <- factor_net(inputs = list(inp), output = layer,
                      config = factor_config(verbose = FALSE))
    result <- fit(net)

    expect_error(runFactorNet(sce, result, layer_name = "wrong_name"),
                 "Available layers")
})
