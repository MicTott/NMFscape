# L1/L2 were accepted by the FactorNet wrappers but never reached nmf_layer(),
# so regularization was silently ignored.

test_that("L1 changes the factorization in FactorNet wrappers", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 120, ncells = 60)
    sce <- logNormCounts(sce)
    sce$batch <- rep(c("A", "B"), length.out = ncol(sce))

    plain <- runConditionedNMF(sce, k = 4, condition_col = "batch",
                               seed = 1, L1 = 0, verbose = FALSE)
    pen <- runConditionedNMF(sce, k = 4, condition_col = "batch",
                             seed = 1, L1 = 0.5, verbose = FALSE)
    expect_false(identical(getBasis(plain, "CondNMF"),
                           getBasis(pen, "CondNMF")))
})

test_that("length-2 penalties are rejected with an explanation", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    sce$batch <- rep(c("A", "B"), length.out = ncol(sce))
    expect_error(runConditionedNMF(sce, k = 3, condition_col = "batch",
                                   L1 = c(0, 0), verbose = FALSE),
                 "single numeric value")
    expect_error(runConditionedNMF(sce, k = 3, condition_col = "batch",
                                   L2 = c(0.1, 0.1), verbose = FALSE),
                 "single numeric value")
})

test_that("runNMFscape still takes length-2 L1/L2", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    plain <- runNMFscape(sce, k = 4, seed = 1, L1 = c(0, 0), verbose = FALSE)
    pen <- runNMFscape(sce, k = 4, seed = 1, L1 = c(0.5, 0.5), verbose = FALSE)
    expect_false(identical(getBasis(plain), getBasis(pen)))
})

test_that("distribution = 'auto' selects a real loss", {
    skip_if_not_installed("scuttle")

    library(scuttle)
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)

    expect_message(
        sce <- runNMFscape(sce, k = 4, distribution = "auto", seed = 1),
        "Selected distribution: [a-z]"
    )
    expect_true("NMF" %in% reducedDimNames(sce))
})
