test_that("selectRank cross-validates a range of ranks", {
    cv <- selectRank(makeSCE(), k = 2:5, seed = 1, verbose = FALSE)

    expect_equal(cv$k, 2:5)
    expect_true(all(cv$test_error > 0))
    expect_s3_class(plotRankSelection(cv), "ggplot")

    reps <- selectRank(makeSCE(), k = 2:4, n_reps = 2, seed = 1,
                       verbose = FALSE)
    expect_true("sd_test_error" %in% names(reps))
})

test_that("tuneNMF picks a rank that transfers back to the fitter", {
    sce <- makeTuneSCE()
    fit <- tuneNMF(sce, params = list(k = c(2, 4)), reps = 2, verbose = FALSE)

    expect_true(fit$best_params$k %in% c(2, 4))
    out <- do.call(runNMFscape, c(list(x = sce, verbose = FALSE),
                                  fit$best_params))
    expect_equal(ncol(reducedDim(out, "NMF")), fit$best_params$k)
    expect_s3_class(plotTuning(fit), "ggplot")
})

test_that("the diagnostics report on a fitted model", {
    sce <- makeSCE()
    fit <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)

    expect_type(evaluateNMF(fit), "double")
    expect_equal(diagnoseNMF(fit)$k, 4L)
    expect_type(selectDistribution(sce, k = 4, verbose = FALSE)$best,
                "character")
    expect_type(assessNMF(fit, label_col = "celltype",
                          min_class_size = 5)$metrics$ari, "double")
    expect_type(gpuInfo()$available, "logical")
})

test_that("alignPrograms matches programs between independent fits", {
    sce <- makeSCE()
    a <- runNMFscape(sce, k = 5, seed = 1, verbose = FALSE)

    # a known permutation of a model's own programs must align perfectly
    b <- a
    basis <- metadata(a)$NMF_basis[, c(4L, 1L, 5L, 2L, 3L)]
    colnames(basis) <- paste0("NMF_", seq_len(5))
    metadata(b)$NMF_basis <- basis

    al <- alignPrograms(a, b)
    expect_true(all(al$mapping$matched))
    expect_equal(al$mapping$similarity, rep(1, 5), tolerance = 1e-8)

    # a smaller k on one side leaves the surplus unmatched rather than forced
    c2 <- runNMFscape(sce, k = 3, seed = 2, name = "NMF2", verbose = FALSE)
    al2 <- alignPrograms(a, c2, name_y = "NMF2")
    expect_equal(sum(al2$mapping$matched), 3L)
    expect_s3_class(plotProgramSimilarity(al2), "ggplot")
})
