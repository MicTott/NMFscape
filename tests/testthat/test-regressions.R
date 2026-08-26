# One block per defect that actually shipped. Each names the bug it pins so a
# future change that reintroduces it fails with an explanation.

test_that("character subset_row names the basis rows (was: NA rownames)", {
    sce <- makeSCE()
    hvg <- rownames(sce)[seq_len(40)]

    for (fit in list(runNMFscape(sce, k = 3, subset_row = hvg, verbose = FALSE),
                     consensusNMF(sce, k = 3, subset_row = hvg, n_runs = 3,
                                  verbose = FALSE),
                     runConditionedNMF(sce, k = 3, condition_col = "batch",
                                       subset_row = hvg, verbose = FALSE))) {
        nm <- setdiff(names(metadata(fit)), "")
        basis <- metadata(fit)[[grep("_basis$", nm, value = TRUE)[1]]]
        expect_equal(rownames(basis), hvg)
    }
    # integer indexing must keep working, and unknown names must error
    expect_equal(rownames(getBasis(runNMFscape(sce, k = 3,
                                               subset_row = seq_len(40),
                                               verbose = FALSE))), hvg)
    expect_error(runNMFscape(sce, k = 3, subset_row = "nope",
                             verbose = FALSE), "not found in rownames")
})

test_that("distribution = 'auto' selects a loss (was: silent no-op)", {
    expect_message(runNMFscape(makeSCE(), k = 4, distribution = "auto",
                               seed = 1),
                   "Selected distribution: [a-z]")
})

test_that("L1 reaches nmf_layer (was: silently dropped by the recipes)", {
    sce <- makeSCE()
    plain <- runConditionedNMF(sce, k = 4, condition_col = "batch",
                               seed = 1, L1 = 0, verbose = FALSE)
    pen <- runConditionedNMF(sce, k = 4, condition_col = "batch",
                             seed = 1, L1 = 0.5, verbose = FALSE)
    expect_false(identical(getBasis(plain, "CondNMF"),
                           getBasis(pen, "CondNMF")))

    # the recipes take a scalar, unlike runNMFscape's c(w, h)
    expect_error(runConditionedNMF(sce, k = 3, condition_col = "batch",
                                   L1 = c(0, 0), verbose = FALSE),
                 "single numeric value")
})

test_that("FactorNet results work with the accessors (was: @ on a list)", {
    sce <- makeSCE()
    cond <- runConditionedNMF(sce, k = 4, condition_col = "batch",
                              verbose = FALSE)

    expect_length(getDiagonal(cond, "CondNMF"), 4)
    expect_equal(dim(reconstructNMF(cond, "CondNMF")), dim(sce))
    # and the ones that genuinely cannot must say why
    expect_error(evaluateNMF(cond, name = "CondNMF"), "total_loss")
})

test_that("refineNMF cycles controls whether W moves (was: undocumented)", {
    fit <- runNMFscape(makeSCE(), k = 4, seed = 1, verbose = FALSE)
    basis0 <- getBasis(fit)

    kept <- refineNMF(fit, label_col = "celltype", verbose = FALSE)
    expect_identical(unname(getBasis(kept, "NMF_refined")), unname(basis0))

    moved <- refineNMF(fit, label_col = "celltype", cycles = 5,
                       refined_name = "c5", verbose = FALSE)
    expect_gt(max(abs(getBasis(moved, "c5") - basis0)), 0)

    # a Matrix-backed assay used to fail the S4 check on cycles > 0
    m <- makeSCE(ngenes = 100, ncells = 80)
    assay(m, "logcounts") <- Matrix::Matrix(as.matrix(assay(m, "logcounts")))
    expect_no_error(refineNMF(runNMFscape(m, k = 3, seed = 1, verbose = FALSE),
                              label_col = "celltype", cycles = 5,
                              verbose = FALSE))
})

test_that("guided removal is refused (was: degenerate factorization)", {
    expect_error(runGuidedNMF(makeSCE(), k = 4, label_col = "celltype",
                              mode = "remove", verbose = FALSE),
                 "runConditionedNMF")
})
