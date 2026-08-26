setup_align_sce <- function(ngenes = 200, ncells = 100) {
    sce <- scuttle::mockSCE(ngenes = ngenes, ncells = ncells)
    scuttle::logNormCounts(sce)
}

test_that("alignPrograms recovers a permutation of a model's own programs", {
    skip_if_not_installed("scuttle")
    sce <- setup_align_sce()
    fit <- runNMFscape(sce, k = 5, seed = 1, verbose = FALSE)

    basis <- getBasis(fit)
    perm <- c(3L, 1L, 5L, 4L, 2L)
    permuted <- basis[, perm, drop = FALSE]
    colnames(permuted) <- paste0("Y_", seq_len(ncol(permuted)))

    res <- alignPrograms(basis, permuted)
    mapping <- res$mapping

    # program i of x sits at position match(i, perm) of the permuted copy
    expected <- paste0("Y_", match(seq_len(5), perm))
    names(expected) <- colnames(basis)

    expect_true(all(mapping$matched))
    expect_equal(unname(mapping$program_y), unname(expected[mapping$program_x]))
    expect_equal(mapping$similarity, rep(1, 5), tolerance = 1e-8)
    expect_length(res$unmatched_y, 0)
})

test_that("alignPrograms returns a mapping and a full similarity matrix", {
    skip_if_not_installed("scuttle")
    sce <- setup_align_sce()
    fit_a <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)
    fit_b <- runNMFscape(sce, k = 4, seed = 2, verbose = FALSE)

    res <- alignPrograms(fit_a, fit_b)

    expect_s3_class(res, "programAlignment")
    expect_s3_class(res$mapping, "data.frame")
    expect_equal(colnames(res$mapping),
                 c("program_x", "program_y", "similarity", "matched"))
    expect_equal(nrow(res$mapping), 4)
    expect_equal(dim(res$similarity), c(4, 4))
    expect_equal(rownames(res$similarity), colnames(getBasis(fit_a)))
    expect_equal(colnames(res$similarity), colnames(getBasis(fit_b)))
    expect_equal(res$features, nrow(sce))
    expect_identical(res$method, "cosine")

    # every program in x gets a distinct partner: a greedy argmax need not
    expect_false(anyDuplicated(res$mapping$program_y) > 0)
    # the reported similarity is the matrix entry for the chosen pair
    for (i in seq_len(nrow(res$mapping))) {
        expect_equal(res$mapping$similarity[i],
                     res$similarity[res$mapping$program_x[i],
                                    res$mapping$program_y[i]])
    }
    # sorted by decreasing similarity
    expect_false(is.unsorted(rev(res$mapping$similarity)))
})

test_that("alignPrograms handles a smaller k in x than in y", {
    skip_if_not_installed("scuttle")
    sce <- setup_align_sce()
    fit_a <- runNMFscape(sce, k = 3, seed = 1, verbose = FALSE)
    fit_b <- runNMFscape(sce, k = 6, seed = 2, verbose = FALSE)

    res <- alignPrograms(fit_a, fit_b)

    expect_equal(dim(res$similarity), c(3, 6))
    expect_equal(nrow(res$mapping), 3)
    expect_true(all(res$mapping$matched))
    expect_length(res$unmatched_y, 3)
    expect_false(any(res$unmatched_y %in% res$mapping$program_y))
})

test_that("alignPrograms reports unmatched programs when k in x is larger", {
    skip_if_not_installed("scuttle")
    sce <- setup_align_sce()
    fit_a <- runNMFscape(sce, k = 6, seed = 1, verbose = FALSE)
    fit_b <- runNMFscape(sce, k = 3, seed = 2, verbose = FALSE)

    res <- alignPrograms(fit_a, fit_b)

    expect_equal(dim(res$similarity), c(6, 3))
    expect_equal(nrow(res$mapping), 6)
    expect_equal(sum(res$mapping$matched), 3)
    expect_true(all(is.na(res$mapping$program_y[!res$mapping$matched])))
    expect_true(all(is.na(res$mapping$similarity[!res$mapping$matched])))
    expect_length(res$unmatched_y, 0)
})

test_that("alignPrograms intersects differing gene sets", {
    skip_if_not_installed("scuttle")
    sce <- setup_align_sce()
    fit_a <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)
    fit_b <- runNMFscape(sce[seq_len(150), ], k = 4, seed = 2, verbose = FALSE)

    expect_message(res <- alignPrograms(fit_a, fit_b), "150 shared features")
    expect_equal(res$features, 150)
    expect_equal(nrow(res$mapping), 4)
})

test_that("alignPrograms accepts nmf models and bare matrices", {
    skip_if_not_installed("scuttle")
    sce <- setup_align_sce()
    fit_a <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)
    fit_b <- runNMFscape(sce, k = 4, seed = 2, verbose = FALSE)

    from_models <- alignPrograms(getModel(fit_a), getModel(fit_b))
    expect_equal(dim(from_models$similarity), c(4, 4))
    expect_equal(nrow(from_models$mapping), 4)

    from_matrices <- alignPrograms(getBasis(fit_a), getBasis(fit_b))
    expect_equal(dim(from_matrices$similarity), c(4, 4))
})

test_that("alignPrograms supports correlation as well as cosine", {
    skip_if_not_installed("scuttle")
    sce <- setup_align_sce()
    fit_a <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)
    fit_b <- runNMFscape(sce, k = 4, seed = 2, verbose = FALSE)

    res <- alignPrograms(fit_a, fit_b, method = "cor")
    expect_identical(res$method, "cor")
    expect_equal(res$similarity,
                 stats::cor(getBasis(fit_a), getBasis(fit_b)))
    expect_error(alignPrograms(fit_a, fit_b, method = "spearman"))
})

test_that("alignPrograms errors informatively on bad input", {
    skip_if_not_installed("scuttle")
    sce <- setup_align_sce()
    fit <- runNMFscape(sce, k = 3, seed = 1, verbose = FALSE)
    basis <- getBasis(fit)

    expect_error(alignPrograms(fit, "not a model"),
                 "SingleCellExperiment, an S4 nmf model")

    disjoint <- basis
    rownames(disjoint) <- paste0("other_", seq_len(nrow(disjoint)))
    expect_error(alignPrograms(basis, disjoint), "share 0 features")

    # no rownames and different feature counts cannot be reconciled
    unnamed_short <- basis[seq_len(50), , drop = FALSE]
    rownames(unnamed_short) <- NULL
    expect_error(alignPrograms(basis, unnamed_short), "no rownames")
})

test_that("alignPrograms warns when feature names are missing", {
    skip_if_not_installed("scuttle")
    sce <- setup_align_sce()
    fit <- runNMFscape(sce, k = 3, seed = 1, verbose = FALSE)
    basis <- getBasis(fit)
    unnamed <- basis
    rownames(unnamed) <- NULL

    expect_warning(res <- alignPrograms(basis, unnamed), "same order")
    expect_equal(res$mapping$similarity, rep(1, 3), tolerance = 1e-8)
})

test_that("alignPrograms works with k = 1", {
    skip_if_not_installed("scuttle")
    sce <- setup_align_sce()
    fit <- runNMFscape(sce, k = 1, seed = 1, verbose = FALSE)

    res <- alignPrograms(fit, fit)
    expect_equal(dim(res$similarity), c(1, 1))
    expect_equal(res$mapping$similarity, 1, tolerance = 1e-8)
})

test_that("plotProgramSimilarity puts matched pairs on the diagonal", {
    skip_if_not_installed("scuttle")
    sce <- setup_align_sce()
    fit_a <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)
    fit_b <- runNMFscape(sce, k = 4, seed = 2, verbose = FALSE)
    res <- alignPrograms(fit_a, fit_b)

    p <- plotProgramSimilarity(res)
    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), 16)
    expect_true(all(c("program_x", "program_y", "similarity", "is_match") %in%
                        colnames(p$data)))
    expect_equal(sum(p$data$is_match), 4)

    # x runs left to right, y top to bottom, so matched pairs land where the
    # index of x equals the reversed index of y
    matched <- p$data[p$data$is_match, ]
    expect_true(all(as.integer(matched$program_x) ==
                        nlevels(matched$program_y) -
                        as.integer(matched$program_y) + 1L))

    unordered <- plotProgramSimilarity(res, order_by_match = FALSE)
    expect_equal(levels(unordered$data$program_y),
                 rev(colnames(res$similarity)))
})

test_that("plotProgramSimilarity handles non-square and rejects bad input", {
    skip_if_not_installed("scuttle")
    sce <- setup_align_sce()
    fit_a <- runNMFscape(sce, k = 3, seed = 1, verbose = FALSE)
    fit_b <- runNMFscape(sce, k = 5, seed = 2, verbose = FALSE)

    p <- plotProgramSimilarity(alignPrograms(fit_a, fit_b), show_values = FALSE)
    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), 15)

    expect_error(plotProgramSimilarity(getBasis(fit_a)),
                 "programAlignment object")
})
