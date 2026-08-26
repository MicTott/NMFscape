setup_align_sce <- function(ngenes = 200, ncells = 100) {
    sce <- scuttle::mockSCE(ngenes = ngenes, ncells = ncells)
    scuttle::logNormCounts(sce)
}

test_that("alignPrograms recovers a known permutation exactly", {
    sce <- setup_align_sce()
    a <- runNMFscape(sce, k = 5, seed = 1, verbose = FALSE)

    perm <- c(4L, 1L, 5L, 2L, 3L)
    b <- a
    basis <- metadata(a)$NMF_basis[, perm]
    colnames(basis) <- paste0("NMF_", seq_len(5))
    metadata(b)$NMF_basis <- basis

    al <- alignPrograms(a, b)
    expect_equal(nrow(al$mapping), 5)
    expect_true(all(al$mapping$matched))
    expect_equal(al$mapping$similarity, rep(1, 5), tolerance = 1e-8)
})

test_that("alignPrograms handles differing k and gene sets", {
    sce <- setup_align_sce()
    a <- runNMFscape(sce, k = 5, seed = 1, verbose = FALSE)
    b <- runNMFscape(sce, k = 3, seed = 2, name = "NMF2", verbose = FALSE)

    al <- alignPrograms(a, b, name_y = "NMF2")
    expect_equal(nrow(al$mapping), 5)
    expect_equal(sum(al$mapping$matched), 3)
    expect_equal(dim(al$similarity), c(5, 3))

    # differing gene sets are intersected
    c2 <- runNMFscape(sce[seq_len(120), ], k = 5, seed = 1, verbose = FALSE)
    expect_equal(nrow(alignPrograms(a, c2)$mapping), 5)
})

test_that("plotProgramSimilarity returns a ggplot", {
    sce <- setup_align_sce()
    a <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)
    b <- runNMFscape(sce, k = 4, seed = 2, name = "NMF2", verbose = FALSE)
    expect_s3_class(plotProgramSimilarity(alignPrograms(a, b, name_y = "NMF2")),
                    "ggplot")
})
