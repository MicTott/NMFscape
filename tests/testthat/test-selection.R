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

test_that("programStability separates real programs from noise splits", {
    # data with exactly three real programs
    set.seed(4)
    ng <- 300; nc <- 200; k_true <- 3
    w <- matrix(0.3, ng, k_true)
    for (j in seq_len(k_true)) {
        w[((j - 1) * 100 + 1):(j * 100), j] <- runif(100, 1, 2)
    }
    grp <- sample(k_true, nc, replace = TRUE)
    h <- matrix(0.3, k_true, nc)
    h[cbind(grp, seq_len(nc))] <- runif(nc, 1, 2)
    mat <- w %*% h
    mat <- log1p(mat * exp(rnorm(length(mat), 0, 0.6)) / max(mat) * 60)
    dimnames(mat) <- list(paste0("g", seq_len(ng)), paste0("c", seq_len(nc)))
    sce <- SingleCellExperiment::SingleCellExperiment(list(logcounts = mat))

    # fitting past the true rank must leave the surplus programs unstable
    over <- consensusNMF(sce, k = 6, n_runs = 10, seed = 1, verbose = FALSE)
    st <- programStability(over)

    expect_equal(nrow(st), 6L)
    # the real programs stay perfectly reproducible; the surplus ones do not
    expect_equal(max(st$frequency), 1)
    expect_lt(min(st$frequency), 1)
    expect_lt(min(st$mean_similarity), max(st$mean_similarity))
})

test_that("runNMFscape chooses k when none is supplied", {
    sce <- makeSCE()
    expect_message(fit <- runNMFscape(sce, k_range = c(2, 4)),
                   "Selected k =")
    expect_true(ncol(reducedDim(fit, "NMF")) %in% c(2, 4))
})

test_that("enrichPrograms recovers the gene set a program was built from", {
    sce <- .testGuidedSCE(n_genes = 150, n_cells = 200, n_prog = 3, noise = 1)
    fit <- runNMFscape(sce, k = 3, seed = 1, verbose = FALSE)
    gene_sets <- split(rownames(fit),
                       rep(c("SET_1", "SET_2", "SET_3"), each = 50))

    res <- enrichPrograms(fit, gene_sets)
    expect_true(all(c("program", "pathway", "pval", "padj", "NES", "size") %in%
                        names(res)))

    # results are sorted by program then p-value, so the first row of each
    # program is its top hit: the gene block that program was simulated from
    top <- res[!duplicated(res$program), ]
    expect_setequal(top$pathway, names(gene_sets))
    expect_true(all(top$padj < 0.05))

    ora <- enrichPrograms(fit, gene_sets, method = "ora", n_top = 50)
    expect_setequal(ora[!duplicated(ora$program), "pathway"], names(gene_sets))
})
