# A 20x20 grid of spots split into four contiguous quadrant domains, each
# driven by its own 30-gene program on top of a shared background. Counts are
# Poisson and deliberately shallow (mean well under one per gene) so that a
# single spot is uninformative and plain NMF produces salt-and-pepper
# programs -- the regime where a spatial prior should help.
makeSpatialGrid <- function(side = 20, n_genes = 120, seed = 1,
                            signal = 0.3, background = 0.3) {
    set.seed(seed)
    grid <- expand.grid(x = seq_len(side), y = seq_len(side))
    half <- side / 2
    domain <- ifelse(grid$x <= half & grid$y <= half, 1L,
              ifelse(grid$x > half & grid$y <= half, 2L,
              ifelse(grid$x <= half & grid$y > half, 3L, 4L)))

    n_spots <- nrow(grid)
    per_domain <- n_genes / 4
    mu <- matrix(background, n_genes, n_spots)
    for (d in seq_len(4)) {
        rows <- ((d - 1) * per_domain + 1):(d * per_domain)
        mu[rows, domain == d] <- background + signal
    }

    counts <- matrix(rpois(n_genes * n_spots, mu), n_genes, n_spots)
    rownames(counts) <- paste0("gene", seq_len(n_genes))
    colnames(counts) <- paste0("spot", seq_len(n_spots))

    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(logcounts = log1p(counts)),
        spatialCoords = as.matrix(grid)
    )
    spe$domain <- factor(domain)
    spe
}

adjustedRandIndex <- function(a, b) {
    tab <- table(a, b)
    choose2 <- function(z) z * (z - 1) / 2
    sum_ij <- sum(choose2(tab))
    sum_i <- sum(choose2(rowSums(tab)))
    sum_j <- sum(choose2(colSums(tab)))
    expected <- sum_i * sum_j / choose2(sum(tab))
    (sum_ij - expected) / ((sum_i + sum_j) / 2 - expected)
}

test_that("runSpatialNMF stores a usable result", {
    spe <- makeSpatialGrid()
    fit <- runSpatialNMF(spe, k = 4, seed = 1, verbose = FALSE)

    expect_equal(dim(reducedDim(fit, "SpatialNMF")), c(ncol(spe), 4))
    expect_equal(rownames(getBasis(fit, "SpatialNMF")), rownames(spe))
    expect_true(all(getBasis(fit, "SpatialNMF") >= 0))
})

test_that("graph_lambda = 0 reproduces plain NMF exactly", {
    spe <- makeSpatialGrid()
    plain <- runNMFscape(spe, k = 4, seed = 1, verbose = FALSE)
    off <- runSpatialNMF(spe, k = 4, graph_lambda = 0, seed = 1,
                         verbose = FALSE)
    expect_equal(unname(reducedDim(off, "SpatialNMF")),
                 unname(reducedDim(plain, "NMF")), tolerance = 1e-8)
})

test_that("spatial regularization recovers the true domains better", {
    spe <- makeSpatialGrid()
    plain <- runNMFscape(spe, k = 4, seed = 1, verbose = FALSE)
    spatial <- runSpatialNMF(spe, k = 4, seed = 1, verbose = FALSE)

    cl <- function(m) apply(m, 1, which.max)
    ari_plain <- adjustedRandIndex(cl(reducedDim(plain, "NMF")),
                                   spe$domain)
    ari_spatial <- adjustedRandIndex(cl(reducedDim(spatial, "SpatialNMF")),
                                     spe$domain)
    expect_gt(ari_spatial, ari_plain)
})

test_that("knn, distance and delaunay graphs all work", {
    spe <- makeSpatialGrid(side = 10, n_genes = 40)
    for (g in c("knn", "distance", "delaunay")) {
        if (g == "delaunay") skip_if_not_installed("deldir")
        fit <- runSpatialNMF(spe, k = 3, graph = g, seed = 1, verbose = FALSE)
        expect_equal(ncol(reducedDim(fit, "SpatialNMF")), 3)
    }
})

test_that("a plain SingleCellExperiment works with a supplied graph", {
    spe <- makeSpatialGrid(side = 10, n_genes = 40)
    ref <- runSpatialNMF(spe, k = 3, seed = 1, verbose = FALSE)
    graph <- metadata(ref)$SpatialNMF_graph

    sce <- SingleCellExperiment::SingleCellExperiment(
        list(logcounts = SummarizedExperiment::assay(spe, "logcounts")))
    fit <- runSpatialNMF(sce, k = 3, graph = graph, seed = 1, verbose = FALSE)
    expect_equal(unname(reducedDim(fit, "SpatialNMF")),
                 unname(reducedDim(ref, "SpatialNMF")), tolerance = 1e-8)

    # a named graph type needs coordinates
    expect_error(runSpatialNMF(sce, k = 3, verbose = FALSE), "SpatialExperiment")
})

test_that("spatialAutocorrelation reports Moran's I per program", {
    spe <- makeSpatialGrid()
    fit <- runSpatialNMF(spe, k = 4, seed = 1, verbose = FALSE)
    mi <- spatialAutocorrelation(fit, "SpatialNMF")

    expect_equal(nrow(mi), 4)
    expect_true(all(mi$morans_i > 0))
})
