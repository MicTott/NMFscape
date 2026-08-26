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


test_that("runSpatialNMF stores results the standard accessors understand", {
    spe <- makeSpatialGrid()
    spe <- runSpatialNMF(spe, k = 4, seed = 7, verbose = FALSE)

    expect_true("SpatialNMF" %in% reducedDimNames(spe))
    expect_equal(dim(reducedDim(spe, "SpatialNMF")), c(ncol(spe), 4L))
    expect_equal(rownames(reducedDim(spe, "SpatialNMF")), colnames(spe))
    expect_equal(colnames(reducedDim(spe, "SpatialNMF")),
                 paste0("NMF_", seq_len(4)))

    basis <- getBasis(spe, "SpatialNMF")
    expect_equal(dim(basis), c(nrow(spe), 4L))
    expect_equal(rownames(basis), rownames(spe))

    expect_equal(dim(getCoefficients(spe, "SpatialNMF")), c(ncol(spe), 4L))
    expect_length(getTopFeatures(spe, "SpatialNMF", n = 5), 4L)
    expect_length(getTopFeatures(spe, "SpatialNMF", n = 5)[[1]], 5L)
    expect_length(getDiagonal(spe, "SpatialNMF"), 4L)
    expect_s4_class(getModel(spe, "SpatialNMF"), "nmf")
    expect_equal(dim(reconstructNMF(spe, "SpatialNMF")), dim(spe))

    # Factors stay in the feasible region.
    expect_true(all(reducedDim(spe, "SpatialNMF") >= 0))
    expect_true(all(basis >= 0))
})

test_that("runSpatialNMF records the graph and its parameters", {
    spe <- makeSpatialGrid()
    spe <- runSpatialNMF(spe, k = 4, n_neighbors = 6, graph_lambda = 2,
                         seed = 7, verbose = FALSE)

    graph <- metadata(spe)$SpatialNMF_graph
    expect_s4_class(graph, "dgCMatrix")
    expect_equal(dim(graph), c(ncol(spe), ncol(spe)))
    expect_true(Matrix::nnzero(graph) > 0)
    # Symmetric, with no self-loops.
    expect_equal(as.matrix(graph), t(as.matrix(graph)))
    expect_true(all(Matrix::diag(graph) == 0))

    params <- metadata(spe)$SpatialNMF_graph_params
    expect_equal(params$graph, "knn")
    expect_equal(params$n_neighbors, 6)
    expect_equal(params$graph_lambda, 2)
    expect_true(params$normalize_laplacian)
    expect_equal(params$seed, 7)
})

test_that("increasing graph_lambda increases spatial smoothness", {
    spe <- makeSpatialGrid()
    lambdas <- c(0.03, 0.1, 0.3, 1, 3)

    fits <- lapply(lambdas, function(lam) {
        runSpatialNMF(spe, k = 4, graph_lambda = lam, seed = 7,
                      verbose = FALSE)
    })

    # Score every fit against one fixed graph so the comparison is fair.
    graph <- metadata(fits[[1]])$SpatialNMF_graph
    moran <- vapply(fits, function(fit) {
        mean(spatialAutocorrelation(fit, graph = graph)$moran_I)
    }, FUN.VALUE = numeric(1))

    expect_true(all(diff(moran) > 0))

    plain <- runNMFscape(spe, k = 4, seed = 7, verbose = FALSE)
    plain_moran <- mean(
        spatialAutocorrelation(plain, name = "NMF", graph = graph)$moran_I
    )
    # The unregularized fit is the least spatially coherent of all.
    expect_true(plain_moran < min(moran))
    # And the penalty makes a large difference, not a rounding one.
    expect_true(max(moran) - plain_moran > 0.2)
})

test_that("spatial regularization recovers the true domains", {
    spe <- makeSpatialGrid()

    spatial <- runSpatialNMF(spe, k = 4, graph_lambda = 1, seed = 7,
                             verbose = FALSE)
    plain <- runNMFscape(spe, k = 4, seed = 7, verbose = FALSE)

    spatial_ari <- adjustedRandIndex(
        max.col(reducedDim(spatial, "SpatialNMF")), spe$domain
    )
    plain_ari <- adjustedRandIndex(
        max.col(reducedDim(plain, "NMF")), spe$domain
    )

    expect_true(spatial_ari > 0.9)
    expect_true(spatial_ari > plain_ari)
})

test_that("graph_lambda = 0 reproduces plain NMF exactly", {
    spe <- makeSpatialGrid()

    spatial <- runSpatialNMF(spe, k = 4, graph_lambda = 0, seed = 7,
                             verbose = FALSE)
    plain <- runNMFscape(spe, k = 4, seed = 7, verbose = FALSE)

    expect_equal(unname(reducedDim(spatial, "SpatialNMF")),
                 unname(reducedDim(plain, "NMF")))
    expect_equal(unname(getBasis(spatial, "SpatialNMF")),
                 unname(getBasis(plain, "NMF")))
})

test_that("knn, distance and delaunay graphs all work", {
    spe <- makeSpatialGrid()

    for (graph_type in c("knn", "distance")) {
        fit <- runSpatialNMF(spe, k = 4, graph = graph_type,
                             graph_lambda = 1, seed = 7, verbose = FALSE)
        expect_equal(metadata(fit)$SpatialNMF_graph_params$graph, graph_type)
        expect_true(Matrix::nnzero(metadata(fit)$SpatialNMF_graph) > 0)
        expect_true(mean(spatialAutocorrelation(fit)$moran_I) > 0.7)
    }

    skip_if_not_installed("deldir")
    fit <- runSpatialNMF(spe, k = 4, graph = "delaunay", graph_lambda = 1,
                         seed = 7, verbose = FALSE)
    expect_equal(metadata(fit)$SpatialNMF_graph_params$graph, "delaunay")
    expect_true(mean(spatialAutocorrelation(fit)$moran_I) > 0.7)
})

test_that("an explicit radius changes the distance graph", {
    spe <- makeSpatialGrid()

    tight <- runSpatialNMF(spe, k = 4, graph = "distance", radius = 1.1,
                           seed = 7, verbose = FALSE)
    wide <- runSpatialNMF(spe, k = 4, graph = "distance", radius = 3,
                          seed = 7, verbose = FALSE)

    expect_true(Matrix::nnzero(metadata(wide)$SpatialNMF_graph) >
                Matrix::nnzero(metadata(tight)$SpatialNMF_graph))
    expect_equal(metadata(tight)$SpatialNMF_graph_params$radius, 1.1)
})

test_that("a plain SingleCellExperiment works with a supplied graph", {
    spe <- makeSpatialGrid()
    reference <- runSpatialNMF(spe, k = 4, graph_lambda = 1, seed = 7,
                               verbose = FALSE)
    graph <- metadata(reference)$SpatialNMF_graph

    sce <- as(spe, "SingleCellExperiment")
    expect_false(is(sce, "SpatialExperiment"))

    fit <- runSpatialNMF(sce, k = 4, graph = graph, graph_lambda = 1,
                         seed = 7, verbose = FALSE)

    expect_true("SpatialNMF" %in% reducedDimNames(fit))
    expect_equal(metadata(fit)$SpatialNMF_graph_params$graph, "custom")
    # Same graph, same penalty, same seed: the same answer.
    expect_equal(unname(reducedDim(fit, "SpatialNMF")),
                 unname(reducedDim(reference, "SpatialNMF")))

    # Without a graph there are no coordinates to fall back on.
    expect_error(runSpatialNMF(sce, k = 4, verbose = FALSE),
                 "SpatialExperiment")
})

test_that("a dense adjacency matrix is accepted and symmetrized", {
    spe <- makeSpatialGrid(side = 6, n_genes = 20)
    n_spots <- ncol(spe)

    # Deliberately asymmetric chain.
    adjacency <- matrix(0, n_spots, n_spots)
    adjacency[cbind(seq_len(n_spots - 1), 2:n_spots)] <- 1

    fit <- runSpatialNMF(spe, k = 3, graph = adjacency, graph_lambda = 1,
                         seed = 7, verbose = FALSE)

    stored <- as.matrix(metadata(fit)$SpatialNMF_graph)
    expect_equal(stored, t(stored))
    expect_true(all(diag(stored) == 0))
})

test_that("the feature graph smooths gene loadings", {
    spe <- makeSpatialGrid()
    n_genes <- nrow(spe)

    # Chain over genes in index order: adjacent genes should get similar
    # loadings when the penalty is switched on.
    feature_graph <- Matrix::sparseMatrix(
        i = seq_len(n_genes - 1), j = 2:n_genes, x = 1,
        dims = c(n_genes, n_genes)
    )

    without <- runSpatialNMF(spe, k = 4, graph_lambda = 1, seed = 7,
                             verbose = FALSE)
    with_graph <- runSpatialNMF(spe, k = 4, graph_lambda = 1,
                                feature_graph = feature_graph,
                                feature_lambda = 1, seed = 7, verbose = FALSE)

    roughness <- function(basis) {
        z <- scale(basis)
        mean(abs(z[-1, ] - z[-nrow(z), ]))
    }

    expect_true(roughness(getBasis(with_graph, "SpatialNMF")) <
                roughness(getBasis(without, "SpatialNMF")))
    expect_true(metadata(with_graph)$SpatialNMF_graph_params$feature_graph)
})

test_that("a feature graph with zero lambda warns", {
    spe <- makeSpatialGrid(side = 6, n_genes = 20)
    feature_graph <- diag(0, nrow(spe))
    feature_graph[1, 2] <- feature_graph[2, 1] <- 1

    expect_warning(
        runSpatialNMF(spe, k = 3, feature_graph = feature_graph,
                      feature_lambda = 0, seed = 7, verbose = FALSE),
        "feature_lambda is 0"
    )
})

test_that("results are reproducible with a seed", {
    spe <- makeSpatialGrid()

    first <- runSpatialNMF(spe, k = 4, graph_lambda = 1, seed = 42,
                           verbose = FALSE)
    second <- runSpatialNMF(spe, k = 4, graph_lambda = 1, seed = 42,
                            verbose = FALSE)

    expect_equal(reducedDim(first, "SpatialNMF"),
                 reducedDim(second, "SpatialNMF"))
    expect_equal(getBasis(first, "SpatialNMF"), getBasis(second, "SpatialNMF"))
})

test_that("subset_row is honored by name", {
    spe <- makeSpatialGrid()
    chosen <- paste0("gene", seq(1, 120, by = 2))

    fit <- runSpatialNMF(spe, k = 4, subset_row = chosen, graph_lambda = 1,
                         seed = 7, verbose = FALSE)

    basis <- getBasis(fit, "SpatialNMF")
    expect_equal(nrow(basis), length(chosen))
    expect_equal(rownames(basis), chosen)
})

test_that("normalize_laplacian changes the fit", {
    spe <- makeSpatialGrid()

    normalized <- runSpatialNMF(spe, k = 4, graph_lambda = 1,
                                normalize_laplacian = TRUE, seed = 7,
                                verbose = FALSE)
    combinatorial <- runSpatialNMF(spe, k = 4, graph_lambda = 1,
                                   normalize_laplacian = FALSE, seed = 7,
                                   verbose = FALSE)

    expect_false(isTRUE(all.equal(
        unname(reducedDim(normalized, "SpatialNMF")),
        unname(reducedDim(combinatorial, "SpatialNMF"))
    )))
    expect_false(
        metadata(combinatorial)$SpatialNMF_graph_params$normalize_laplacian
    )
})

test_that("runSpatialNMF validates its inputs", {
    spe <- makeSpatialGrid(side = 6, n_genes = 20)

    expect_error(runSpatialNMF(matrix(1:10, nrow = 2), k = 2),
                 "SingleCellExperiment")
    expect_error(runSpatialNMF(spe, k = 3, assay = "missing"),
                 "not found")
    expect_error(runSpatialNMF(spe, k = 3, graph_lambda = -1),
                 "non-negative")
    expect_error(runSpatialNMF(spe, k = 3, feature_lambda = -1),
                 "non-negative")
    expect_error(runSpatialNMF(spe, k = 3, graph = "voronoi"),
                 "should be one of")
    expect_error(runSpatialNMF(spe, k = 3, graph = diag(5)),
                 "must be a 36 x 36 matrix")
    expect_error(runSpatialNMF(spe, k = 3, graph = -matrix(1, 36, 36)),
                 "non-negative edge weights")
    expect_error(runSpatialNMF(spe, k = 3, n_neighbors = 100),
                 "must be smaller than")
    expect_error(runSpatialNMF(spe, k = 3, graph = "distance", radius = -1),
                 "positive number")
    # RcppML's own graph arguments are not a supported passthrough.
    expect_error(runSpatialNMF(spe, k = 3, graph_H = diag(36)),
                 "builds graph_H")
})

test_that("spatialAutocorrelation reports Moran's I per program", {
    spe <- makeSpatialGrid()
    spe <- runSpatialNMF(spe, k = 4, graph_lambda = 1, seed = 7,
                         verbose = FALSE)

    result <- spatialAutocorrelation(spe)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 4L)
    expect_equal(colnames(result), c("program", "moran_I", "expected_I"))
    expect_setequal(result$program, paste0("NMF_", seq_len(4)))
    # Sorted most to least spatially structured.
    expect_equal(result$moran_I, sort(result$moran_I, decreasing = TRUE))
    expect_equal(unique(result$expected_I), -1 / (ncol(spe) - 1))
    expect_true(all(result$moran_I > 0.7))
})

test_that("spatialAutocorrelation gives a spatially random program I near 0", {
    spe <- makeSpatialGrid()
    spe <- runSpatialNMF(spe, k = 4, graph_lambda = 1, seed = 7,
                         verbose = FALSE)

    # Shuffling spot loadings destroys the spatial structure but not the
    # values, so Moran's I should collapse to its expectation.
    scrambled <- spe
    set.seed(99)
    shuffled <- reducedDim(spe, "SpatialNMF")[sample(ncol(spe)), , drop = FALSE]
    rownames(shuffled) <- colnames(spe)
    reducedDim(scrambled, "SpatialNMF") <- shuffled

    result <- spatialAutocorrelation(scrambled)
    expect_true(all(abs(result$moran_I) < 0.1))
})

test_that("spatialAutocorrelation permutation p-values work", {
    spe <- makeSpatialGrid()
    spe <- runSpatialNMF(spe, k = 4, graph_lambda = 1, seed = 7,
                         verbose = FALSE)

    result <- spatialAutocorrelation(spe, n_permutations = 49, seed = 1)

    expect_true("p_value" %in% colnames(result))
    expect_true(all(result$p_value >= 1 / 50 & result$p_value <= 1))
    # Strongly structured programs are significant at this resolution.
    expect_true(all(result$p_value < 0.05))

    # Same seed, same answer.
    repeated <- spatialAutocorrelation(spe, n_permutations = 49, seed = 1)
    expect_equal(result, repeated)
})

test_that("spatialAutocorrelation falls back to coordinates for plain NMF", {
    spe <- makeSpatialGrid()
    plain <- runNMFscape(spe, k = 4, seed = 7, verbose = FALSE)

    # No graph is stored by runNMFscape, so one is built from spatialCoords.
    expect_false("NMF_graph" %in% names(metadata(plain)))
    result <- spatialAutocorrelation(plain, name = "NMF")
    expect_equal(nrow(result), 4L)
    expect_true(all(is.finite(result$moran_I)))
})

test_that("spatialAutocorrelation validates its inputs", {
    spe <- makeSpatialGrid(side = 6, n_genes = 20)
    spe <- runSpatialNMF(spe, k = 3, seed = 7, verbose = FALSE)

    expect_error(spatialAutocorrelation(spe, name = "Missing"), "not found")
    expect_error(spatialAutocorrelation(spe, n_permutations = -1),
                 "non-negative")
    expect_error(spatialAutocorrelation(spe, graph = diag(5)),
                 "must be a 36 x 36 matrix")
    expect_error(spatialAutocorrelation(matrix(1:10, nrow = 2)),
                 "SingleCellExperiment")
})
