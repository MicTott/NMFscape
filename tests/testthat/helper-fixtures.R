# Shared fixtures. testthat sources helper-*.R before the test files.

# A small SingleCellExperiment with log-normalized counts and a grouping.
makeSCE <- function(ngenes = 150, ncells = 120, seed = 1) {
    set.seed(seed)
    sce <- scuttle::logNormCounts(scuttle::mockSCE(ngenes = ngenes,
                                                   ncells = ncells))
    sce$celltype <- rep(c("A", "B", "C"), length.out = ncol(sce))
    sce$batch <- rep(c("b1", "b2"), length.out = ncol(sce))
    sce
}

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

# A noisy k-program simulation where the unguided fit is good but imperfect,
# leaving room for guidance to help or hurt measurably.
.testGuidedSCE <- function(seed = 1, n_genes = 200, n_cells = 400,
                           n_prog = 5, noise = 5, overlap = 0.5) {
    set.seed(seed)
    genes_per <- n_genes %/% n_prog
    group <- rep(seq_len(n_prog), length.out = n_cells)

    w_mat <- matrix(runif(n_genes * n_prog, 0, overlap), n_genes, n_prog)
    for (i in seq_len(n_prog)) {
        idx <- seq((i - 1) * genes_per + 1, i * genes_per)
        w_mat[idx, i] <- w_mat[idx, i] + runif(genes_per, 1, 3)
    }
    h_mat <- matrix(runif(n_prog * n_cells, 0, overlap), n_prog, n_cells)
    for (j in seq_len(n_cells)) {
        h_mat[group[j], j] <- h_mat[group[j], j] + runif(1, 1, 2)
    }
    mat <- w_mat %*% h_mat +
        matrix(abs(rnorm(n_genes * n_cells, 0, noise)), n_genes, n_cells)

    rownames(mat) <- sprintf("Gene_%03d", seq_len(n_genes))
    colnames(mat) <- sprintf("Cell_%03d", seq_len(n_cells))

    sce <- SingleCellExperiment::SingleCellExperiment(list(logcounts = mat))
    sce$program <- paste0("P", group)
    sce
}

.argmaxLabels <- function(sce, name) {
    apply(reducedDim(sce, name), 1, which.max)
}

makeTuneSCE <- function(ngenes = 60, ncells = 40, seed = 42) {
    set.seed(seed)
    sce <- scuttle::mockSCE(ngenes = ngenes, ncells = ncells)
    scuttle::logNormCounts(sce)
}

makeMultiSCE <- function() {
    sce <- makeTuneSCE()
    set.seed(7)
    atac <- SingleCellExperiment(list(
        counts = abs(matrix(rnorm(20 * ncol(sce)), 20, ncol(sce)))
    ))
    rownames(atac) <- paste0("peak_", seq_len(20))
    colnames(atac) <- colnames(sce)
    altExp(sce, "ATAC") <- atac
    sce
}
