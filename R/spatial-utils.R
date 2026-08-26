# Internal helpers for spatially-aware NMF
# @keywords internal
# @noRd

# Pull spatial coordinates off a SpatialExperiment, with a message that points
# at the "supply your own graph" escape hatch when they are missing.
.spatialCoordsMatrix <- function(x) {
    if (!is(x, "SpatialExperiment")) {
        stop("Building a spatial graph from coordinates requires a ",
             "SpatialExperiment. For a SingleCellExperiment, pass a ",
             "precomputed cell-by-cell adjacency matrix as 'graph' instead.")
    }
    coords <- SpatialExperiment::spatialCoords(x)
    if (is.null(coords) || nrow(coords) == 0 || ncol(coords) == 0) {
        stop("spatialCoords(x) is empty; cannot build a spatial graph. ",
             "Pass a precomputed adjacency matrix as 'graph' instead.")
    }
    coords <- as.matrix(coords)
    if (anyNA(coords)) {
        stop("spatialCoords(x) contains NA values")
    }
    storage.mode(coords) <- "double"
    coords
}

.asSparseAdjacency <- function(m) {
    m <- methods::as(methods::as(methods::as(m, "dMatrix"),
                                 "generalMatrix"), "CsparseMatrix")
    Matrix::diag(m) <- 0
    Matrix::drop0(m)
}

# kNN adjacency is directed; a spatial smoothness penalty needs a symmetric
# graph, so average the two directions (equivalent to an OR for 0/1 weights
# up to a factor of two, which is absorbed by graph_lambda).
.symmetrizeAdjacency <- function(a_mat) {
    .asSparseAdjacency((a_mat + Matrix::t(a_mat)) / 2)
}

# Radius that gives the median spot about n_neighbors neighbours, used when
# graph = "distance" is requested without an explicit radius.
.autoRadius <- function(coords, n_neighbors) {
    knn <- BiocNeighbors::findKNN(coords, k = n_neighbors, get.index = FALSE)
    stats::median(knn$distance[, n_neighbors])
}

.knnAdjacency <- function(coords, n_neighbors) {
    n_spots <- nrow(coords)
    knn <- BiocNeighbors::findKNN(coords, k = n_neighbors, get.distance = FALSE)
    Matrix::sparseMatrix(
        i = rep(seq_len(n_spots), times = n_neighbors),
        j = as.vector(knn$index),
        x = 1,
        dims = c(n_spots, n_spots)
    )
}

.distanceAdjacency <- function(coords, radius) {
    n_spots <- nrow(coords)
    nbr <- BiocNeighbors::findNeighbors(coords, threshold = radius,
                                        get.distance = FALSE)
    counts <- lengths(nbr$index)
    if (sum(counts) == 0) {
        stop("radius = ", signif(radius, 4), " yields no neighbours; ",
             "increase 'radius'")
    }
    Matrix::sparseMatrix(
        i = rep(seq_len(n_spots), times = counts),
        j = unlist(nbr$index, use.names = FALSE),
        x = 1,
        dims = c(n_spots, n_spots)
    )
}

.delaunayAdjacency <- function(coords) {
    if (!requireNamespace("deldir", quietly = TRUE)) {
        stop("graph = \"delaunay\" requires the 'deldir' package. ",
             "Install it, or use graph = \"knn\" or graph = \"distance\".")
    }
    if (ncol(coords) != 2) {
        stop("graph = \"delaunay\" requires 2D coordinates, but ",
             "spatialCoords(x) has ", ncol(coords), " columns")
    }
    tri <- deldir::deldir(coords[, 1], coords[, 2], suppressMsge = TRUE)
    edges <- tri$delsgs
    if (nrow(edges) == 0) {
        stop("Delaunay triangulation produced no edges")
    }
    Matrix::sparseMatrix(
        i = edges$ind1,
        j = edges$ind2,
        x = 1,
        dims = c(nrow(coords), nrow(coords))
    )
}

.buildSpatialAdjacency <- function(coords, graph_type, n_neighbors, radius,
                                   verbose = FALSE) {
    n_spots <- nrow(coords)

    if (graph_type %in% c("knn", "distance")) {
        if (length(n_neighbors) != 1 || !is.numeric(n_neighbors) ||
            is.na(n_neighbors) || n_neighbors < 1) {
            stop("n_neighbors must be a single positive integer")
        }
        if (n_neighbors >= n_spots) {
            stop("n_neighbors (", n_neighbors, ") must be smaller than the ",
                 "number of spots (", n_spots, ")")
        }
    }

    a_mat <- switch(
        graph_type,
        knn = .knnAdjacency(coords, n_neighbors),
        distance = {
            if (is.null(radius)) {
                radius <- .autoRadius(coords, n_neighbors)
                if (verbose) {
                    message("Using auto radius = ", signif(radius, 4),
                            " (median distance to neighbour ", n_neighbors, ")")
                }
            }
            if (length(radius) != 1 || !is.numeric(radius) || is.na(radius) ||
                radius <= 0) {
                stop("radius must be a single positive number")
            }
            .distanceAdjacency(coords, radius)
        },
        delaunay = .delaunayAdjacency(coords)
    )

    a_mat <- .symmetrizeAdjacency(a_mat)

    if (verbose) {
        message("Built ", graph_type, " graph: ",
                format(Matrix::nnzero(a_mat) / 2), " edges, mean degree ",
                signif(mean(Matrix::rowSums(a_mat)), 3))
    }
    a_mat
}

# Coerce and check a user-supplied adjacency matrix.
.userAdjacency <- function(g, n_expected, label) {
    if (is.data.frame(g)) {
        g <- as.matrix(g)
    }
    if (!(is.matrix(g) || is(g, "Matrix"))) {
        stop("'", label, "' must be a matrix or Matrix object")
    }
    if (nrow(g) != n_expected || ncol(g) != n_expected) {
        stop("'", label, "' must be a ", n_expected, " x ", n_expected,
             " matrix, but is ", nrow(g), " x ", ncol(g))
    }
    a_mat <- .asSparseAdjacency(g)
    if (anyNA(a_mat@x)) {
        stop("'", label, "' contains NA values")
    }
    if (any(a_mat@x < 0)) {
        stop("'", label, "' must have non-negative edge weights")
    }
    .symmetrizeAdjacency(a_mat)
}

# Split the Laplacian into the non-negative pieces L = D - S that the
# multiplicative updates need: S is an affinity matrix and D a diagonal.
# Combinatorial: S = A, D = degree. Symmetric normalized: S = D^(-1/2) A
# D^(-1/2), D = I. Isolated nodes get zero on both sides so they contribute
# no penalty rather than an undefined one.
.laplacianParts <- function(a_mat, normalize) {
    n_nodes <- nrow(a_mat)
    degree <- Matrix::rowSums(a_mat)
    connected <- degree > 0

    if (normalize) {
        inv_sqrt <- numeric(n_nodes)
        inv_sqrt[connected] <- 1 / sqrt(degree[connected])
        d_mat <- Matrix::Diagonal(n_nodes, x = inv_sqrt)
        affinity <- d_mat %*% a_mat %*% d_mat
        diagonal <- as.numeric(connected)
    } else {
        affinity <- a_mat
        diagonal <- degree
    }

    list(
        affinity = methods::as(methods::as(methods::as(affinity, "dMatrix"),
                                           "generalMatrix"), "CsparseMatrix"),
        diagonal = diagonal
    )
}

# tr(M L M^T) for the row-wise form used by H (k x n), computed from the
# split parts so no dense Laplacian is ever formed.
.graphPenalty <- function(m_mat, parts) {
    sum(m_mat * (m_mat %*% Matrix::Diagonal(x = parts$diagonal))) -
        sum(m_mat * (m_mat %*% parts$affinity))
}

# Graph-regularized NMF by multiplicative updates (Cai et al. 2011, extended
# to penalize both factors). Minimizes
#   ||A - W H||_F^2 + lambda_h tr(H L_h H^T) + lambda_w tr(W^T L_w W)
# Splitting each Laplacian as L = D - S with both parts non-negative keeps
# every term of the updates non-negative, so W and H never leave the feasible
# region and the objective is non-increasing.
#
# RcppML::nmf() supplies the initialization. Its own graph_W/graph_H/
# graph_lambda arguments are deliberately not used: in RcppML 1.0.0 they
# perturb the fit without making the factors smoother over the supplied graph
# (a randomly rewired graph changes the result as much as the true one), so
# they cannot back this feature.
#
# graph_lambda and feature_lambda arrive dimensionless. The raw penalty and
# the raw residual differ by many orders of magnitude and by dataset (the
# residual scales with genes x spots, the penalty with the size of H), so a
# fixed lambda would be a no-op on one dataset and catastrophic on the next.
# Each is therefore rescaled against the unregularized initialization, so that
# lambda = 1 makes the graph penalty one tenth of the initial reconstruction
# error. That anchor puts the default at the top of the accuracy curve on
# synthetic tissue, with useful values spanning roughly 0.01 to 10.
.gnmfRefine <- function(a_mat, w_mat, h_mat, h_parts, w_parts,
                        graph_lambda, feature_lambda, tol, maxit,
                        verbose = FALSE) {
    eps <- .Machine$double.eps

    # W H is invariant under W -> Wc, H -> H/c, but tr(H L H^T) scales as
    # 1/c^2, so an unconstrained optimizer drives the penalty to zero by
    # inflating W and shrinking H rather than by making H smooth. Pinning the
    # columns of W to unit norm after every update removes that escape route.
    normalized <- .unitColumns(w_mat, h_mat)
    w_mat <- normalized$w
    h_mat <- normalized$h

    residual0 <- sum((a_mat - w_mat %*% h_mat)^2)
    anchor <- residual0 / 10
    if (!is.null(h_parts)) {
        graph_lambda <- graph_lambda * anchor /
            max(.graphPenalty(h_mat, h_parts), eps)
    }
    if (!is.null(w_parts)) {
        feature_lambda <- feature_lambda * anchor /
            max(.graphPenalty(t(w_mat), w_parts), eps)
    }

    objective <- function(w, h) {
        value <- sum((a_mat - w %*% h)^2)
        if (!is.null(h_parts)) {
            value <- value + graph_lambda * .graphPenalty(h, h_parts)
        }
        if (!is.null(w_parts)) {
            value <- value + feature_lambda * .graphPenalty(t(w), w_parts)
        }
        value
    }

    previous <- objective(w_mat, h_mat)
    iter <- 0L

    for (iter in seq_len(maxit)) {
        # H update: H * (W'A + lambda H S) / (W'W H + lambda H D)
        numer <- as.matrix(crossprod(w_mat, a_mat))
        denom <- crossprod(w_mat) %*% h_mat
        if (!is.null(h_parts)) {
            numer <- numer +
                graph_lambda * as.matrix(h_mat %*% h_parts$affinity)
            denom <- denom +
                graph_lambda * sweep(h_mat, 2, h_parts$diagonal, "*")
        }
        h_mat <- h_mat * (numer / (denom + eps))

        # W update: W * (A H' + lambda S W) / (W H H' + lambda D W)
        numer <- as.matrix(tcrossprod(a_mat, h_mat))
        denom <- w_mat %*% tcrossprod(h_mat)
        if (!is.null(w_parts)) {
            numer <- numer +
                feature_lambda * as.matrix(w_parts$affinity %*% w_mat)
            denom <- denom +
                feature_lambda * sweep(w_mat, 1, w_parts$diagonal, "*")
        }
        w_mat <- w_mat * (numer / (denom + eps))

        normalized <- .unitColumns(w_mat, h_mat)
        w_mat <- normalized$w
        h_mat <- normalized$h

        current <- objective(w_mat, h_mat)
        change <- abs(previous - current) / (abs(previous) + eps)
        previous <- current
        if (change < tol) {
            break
        }
    }

    if (verbose) {
        message("Graph-regularized refinement: ", iter,
                " iterations, objective ", signif(previous, 6))
    }

    list(w = w_mat, h = h_mat, iter = iter, objective = previous,
         scaled_graph_lambda = graph_lambda,
         scaled_feature_lambda = feature_lambda)
}

# Push all scale into H so that the columns of W have unit L2 norm, leaving
# the product W H unchanged.
.unitColumns <- function(w_mat, h_mat) {
    col_norm <- sqrt(colSums(w_mat^2))
    keep <- col_norm > 0
    w_mat[, keep] <- sweep(w_mat[, keep, drop = FALSE], 2, col_norm[keep], "/")
    h_mat[keep, ] <- sweep(h_mat[keep, , drop = FALSE], 1, col_norm[keep], "*")
    list(w = w_mat, h = h_mat)
}

# Restore the RcppML convention A ~ W diag(d) H, with columns of W and rows of
# H summing to one.
.rescaleFactors <- function(w_mat, h_mat) {
    w_scale <- colSums(w_mat)
    h_scale <- rowSums(h_mat)
    d_vec <- w_scale * h_scale
    keep <- w_scale > 0 & h_scale > 0
    w_mat[, keep] <- sweep(w_mat[, keep, drop = FALSE], 2, w_scale[keep], "/")
    h_mat[keep, ] <- sweep(h_mat[keep, , drop = FALSE], 1, h_scale[keep], "/")
    list(w = w_mat, d = d_vec, h = h_mat)
}

.checkLambda <- function(value, label) {
    if (length(value) != 1 || !is.numeric(value) || is.na(value) ||
        value < 0) {
        stop("'", label, "' must be a single non-negative number")
    }
    invisible(value)
}

# Moran's I for one vector against a sparse weight matrix.
.moransI <- function(v, a_mat, s0) {
    z <- v - mean(v)
    den <- sum(z^2)
    if (den == 0) {
        return(NA_real_)
    }
    (length(z) / s0) * sum(z * as.numeric(a_mat %*% z)) / den
}
