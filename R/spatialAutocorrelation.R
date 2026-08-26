#' Measure spatial autocorrelation of NMF programs
#'
#' Computes Moran's I for each program's spot loadings against a spatial
#' neighbour graph. Moran's I is near \eqn{-1/(n-1)} for a spatially random
#' (salt-and-pepper) program and approaches 1 for one that forms contiguous
#' domains, so it is the natural way to ask whether programs are spatially
#' structured and to choose \code{graph_lambda} for
#' \code{\link{runSpatialNMF}}: sweep lambda, and take the value past which
#' Moran's I stops rising appreciably but reconstruction error keeps climbing.
#'
#' By default the graph stored by \code{\link{runSpatialNMF}} is reused, so the
#' statistic is computed against exactly the graph the model was regularized
#' with. For results from \code{\link{runNMFscape}}, which store no graph, a
#' kNN graph is built from \code{spatialCoords(x)}. To compare a plain fit and
#' a spatial fit fairly, pass the same \code{graph} to both calls.
#'
#' @param x A SpatialExperiment or SingleCellExperiment with NMF results
#' @param name Character, name of the NMF result to score (default
#'   "SpatialNMF")
#' @param graph Optional square spot-by-spot adjacency matrix to score
#'   against. If NULL (default), uses the graph stored by
#'   \code{\link{runSpatialNMF}} under \code{name}, or builds a kNN graph from
#'   \code{spatialCoords(x)} when no graph was stored.
#' @param n_neighbors Integer, neighbours to use when a kNN graph has to be
#'   built (default 6)
#' @param n_permutations Integer, number of random permutations used to get a
#'   one-sided p-value for positive autocorrelation (default 0, no test)
#' @param seed Integer, random seed for the permutations
#'
#' @return A data.frame with one row per program and columns:
#'   \itemize{
#'     \item \code{program}: program (factor) name
#'     \item \code{moran_I}: observed Moran's I
#'     \item \code{expected_I}: expectation under spatial randomness,
#'       \eqn{-1/(n-1)}
#'   }
#'   plus \code{p_value} when \code{n_permutations > 0}. Rows are ordered by
#'   decreasing \code{moran_I}.
#'
#' @seealso \code{\link{runSpatialNMF}}
#'
#' @examples
#' library(SpatialExperiment)
#'
#' grid <- expand.grid(x = seq_len(10), y = seq_len(10))
#' domain <- ifelse(grid$x <= 5, 1, 2)
#' set.seed(42)
#' counts <- matrix(rpois(40 * nrow(grid), 1), nrow = 40)
#' counts[1:20, domain == 1] <- counts[1:20, domain == 1] + 5
#' counts[21:40, domain == 2] <- counts[21:40, domain == 2] + 5
#' rownames(counts) <- paste0("gene", seq_len(40))
#' colnames(counts) <- paste0("spot", seq_len(nrow(grid)))
#'
#' spe <- SpatialExperiment(
#'     assays = list(logcounts = log1p(counts)),
#'     spatialCoords = as.matrix(grid)
#' )
#' spe <- runSpatialNMF(spe, k = 2, seed = 1, verbose = FALSE)
#'
#' spatialAutocorrelation(spe)
#' spatialAutocorrelation(spe, n_permutations = 99, seed = 1)
#'
#' @export
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom stats median
spatialAutocorrelation <- function(x, name = "SpatialNMF", graph = NULL,
                                   n_neighbors = 6, n_permutations = 0,
                                   seed = NULL) {

    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment or SpatialExperiment object")
    }
    if (!name %in% reducedDimNames(x)) {
        stop("reducedDim '", name, "' not found. Run runSpatialNMF() first.")
    }
    if (length(n_permutations) != 1 || !is.numeric(n_permutations) ||
        is.na(n_permutations) || n_permutations < 0) {
        stop("n_permutations must be a single non-negative integer")
    }

    graph_key <- paste0(name, "_graph")
    if (!is.null(graph)) {
        a_mat <- .userAdjacency(graph, ncol(x), "graph")
    } else if (graph_key %in% names(metadata(x))) {
        a_mat <- metadata(x)[[graph_key]]
    } else {
        coords <- .spatialCoordsMatrix(x)
        a_mat <- .buildSpatialAdjacency(coords, "knn", n_neighbors, NULL,
                                        verbose = FALSE)
    }

    coeff <- as.matrix(reducedDim(x, name))
    n_spots <- nrow(coeff)
    s0 <- sum(a_mat)
    if (s0 == 0) {
        stop("The spatial graph has no edges; Moran's I is undefined")
    }

    observed <- vapply(seq_len(ncol(coeff)),
                       function(i) .moransI(coeff[, i], a_mat, s0),
                       FUN.VALUE = numeric(1))

    program_names <- colnames(coeff)
    if (is.null(program_names)) {
        program_names <- paste0("NMF_", seq_len(ncol(coeff)))
    }

    result <- data.frame(
        program = program_names,
        moran_I = observed,
        expected_I = -1 / (n_spots - 1),
        stringsAsFactors = FALSE
    )

    if (n_permutations > 0) {
        if (!is.null(seed)) {
            set.seed(seed)
        }
        n_perm <- as.integer(n_permutations)
        # Under the permutation null only the cross-product term changes, but
        # recomputing the whole statistic keeps this readable and it is a
        # single sparse multiply per draw.
        p_values <- vapply(seq_len(ncol(coeff)), function(i) {
            if (is.na(observed[i])) {
                return(NA_real_)
            }
            null_i <- vapply(seq_len(n_perm), function(b) {
                .moransI(coeff[sample.int(n_spots), i], a_mat, s0)
            }, FUN.VALUE = numeric(1))
            (1 + sum(null_i >= observed[i])) / (n_perm + 1)
        }, FUN.VALUE = numeric(1))
        result$p_value <- p_values
    }

    result <- result[order(result$moran_I, decreasing = TRUE), , drop = FALSE]
    rownames(result) <- NULL
    result
}
