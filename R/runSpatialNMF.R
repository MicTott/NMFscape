#' Run spatially-aware (graph-regularized) NMF
#'
#' Performs NMF with a graph-Laplacian smoothness penalty on the cell/spot
#' factor loadings, so that spots which neighbour each other in tissue are
#' encouraged to use similar programs. This adds
#' \eqn{\lambda \cdot tr(H L H^T)} to the NMF objective, where \eqn{L} is the
#' Laplacian of a spatial neighbour graph. Because
#' \eqn{tr(H L H^T) = \frac{1}{2}\sum_{ij} A_{ij} \|h_i - h_j\|^2}, the penalty
#' is simply the total squared difference in program usage across graph edges.
#'
#' The constraint is soft: \code{graph_lambda = 0} recovers plain
#' \code{\link{runNMFscape}}, small values remove salt-and-pepper noise while
#' preserving boundaries, and large values oversmooth until programs blur into
#' each other. Use \code{\link{spatialAutocorrelation}} to sweep
#' \code{graph_lambda} and see the bias-variance trade-off directly: Moran's I
#' climbs monotonically with lambda, while agreement with the true domains
#' peaks and then falls away.
#'
#' @details
#' The factorization is initialized with \code{\link[RcppML]{nmf}} and then
#' refined with graph-regularized multiplicative updates (Cai et al. 2011),
#' which minimize
#' \eqn{\|A - WH\|_F^2 + \lambda_h tr(H L_h H^T) + \lambda_w tr(W^T L_w W)}.
#' Splitting each Laplacian into its non-negative parts \eqn{L = D - S} keeps
#' every factor non-negative and the objective non-increasing. The columns of
#' W are held at unit norm throughout: \eqn{WH} is invariant under
#' \eqn{W \to Wc, H \to H/c} while the penalty scales as \eqn{1/c^2}, so
#' without that constraint the optimizer would zero out the penalty by
#' inflating W instead of by making H spatially smooth.
#'
#' RcppML 1.0.0 exposes its own \code{graph_W}, \code{graph_H} and
#' \code{graph_lambda} arguments, but they are not used here and are rejected
#' if passed through \code{...}. In testing they change the fit without making
#' the factors any smoother over the supplied graph: a randomly rewired graph
#' with the same degree sequence perturbs the result as much as the true one,
#' and relative smoothness over the graph's own edges gets worse, not better,
#' as \code{graph_lambda} rises. What looks like smoothing there is uniform
#' shrinkage of H.
#'
#' @section Supplying your own graph:
#' Set \code{graph} to a square spot-by-spot adjacency matrix (dense or
#' \pkg{Matrix} sparse) to use a precomputed neighbour graph instead of
#' building one from coordinates. This is the way to run spatially-aware NMF on
#' a plain \code{SingleCellExperiment}, or on any object whose neighbour
#' structure is not Euclidean distance in \code{spatialCoords()} --- for
#' example a graph built across tissue sections, a histology-aware graph, or a
#' kNN graph in expression space. The matrix is symmetrized and converted to a
#' Laplacian exactly as a coordinate-derived graph would be, so it should hold
#' non-negative edge weights, not a Laplacian. Character values of \code{graph}
#' require a \code{SpatialExperiment}.
#'
#' @param x A SpatialExperiment object, or a SingleCellExperiment when
#'   \code{graph} is a precomputed adjacency matrix
#' @param k Integer, number of factors for NMF (rank)
#' @param graph Either a character string naming how to build the neighbour
#'   graph from \code{spatialCoords(x)} --- \code{"knn"} (default),
#'   \code{"delaunay"} or \code{"distance"} --- or a square spot-by-spot
#'   adjacency matrix of non-negative weights. \code{"delaunay"} requires the
#'   \pkg{deldir} package and 2D coordinates.
#' @param n_neighbors Integer, number of nearest neighbours for
#'   \code{graph = "knn"} (default 6, the Visium lattice degree). Also sets the
#'   auto-selected radius for \code{graph = "distance"}.
#' @param radius Numeric, neighbourhood radius for \code{graph = "distance"}.
#'   If NULL (default), uses the median distance to the \code{n_neighbors}-th
#'   nearest neighbour, so the median spot gets about \code{n_neighbors}
#'   neighbours.
#' @param graph_lambda Numeric, strength of the spatial smoothness penalty on
#'   the cell/spot factors (default 1). Zero disables it, giving plain NMF.
#'   The value is dimensionless: it is rescaled internally against the
#'   unregularized fit, so 1 means the graph penalty is one tenth of the
#'   initial reconstruction error regardless of dataset size or scale. Useful
#'   values span roughly 0.01 (barely any smoothing) to 10 (heavy).
#' @param feature_graph Optional square gene-by-gene adjacency matrix (a
#'   co-expression or protein-protein interaction prior) used to smooth the
#'   gene loadings. Must match the features actually factorized, i.e. after
#'   \code{subset_row}.
#' @param feature_lambda Numeric, strength of the penalty applied to
#'   \code{feature_graph} (default 0), on the same dimensionless scale as
#'   \code{graph_lambda}
#' @param assay Character or integer, which assay to use (default "logcounts")
#' @param name Character, name for the reducedDim slot (default "SpatialNMF")
#' @param subset_row Vector specifying which features to use
#' @param normalize_laplacian Logical, whether to use the symmetric normalized
#'   Laplacian \eqn{I - D^{-1/2} A D^{-1/2}} rather than the combinatorial
#'   \eqn{D - A} (default TRUE). Normalization keeps low-degree spots at tissue
#'   edges from being penalized less than interior spots.
#' @param tol Numeric, relative convergence tolerance, used both for the
#'   \code{\link[RcppML]{nmf}} initialization and for the graph-regularized
#'   updates (default 1e-5)
#' @param maxit Integer, maximum iterations for each of those two stages
#'   (default 100)
#' @param absorb_d Logical, whether to absorb the diagonal scaling factor into
#'   W and H symmetrically (default TRUE), as in \code{\link{runNMFscape}}
#' @param seed Integer, random seed for reproducibility
#' @param verbose Logical, whether to print progress (default TRUE)
#' @param ... Additional arguments passed to \code{\link[RcppML]{nmf}} for the
#'   initialization, including \code{L1}, \code{L2} and \code{loss}. RcppML's
#'   own \code{graph_W}, \code{graph_H} and \code{graph_lambda} arguments are
#'   rejected; see Details.
#'
#' @return The input object with spatially-aware NMF results stored in:
#'   \itemize{
#'     \item \code{reducedDim(x, name)}: spots x k coefficient matrix
#'     \item \code{metadata(x)[[paste0(name, "_basis")]]}: genes x k basis
#'     \item \code{metadata(x)[[paste0(name, "_model")]]}: raw RcppML S4 nmf
#'     \item \code{metadata(x)[[paste0(name, "_graph")]]}: the symmetrized
#'       sparse adjacency matrix used to build the Laplacian
#'     \item \code{metadata(x)[[paste0(name, "_graph_params")]]}: list recording
#'       the graph type and penalty settings
#'   }
#'
#' @seealso \code{\link{spatialAutocorrelation}} to measure how spatially
#'   structured the resulting programs are, and \code{\link{runNMFscape}} for
#'   the unregularized version.
#'
#' @examples
#' library(SpatialExperiment)
#'
#' # A small grid of spots with two contiguous domains
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
#'
#' spe <- runSpatialNMF(spe, k = 2, graph_lambda = 1, seed = 1, verbose = FALSE)
#' dim(reducedDim(spe, "SpatialNMF"))
#' spatialAutocorrelation(spe)
#'
#' @export
#' @importFrom RcppML nmf
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom BiocNeighbors findKNN findNeighbors
#' @importFrom Matrix sparseMatrix Diagonal rowSums drop0 nnzero
runSpatialNMF <- function(x, k, graph = c("knn", "delaunay", "distance"),
                          n_neighbors = 6, radius = NULL,
                          graph_lambda = 1,
                          feature_graph = NULL, feature_lambda = 0,
                          assay = "logcounts", name = "SpatialNMF",
                          subset_row = NULL, normalize_laplacian = TRUE,
                          tol = 1e-5, maxit = 100,
                          absorb_d = TRUE, seed = NULL, verbose = TRUE, ...) {

    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment or SpatialExperiment object")
    }
    if (!assay %in% assayNames(x)) {
        stop("assay '", assay, "' not found in x")
    }

    .checkLambda(graph_lambda, "graph_lambda")
    .checkLambda(feature_lambda, "feature_lambda")

    dot_names <- names(list(...))
    reserved <- intersect(dot_names, c("graph_W", "graph_H", "graph_lambda"))
    if (length(reserved) > 0) {
        stop("runSpatialNMF() builds ", paste(reserved, collapse = ", "),
             " itself; use the 'graph', 'feature_graph', 'graph_lambda' and ",
             "'feature_lambda' arguments instead")
    }

    feature_names <- .subsetRowNames(x, subset_row)
    mat <- assay(x, assay)
    if (!is.null(subset_row)) {
        mat <- mat[subset_row, , drop = FALSE]
    }
    if (any(mat < 0)) {
        warning("Negative values detected. ",
                "Consider using log-transformed data.")
        mat[mat < 0] <- 0
    }

    # Resolve the spot graph: either build it from coordinates or take the
    # user's own adjacency matrix.
    if (is.character(graph)) {
        graph_type <- match.arg(graph)
        coords <- .spatialCoordsMatrix(x)
        if (nrow(coords) != ncol(x)) {
            stop("spatialCoords(x) has ", nrow(coords), " rows but x has ",
                 ncol(x), " columns")
        }
        adjacency <- .buildSpatialAdjacency(coords, graph_type, n_neighbors,
                                            radius, verbose)
    } else {
        graph_type <- "custom"
        adjacency <- .userAdjacency(graph, ncol(x), "graph")
        if (verbose) {
            message("Using supplied spot graph: ",
                    format(Matrix::nnzero(adjacency) / 2), " edges")
        }
    }

    h_parts <- NULL
    if (graph_lambda > 0) {
        h_parts <- .laplacianParts(adjacency, normalize_laplacian)
    }

    w_parts <- NULL
    feature_adjacency <- NULL
    if (!is.null(feature_graph)) {
        feature_adjacency <- .userAdjacency(feature_graph, nrow(mat),
                                            "feature_graph")
        if (feature_lambda > 0) {
            w_parts <- .laplacianParts(feature_adjacency, normalize_laplacian)
        } else {
            warning("feature_graph supplied but feature_lambda is 0; ",
                    "the gene network will not be applied")
        }
    }

    if (verbose) {
        if (is.null(h_parts) && is.null(w_parts)) {
            message("Running NMF with k=", k,
                    " factors (graph_lambda=0, no graph penalty)...")
        } else {
            message("Running spatially-aware NMF with k=", k,
                    " factors (graph_lambda=", graph_lambda, ")...")
        }
    }

    nmf_result <- RcppML::nmf(data = mat, k = k, tol = tol, maxit = maxit,
                              seed = seed, verbose = FALSE, ...)

    refinement <- NULL
    if (!is.null(h_parts) || !is.null(w_parts)) {
        # Start the graph-regularized updates from the unregularized fit,
        # splitting the diagonal symmetrically so that A ~ W H with both
        # factors on comparable scales.
        sqrt_d <- sqrt(nmf_result@d)
        refinement <- .gnmfRefine(
            a_mat = as.matrix(mat),
            w_mat = nmf_result@w %*% diag(sqrt_d, nrow = k),
            h_mat = diag(sqrt_d, nrow = k) %*% nmf_result@h,
            h_parts = h_parts, w_parts = w_parts,
            graph_lambda = graph_lambda, feature_lambda = feature_lambda,
            tol = tol, maxit = maxit, verbose = verbose
        )
        scaled <- .rescaleFactors(refinement$w, refinement$h)
        nmf_result@w <- scaled$w
        nmf_result@d <- scaled$d
        nmf_result@h <- scaled$h
        nmf_result@misc$graph_iter <- refinement$iter
        nmf_result@misc$graph_objective <- refinement$objective
        nmf_result@misc$scaled_graph_lambda <- refinement$scaled_graph_lambda
        nmf_result@misc$scaled_feature_lambda <-
            refinement$scaled_feature_lambda
    }

    parts <- .absorbDiagonal(nmf_result@w, nmf_result@h, nmf_result@d,
                             absorb_d)
    parts <- .setFactorNames(parts$basis, parts$coeff, feature_names,
                             colnames(x), k)

    reducedDim(x, name) <- parts$coeff
    metadata(x)[[paste0(name, "_basis")]] <- parts$basis
    metadata(x)[[paste0(name, "_model")]] <- nmf_result
    metadata(x)[[paste0(name, "_graph")]] <- adjacency
    metadata(x)[[paste0(name, "_graph_params")]] <- list(
        graph = graph_type,
        n_neighbors = if (graph_type %in% c("knn", "distance")) {
            n_neighbors
        } else {
            NA_integer_
        },
        radius = if (identical(graph_type, "distance")) radius else NULL,
        graph_lambda = graph_lambda,
        feature_lambda = feature_lambda,
        feature_graph = !is.null(feature_adjacency),
        normalize_laplacian = normalize_laplacian,
        seed = seed
    )

    if (verbose) {
        message("Spatial NMF completed. Results stored in reducedDim(x, '",
                name, "')")
    }

    return(x)
}
