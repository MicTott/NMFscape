#' Run deep (hierarchical) NMF on a SingleCellExperiment
#'
#' Performs multi-layer hierarchical NMF where each layer further compresses
#' the representation. Layer 1 finds fine-grained gene programs, subsequent
#' layers group them into meta-programs. Uses the RcppML FactorNet graph DSL.
#'
#' Layer 1 factorizes the assay matrix; each subsequent layer factorizes the
#' transpose of the previous layer's H, so the orientation alternates and cells
#' appear in the even layers' W while features appear in the odd layers' W.
#' The effective basis (genes x k_final) and cell embedding (cells x k_final)
#' are formed by propagating through all layers, so both are expressed in the
#' deepest layer's meta-program space regardless of depth. Per-layer W, H and d
#' are stored separately for exploration.
#'
#' @param x A SingleCellExperiment or SpatialExperiment object
#' @param k Integer vector of length >= 2, e.g. \code{c(20, 10)}. Each element
#'   is the rank for that layer (first layer operates on features, subsequent
#'   layers compress the previous layer's output).
#' @param assay Character, which assay to use (default "logcounts")
#' @param name Character, name for the reducedDim slot (default "DeepNMF")
#' @param subset_row Vector specifying which features to use
#' @param tol Numeric, tolerance for convergence (default 1e-5)
#' @param maxit Integer, maximum iterations (default 100)
#' @param L1 Numeric scalar, L1 regularization applied to each NMF layer
#'   (default 0). Note this differs from \code{\link{runNMFscape}}, which takes a
#'   length-2 c(w, h) vector; \code{\link[RcppML]{nmf_layer}} uses one penalty
#'   per layer.
#' @param L2 Numeric scalar, L2 regularization applied to each NMF layer
#'   (default 0)
#' @param distribution Character, loss function (default "mse")
#' @param absorb_d Logical, whether to absorb diagonal scaling (default TRUE)
#' @param seed Integer, random seed for reproducibility
#' @param verbose Logical, whether to print progress (default TRUE)
#' @param ... Additional arguments passed to \code{\link[RcppML]{factor_config}}
#'
#' @return The input object with deep NMF results stored in:
#'   \itemize{
#'     \item \code{reducedDim(x, name)}: cells x k[last] coefficient matrix,
#'       propagated through every layer.
#'     \item \code{metadata(x)[[paste0(name, "_basis")]]}: genes x k[last]
#'       effective basis matrix (product through all layers), mapping original
#'       features directly to meta-programs
#'     \item \code{metadata(x)[[paste0(name, "_layers")]]}: named list with
#'       per-layer W, H, d for exploration
#'     \item \code{metadata(x)[[paste0(name, "_model")]]}: raw FactorNet result
#'   }
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#' sce <- runDeepNMF(sce, k = c(15, 5))
#'
#' # Effective basis maps genes to 5 meta-programs
#' dim(getBasis(sce, "DeepNMF"))
#'
#' # Per-layer details
#' layers <- metadata(sce)$DeepNMF_layers
#' dim(layers$layer_1$W)  # genes x 15
#' dim(layers$layer_2$W)  # cells x 5
#'
#' # Three or more layers are supported; the embedding is always k[last]
#' sce <- runDeepNMF(sce, k = c(20, 10, 5), name = "DeepNMF3")
#' dim(reducedDim(sce, "DeepNMF3"))
#'
#' @export
#' @importFrom RcppML factor_input nmf_layer factor_net factor_config fit
runDeepNMF <- function(x, k, assay = "logcounts", name = "DeepNMF",
                       subset_row = NULL, tol = 1e-5, maxit = 100,
                       L1 = 0, L2 = 0,
                       distribution = "mse", absorb_d = TRUE,
                       seed = NULL, verbose = TRUE, ...) {

    .validateSCE(x, assay)

    .checkLayerPenalty(L1, "L1")
    .checkLayerPenalty(L2, "L2")

    k <- as.integer(k)
    n_layers <- length(k)
    if (n_layers < 2) {
        stop("k must have length >= 2 for deep NMF. ",
             "Use runNMFscape() for single-layer NMF.")
    }

    if (!all(diff(k) < 0)) {
        warning("k values are not strictly decreasing; typically k should ",
                "decrease with depth for hierarchical factorization")
    }

    feature_names <- .subsetRowNames(x, subset_row)
    mat <- .extractAssayMatrix(x, assay, subset_row)

    if (verbose) {
        message("Running deep NMF with ", n_layers, " layers (k=",
                paste(k, collapse = " -> "), ")...")
    }

    # Build FactorNet graph: chain of nmf_layers
    input <- RcppML::factor_input(mat, name = "input")
    current <- input
    for (i in seq_along(k)) {
        layer_name <- paste0("layer_", i)
        current <- RcppML::nmf_layer(current, k = k[i], L1 = L1, L2 = L2,
                                     name = layer_name)
    }

    config <- RcppML::factor_config(
        tol = tol, maxit = maxit, loss = distribution,
        seed = seed, verbose = FALSE, ...
    )
    net <- RcppML::factor_net(inputs = list(input), output = current,
                              config = config)
    result <- RcppML::fit(net)

    # Extract per-layer results
    layer_data <- lapply(seq_along(k), function(i) {
        ln <- paste0("layer_", i)
        l <- result$layers[[ln]]
        list(W = l$W, H = l$H, d = l$d)
    })
    names(layer_data) <- paste0("layer_", seq_along(k))

    last_layer <- result$layers[[paste0("layer_", n_layers)]]

    # Compute the effective decomposition A ~ Basis %*% diag(d_n) %*% t(Coeff).
    #
    # Layer 1 factorizes A; layer i > 1 factorizes t(H_{i-1}), so the
    # orientation alternates: odd layers carry the feature/program side,
    # even layers carry the cell side. Unrolling
    #   A         ~ W1 D1 H1
    #   H_{i-1}   ~ t(H_i) D_i t(W_i)
    # gives an alternating product, e.g. for n = 3
    #   A ~ (W1 D1 W3) D3 (H3 D2 t(W2))
    # so the feature side is the ascending product over odd layers and the
    # cell side the ascending product over even layers, with the deepest
    # layer's H joining whichever side is one step short.
    .chainProduct <- function(idx, hold_last_d) {
        acc <- NULL
        for (j in seq_along(idx)) {
            li <- result$layers[[paste0("layer_", idx[j])]]
            term <- if (hold_last_d && j == length(idx)) {
                li$W
            } else {
                li$W %*% diag(li$d, nrow = length(li$d))
            }
            acc <- if (is.null(acc)) term else acc %*% term
        }
        acc
    }

    odd <- seq(1, n_layers, by = 2)
    even <- seq(2, n_layers, by = 2)
    n_is_even <- n_layers %% 2 == 0

    # The deepest layer's d is held back so it can be split symmetrically.
    feature_side <- .chainProduct(odd, hold_last_d = !n_is_even)
    cell_side <- .chainProduct(even, hold_last_d = n_is_even)

    if (n_is_even) {
        eff_basis <- feature_side %*% t(last_layer$H)
        coeff_matrix <- cell_side
    } else {
        eff_basis <- feature_side
        coeff_matrix <- cell_side %*% t(last_layer$H)
    }

    d_last <- last_layer$d
    k_last <- length(d_last)
    if (absorb_d) {
        sqrt_d <- sqrt(d_last)
        eff_basis <- eff_basis %*% diag(sqrt_d, nrow = k_last)
        coeff_matrix <- coeff_matrix %*% diag(sqrt_d, nrow = k_last)
    } else {
        eff_basis <- eff_basis %*% diag(d_last, nrow = k_last)
    }

    k_final <- k[n_layers]
    named <- .setFactorNames(eff_basis, coeff_matrix,
                             feature_names, colnames(x), k_final)

    # Store
    reducedDim(x, name) <- named$coeff
    metadata(x)[[paste0(name, "_basis")]] <- named$basis
    metadata(x)[[paste0(name, "_layers")]] <- layer_data
    metadata(x)[[paste0(name, "_layer")]] <- paste0("layer_", n_layers)
    metadata(x)[[paste0(name, "_model")]] <- result

    if (verbose) {
        message("Deep NMF completed (", n_layers, " layers). ",
                "Results stored in reducedDim(x, '", name, "')")
    }

    return(x)
}
