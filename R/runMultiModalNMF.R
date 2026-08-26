#' Run multi-modal joint NMF on a SingleCellExperiment
#'
#' Performs joint NMF factorization across multiple modalities from the same
#' cells. Each modality gets its own feature loadings (W) while sharing a
#' common cell embedding (H). Uses \code{\link[RcppML]{factor_shared}}
#' internally.
#'
#' Modalities can be specified as either:
#' \itemize{
#'   \item Multiple assays in the main SCE (same feature space, e.g., different
#'     normalization layers) via the \code{assays} parameter
#'   \item Alternative experiments (different feature spaces, e.g., RNA + ATAC)
#'     via the \code{alt_exps} parameter, where each altExp is a separate SCE
#'     stored via \code{\link[SingleCellExperiment]{altExp}}
#' }
#'
#' @param x A SingleCellExperiment or SpatialExperiment object
#' @param k Integer, number of shared factors
#' @param assays Character vector of >= 2 assay names from the main SCE.
#'   Use this when modalities share the same feature space. Mutually exclusive
#'   with \code{alt_exps}.
#' @param alt_exps Character vector of altExp names plus optionally the main
#'   experiment. E.g., \code{c("main", "ATAC")} uses the main SCE's assay
#'   as one modality and \code{altExp(x, "ATAC")} as another.
#' @param main_assay Character, assay to use from the main SCE when using
#'   \code{alt_exps} (default "logcounts")
#' @param alt_assay Character, assay to use from each altExp (default "counts")
#' @param modality_names Character vector of short names for each modality.
#'   Defaults to assay or altExp names.
#' @param name Character, name for the reducedDim slot (default "MultiNMF")
#' @param subset_rows Named list of per-modality feature subsets
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
#' @return The input object with multi-modal NMF results stored in:
#'   \itemize{
#'     \item \code{reducedDim(x, name)}: cells x k shared coefficient matrix
#'     \item \code{metadata(x)[[paste0(name, "_basis_", modality)]]}: per-modality
#'       basis matrices (features_i x k)
#'     \item \code{metadata(x)[[paste0(name, "_modalities")]]}: character vector
#'       of modality names
#'     \item \code{metadata(x)[[paste0(name, "_model")]]}: raw FactorNet result
#'   }
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#'
#' # Create altExp for second modality
#' atac <- SingleCellExperiment(list(
#'     counts = abs(matrix(rnorm(50 * ncol(sce)), 50, ncol(sce)))
#' ))
#' rownames(atac) <- paste0("peak_", seq_len(50))
#' colnames(atac) <- colnames(sce)
#' altExp(sce, "ATAC") <- atac
#'
#' sce <- runMultiModalNMF(sce, k = 5,
#'                         alt_exps = c("main", "ATAC"),
#'                         modality_names = c("RNA", "ATAC"))
#' getBasis(sce, "MultiNMF", modality = "RNA")
#'
#' @export
#' @importFrom RcppML factor_input factor_shared nmf_layer factor_net
#'   factor_config fit
#' @importFrom SingleCellExperiment altExp altExpNames
runMultiModalNMF <- function(x, k, assays = NULL, alt_exps = NULL,
                             main_assay = "logcounts", alt_assay = "counts",
                             modality_names = NULL,
                             name = "MultiNMF", subset_rows = NULL,
                             tol = 1e-5, maxit = 100,
                             L1 = 0, L2 = 0,
                             distribution = "mse", absorb_d = TRUE,
                             seed = NULL, verbose = TRUE, ...) {

    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment or SpatialExperiment object")
    }

    .checkLayerPenalty(L1, "L1")
    .checkLayerPenalty(L2, "L2")

    if (is.null(assays) && is.null(alt_exps)) {
        stop("Specify either 'assays' or 'alt_exps'")
    }
    if (!is.null(assays) && !is.null(alt_exps)) {
        stop("Specify either 'assays' or 'alt_exps', not both")
    }

    # Extract matrices from either assays or altExps
    if (!is.null(assays)) {
        if (length(assays) < 2) {
            stop("runMultiModalNMF requires at least 2 modalities. ",
                 "Use runNMFscape() for single-modality NMF.")
        }
        if (is.null(modality_names)) modality_names <- assays

        mats <- lapply(seq_along(assays), function(i) {
            sub <- if (!is.null(subset_rows)) {
                subset_rows[[modality_names[i]]]
            } else {
                NULL
            }
            .extractAssayMatrix(x, assays[i], sub)
        })
        names(mats) <- modality_names

    } else {
        if (length(alt_exps) < 2) {
            stop("runMultiModalNMF requires at least 2 modalities. ",
                 "Use runNMFscape() for single-modality NMF.")
        }
        if (is.null(modality_names)) modality_names <- alt_exps

        mats <- lapply(seq_along(alt_exps), function(i) {
            ae <- alt_exps[i]
            mod <- modality_names[i]
            sub <- if (!is.null(subset_rows)) subset_rows[[mod]] else NULL

            if (ae == "main") {
                .extractAssayMatrix(x, main_assay, sub)
            } else {
                if (!ae %in% altExpNames(x)) {
                    stop("altExp '", ae, "' not found in x")
                }
                ae_sce <- altExp(x, ae)
                mat <- assay(ae_sce, alt_assay)
                if (!is.null(sub)) mat <- mat[sub, , drop = FALSE]
                if (any(mat < 0)) {
                    warning("Negative values in altExp '", ae,
                            "'. Setting to 0.")
                    mat[mat < 0] <- 0
                }
                mat
            }
        })
        names(mats) <- modality_names
    }

    if (length(modality_names) != length(mats)) {
        stop("modality_names must have the same length as the number of modalities")
    }

    # Validate matching cell counts
    n_cells <- vapply(mats, ncol, integer(1))
    if (length(unique(n_cells)) != 1) {
        stop("All modalities must have the same number of cells. Found: ",
             paste(paste0(modality_names, "=", n_cells), collapse = ", "))
    }

    if (verbose) {
        message("Running multi-modal NMF with k=", k, " across ",
                length(mats), " modalities (",
                paste(modality_names, collapse = ", "), ")...")
    }

    # Build FactorNet graph
    inputs <- lapply(seq_along(mats), function(i) {
        RcppML::factor_input(mats[[i]], name = modality_names[i])
    })
    shared_node <- do.call(RcppML::factor_shared, inputs)
    output <- RcppML::nmf_layer(shared_node, k = k, L1 = L1, L2 = L2,
                                name = "shared_nmf")

    config <- RcppML::factor_config(
        tol = tol, maxit = maxit, loss = distribution,
        seed = seed, verbose = FALSE, ...
    )
    net <- RcppML::factor_net(inputs = inputs, output = output, config = config)
    result <- RcppML::fit(net)

    # Extract shared results
    layer <- result$layers$shared_nmf
    d_vec <- layer$d
    h_mat <- layer$H
    w_list <- layer$W

    if (absorb_d) {
        sqrt_d <- sqrt(d_vec)
        k_len <- length(sqrt_d)
        coeff_matrix <- t(diag(sqrt_d, nrow = k_len) %*% h_mat)
    } else {
        coeff_matrix <- t(h_mat)
    }

    factor_names <- paste0("NMF_", seq_len(k))
    rownames(coeff_matrix) <- colnames(x)
    colnames(coeff_matrix) <- factor_names

    reducedDim(x, name) <- coeff_matrix

    # Store per-modality basis matrices
    for (i in seq_along(modality_names)) {
        mod <- modality_names[i]
        w_mat <- w_list[[mod]]

        if (absorb_d) {
            basis <- w_mat %*% diag(sqrt_d, nrow = k_len)
        } else {
            basis <- w_mat
        }

        rownames(basis) <- rownames(mats[[mod]])
        colnames(basis) <- factor_names
        metadata(x)[[paste0(name, "_basis_", mod)]] <- basis
    }

    metadata(x)[[paste0(name, "_modalities")]] <- modality_names
    metadata(x)[[paste0(name, "_model")]] <- result

    if (verbose) {
        message("Multi-modal NMF completed. Shared embedding in ",
                "reducedDim(x, '", name, "')")
    }

    return(x)
}
