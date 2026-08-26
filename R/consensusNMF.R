#' Run consensus NMF on a SingleCellExperiment object
#'
#' Runs multiple NMF replicates and computes a consensus clustering to
#' identify stable sample groupings. Wraps
#' \code{\link[RcppML]{consensus_nmf}} with a SingleCellExperiment interface.
#'
#' @param x A SingleCellExperiment or SpatialExperiment object
#' @param k Integer, number of factors (rank)
#' @param assay Character, which assay to use (default "logcounts")
#' @param name Character, name for the reducedDim slot (default "cNMF")
#' @param subset_row Vector specifying which features to use
#' @param n_runs Integer, number of NMF replicates (default 50)
#' @param absorb_d Logical, whether to absorb diagonal scaling into W/H
#'   (default TRUE)
#' @param seed Integer, random seed for reproducibility
#' @param verbose Logical, whether to print progress (default TRUE)
#' @param ... Additional arguments passed to
#'   \code{\link[RcppML]{consensus_nmf}}
#'
#' @return The input object with consensus NMF results stored in:
#'   \itemize{
#'     \item \code{reducedDim(x, name)}: cells x k coefficient matrix
#'       (from first replicate)
#'     \item \code{metadata(x)[[paste0(name, "_basis")]]}: genes x k basis
#'       matrix
#'     \item \code{metadata(x)[[paste0(name, "_model")]]}: the representative
#'       replicate as a raw RcppML S4 nmf object, so that \code{\link{getModel}},
#'       \code{\link{getDiagonal}}, \code{\link{evaluateNMF}} and
#'       \code{\link{diagnoseNMF}} work on consensus results
#'     \item \code{metadata(x)[[paste0(name, "_consensus")]]}: full consensus
#'       result list (consensus matrix, clusters, cophenetic, hclust)
#'   }
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#' sce <- consensusNMF(sce, k = 3, n_runs = 10)
#'
#' @export
#' @importFrom RcppML consensus_nmf
consensusNMF <- function(x, k, assay = "logcounts", name = "cNMF",
                         subset_row = NULL, n_runs = 50,
                         absorb_d = TRUE, seed = NULL,
                         verbose = TRUE, ...) {

    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment or SpatialExperiment object")
    }

    if (!assay %in% assayNames(x)) {
        stop("assay '", assay, "' not found in x")
    }

    feature_names <- .subsetRowNames(x, subset_row)

    mat <- assay(x, assay)

    if (!is.null(subset_row)) {
        mat <- mat[subset_row, , drop = FALSE]
    }

    if (any(mat < 0)) {
        warning("Negative values detected. Setting to 0 for NMF.")
        mat[mat < 0] <- 0
    }

    if (verbose) {
        message("Running consensus NMF with k=", k,
                " (", n_runs, " replicates)...")
    }

    cres <- RcppML::consensus_nmf(
        data = mat, k = k, reps = n_runs,
        seed = seed, verbose = FALSE, ...
    )

    # Use the first model for representative W/H
    model <- cres$models[[1]]
    w_mat <- model@w
    h_mat <- model@h
    d_vec <- model@d

    if (absorb_d) {
        sqrt_d <- sqrt(d_vec)
        basis_matrix <- w_mat %*% diag(sqrt_d, nrow = length(sqrt_d))
        coeff_matrix <- t(diag(sqrt_d, nrow = length(sqrt_d)) %*% h_mat)
    } else {
        basis_matrix <- w_mat
        coeff_matrix <- t(h_mat)
    }

    # Set names
    factor_names <- paste0("NMF_", seq_len(k))
    rownames(basis_matrix) <- feature_names
    colnames(basis_matrix) <- factor_names
    rownames(coeff_matrix) <- colnames(x)
    colnames(coeff_matrix) <- factor_names

    # Store results
    reducedDim(x, name) <- coeff_matrix
    metadata(x)[[paste0(name, "_basis")]] <- basis_matrix
    metadata(x)[[paste0(name, "_model")]] <- model
    metadata(x)[[paste0(name, "_consensus")]] <- cres

    if (verbose) {
        message("Consensus NMF completed (cophenetic = ",
                round(cres$cophenetic, 4), ")")
    }

    return(x)
}
