#' Run NMFscape Non-negative Matrix Factorization on SingleCellExperiment objects
#'
#' Performs NMF using RcppML on SingleCellExperiment or SpatialExperiment objects
#' and stores results in reducedDims slot with metadata. The raw RcppML model
#' object is also stored for use by downstream functions like
#' \code{\link{predictNMF}} and \code{\link{refineNMF}}.
#'
#' @param x A SingleCellExperiment or SpatialExperiment object
#' @param k Integer, number of factors for NMF (rank)
#' @param assay Character or integer, which assay to use (default "logcounts")
#' @param name Character, name for the reducedDim slot (default "NMF")
#' @param subset_row Vector specifying which features to use
#' @param tol Numeric, tolerance for convergence (default 1e-5)
#' @param maxit Integer, maximum iterations (default 100)
#' @param L1 Numeric vector of length 2, L1 regularization for [w, h] (default c(0,0))
#' @param L2 Numeric vector of length 2, L2 regularization for [w, h] (default c(0,0))
#' @param distribution Character, loss function distribution. One of "mse",
#'   "gp" (generalized Poisson), "nb" (negative binomial), "gamma",
#'   "inverse_gaussian", "tweedie", or "auto" to select via BIC using
#'   \code{\link[RcppML]{auto_nmf_distribution}} (default "mse")
#' @param absorb_d Logical, whether to absorb the diagonal scaling factor into
#'   W and H symmetrically (default TRUE). When TRUE, the stored basis and
#'   coefficient matrices satisfy A ~ W_eff \%*\% t(H_eff). When FALSE, raw W
#'   and H are stored and the diagonal is available via \code{\link{getDiagonal}}.
#' @param seed Integer, random seed for reproducibility
#' @param verbose Logical, whether to print progress (default TRUE)
#' @param ... Additional arguments passed to \code{\link[RcppML]{nmf}}, including
#'   \code{nonneg}, \code{robust}, \code{zi}, \code{projective}, \code{symmetric}
#'
#' @return The input object with NMF results stored in:
#'   \itemize{
#'     \item \code{reducedDim(x, name)}: cells x k coefficient matrix
#'     \item \code{metadata(x)[[paste0(name, "_basis")]]}: genes x k basis matrix
#'     \item \code{metadata(x)[[paste0(name, "_model")]]}: raw RcppML S4 nmf object
#'   }
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' sce <- runNMFscape(sce, k = 10)
#'
#' # Access results
#' nmf_coords <- reducedDim(sce, "NMF")
#' basis <- metadata(sce)$NMF_basis
#' model <- getModel(sce)
#'
#' @export
#' @importFrom RcppML nmf
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom S4Vectors metadata<-
#' @importFrom Matrix Matrix
runNMFscape <- function(x, k, assay = "logcounts", name = "NMF",
                        subset_row = NULL, tol = 1e-5, maxit = 100,
                        L1 = c(0, 0), L2 = c(0, 0),
                        distribution = "mse", absorb_d = TRUE,
                        seed = NULL, verbose = TRUE, ...) {

    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment or SpatialExperiment object")
    }

    if (!assay %in% assayNames(x)) {
        stop("assay '", assay, "' not found in x")
    }

    mat <- assay(x, assay)

    if (!is.null(subset_row)) {
        mat <- mat[subset_row, , drop = FALSE]
    }

    if (any(mat < 0)) {
        warning("Negative values detected. Consider using log-transformed data.")
        mat[mat < 0] <- 0
    }

    # Auto-select distribution if requested
    loss <- distribution
    if (identical(loss, "auto")) {
        if (verbose) message("Selecting optimal distribution via BIC...")
        dist_result <- RcppML::auto_nmf_distribution(
            data = mat, k = k, seed = seed, verbose = FALSE, ...
        )
        loss <- dist_result$best
        if (verbose) message("Selected distribution: ", loss)
    }

    if (verbose) {
        message("Running NMF with k=", k, " factors (loss=", loss, ")...")
    }

    nmf_result <- RcppML::nmf(data = mat, k = k, tol = tol, maxit = maxit,
                              L1 = L1, L2 = L2, loss = loss,
                              seed = seed, verbose = FALSE, ...)

    # Extract and optionally absorb diagonal scaling
    w_mat <- nmf_result@w
    h_mat <- nmf_result@h
    d_vec <- nmf_result@d

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
    if (!is.null(subset_row)) {
        rownames(basis_matrix) <- rownames(x)[subset_row]
    } else {
        rownames(basis_matrix) <- rownames(x)
    }
    colnames(basis_matrix) <- factor_names
    rownames(coeff_matrix) <- colnames(x)
    colnames(coeff_matrix) <- factor_names

    # Store results
    reducedDim(x, name) <- coeff_matrix
    metadata(x)[[paste0(name, "_basis")]] <- basis_matrix
    metadata(x)[[paste0(name, "_model")]] <- nmf_result

    if (verbose) {
        message("NMF completed. Results stored in reducedDim(x, '", name, "')")
    }

    return(x)
}