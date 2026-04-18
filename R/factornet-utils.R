# Internal helpers for FactorNet wrappers
# @keywords internal
# @noRd

.validateSCE <- function(x, assay) {
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment or SpatialExperiment object")
    }
    if (!assay %in% assayNames(x)) {
        stop("assay '", assay, "' not found in x")
    }
}

.extractAssayMatrix <- function(x, assay_name, subset_row = NULL) {
    mat <- assay(x, assay_name)
    if (!is.null(subset_row)) {
        mat <- mat[subset_row, , drop = FALSE]
    }
    if (any(mat < 0)) {
        warning("Negative values detected. Setting to 0 for NMF.")
        mat[mat < 0] <- 0
    }
    mat
}

.absorbDiagonal <- function(w_mat, h_mat, d_vec, absorb_d) {
    if (absorb_d) {
        sqrt_d <- sqrt(d_vec)
        k <- length(sqrt_d)
        basis <- w_mat %*% diag(sqrt_d, nrow = k)
        coeff <- t(diag(sqrt_d, nrow = k) %*% h_mat)
    } else {
        basis <- w_mat
        coeff <- t(h_mat)
    }
    list(basis = basis, coeff = coeff)
}

.setFactorNames <- function(basis, coeff, feature_names, cell_names, k) {
    factor_names <- paste0("NMF_", seq_len(k))
    rownames(basis) <- feature_names
    colnames(basis) <- factor_names
    rownames(coeff) <- cell_names
    colnames(coeff) <- factor_names
    list(basis = basis, coeff = coeff)
}
