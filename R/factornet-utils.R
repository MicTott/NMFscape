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

.subsetRowNames <- function(x, subset_row) {
    if (is.null(subset_row)) {
        return(rownames(x))
    }
    if (is.character(subset_row)) {
        missing <- setdiff(subset_row, rownames(x))
        if (length(missing) > 0) {
            stop("subset_row features not found in rownames(x): ",
                 paste(utils::head(missing, 5), collapse = ", "),
                 if (length(missing) > 5) ", ..." else "")
        }
        return(subset_row)
    }
    rownames(x)[subset_row]
}

.checkLayerPenalty <- function(value, name) {
    if (length(value) != 1 || !is.numeric(value) || is.na(value)) {
        stop(name, " must be a single numeric value for FactorNet-based ",
             "factorizations. Unlike runNMFscape(), which takes a length-2 ",
             "c(w, h) vector, RcppML::nmf_layer() applies one penalty per ",
             "layer.")
    }
    invisible(value)
}

.isFactorNet <- function(model) {
    is(model, "factor_net_result") ||
        (is.list(model) && !is.null(model$layers))
}

# The layer a FactorNet wrapper treats as its primary output, recorded at
# write time so accessors can find W/d/H without knowing which recipe ran.
.primaryLayer <- function(x, name) {
    key <- paste0(name, "_layer")
    if (!key %in% names(metadata(x))) {
        return(NULL)
    }
    metadata(x)[[key]]
}
