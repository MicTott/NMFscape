#' Extract NMF basis matrix from SingleCellExperiment object
#'
#' For multi-modal results from \code{\link{runMultiModalNMF}}, use the
#' \code{modality} parameter to retrieve the basis for a specific modality.
#'
#' @param x A SingleCellExperiment object with NMF results
#' @param name Character, name of the NMF result to extract (default "NMF")
#' @param modality Character, optional modality name for multi-modal results
#'   (e.g., "RNA", "ATAC"). If NULL, returns the standard basis matrix.
#'
#' @return Matrix with features x factors (basis matrix)
#' @export
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#' sce <- runNMFscape(sce, k = 3, verbose = FALSE)
#' basis <- getBasis(sce)
#' dim(basis)
getBasis <- function(x, name = "NMF", modality = NULL) {
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment object")
    }

    if (!is.null(modality)) {
        basis_name <- paste0(name, "_basis_", modality)
    } else {
        basis_name <- paste0(name, "_basis")
    }

    if (!basis_name %in% names(metadata(x))) {
        # Check if this is a multi-modal result and provide helpful error
        modalities_key <- paste0(name, "_modalities")
        if (is.null(modality) && modalities_key %in% names(metadata(x))) {
            mods <- metadata(x)[[modalities_key]]
            stop("Multi-modal NMF result '", name, "' requires a modality. ",
                 "Available: ", paste(mods, collapse = ", "))
        }
        stop("Basis matrix '", basis_name, "' not found in metadata. ",
             "Run runNMFscape() first.")
    }

    return(metadata(x)[[basis_name]])
}

#' Extract NMF coefficient matrix from SingleCellExperiment object
#'
#' @param x A SingleCellExperiment object with NMF results  
#' @param name Character, name of the NMF result to extract (default "NMF")
#'
#' @return Matrix with cells x factors (coefficient matrix, same as reducedDim)
#' @export
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#' sce <- runNMFscape(sce, k = 3, verbose = FALSE)
#' coeffs <- getCoefficients(sce)
#' dim(coeffs)
getCoefficients <- function(x, name = "NMF") {
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment object")
    }
    
    if (!name %in% reducedDimNames(x)) {
        stop("reducedDim '", name, "' not found. Run runNMFscape() first.")
    }
    
    return(reducedDim(x, name))
}

#' Calculate feature loadings (which genes contribute most to each factor)
#'
#' @param x A SingleCellExperiment object with NMF results
#' @param name Character, name of the NMF result to use (default "NMF")
#' @param n Integer, number of top features to return per factor (default 10)
#' @param modality Character, optional modality name for multi-modal results
#'
#' @return List of character vectors, each containing top features for a factor
#' @export
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#' sce <- runNMFscape(sce, k = 3, verbose = FALSE)
#' top_genes <- getTopFeatures(sce, n = 10)
#' head(top_genes[[1]])
getTopFeatures <- function(x, name = "NMF", n = 10, modality = NULL) {
    basis <- getBasis(x, name, modality = modality)
    
    top_features <- apply(basis, 2, function(col) {
        idx <- order(col, decreasing = TRUE)[seq_len(min(n, length(col)))]
        rownames(basis)[idx]
    })
    
    # Convert to named list
    factor_names <- colnames(basis)
    if (is.null(factor_names)) {
        factor_names <- paste0("Factor_", seq_len(ncol(basis)))
    }
    
    if (is.matrix(top_features)) {
        result <- as.list(as.data.frame(top_features, stringsAsFactors = FALSE))
    } else {
        result <- list(top_features)
    }
    names(result) <- factor_names
    
    return(result)
}

#' Reconstruct original matrix from NMF factors
#'
#' Uses the raw model (W * diag(d) * H) when available for exact reconstruction,
#' otherwise falls back to basis \%*\% t(coefficients) which is correct when
#' d was absorbed.
#'
#' @param x A SingleCellExperiment object with NMF results
#' @param name Character, name of the NMF result to use (default "NMF")
#'
#' @return Reconstructed matrix (features x cells)
#' @export
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#' sce <- runNMFscape(sce, k = 3, verbose = FALSE)
#' reconstructed <- reconstructNMF(sce)
#' dim(reconstructed)
reconstructNMF <- function(x, name = "NMF") {
    model_name <- paste0(name, "_model")
    if (model_name %in% names(metadata(x))) {
        model <- metadata(x)[[model_name]]
        if (!.isFactorNet(model)) {
            return(model@w %*% diag(model@d, nrow = length(model@d)) %*% model@h)
        }
        # FactorNet layers do not reconstruct the input on their own; the
        # effective basis and embedding do, exactly when absorb_d was TRUE.
    }
    basis <- getBasis(x, name)
    coeffs <- getCoefficients(x, name)
    return(basis %*% t(coeffs))
}

#' Extract raw RcppML NMF model from SingleCellExperiment object
#'
#' Returns the raw model object stored by \code{\link{runNMFscape}} and friends.
#' Useful for passing to RcppML functions like \code{predict},
#' \code{refine}, \code{align}, etc.
#'
#' The class depends on which function produced the result.
#' \code{\link{runNMFscape}}, \code{\link{consensusNMF}} and
#' \code{\link{refineNMF}} store an S4 \code{nmf}; the FactorNet recipes
#' (\code{\link{runMultiModalNMF}},
#' \code{\link{runConditionedNMF}}, \code{\link{runFactorNet}}) store a
#' \code{factor_net_result}, which is a list of layers rather than an S4
#' object and so does not support \code{@@w} / \code{@@h} access.
#'
#' @param x A SingleCellExperiment object with NMF results
#' @param name Character, name of the NMF result (default "NMF")
#'
#' @return Either an S4 \code{nmf} object with slots \code{@@w}, \code{@@d},
#'   \code{@@h}, \code{@@misc}, or a \code{factor_net_result} list whose
#'   \code{$layers} element holds per-layer \code{W}, \code{d} and \code{H}
#' @export
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#' sce <- runNMFscape(sce, k = 3, verbose = FALSE)
#' model <- getModel(sce)
#' dim(model@@w)
getModel <- function(x, name = "NMF") {
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment object")
    }
    model_name <- paste0(name, "_model")
    if (!model_name %in% names(metadata(x))) {
        stop("NMF model '", model_name, "' not found in metadata. ",
             "Run runNMFscape() first.")
    }
    metadata(x)[[model_name]]
}

#' Extract diagonal scaling vector from NMF decomposition
#'
#' Returns the diagonal scaling vector d from the decomposition A = W * diag(d) * H.
#'
#' @param x A SingleCellExperiment object with NMF results
#' @param name Character, name of the NMF result (default "NMF")
#'
#' @return Numeric vector of length k (diagonal scaling factors)
#' @export
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#' sce <- runNMFscape(sce, k = 3, verbose = FALSE)
#' getDiagonal(sce)
getDiagonal <- function(x, name = "NMF") {
    model <- getModel(x, name)
    if (.isFactorNet(model)) {
        layer_name <- .primaryLayer(x, name)
        if (is.null(layer_name) || !layer_name %in% names(model$layers)) {
            stop("Cannot locate the primary layer for FactorNet result '",
                 name, "'. Available layers: ",
                 paste(names(model$layers), collapse = ", "))
        }
        return(model$layers[[layer_name]]$d)
    }
    model@d
}