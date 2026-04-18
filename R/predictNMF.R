#' Project new data onto an existing NMF model
#'
#' Projects a new SingleCellExperiment onto an NMF model learned from a
#' reference dataset. Handles feature alignment between the new data and
#' the model's basis matrix. Wraps \code{RcppML::predict} with
#' SCE integration.
#'
#' @param new_data A SingleCellExperiment or SpatialExperiment with new samples
#' @param reference Either a SingleCellExperiment containing a stored NMF model
#'   (from \code{\link{runNMFscape}}), or a raw RcppML \code{nmf} S4 object
#' @param ref_name Character, name of the NMF result in the reference SCE
#'   (default "NMF"). Ignored if \code{reference} is a raw model.
#' @param name Character, name for the reducedDim slot in new_data
#'   (default "NMF_projected")
#' @param assay Character, which assay to use from new_data (default "logcounts")
#' @param absorb_d Logical, whether to absorb diagonal scaling into the
#'   projected coefficients (default TRUE)
#' @param verbose Logical, whether to print progress (default TRUE)
#' @param ... Additional arguments passed to \code{RcppML::predict}
#'
#' @return The \code{new_data} object with projected NMF coefficients stored
#'   in \code{reducedDim(new_data, name)}
#'
#' @examples
#' library(scuttle)
#' # Create reference and query datasets
#' ref_sce <- mockSCE(ngenes = 100, ncells = 50)
#' ref_sce <- logNormCounts(ref_sce)
#' ref_sce <- runNMFscape(ref_sce, k = 5, verbose = FALSE)
#'
#' query_sce <- mockSCE(ngenes = 100, ncells = 30)
#' query_sce <- logNormCounts(query_sce)
#' query_sce <- predictNMF(query_sce, reference = ref_sce, verbose = FALSE)
#'
#' @export
predictNMF <- function(new_data, reference, ref_name = "NMF",
                       name = "NMF_projected", assay = "logcounts",
                       absorb_d = TRUE, verbose = TRUE, ...) {

    if (!is(new_data, "SingleCellExperiment")) {
        stop("new_data must be a SingleCellExperiment or SpatialExperiment object")
    }

    # Extract model
    if (is(reference, "SingleCellExperiment")) {
        model <- getModel(reference, ref_name)
    } else if (is(reference, "nmf")) {
        model <- reference
    } else {
        stop("reference must be a SingleCellExperiment with stored NMF results ",
             "or a raw RcppML nmf model object")
    }

    if (!assay %in% assayNames(new_data)) {
        stop("assay '", assay, "' not found in new_data")
    }

    mat <- assay(new_data, assay)

    # Align features
    model_features <- rownames(model@w)
    if (!is.null(model_features)) {
        shared <- intersect(rownames(mat), model_features)
        if (length(shared) == 0) {
            stop("No shared features between new_data and the NMF model")
        }
        if (length(shared) < length(model_features)) {
            if (verbose) {
                message("Using ", length(shared), " of ",
                        length(model_features), " model features")
            }
            # Subset both to shared features in model order
            mat <- mat[shared, , drop = FALSE]
            # Build a subsetted model with matching w rows
            w_sub <- model@w[shared, , drop = FALSE]
            model <- new("nmf",
                         w = w_sub, d = model@d, h = model@h,
                         misc = model@misc)
        } else {
            mat <- mat[model_features, , drop = FALSE]
        }
    }

    if (any(mat < 0)) {
        mat[mat < 0] <- 0
    }

    if (verbose) message("Projecting new data onto NMF model...")

    pred <- predict(model, data = mat, ...)

    # Extract projected coefficients
    h_proj <- pred@h
    d_vec <- pred@d

    if (absorb_d) {
        sqrt_d <- sqrt(d_vec)
        coeff_matrix <- t(diag(sqrt_d, nrow = length(sqrt_d)) %*% h_proj)
    } else {
        coeff_matrix <- t(h_proj)
    }

    k <- ncol(coeff_matrix)
    factor_names <- paste0("NMF_", seq_len(k))
    rownames(coeff_matrix) <- colnames(new_data)
    colnames(coeff_matrix) <- factor_names

    reducedDim(new_data, name) <- coeff_matrix

    if (verbose) {
        message("Projection stored in reducedDim(new_data, '", name, "')")
    }

    return(new_data)
}
