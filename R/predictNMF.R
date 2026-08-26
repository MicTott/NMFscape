#' Project new data onto an existing NMF model
#'
#' Projects a new SingleCellExperiment onto an NMF model learned from a
#' reference dataset. Handles feature alignment between the new data and
#' the model's basis matrix. Wraps \code{RcppML::predict} with
#' SCE integration.
#'
#' @param new_data A SingleCellExperiment or SpatialExperiment with new samples
#' @param reference Either a SingleCellExperiment containing a stored NMF model
#'   (from \code{\link{runNMFscape}} or one of the FactorNet recipes), a raw
#'   RcppML \code{nmf} S4 object, or a fitted \code{factor_net_result}.
#'   FactorNet references are projected through the whole graph and require
#'   \code{new_data} to contain every feature used to fit the reference.
#' @param ref_name Character, name of the NMF result in the reference SCE
#'   (default "NMF"). Ignored if \code{reference} is a raw model.
#' @param name Character, name for the reducedDim slot in new_data
#'   (default "NMF_projected")
#' @param assay Character, which assay to use from new_data (default "logcounts")
#' @param absorb_d Logical, whether to absorb diagonal scaling into the
#'   projected coefficients (default TRUE). Ignored for FactorNet references,
#'   which return coefficients with scaling already applied.
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
    } else if (is(reference, "nmf") || .isFactorNet(reference)) {
        model <- reference
    } else {
        stop("reference must be a SingleCellExperiment with stored NMF results, ",
             "a raw RcppML nmf model, or a fitted factor_net_result")
    }

    if (!assay %in% assayNames(new_data)) {
        stop("assay '", assay, "' not found in new_data")
    }

    mat <- assay(new_data, assay)

    # FactorNet results are lists of layers, not S4 nmf objects, and RcppML
    # projects them through the whole graph rather than a single basis.
    if (.isFactorNet(model)) {
        # RcppML 1.0.0 only projects single-layer, single-modal graphs; deeper
        # graphs error inside predict() and multi-modal ones need new data for
        # every modality row-bound in training order.
        if (isTRUE(model$multi_modal)) {
            stop("predictNMF() cannot project multi-modal FactorNet results. ",
                 "RcppML::predict() requires new data for every modality ",
                 "concatenated in training order; project each modality ",
                 "separately with runNMFscape() instead.")
        }
        if (!is.null(model$n_layers) && model$n_layers > 1) {
            stop("predictNMF() cannot project multi-layer FactorNet results ",
                 "such as those from runDeepNMF(). RcppML 1.0.0 supports ",
                 "predict() only for single-layer graphs.")
        }

        ref_features <- NULL
        if (is(reference, "SingleCellExperiment")) {
            basis_name <- paste0(ref_name, "_basis")
            if (basis_name %in% names(metadata(reference))) {
                ref_features <- rownames(metadata(reference)[[basis_name]])
            }
        }
        if (!is.null(ref_features)) {
            missing <- setdiff(ref_features, rownames(mat))
            if (length(missing) > 0) {
                stop("new_data is missing ", length(missing), " of ",
                     length(ref_features), " features used to fit '",
                     ref_name, "'. FactorNet projection requires the full ",
                     "feature set; subset the reference instead.")
            }
            mat <- mat[ref_features, , drop = FALSE]
        }
        if (any(mat < 0)) {
            mat[mat < 0] <- 0
        }
        if (verbose) message("Projecting new data through FactorNet graph...")

        h_proj <- stats::predict(model, newdata = as.matrix(mat), ...)
        coeff_matrix <- t(h_proj)
        rownames(coeff_matrix) <- colnames(new_data)
        colnames(coeff_matrix) <- paste0("NMF_", seq_len(ncol(coeff_matrix)))
        reducedDim(new_data, name) <- coeff_matrix

        if (verbose) {
            message("Projection stored in reducedDim(new_data, '", name, "')")
        }
        return(new_data)
    }

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
