#' Refine NMF model using cell type labels
#'
#' Performs label-guided refinement of an NMF model by shifting embeddings
#' toward class centroids. Optionally includes batch correction. Wraps
#' \code{\link[RcppML]{refine}} with a SingleCellExperiment interface.
#'
#' @param x A SingleCellExperiment or SpatialExperiment object with NMF results
#' @param label_col Character, column name in \code{colData(x)} containing
#'   cell type or group labels
#' @param name Character, name of the NMF result to refine (default "NMF")
#' @param refined_name Character, name for the refined NMF result
#'   (default "NMF_refined")
#' @param assay Character, which assay to use (default "logcounts")
#' @param batch_col Character, optional column name in \code{colData(x)} for
#'   batch correction (default NULL)
#' @param absorb_d Logical, whether to absorb diagonal scaling into the
#'   refined W/H (default TRUE)
#' @param verbose Logical, whether to print progress (default TRUE)
#' @param ... Additional arguments passed to \code{\link[RcppML]{refine}},
#'   including \code{lambda}, \code{cycles}, \code{nonneg}, \code{whiten}
#'
#' @return The input object with refined NMF results stored in:
#'   \itemize{
#'     \item \code{reducedDim(x, refined_name)}: cells x k refined coefficients
#'     \item \code{metadata(x)[[paste0(refined_name, "_basis")]]}: genes x k
#'       refined basis matrix
#'     \item \code{metadata(x)[[paste0(refined_name, "_model")]]}: raw refined
#'       RcppML nmf object
#'   }
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#' sce <- runNMFscape(sce, k = 5, verbose = FALSE)
#'
#' # Add fake labels for demonstration
#' sce$celltype <- sample(c("TypeA", "TypeB"), ncol(sce), replace = TRUE)
#' sce <- refineNMF(sce, label_col = "celltype", verbose = FALSE)
#'
#' @export
#' @importFrom RcppML refine
refineNMF <- function(x, label_col, name = "NMF",
                      refined_name = "NMF_refined",
                      assay = "logcounts", batch_col = NULL,
                      absorb_d = TRUE, verbose = TRUE, ...) {

    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment or SpatialExperiment object")
    }

    if (!label_col %in% colnames(colData(x))) {
        stop("label_col '", label_col, "' not found in colData(x)")
    }

    model <- getModel(x, name)
    labels <- colData(x)[[label_col]]

    if (!assay %in% assayNames(x)) {
        stop("assay '", assay, "' not found in x")
    }

    mat <- assay(x, assay)

    if (any(mat < 0)) {
        mat[mat < 0] <- 0
    }

    batch <- NULL
    if (!is.null(batch_col)) {
        if (!batch_col %in% colnames(colData(x))) {
            stop("batch_col '", batch_col, "' not found in colData(x)")
        }
        batch <- colData(x)[[batch_col]]
    }

    if (verbose) message("Refining NMF with labels from '", label_col, "'...")

    refined <- RcppML::refine(model, data = mat, labels = labels,
                              batch = batch, ...)

    # Extract and store refined results
    w_mat <- refined@w
    h_mat <- refined@h
    d_vec <- refined@d

    if (absorb_d) {
        sqrt_d <- sqrt(d_vec)
        k <- length(sqrt_d)
        basis_matrix <- w_mat %*% diag(sqrt_d, nrow = k)
        coeff_matrix <- t(diag(sqrt_d, nrow = k) %*% h_mat)
    } else {
        basis_matrix <- w_mat
        coeff_matrix <- t(h_mat)
    }

    k <- ncol(basis_matrix)
    factor_names <- paste0("NMF_", seq_len(k))
    rownames(basis_matrix) <- rownames(x)
    colnames(basis_matrix) <- factor_names
    rownames(coeff_matrix) <- colnames(x)
    colnames(coeff_matrix) <- factor_names

    reducedDim(x, refined_name) <- coeff_matrix
    metadata(x)[[paste0(refined_name, "_basis")]] <- basis_matrix
    metadata(x)[[paste0(refined_name, "_model")]] <- refined

    if (verbose) {
        message("Refined NMF stored in reducedDim(x, '", refined_name, "')")
    }

    return(x)
}
