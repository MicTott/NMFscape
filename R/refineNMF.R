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
#' @param cycles Integer, number of refit cycles (default 0). This argument
#'   decides what the function actually does, so it is surfaced rather than
#'   left to \code{...}. At \code{cycles = 0} refinement is purely an
#'   embedding operation: the coefficients move toward class centroids and the
#'   basis matrix W is returned \emph{bit-identical} to the input model's, so
#'   the gene programs themselves are unchanged. Any positive value refits W
#'   against the shifted embedding, which is what changes the programs. See
#'   "Choosing cycles".
#' @param ... Additional arguments passed to \code{\link[RcppML]{refine}},
#'   including \code{lambda}, \code{nonneg} and \code{whiten}
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
#' @section Choosing cycles:
#' The default of 0 is inherited from \code{\link[RcppML]{refine}} and is
#' deliberately conservative, but it is easy to misread: a function called
#' "refine" that leaves the gene programs untouched surprises most users.
#'
#' The difference is large where it matters. On simulated data containing a
#' rare cell subtype, F1 for recovering that subtype at 0.5\% abundance was
#' 0.50 at \code{cycles = 0} and 1.00 at \code{cycles = 20}, with the higher
#' setting also being effectively deterministic across seeds. If the rare
#' population was never given its own factor by the unguided fit, no amount of
#' embedding-only refinement can recover it, because there is no program to
#' move. Raising \code{cycles} lets the basis be rewritten so that one can
#' appear.
#'
#' Set \code{cycles = 0} to keep the discovered programs fixed and only adjust
#' the embedding. Raise it when you want the labels to reshape the programs
#' themselves, and compare the two with \code{\link{alignPrograms}}, which
#' reports a cosine of exactly 1 when W is untouched.
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
                      cycles = 0L, absorb_d = TRUE, verbose = TRUE, ...) {

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

    # RcppML's nmf validity requires a base matrix in @w; a dgeMatrix input
    # propagates through refine() and fails the S4 check on cycles > 0.
    if (!is.matrix(mat)) {
        mat <- as.matrix(mat)
    }

    refined <- RcppML::refine(model, data = mat, labels = labels,
                              batch = batch, cycles = cycles, ...)

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
