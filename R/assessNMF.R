#' Assess NMF embedding quality against known labels
#'
#' Evaluates how well NMF factors separate known cell types or groups using
#' clustering metrics (ARI, NMI, silhouette) and cross-validated
#' classification. Wraps \code{\link[RcppML]{assess}} with a
#' SingleCellExperiment interface.
#'
#' @param x A SingleCellExperiment or SpatialExperiment object with NMF results
#' @param name Character, name of the NMF result to assess (default "NMF")
#' @param label_col Character, column name in \code{colData(x)} containing
#'   group labels (e.g., cell type annotations)
#' @param batch_col Character, optional column name for batch variable
#'   (default NULL)
#' @param min_class_size Integer, minimum samples per class to include
#'   (default 10)
#' @param ... Additional arguments passed to \code{\link[RcppML]{assess}},
#'   including \code{metrics}, \code{n_folds}, \code{classifiers},
#'   \code{k_nn}, \code{seed}
#'
#' @return A list of class \code{nmf_assessment} with components:
#'   \itemize{
#'     \item \code{metrics}: named list of metric values (ari, nmi, silhouette,
#'       accuracy, f1, precision, recall, auroc per classifier)
#'     \item \code{classification}: data.frame of per-fold classification results
#'     \item \code{params}: list of assessment parameters
#'   }
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 200)
#' sce <- logNormCounts(sce)
#' sce <- runNMFscape(sce, k = 5, verbose = FALSE)
#' sce$celltype <- sample(c("TypeA", "TypeB", "TypeC"), ncol(sce),
#'                        replace = TRUE)
#' result <- assessNMF(sce, label_col = "celltype", min_class_size = 5)
#' result$metrics$ari
#'
#' @export
#' @importFrom RcppML assess
assessNMF <- function(x, name = "NMF", label_col,
                      batch_col = NULL, min_class_size = 10, ...) {

    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment or SpatialExperiment object")
    }

    if (!label_col %in% colnames(colData(x))) {
        stop("label_col '", label_col, "' not found in colData(x)")
    }

    model <- getModel(x, name)
    labels <- colData(x)[[label_col]]

    batch <- NULL
    if (!is.null(batch_col)) {
        if (!batch_col %in% colnames(colData(x))) {
            stop("batch_col '", batch_col, "' not found in colData(x)")
        }
        batch <- colData(x)[[batch_col]]
    }

    RcppML::assess(model, labels = labels, batch = batch,
                   min_class_size = min_class_size, ...)
}
