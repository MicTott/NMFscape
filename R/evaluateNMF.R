#' Evaluate NMF model reconstruction quality
#'
#' Computes the reconstruction loss for an NMF model against its original
#' data matrix. Wraps \code{\link[RcppML]{evaluate}} with a
#' SingleCellExperiment interface.
#'
#' @param x A SingleCellExperiment or SpatialExperiment object with NMF results
#' @param name Character, name of the NMF result to evaluate (default "NMF")
#' @param assay Character, which assay was used for NMF (default "logcounts")
#'
#' @return Numeric scalar, the reconstruction loss (lower is better)
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#' sce <- runNMFscape(sce, k = 5, verbose = FALSE)
#' evaluateNMF(sce)
#'
#' @export
#' @importFrom RcppML evaluate
evaluateNMF <- function(x, name = "NMF", assay = "logcounts") {

    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment or SpatialExperiment object")
    }

    if (!assay %in% assayNames(x)) {
        stop("assay '", assay, "' not found in x")
    }

    model <- getModel(x, name)
    mat <- assay(x, assay)

    if (any(mat < 0)) {
        mat[mat < 0] <- 0
    }

    RcppML::evaluate(model, mat)
}
