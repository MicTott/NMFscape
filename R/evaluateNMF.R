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
#' @seealso \code{\link{getModel}} for which functions store an S4
#'   \code{nmf} model. FactorNet-based results are not supported here;
#'   use their stored \code{$total_loss} instead.
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

    if (.isFactorNet(model)) {
        stop("evaluateNMF() requires an S4 nmf model, but '", name,
             "' was produced by a FactorNet recipe (",
             "runMultiModalNMF, runConditionedNMF or runFactorNet). ",
             "The fitted loss for such models is available as ",
             "metadata(x)[['", name, "_model']]$total_loss, and per-layer ",
             "losses via $layers.")
    }

    mat <- assay(x, assay)

    if (any(mat < 0)) {
        mat[mat < 0] <- 0
    }

    RcppML::evaluate(model, mat)
}
