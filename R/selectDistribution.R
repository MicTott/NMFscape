#' Select optimal NMF distribution for a dataset
#'
#' Evaluates multiple statistical distributions (loss functions) for NMF and
#' recommends the best one based on BIC or AIC. Wraps
#' \code{\link[RcppML]{auto_nmf_distribution}} with a SingleCellExperiment
#' interface.
#'
#' @param x A SingleCellExperiment or SpatialExperiment object
#' @param k Integer, rank to use for comparison (default 10)
#' @param assay Character, which assay to use (default "logcounts")
#' @param subset_row Vector specifying which features to use
#' @param distributions Character vector of distributions to compare
#'   (default c("mse", "gp", "nb"))
#' @param criterion Character, model selection criterion: "bic" or "aic"
#'   (default "bic")
#' @param seed Integer, random seed for reproducibility
#' @param verbose Logical, whether to print progress (default TRUE)
#' @param ... Additional arguments passed to
#'   \code{\link[RcppML]{auto_nmf_distribution}}
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{best}: character, name of the best distribution
#'     \item \code{comparison}: data.frame with columns \code{distribution},
#'       \code{nll}, \code{df}, \code{aic}, \code{bic}, \code{selected}
#'   }
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#' result <- selectDistribution(sce, k = 5)
#' result$best
#' result$comparison
#'
#' @export
#' @importFrom RcppML auto_nmf_distribution
selectDistribution <- function(x, k = 10, assay = "logcounts",
                               subset_row = NULL,
                               distributions = c("mse", "gp", "nb"),
                               criterion = "bic",
                               seed = NULL, verbose = TRUE, ...) {

    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment or SpatialExperiment object")
    }

    if (!assay %in% assayNames(x)) {
        stop("assay '", assay, "' not found in x")
    }

    mat <- assay(x, assay)

    if (!is.null(subset_row)) {
        mat <- mat[subset_row, , drop = FALSE]
    }

    if (any(mat < 0)) {
        warning("Negative values detected. Setting to 0 for NMF.")
        mat[mat < 0] <- 0
    }

    if (verbose) {
        message("Comparing distributions: ",
                paste(distributions, collapse = ", "), "...")
    }

    result <- RcppML::auto_nmf_distribution(
        data = mat, k = k, distributions = distributions,
        criterion = criterion, seed = seed, verbose = FALSE, ...
    )

    if (verbose) {
        message("Best distribution: ", result$loss)
    }

    list(
        best = result$loss,
        comparison = result$comparison
    )
}
