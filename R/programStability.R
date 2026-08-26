#' Score how reproducible each NMF program is
#'
#' \code{\link{consensusNMF}} reports a single cophenetic coefficient for the
#' whole factorization, which answers "is this rank stable?" but not "which of
#' my programs can I trust?". This function answers the second question by
#' matching the programs of every consensus replicate back onto the
#' representative fit with \code{\link{alignPrograms}}, and reporting, per
#' program, how often and how closely it reappeared.
#'
#' A program that is recovered in nearly every replicate at high similarity is
#' a robust feature of the data. One that appears sporadically, or only as a
#' weak match, is more likely to be a partition of the noise and should be
#' treated with caution regardless of how interpretable its gene loadings look.
#'
#' @param x A SingleCellExperiment or SpatialExperiment processed by
#'   \code{\link{consensusNMF}}
#' @param name Character, name of the consensus result (default "cNMF")
#' @param threshold Numeric, the similarity above which a replicate counts as
#'   having recovered a program (default 0.8)
#' @param method Character, similarity measure passed to
#'   \code{\link{alignPrograms}}: "cosine" (default) or "cor"
#'
#' @return A data.frame with one row per program, ordered by decreasing
#'   \code{frequency}, containing:
#'   \itemize{
#'     \item \code{program}: program name
#'     \item \code{frequency}: fraction of replicates recovering it above
#'       \code{threshold}
#'     \item \code{mean_similarity}, \code{min_similarity}: similarity to the
#'       best match in each replicate
#'     \item \code{n_replicates}: number of replicates compared
#'   }
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 200, ncells = 100)
#' sce <- logNormCounts(sce)
#' sce <- consensusNMF(sce, k = 3, n_runs = 5, seed = 1, verbose = FALSE)
#' programStability(sce)
#'
#' @seealso \code{\link{consensusNMF}} for the factorization,
#'   \code{\link{alignPrograms}} for the underlying matching
#' @export
programStability <- function(x, name = "cNMF", threshold = 0.8,
                             method = c("cosine", "cor")) {
    method <- match.arg(method)

    key <- paste0(name, "_consensus")
    if (!key %in% names(metadata(x))) {
        stop("Consensus result '", name, "' not found in metadata. ",
             "Run consensusNMF() first.")
    }
    models <- metadata(x)[[key]]$models
    if (!is.list(models) || length(models) < 2) {
        stop("Consensus result '", name, "' stores fewer than two replicate ",
             "models, so stability cannot be assessed. Re-run consensusNMF() ",
             "with a larger n_runs.")
    }

    reference <- getBasis(x, name)
    others <- models[-1]

    sims <- vapply(others, function(model) {
        basis <- model@w
        rownames(basis) <- rownames(reference)
        colnames(basis) <- colnames(reference)
        mapping <- alignPrograms(reference, basis, method = method)$mapping
        # restore the reference program order; unmatched programs score 0
        out <- rep(0, ncol(reference))
        idx <- match(colnames(reference), mapping$program_x)
        matched <- !is.na(idx) & !is.na(mapping$similarity[idx])
        out[matched] <- mapping$similarity[idx][matched]
        out
    }, numeric(ncol(reference)))

    sims <- matrix(sims, nrow = ncol(reference))
    result <- data.frame(
        program = colnames(reference),
        frequency = rowMeans(sims >= threshold),
        mean_similarity = rowMeans(sims),
        min_similarity = apply(sims, 1, min),
        n_replicates = ncol(sims),
        stringsAsFactors = FALSE
    )
    result <- result[order(-result$frequency, -result$mean_similarity), ]
    rownames(result) <- NULL
    result
}
