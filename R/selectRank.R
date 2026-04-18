#' Select optimal NMF rank via cross-validation
#'
#' Runs NMF across a range of ranks using held-out data to estimate
#' reconstruction error. Wraps \code{\link[RcppML]{nmf}} with
#' \code{test_fraction > 0} and k as a vector.
#'
#' @param x A SingleCellExperiment or SpatialExperiment object
#' @param k Integer vector, ranks to evaluate (default 2:20)
#' @param assay Character, which assay to use (default "logcounts")
#' @param subset_row Vector specifying which features to use
#' @param test_fraction Numeric, fraction of data to hold out (default 0.2)
#' @param distribution Character, loss function (default "mse")
#' @param n_reps Integer, number of replicate runs with different seeds
#'   (default 1). If > 1, results are aggregated with mean and sd.
#' @param seed Integer, random seed for reproducibility
#' @param verbose Logical, whether to print progress (default TRUE)
#' @param ... Additional arguments passed to \code{\link[RcppML]{nmf}}
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item \code{k}: rank tested
#'     \item \code{train_error}: training reconstruction error
#'     \item \code{test_error}: test (held-out) reconstruction error
#'     \item \code{sd_train_error}: standard deviation across reps (if n_reps > 1)
#'     \item \code{sd_test_error}: standard deviation across reps (if n_reps > 1)
#'   }
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#' cv_results <- selectRank(sce, k = 2:8)
#' cv_results
#'
#' @export
selectRank <- function(x, k = seq(2, 20), assay = "logcounts",
                       subset_row = NULL, test_fraction = 0.2,
                       distribution = "mse", n_reps = 1,
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
        message("Running cross-validation for k = ",
                min(k), ":", max(k), " (", n_reps, " rep(s))...")
    }

    all_results <- vector("list", n_reps)
    for (i in seq_len(n_reps)) {
        rep_seed <- if (!is.null(seed)) seed + i - 1L else NULL
        cv <- RcppML::nmf(data = mat, k = k,
                          test_fraction = test_fraction,
                          loss = distribution,
                          seed = rep_seed, verbose = FALSE, ...)
        all_results[[i]] <- data.frame(
            k = cv$k,
            train_error = cv$train_mse,
            test_error = cv$test_mse,
            rep = i,
            stringsAsFactors = FALSE
        )
    }

    combined <- do.call(rbind, all_results)

    if (n_reps > 1) {
        result <- do.call(rbind, lapply(unique(combined$k), function(ki) {
            sub <- combined[combined$k == ki, , drop = FALSE]
            data.frame(
                k = ki,
                train_error = mean(sub$train_error),
                test_error = mean(sub$test_error),
                sd_train_error = stats::sd(sub$train_error),
                sd_test_error = stats::sd(sub$test_error),
                stringsAsFactors = FALSE
            )
        }))
    } else {
        result <- combined[, c("k", "train_error", "test_error"), drop = FALSE]
    }

    rownames(result) <- NULL

    if (verbose) {
        best_k <- result$k[which.min(result$test_error)]
        message("Best rank by test error: k = ", best_k)
    }

    result
}
