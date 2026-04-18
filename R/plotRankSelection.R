#' Plot cross-validation results for NMF rank selection
#'
#' Visualizes the output of \code{\link{selectRank}} as a line plot of
#' test error vs rank k, with optional error bars and highlight of the
#' best rank.
#'
#' @param cv_results A data.frame from \code{\link{selectRank}}
#' @param show_sd Logical, whether to show error bars when available
#'   (default TRUE)
#' @param highlight_best Logical, whether to highlight the rank with
#'   lowest test error (default TRUE)
#' @param ... Additional arguments (currently unused)
#'
#' @return A ggplot object
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#' cv_results <- selectRank(sce, k = 2:8, verbose = FALSE)
#' plotRankSelection(cv_results)
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_line theme_bw labs
#'   geom_vline geom_errorbar
plotRankSelection <- function(cv_results, show_sd = TRUE,
                              highlight_best = TRUE, ...) {

    has_sd <- "sd_test_error" %in% colnames(cv_results)

    p <- ggplot(cv_results, aes(x = k, y = test_error)) +
        geom_line(color = "steelblue", linewidth = 0.8) +
        geom_point(color = "steelblue", size = 2.5) +
        theme_bw() +
        labs(x = "Rank (k)", y = "Test Error (MSE)")

    if (show_sd && has_sd) {
        p <- p + geom_errorbar(
            aes(ymin = test_error - sd_test_error,
                ymax = test_error + sd_test_error),
            width = 0.3, color = "steelblue", alpha = 0.5
        )
    }

    if (highlight_best) {
        best_idx <- which.min(cv_results$test_error)
        best_k <- cv_results$k[best_idx]
        p <- p + geom_vline(xintercept = best_k, linetype = "dashed",
                            color = "firebrick", alpha = 0.7)
    }

    p
}
