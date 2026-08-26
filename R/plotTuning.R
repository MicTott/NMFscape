#' Plot hyperparameter cross-validation results
#'
#' Visualizes the output of \code{\link{tuneNMF}}. The layout adapts to how
#' many parameters were searched:
#' \itemize{
#'   \item One parameter: mean held-out loss against that parameter, with
#'     standard-error bars and the best value marked, in the style of
#'     \code{\link{plotRankSelection}}.
#'   \item Two parameters: a tile heatmap of mean held-out loss, with the best
#'     cell outlined.
#'   \item Three or more: a ranked dot plot of the best \code{top_n}
#'     combinations.
#' }
#' Deep-recipe results, which cross-validate one layer at a time, are faceted by
#' layer.
#'
#' @param tuning An \code{nmfscape_tuning} object from \code{\link{tuneNMF}}
#' @param show_se Logical, whether to draw standard-error bars where they apply
#'   (default TRUE)
#' @param highlight_best Logical, whether to mark the best combination
#'   (default TRUE)
#' @param top_n Integer, how many combinations to show in the ranked dot plot
#'   used for three or more parameters (default 20)
#' @param ... Additional arguments (currently unused)
#'
#' @return A ggplot object
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 80, ncells = 60)
#' sce <- logNormCounts(sce)
#'
#' tuned <- tuneNMF(sce, params = list(k = c(2, 4, 6)), reps = 2,
#'                  verbose = FALSE)
#' plotTuning(tuned)
#'
#' @seealso \code{\link{tuneNMF}}, \code{\link{plotRankSelection}}
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_errorbar geom_tile
#'   geom_vline facet_wrap labs theme theme_bw element_text
#'   scale_fill_viridis_c
plotTuning <- function(tuning, show_se = TRUE, highlight_best = TRUE,
                       top_n = 20, ...) {

    if (!is(tuning, "nmfscape_tuning")) {
        stop("tuning must be an nmfscape_tuning object from tuneNMF()")
    }

    summary_df <- tuning$summary
    faceted <- "layer" %in% colnames(summary_df)
    varying <- .varyingParams(summary_df, tuning$param_names)

    if (length(varying) == 0) {
        stop("No parameter varies across the cross-validation results; ",
             "nothing to plot")
    }

    if (length(varying) == 1) {
        return(.plotTuningOne(summary_df, varying, faceted, show_se,
                              highlight_best))
    }
    if (length(varying) == 2 && !faceted) {
        return(.plotTuningTwo(summary_df, varying, highlight_best))
    }
    .plotTuningMany(summary_df, varying, faceted, show_se, top_n)
}


# Parameter columns that actually take more than one value. `k` is dropped from
# the check when it is a list column (deep recipe stores per-layer grids).
# @noRd
.varyingParams <- function(summary_df, param_names) {
    cols <- intersect(param_names, colnames(summary_df))
    keep <- vapply(cols, function(nm) {
        v <- summary_df[[nm]]
        if (is.list(v)) v <- vapply(v, paste, character(1), collapse = ",")
        length(unique(v)) > 1
    }, logical(1))
    cols[keep]
}

# @noRd
.tuningLabel <- function(v) {
    if (is.list(v)) {
        return(vapply(v, paste, character(1), collapse = ","))
    }
    as.character(v)
}

# @noRd
.plotTuningOne <- function(summary_df, param, faceted, show_se,
                           highlight_best) {

    dat <- data.frame(
        xval = summary_df[[param]],
        mean_test_loss = summary_df$mean_test_loss,
        se_test_loss = summary_df$se_test_loss,
        stringsAsFactors = FALSE
    )
    if (faceted) dat$layer <- factor(paste("Layer", summary_df$layer))
    dat <- dat[order(dat$xval), , drop = FALSE]

    p <- ggplot(dat, aes(x = xval, y = mean_test_loss)) +
        geom_line(color = "steelblue", linewidth = 0.8) +
        geom_point(color = "steelblue", size = 2.5) +
        theme_bw() +
        labs(x = param, y = "Mean held-out loss")

    if (show_se && any(!is.na(dat$se_test_loss))) {
        p <- p + geom_errorbar(
            aes(ymin = mean_test_loss - se_test_loss,
                ymax = mean_test_loss + se_test_loss),
            width = 0.3, color = "steelblue", alpha = 0.5
        )
    }

    if (highlight_best) {
        best <- if (faceted) {
            do.call(rbind, lapply(split(dat, dat$layer), function(d) {
                d[which.min(d$mean_test_loss), , drop = FALSE]
            }))
        } else {
            dat[which.min(dat$mean_test_loss), , drop = FALSE]
        }
        p <- p + geom_vline(data = best, aes(xintercept = xval),
                            linetype = "dashed", color = "firebrick",
                            alpha = 0.7)
    }

    if (faceted) p <- p + facet_wrap(~ layer, scales = "free")

    p
}

# @noRd
.plotTuningTwo <- function(summary_df, params, highlight_best) {

    dat <- data.frame(
        xval = .tuningLabel(summary_df[[params[1]]]),
        yval = .tuningLabel(summary_df[[params[2]]]),
        mean_test_loss = summary_df$mean_test_loss,
        stringsAsFactors = FALSE
    )
    dat$xval <- factor(dat$xval,
                       levels = .sortLevels(summary_df[[params[1]]]))
    dat$yval <- factor(dat$yval,
                       levels = .sortLevels(summary_df[[params[2]]]))

    p <- ggplot(dat, aes(x = xval, y = yval, fill = mean_test_loss)) +
        geom_tile(color = "white") +
        scale_fill_viridis_c(name = "Mean\nheld-out\nloss", direction = -1) +
        theme_bw() +
        labs(x = params[1], y = params[2])

    if (highlight_best) {
        best <- dat[which.min(dat$mean_test_loss), , drop = FALSE]
        p <- p + geom_tile(data = best, aes(x = xval, y = yval),
                           fill = NA, color = "firebrick", linewidth = 1.1)
    }

    p
}

# @noRd
.plotTuningMany <- function(summary_df, params, faceted, show_se, top_n) {

    labels <- do.call(paste, c(
        lapply(params, function(nm) {
            paste0(nm, "=", .tuningLabel(summary_df[[nm]]))
        }),
        list(sep = ", ")
    ))

    dat <- data.frame(
        combo_label = labels,
        mean_test_loss = summary_df$mean_test_loss,
        se_test_loss = summary_df$se_test_loss,
        stringsAsFactors = FALSE
    )
    if (faceted) {
        dat$layer <- factor(paste("Layer", summary_df$layer))
        dat$combo_label <- paste0(dat$combo_label, " [L",
                                  summary_df$layer, "]")
    }

    dat <- dat[order(dat$mean_test_loss), , drop = FALSE]
    dat <- utils::head(dat, top_n)
    dat$combo_label <- factor(dat$combo_label, levels = dat$combo_label)

    p <- ggplot(dat, aes(x = combo_label, y = mean_test_loss)) +
        geom_point(color = "steelblue", size = 2.5) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = NULL, y = "Mean held-out loss")

    if (show_se && any(!is.na(dat$se_test_loss))) {
        p <- p + geom_errorbar(
            aes(ymin = mean_test_loss - se_test_loss,
                ymax = mean_test_loss + se_test_loss),
            width = 0.3, color = "steelblue", alpha = 0.5
        )
    }

    if (faceted) p <- p + facet_wrap(~ layer, scales = "free_x")

    p
}

# Order tile axes numerically when the values are numbers, alphabetically
# otherwise.
# @noRd
.sortLevels <- function(v) {
    labels <- .tuningLabel(v)
    if (is.numeric(v)) {
        keep <- !duplicated(labels)
        return(labels[keep][order(v[keep])])
    }
    sort(unique(labels))
}
