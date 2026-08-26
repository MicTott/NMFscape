#' Plot the program similarity matrix from an alignment
#'
#' Draws the program-by-program similarity matrix produced by
#' \code{\link{alignPrograms}} as a heatmap. Programs in \code{y} are reordered
#' so that the pairs chosen by the bipartite matching lie on the diagonal: a
#' pair of runs that recovered the same programs shows a bright diagonal and a
#' dim off-diagonal, while off-diagonal brightness marks programs that split,
#' merged, or are specific to one run.
#'
#' @param alignment A \code{programAlignment} object, the output of
#'   \code{\link{alignPrograms}}.
#' @param order_by_match Logical, whether to reorder the programs of \code{y}
#'   so matched pairs fall on the diagonal (default TRUE). Set FALSE to keep
#'   the original program order.
#' @param show_values Logical, whether to print the similarity in each tile
#'   (default TRUE).
#' @param digits Integer, number of decimal places for the printed values
#'   (default 2).
#' @param highlight_matches Logical, whether to outline the matched pairs
#'   (default TRUE).
#' @param highlight_color Character, outline color for matched pairs
#'   (default "black").
#' @param color_palette Character vector of colors defining the fill gradient
#'   from low to high similarity
#'   (default c("#f7fbff", "#6baed6", "#08306b")).
#' @param x_label Character, axis label for the programs of \code{x}
#'   (default "Programs in x").
#' @param y_label Character, axis label for the programs of \code{y}
#'   (default "Programs in y").
#'
#' @return A ggplot object.
#'
#' @seealso \code{\link{alignPrograms}}
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 200, ncells = 100)
#' sce <- logNormCounts(sce)
#' fit_a <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)
#' fit_b <- runNMFscape(sce, k = 4, seed = 2, verbose = FALSE)
#'
#' res <- alignPrograms(fit_a, fit_b)
#' plotProgramSimilarity(res)
#'
#' # keep the original program order and drop the printed values
#' plotProgramSimilarity(res, order_by_match = FALSE, show_values = FALSE)
#'
#' @export
#' @importFrom ggplot2 geom_tile geom_text scale_fill_gradientn coord_fixed
plotProgramSimilarity <- function(alignment, order_by_match = TRUE,
                                  show_values = TRUE, digits = 2,
                                  highlight_matches = TRUE,
                                  highlight_color = "black",
                                  color_palette = c("#f7fbff", "#6baed6",
                                                    "#08306b"),
                                  x_label = "Programs in x",
                                  y_label = "Programs in y") {

    if (!is(alignment, "programAlignment")) {
        stop("'alignment' must be a programAlignment object returned by ",
             "alignPrograms(), not an object of class '",
             class(alignment)[1], "'.")
    }

    sim <- alignment$similarity
    mapping <- alignment$mapping

    x_levels <- rownames(sim)
    y_levels <- colnames(sim)

    if (order_by_match) {
        # walk x in the original program order so the diagonal is readable,
        # then append the programs of y that nothing was matched to
        matched <- mapping[mapping$matched, , drop = FALSE]
        matched <- matched[match(x_levels, matched$program_x), , drop = FALSE]
        matched <- matched[!is.na(matched$program_x), , drop = FALSE]
        y_levels <- c(matched$program_y,
                      setdiff(y_levels, matched$program_y))
    }

    plot_data <- expand.grid(
        program_x = x_levels,
        program_y = y_levels,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    plot_data$similarity <- sim[cbind(plot_data$program_x,
                                      plot_data$program_y)]
    plot_data$is_match <- paste(plot_data$program_x, plot_data$program_y) %in%
        paste(mapping$program_x[mapping$matched],
              mapping$program_y[mapping$matched])

    plot_data$program_x <- factor(plot_data$program_x, levels = x_levels)
    # reverse so the first program of y sits at the top of the panel
    plot_data$program_y <- factor(plot_data$program_y, levels = rev(y_levels))

    p <- ggplot(plot_data, aes(x = program_x, y = program_y)) +
        geom_tile(aes(fill = similarity), color = "white", linewidth = 0.3) +
        scale_fill_gradientn(colours = color_palette,
                             name = paste0(alignment$method, "\nsimilarity"))

    if (highlight_matches && any(plot_data$is_match)) {
        p <- p + geom_tile(
            data = plot_data[plot_data$is_match, , drop = FALSE],
            fill = NA, color = highlight_color, linewidth = 0.8
        )
    }

    if (show_values) {
        p <- p + geom_text(
            aes(label = format(round(similarity, digits), nsmall = digits)),
            size = 3
        )
    }

    p +
        coord_fixed() +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        labs(
            x = x_label,
            y = y_label,
            title = "NMF program alignment"
        )
}
