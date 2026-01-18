#' Plot dot plot of NMF program weights by group
#'
#' Creates a dot plot showing NMF program usage across cell groups (e.g., cell types).
#' Color represents the sum program weight in each group, and dot size represents
#' the percentage of cells with non-zero weights.
#'
#' @param x A SingleCellExperiment object with NMF results
#' @param group Character, column name in colData(x) defining cell groups (e.g., "celltype")
#' @param nmf_name Character, name of NMF result to use (default "NMF")
#' @param programs Integer vector of programs to plot. If NULL (default), plots all programs
#' @param scale Logical, whether to scale program weights to 0-1 per program (default TRUE)
#' @param min_pct Numeric, minimum percentage of cells (0-100) required to show a dot (default 0)
#' @param color_palette Character vector of length 2 for low and high colors
#'   (default c("lightgrey", "blue"))
#' @param dot_scale Numeric, scaling factor for dot sizes (default 6)
#'
#' @return A ggplot object
#' @export
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 500, ncells = 200)
#' sce <- logNormCounts(sce)
#' sce <- runNMFscape(sce, k = 5, verbose = FALSE)
#'
#' # Add mock cell type labels
#' sce$celltype <- sample(c("TypeA", "TypeB", "TypeC"), ncol(sce), replace = TRUE)
#'
#' # Plot all programs across cell types
#' plotProgramDots(sce, group = "celltype")
#'
#' # Plot specific programs
#' plotProgramDots(sce, group = "celltype", programs = c(1, 3, 5))
#'
#' # Without scaling
#' plotProgramDots(sce, group = "celltype", scale = FALSE)
plotProgramDots <- function(x, group, nmf_name = "NMF", programs = NULL,
                            scale = FALSE, min_pct = 0,
                            color_palette = c("lightgrey", "blue"),
                            dot_scale = 6) {

    # Check if group exists in colData
    if (!group %in% names(colData(x))) {
        stop("Group '", group, "' not found in colData(x)")
    }

    # Get NMF coefficients from reducedDim (cells x programs)
    if (!nmf_name %in% reducedDimNames(x)) {
        stop("NMF result '", nmf_name, "' not found. ",
             "Run runNMFscape() first.")
    }

    coeffs <- reducedDim(x, nmf_name)  # cells x programs
    k <- ncol(coeffs)

    # Select programs to plot
    if (is.null(programs)) {
        programs <- seq_len(k)
    } else {
        if (is.character(programs)) {
            # Convert program names to indices
            programs <- match(programs, colnames(coeffs))
        }
        if (any(is.na(programs)) || any(programs < 1) || any(programs > k)) {
            stop("Invalid program selection")
        }
    }

    # Get group labels
    groups <- as.character(colData(x)[[group]])

    # Calculate statistics for each group and program
    results_list <- list()

    for (prog_idx in programs) {
        prog_weights <- coeffs[, prog_idx]  # Extract column (all cells for this program)

        # Scale if requested
        if (scale) {
            prog_weights <- (prog_weights - min(prog_weights)) /
                           (max(prog_weights) - min(prog_weights))
        }

        # Calculate sum and percent non-zero for each group
        for (grp in unique(groups)) {
            grp_cells <- groups == grp
            grp_weights <- prog_weights[grp_cells]

            sum_weight <- sum(grp_weights, na.rm = TRUE)
            pct_nonzero <- sum(grp_weights > 0) / length(grp_weights) * 100

            results_list[[length(results_list) + 1]] <- data.frame(
                Program = colnames(coeffs)[prog_idx],
                Group = grp,
                sumWeight = sum_weight,
                PctNonzero = pct_nonzero,
                stringsAsFactors = FALSE
            )
        }
    }

    # Combine into single data frame
    plot_data <- do.call(rbind, results_list)

    # Filter by minimum percentage
    plot_data <- plot_data[plot_data$PctNonzero >= min_pct, ]

    # Convert to factor to control ordering
    plot_data$Program <- factor(plot_data$Program,
                                levels = colnames(coeffs)[programs])
    plot_data$Group <- factor(plot_data$Group)

    # Create dot plot
    p <- ggplot(plot_data, aes(x = Program, y = Group)) +
        geom_point(aes(size = PctNonzero, color = sumWeight)) +
        scale_size_continuous(name = "Percent Non-zero", range = c(0, dot_scale)) +
        scale_color_gradient(
            low = "white",
            high = "black",
            name = if (scale) "Scaled\nWeight" else "Sum \nWeight"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_blank()
        ) +
        labs(
            x = "NMF Program",
            y = group,
            title = "NMF Program Usage by Group"
        )

    return(p)
}
