#' Plot heatmap of differentially expressed programs across groups
#'
#' Creates a heatmap showing log fold change of NMF program usage across cell groups
#' (e.g., cell types) using results from FindAllDEPs. Stars indicate programs with
#' log fold change above a threshold. Log fold change is calculated as log2(fold_enrichment)
#' from the differential expression analysis.
#'
#' @param deps_results SimpleList or list of DataFrames from FindAllDEPs()
#' @param programs Character vector of programs to plot. If NULL (default), plots all programs
#' @param logfc_threshold Numeric, log2 fold change threshold for marking with stars (default 1)
#' @param color_palette Character vector of length 3 for low, mid, high colors
#'   (default c("blue", "white", "red"))
#' @param cluster_rows Logical, whether to cluster programs (default FALSE)
#' @param cluster_cols Logical, whether to cluster groups (default FALSE)
#' @param show_rownames Logical, whether to show program names (default TRUE)
#' @param show_colnames Logical, whether to show group names (default TRUE)
#' @param star_size Numeric, size of star markers (default 8)
#' @param ... Additional arguments passed to pheatmap::pheatmap
#'
#' @return A pheatmap object
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
#' # Find differentially expressed programs
#' deps <- FindAllDEPs(sce, cell_type_col = "celltype")
#'
#' # Plot enrichment heatmap
#' plotDEPsHeatmap(deps)
#'
#' # Plot with custom threshold
#' plotDEPsHeatmap(deps, logfc_threshold = 0.5)
#'
#' # Plot specific programs
#' plotDEPsHeatmap(deps, programs = c("NMF_1", "NMF_2", "NMF_3"))
plotDEPsHeatmap <- function(deps_results, programs = NULL,
                            logfc_threshold = 1,
                            color_palette = c("dodgerblue", "white", "red"),
                            cluster_rows = FALSE, cluster_cols = FALSE,
                            show_rownames = TRUE, show_colnames = TRUE,
                            star_size = 14, ...) {

    # Check if pheatmap is available
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required for this function. ",
             "Install it with: install.packages('pheatmap')")
    }

    # Validate input - accept both SimpleList and regular list
    if (length(deps_results) == 0) {
        stop("deps_results must be a non-empty output from FindAllDEPs()")
    }

    # Get cell types and all programs
    cell_types <- names(deps_results)
    all_programs <- rownames(deps_results[[1]])

    # Select programs to plot
    if (is.null(programs)) {
        programs <- all_programs
    } else {
        if (!all(programs %in% all_programs)) {
            stop("Some programs not found in deps_results")
        }
    }

    # Build log fold change matrix from deps_results
    logfc_matrix <- matrix(NA, nrow = length(programs), ncol = length(cell_types))
    rownames(logfc_matrix) <- programs
    colnames(logfc_matrix) <- cell_types

    for (ct in cell_types) {
        deps_df <- deps_results[[ct]]

        for (prog in programs) {
            # Calculate log2 fold change from fold_enrichment
            fold_enrich <- deps_df[prog, "fold_enrichment"]

            # Handle infinite fold enrichment
            if (is.infinite(fold_enrich)) {
                logfc_matrix[prog, ct] <- ifelse(fold_enrich > 0, 10, -10)
            } else {
                logfc_matrix[prog, ct] <- log2(fold_enrich)
            }
        }
    }

    # Create color palette
    color_fun <- grDevices::colorRampPalette(color_palette)(100)

    # Determine breaks for symmetric color scale around 0
    max_abs_logfc <- max(abs(logfc_matrix), na.rm = TRUE)
    breaks <- seq(-max_abs_logfc, max_abs_logfc, length.out = 101)

    # Create display labels with stars for high enrichment
    display_matrix <- logfc_matrix
    for (i in seq_len(nrow(logfc_matrix))) {
        for (j in seq_len(ncol(logfc_matrix))) {
            if (!is.na(logfc_matrix[i, j]) && abs(logfc_matrix[i, j]) >= logfc_threshold) {
                display_matrix[i, j] <- logfc_matrix[i, j]
            }
        }
    }

    # Plot heatmap
    p <- pheatmap::pheatmap(
        logfc_matrix,
        color = color_fun,
        breaks = breaks,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        show_rownames = show_rownames,
        show_colnames = show_colnames,
        fontsize = 10,
        fontsize_row = if (length(programs) > 20) 8 else 10,
        fontsize_col = if (length(cell_types) > 20) 8 else 10,
        display_numbers = matrix(
            ifelse(!is.na(logfc_matrix) & logfc_matrix >= logfc_threshold, "*", ""),
            nrow = nrow(logfc_matrix),
            ncol = ncol(logfc_matrix)
        ),
        number_color = "black",
        fontsize_number = star_size,
        legend = TRUE,
        ...
    )

    return(p)
}
