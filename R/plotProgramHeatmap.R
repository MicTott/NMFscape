#' Plot heatmap of top genes per NMF program
#'
#' Creates a heatmap showing the weights of top-ranked genes for each NMF program.
#' Genes are ranked by their loadings in the basis matrix (W) and displayed with
#' their weights across all programs.
#'
#' @param x A SingleCellExperiment object with NMF results
#' @param nmf_name Character, name of NMF result to use (default "NMF")
#' @param programs Integer vector, specific programs to plot (default NULL for all).
#'   Can specify program indices (e.g., c(1, 3, 5)) or names (e.g., c("GEP_1", "GEP_3"))
#' @param n_genes Integer, number of top genes to show per program (default 10)
#' @param make_unique Logical, whether to select unique top genes per program (default TRUE).
#'   When TRUE, genes are iteratively assigned to programs based on highest weight,
#'   ensuring each gene appears only once. When FALSE, genes can appear in multiple programs.
#' @param scale_rows Logical, whether to scale gene weights across programs (default TRUE)
#' @param cluster_rows Logical, whether to cluster genes (default TRUE)
#' @param cluster_cols Logical, whether to cluster programs (default FALSE)
#' @param show_rownames Logical, whether to show gene names (default TRUE)
#' @param color_palette Character vector of length 3 for low, mid, high colors
#'   (default c("#2166ac", "white", "#b2182b") - a diverging blue-white-red palette)
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
#' # Plot top 10 unique genes per program (default)
#' plotProgramHeatmap(sce, n_genes = 10)
#'
#' # Allow shared genes across programs
#' plotProgramHeatmap(sce, n_genes = 10, make_unique = FALSE)
#'
#' # Plot top 5 genes without scaling
#' plotProgramHeatmap(sce, n_genes = 5, scale_rows = FALSE)
#'
#' # Custom colors
#' plotProgramHeatmap(sce, n_genes = 8,
#'                    color_palette = c("navy", "white", "darkred"))
#'
#' # Plot only specific programs
#' plotProgramHeatmap(sce, programs = c(1, 3, 5), n_genes = 10)
plotProgramHeatmap <- function(x, nmf_name = "NMF", programs = NULL, n_genes = 10, make_unique = TRUE,
                               scale_rows = TRUE, cluster_rows = TRUE,
                               cluster_cols = FALSE, show_rownames = TRUE,
                               color_palette = c("dodgerblue", "white", "red"), ...) {

    # Check if pheatmap is available
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required for this function. ",
             "Install it with: install.packages('pheatmap')")
    }

    # Get basis matrix
    basis_name <- paste0(nmf_name, "_basis")
    if (!basis_name %in% names(metadata(x))) {
        stop("NMF result '", nmf_name, "' not found. ",
             "Run runNMFscape() first.")
    }

    basis <- metadata(x)[[basis_name]]
    k <- ncol(basis)

    # Subset programs if requested
    if (!is.null(programs)) {
        # Handle both numeric indices and character names
        if (is.numeric(programs)) {
            if (any(programs < 1 | programs > k)) {
                stop("Program indices must be between 1 and ", k)
            }
            basis <- basis[, programs, drop = FALSE]
        } else if (is.character(programs)) {
            if (!all(programs %in% colnames(basis))) {
                missing <- setdiff(programs, colnames(basis))
                stop("Programs not found: ", paste(missing, collapse = ", "))
            }
            basis <- basis[, programs, drop = FALSE]
        } else {
            stop("'programs' must be numeric or character vector")
        }
        k <- ncol(basis)
    }

    if (make_unique) {
        # Select unique top genes per program
        # Strategy: iteratively assign genes to programs based on highest weight
        all_top_genes <- character(0)
        genes_per_program <- list()
        available_genes <- rownames(basis)

        for (prog_idx in seq_len(k)) {
            # Get weights for this program from available genes
            prog_weights <- basis[available_genes, prog_idx]

            # Select top n_genes
            n_to_select <- min(n_genes, length(available_genes))
            top_idx <- order(prog_weights, decreasing = TRUE)[seq_len(n_to_select)]
            selected_genes <- available_genes[top_idx]

            # Store genes for this program
            genes_per_program[[colnames(basis)[prog_idx]]] <- selected_genes
            all_top_genes <- c(all_top_genes, selected_genes)

            # Remove selected genes from available pool
            available_genes <- setdiff(available_genes, selected_genes)
        }

        top_genes_list <- genes_per_program
    } else {
        # Get top genes for each program (allows overlap)
        # Use the subsetted basis matrix
        top_genes_list <- lapply(seq_len(k), function(prog_idx) {
            prog_weights <- basis[, prog_idx]
            top_idx <- order(prog_weights, decreasing = TRUE)[seq_len(n_genes)]
            rownames(basis)[top_idx]
        })
        names(top_genes_list) <- colnames(basis)

        # Create set of top genes across all programs
        all_top_genes <- unique(unlist(top_genes_list))
    }

    # Extract weights for these genes
    gene_weights <- basis[all_top_genes, , drop = FALSE]

    # Scale rows if requested (z-score across programs)
    if (scale_rows) {
        gene_weights <- t(scale(t(gene_weights)))
    }

    # Create color palette for heatmap
    color_fun <- grDevices::colorRampPalette(color_palette)(100)

    # Generate distinct colors for each program annotation using polychrome palette
    n_programs <- ncol(basis)
    if (!requireNamespace("pals", quietly = TRUE)) {
        stop("Package 'pals' is required for this function. ",
             "Install it with: install.packages('pals')")
    }

    # Use polychrome palette for maximally distinct colors
    if (n_programs <= 36) {
        program_colors <- pals::polychrome(n = n_programs)
    } else {
        # For >36 programs, supplement with alphabet2
        base_colors <- pals::polychrome(n = 36)
        extra_colors <- pals::alphabet2(n = n_programs - 36)
        program_colors <- c(base_colors, extra_colors)
    }
    names(program_colors) <- colnames(basis)

    # Create row annotation for which program each gene is top-ranked in
    gene_annotations <- data.frame(
        Program = factor(
            vapply(rownames(gene_weights), function(gene) {
                # Find which program this gene is top-ranked in
                prog_idx <- which(vapply(top_genes_list, function(x) gene %in% x,
                                        FUN.VALUE = logical(1)))
                if (length(prog_idx) > 1) {
                    # If gene is top in multiple programs, pick the one with highest weight
                    prog_idx <- which.max(basis[gene, ])
                }
                colnames(basis)[prog_idx]
            }, FUN.VALUE = character(1)),
            levels = colnames(basis)
        ),
        row.names = rownames(gene_weights)
    )

    # Color list for row annotation
    annotation_colors <- list(
        Program = program_colors
    )

    # Plot heatmap
    pheatmap::pheatmap(
        gene_weights,
        color = color_fun,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        show_rownames = show_rownames,
        show_colnames = TRUE,
        annotation_row = gene_annotations,
        annotation_colors = annotation_colors,
        fontsize_row = if (length(all_top_genes) > 50) 6 else 8,
        ...
    )
}
