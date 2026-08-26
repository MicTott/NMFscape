#' Gene-set enrichment of NMF programs
#'
#' Tests each NMF program against a collection of gene sets. A program's
#' column of the basis matrix is already a ranking of every measured gene, so
#' the default method runs preranked GSEA (\code{fgsea}) on that ranking and no
#' arbitrary cutoff is needed. The cruder alternative, \code{method = "ora"},
#' runs a hypergeometric test on the program's top \code{n_top} genes against
#' the other measured genes as background.
#'
#' The background for over-representation is the set of genes present in the
#' basis matrix, not the full gene-set collection. Testing against genes that
#' were never measured inflates significance.
#'
#' @param x A SingleCellExperiment object with NMF results
#' @param gene_sets A named list of character vectors of gene identifiers, or a
#'   data.frame in the \pkg{msigdbr} shape (a gene-set name column such as
#'   \code{gs_name} and a gene column such as \code{gene_symbol})
#' @param name Character, name of the NMF result to use (default "NMF")
#' @param method Character, "gsea" (preranked, default) or "ora"
#'   (over-representation on the top \code{n_top} genes)
#' @param modality Character, modality of a multi-modal result (default NULL)
#' @param programs Integer or character vector, programs to test
#'   (default NULL for all)
#' @param n_top Integer, number of top genes per program used by
#'   \code{method = "ora"} (default 200, ignored by "gsea")
#' @param min_size Integer, minimum gene-set size after intersecting with the
#'   measured genes (default 10)
#' @param max_size Integer, maximum gene-set size after intersecting with the
#'   measured genes (default 500)
#' @param ... Additional arguments passed to \code{fgsea::fgsea}, for example
#'   \code{BPPARAM} or \code{eps}. \code{scoreType} defaults to "pos"
#'   because NMF loadings are non-negative and the ranking is therefore
#'   one-sided
#'
#' @return A data.frame with one row per program and gene set, sorted by
#'   program and then by p-value. Columns are \code{program}, \code{pathway},
#'   \code{pval}, \code{padj} and \code{size}, plus \code{ES}, \code{NES} and
#'   \code{leading_edge} for \code{method = "gsea"}, or \code{overlap},
#'   \code{expected}, \code{odds_ratio} and \code{overlap_genes} for
#'   \code{method = "ora"}.
#' @export
#' @examples
#' library(scuttle)
#' set.seed(1)
#' sce <- logNormCounts(mockSCE(ngenes = 200, ncells = 100))
#' sce <- runNMFscape(sce, k = 3, verbose = FALSE)
#'
#' # gene sets are supplied by the user; here two hand-made ones
#' genes <- rownames(sce)
#' gene_sets <- list(SET_A = genes[1:30], SET_B = genes[31:60])
#'
#' res <- enrichPrograms(sce, gene_sets, min_size = 5)
#' head(res)
#'
#' # over-representation on the top 50 genes instead
#' ora <- enrichPrograms(sce, gene_sets, method = "ora", n_top = 50,
#'                       min_size = 5)
#' head(ora)
enrichPrograms <- function(x, gene_sets, name = "NMF",
                           method = c("gsea", "ora"), modality = NULL,
                           programs = NULL, n_top = 200, min_size = 10,
                           max_size = 500, ...) {

    method <- match.arg(method)

    basis <- getBasis(x, name, modality = modality)
    if (is.null(rownames(basis))) {
        stop("The basis matrix has no rownames, so its genes cannot be ",
             "matched to 'gene_sets'.")
    }
    basis <- .selectPrograms(basis, programs)

    gene_sets <- .asGeneSetList(gene_sets)
    gene_sets <- .filterGeneSets(gene_sets, rownames(basis), min_size,
                                 max_size)

    program_names <- colnames(basis)
    if (is.null(program_names)) {
        program_names <- paste0("Program_", seq_len(ncol(basis)))
    }

    if (length(gene_sets) == 0L) {
        return(.emptyEnrichment(method))
    }

    res <- lapply(seq_len(ncol(basis)), function(j) {
        if (method == "gsea") {
            out <- .gseaProgram(basis[, j], gene_sets, min_size, max_size,
                                ...)
        } else {
            out <- .oraProgram(basis[, j], gene_sets, n_top)
        }
        if (nrow(out) == 0L) {
            return(NULL)
        }
        cbind(program = program_names[j], out, stringsAsFactors = FALSE)
    })

    res <- do.call(rbind, Filter(Negate(is.null), res))
    if (is.null(res) || nrow(res) == 0L) {
        return(.emptyEnrichment(method))
    }

    res <- res[order(factor(res$program, levels = program_names), res$pval), ,
               drop = FALSE]
    rownames(res) <- NULL
    res
}

# Subset a basis matrix by program index or program name.
.selectPrograms <- function(basis, programs) {
    if (is.null(programs)) {
        return(basis)
    }
    if (is.numeric(programs)) {
        if (any(programs < 1 | programs > ncol(basis))) {
            stop("Program indices must be between 1 and ", ncol(basis))
        }
    } else if (is.character(programs)) {
        missing_programs <- setdiff(programs, colnames(basis))
        if (length(missing_programs) > 0L) {
            stop("Programs not found: ",
                 paste(missing_programs, collapse = ", "))
        }
    } else {
        stop("'programs' must be a numeric or character vector")
    }
    basis[, programs, drop = FALSE]
}

# Accept either a named list or an msigdbr-style data.frame.
.asGeneSetList <- function(gene_sets) {
    if (is.data.frame(gene_sets)) {
        set_col <- intersect(c("gs_name", "gs_id", "set", "pathway"),
                             colnames(gene_sets))
        gene_col <- intersect(c("gene_symbol", "db_gene_symbol", "entrez_gene",
                                "ensembl_gene", "gene"), colnames(gene_sets))
        if (length(set_col) == 0L || length(gene_col) == 0L) {
            stop("A data.frame 'gene_sets' needs a gene-set name column ",
                 "(e.g. 'gs_name') and a gene column (e.g. 'gene_symbol').")
        }
        gene_sets <- split(as.character(gene_sets[[gene_col[1]]]),
                           as.character(gene_sets[[set_col[1]]]))
    }
    if (!is.list(gene_sets) || is.null(names(gene_sets))) {
        stop("'gene_sets' must be a named list of character vectors or an ",
             "msigdbr-style data.frame.")
    }
    lapply(gene_sets, function(genes) unique(as.character(genes)))
}

# Restrict sets to the measured genes and drop those outside the size range.
.filterGeneSets <- function(gene_sets, measured, min_size, max_size) {
    gene_sets <- lapply(gene_sets, function(genes) intersect(genes, measured))
    sizes <- lengths(gene_sets)
    keep <- sizes >= min_size & sizes <= max_size
    if (sum(keep) < 3L) {
        warning(sum(keep), " of ", length(gene_sets), " gene sets have ",
                min_size, "-", max_size, " genes in common with the ",
                length(measured), " measured genes. Check that 'gene_sets' ",
                "uses the same gene identifiers as the basis matrix, or ",
                "relax 'min_size'/'max_size'.")
    }
    gene_sets[keep]
}

.emptyEnrichment <- function(method) {
    base <- data.frame(program = character(0), pathway = character(0),
                       pval = numeric(0), padj = numeric(0),
                       size = integer(0), stringsAsFactors = FALSE)
    if (method == "gsea") {
        return(cbind(base, ES = numeric(0), NES = numeric(0),
                     leading_edge = character(0), stringsAsFactors = FALSE))
    }
    cbind(base, overlap = integer(0), expected = numeric(0),
          odds_ratio = numeric(0), overlap_genes = character(0),
          stringsAsFactors = FALSE)
}

# Preranked GSEA on one program's gene loadings.
.gseaProgram <- function(loadings, gene_sets, min_size, max_size, ...) {
    if (!requireNamespace("fgsea", quietly = TRUE)) {
        stop("Package 'fgsea' is required for method = 'gsea'. ",
             "Install it with: BiocManager::install('fgsea')")
    }
    args <- list(...)
    # NMF loadings are non-negative, so the ranking is one-sided and fgsea's
    # two-sided default cannot estimate p-values for depleted sets.
    if (is.null(args$scoreType)) {
        args$scoreType <- "pos"
    }
    res <- do.call(fgsea::fgsea,
                   c(list(pathways = gene_sets, stats = loadings,
                          minSize = min_size, maxSize = max_size), args))
    res <- as.data.frame(res, stringsAsFactors = FALSE)
    if (nrow(res) == 0L) {
        return(.emptyEnrichment("gsea")[, -1, drop = FALSE])
    }
    leading_edge <- vapply(res$leadingEdge,
                           function(genes) paste(genes, collapse = ","),
                           character(1))
    data.frame(pathway = as.character(res$pathway), pval = res$pval,
               padj = res$padj, size = as.integer(res$size), ES = res$ES,
               NES = res$NES, leading_edge = leading_edge,
               stringsAsFactors = FALSE)
}

# Hypergeometric over-representation of one program's top genes, using the
# measured genes as the background.
.oraProgram <- function(loadings, gene_sets, n_top) {
    universe <- names(loadings)
    n_universe <- length(universe)
    n_top <- min(n_top, n_universe)
    top_genes <- universe[order(loadings, decreasing = TRUE)[seq_len(n_top)]]

    sizes <- lengths(gene_sets)
    overlap_genes <- lapply(gene_sets, function(genes)
        intersect(genes, top_genes))
    overlap <- lengths(overlap_genes)

    pval <- stats::phyper(overlap - 1L, sizes, n_universe - sizes, n_top,
                          lower.tail = FALSE)
    expected <- n_top * sizes / n_universe
    odds_ratio <- (overlap / (n_top - overlap)) /
        ((sizes - overlap) / (n_universe - sizes - n_top + overlap))

    data.frame(pathway = names(gene_sets), pval = pval,
               padj = stats::p.adjust(pval, method = "BH"),
               size = as.integer(sizes), overlap = as.integer(overlap),
               expected = expected, odds_ratio = odds_ratio,
               overlap_genes = vapply(overlap_genes,
                                      function(genes)
                                          paste(genes, collapse = ","),
                                      character(1)),
               stringsAsFactors = FALSE)
}

#' Dot plot of program enrichment results
#'
#' Plots the top gene sets per program from \code{enrichPrograms()}. Point size
#' is the adjusted p-value on a -log10 scale and colour is the normalized
#' enrichment score (GSEA) or the odds ratio (over-representation).
#'
#' @param enrichment A data.frame returned by \code{enrichPrograms()}
#' @param top_n Integer, number of gene sets to show per program (default 5)
#' @param padj_cutoff Numeric, keep only gene sets with an adjusted p-value at
#'   or below this (default 0.05)
#' @param facet Logical, give each program its own panel (default FALSE)
#'
#' @return A ggplot object
#' @export
#' @examples
#' library(scuttle)
#' set.seed(1)
#' sce <- logNormCounts(mockSCE(ngenes = 200, ncells = 100))
#' sce <- runNMFscape(sce, k = 3, verbose = FALSE)
#'
#' genes <- rownames(sce)
#' gene_sets <- list(SET_A = genes[1:30], SET_B = genes[31:60])
#' res <- enrichPrograms(sce, gene_sets, min_size = 5)
#' plotProgramEnrichment(res, padj_cutoff = 1)
plotProgramEnrichment <- function(enrichment, top_n = 5, padj_cutoff = 0.05,
                                  facet = FALSE) {

    if (!is.data.frame(enrichment) ||
        !all(c("program", "pathway", "padj") %in% colnames(enrichment))) {
        stop("'enrichment' must be a data.frame returned by enrichPrograms()")
    }

    keep <- !is.na(enrichment$padj) & enrichment$padj <= padj_cutoff
    if (!any(keep)) {
        stop("No gene set has an adjusted p-value at or below ", padj_cutoff,
             ". Raise 'padj_cutoff' to plot the results anyway.")
    }
    enrichment <- enrichment[keep, , drop = FALSE]

    enrichment <- do.call(rbind, lapply(
        split(enrichment, enrichment$program), function(df) {
            df[order(df$padj)[seq_len(min(top_n, nrow(df)))], , drop = FALSE]
        }))

    enrichment$neg_log10_padj <- -log10(pmax(enrichment$padj,
                                             .Machine$double.xmin))
    use_nes <- "NES" %in% colnames(enrichment)
    enrichment$colour_by <- if (use_nes) enrichment$NES else
        enrichment$odds_ratio

    p <- ggplot(enrichment, aes(x = program, y = reorder(pathway, colour_by),
                                size = neg_log10_padj, color = colour_by)) +
        geom_point() +
        scale_size_continuous(name = "-log10(padj)") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Program", y = NULL,
             color = if (use_nes) "NES" else "Odds ratio")

    if (use_nes) {
        p <- p + scale_color_gradient2(low = "#2166ac", mid = "grey90",
                                       high = "#b2182b", midpoint = 0)
    } else {
        p <- p + scale_color_viridis_c()
    }

    if (facet) {
        p <- p + facet_wrap(~program, scales = "free_y")
    }

    p
}
