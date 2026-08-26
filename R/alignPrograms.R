#' Align programs between two independently fit NMF models
#'
#' Compares two NMF factorizations that were fit separately and reports which
#' program in one run corresponds to which program in the other. This answers
#' "did these two runs find the same programs?", which is a different question
#' from \code{\link{predictNMF}}: that function projects new data into an
#' existing basis, so both datasets end up sharing one set of programs by
#' construction, whereas \code{alignPrograms} takes two sets of programs that
#' were discovered independently and matches them up after the fact.
#'
#' The two bases are compared on their shared features. Program similarity is
#' the column-wise cosine similarity (or Pearson correlation) of the basis
#' matrices, and the one-to-one assignment is the minimum-cost bipartite
#' matching of \code{1 - similarity} via
#' \code{\link[RcppML]{bipartiteMatch}} (the Hungarian algorithm), not a greedy
#' argmax. A greedy pass can assign the same program in \code{y} to several
#' programs in \code{x}; the bipartite matching cannot.
#'
#' The two models need not share a rank. When \code{y} has more programs than
#' \code{x}, every program in \code{x} is matched and the surplus programs in
#' \code{y} are left over. When \code{x} has more programs than \code{y}, the
#' surplus programs in \code{x} are reported with \code{matched = FALSE} and a
#' missing partner: those are candidate run-specific programs.
#'
#' @param x A SingleCellExperiment or SpatialExperiment object carrying NMF
#'   results, an S4 \code{nmf} model, or a features x programs basis matrix.
#' @param y A second object of any of those three kinds, fit independently of
#'   \code{x}.
#' @param name Character, name of the NMF result to pull from \code{x} when it
#'   is a SingleCellExperiment (default "NMF"). Ignored otherwise.
#' @param name_y Character, name of the NMF result to pull from \code{y} when
#'   it is a SingleCellExperiment. Defaults to \code{name}.
#' @param modality Character, optional modality name used when \code{x} is a
#'   multi-modal result from \code{\link{runMultiModalNMF}}.
#' @param modality_y Character, optional modality for \code{y}. Defaults to
#'   \code{modality}.
#' @param method Character, "cosine" (default) or "cor", the program-to-program
#'   similarity measure.
#' @param ... Further arguments, currently unused.
#'
#' @return An object of class \code{programAlignment}, a list with elements:
#'   \describe{
#'     \item{\code{mapping}}{A data.frame with one row per program in \code{x}
#'       and columns \code{program_x}, \code{program_y} (\code{NA} when
#'       unmatched), \code{similarity} and \code{matched}. Rows are ordered by
#'       decreasing similarity.}
#'     \item{\code{similarity}}{The full k_x by k_y similarity matrix, with
#'       program names as dimnames.}
#'     \item{\code{unmatched_y}}{Character vector of programs in \code{y} that
#'       no program in \code{x} was assigned to.}
#'     \item{\code{method}}{The similarity measure used.}
#'     \item{\code{features}}{Number of shared features the comparison used.}
#'     \item{\code{cost}}{Total cost of the bipartite matching.}
#'   }
#'
#' @seealso \code{\link{plotProgramSimilarity}} to draw the similarity matrix,
#'   \code{\link{predictNMF}} to project data into an existing basis instead.
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 200, ncells = 100)
#' sce <- logNormCounts(sce)
#'
#' # two independent fits of the same data
#' fit_a <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)
#' fit_b <- runNMFscape(sce, k = 4, seed = 2, verbose = FALSE)
#'
#' res <- alignPrograms(fit_a, fit_b)
#' res$mapping
#' dim(res$similarity)
#'
#' @export
#' @importFrom RcppML cosine bipartiteMatch
#' @importFrom stats cor
alignPrograms <- function(x, y, name = "NMF", name_y = NULL,
                          modality = NULL, modality_y = NULL,
                          method = c("cosine", "cor"), ...) {

    method <- match.arg(method)

    if (is.null(name_y)) {
        name_y <- name
    }
    if (is.null(modality_y)) {
        modality_y <- modality
    }

    w_x <- .programBasis(x, name, modality, "x")
    w_y <- .programBasis(y, name_y, modality_y, "y")

    shared <- .sharedFeatures(w_x, w_y)
    w_x <- shared$x
    w_y <- shared$y

    sim <- .programSimilarity(w_x, w_y, method)

    # bipartiteMatch minimizes cost and rejects negative entries, so the
    # similarity is turned into a distance and floored at zero: cosine can
    # return 1 + epsilon for identical programs, which would otherwise trip
    # a malformed error inside RcppML.
    cost <- pmax(1 - sim, 0)
    dim(cost) <- dim(sim)

    match_res <- RcppML::bipartiteMatch(cost)
    assignment <- match_res$assignment

    # assignment is 0-indexed, with -1 marking a program in x that has no
    # partner (only possible when x has more programs than y).
    partner <- ifelse(assignment >= 0, assignment + 1L, NA_integer_)

    # rows of sim index programs in x, columns index programs in y
    program_x <- rownames(sim)
    program_y <- colnames(sim)

    similarity <- vapply(
        seq_along(partner),
        function(i) {
            if (is.na(partner[i])) NA_real_ else sim[i, partner[i]]
        },
        FUN.VALUE = numeric(1)
    )

    mapping <- data.frame(
        program_x = program_x,
        program_y = ifelse(is.na(partner), NA_character_,
                           program_y[partner]),
        similarity = similarity,
        matched = !is.na(partner),
        stringsAsFactors = FALSE
    )
    mapping <- mapping[order(mapping$similarity, decreasing = TRUE,
                             na.last = TRUE), , drop = FALSE]
    rownames(mapping) <- NULL

    result <- list(
        mapping = mapping,
        similarity = sim,
        unmatched_y = setdiff(program_y, stats::na.omit(mapping$program_y)),
        method = method,
        features = nrow(w_x),
        cost = match_res$cost
    )
    class(result) <- "programAlignment"
    result
}

# Pull a features x programs basis out of whichever object the user supplied.
# @noRd
.programBasis <- function(obj, name, modality, arg) {
    if (is(obj, "SingleCellExperiment")) {
        return(as.matrix(getBasis(obj, name, modality)))
    }
    if (is(obj, "nmf")) {
        return(as.matrix(obj@w))
    }
    if (is.matrix(obj) || is(obj, "Matrix")) {
        return(as.matrix(obj))
    }
    stop("'", arg, "' must be a SingleCellExperiment, an S4 nmf model, or a ",
         "features x programs basis matrix, not an object of class '",
         class(obj)[1], "'.")
}

# Restrict two bases to their shared features and give programs stable names.
# @noRd
.sharedFeatures <- function(w_x, w_y) {
    if (ncol(w_x) < 1 || ncol(w_y) < 1) {
        stop("Both bases must contain at least one program.")
    }

    if (is.null(colnames(w_x))) {
        colnames(w_x) <- paste0("x_", seq_len(ncol(w_x)))
    }
    if (is.null(colnames(w_y))) {
        colnames(w_y) <- paste0("y_", seq_len(ncol(w_y)))
    }

    if (is.null(rownames(w_x)) || is.null(rownames(w_y))) {
        if (nrow(w_x) != nrow(w_y)) {
            stop("Cannot align programs: at least one basis has no rownames ",
                 "and the two bases have different numbers of features (",
                 nrow(w_x), " vs ", nrow(w_y), "). Supply bases with feature ",
                 "names so the shared features can be identified.")
        }
        warning("At least one basis has no rownames; assuming both bases ",
                "list the same features in the same order.")
        return(list(x = w_x, y = w_y))
    }

    common <- intersect(rownames(w_x), rownames(w_y))
    if (length(common) < 2) {
        stop("Cannot align programs: the two bases share ", length(common),
             " features. At least 2 shared features are required.")
    }
    if (length(common) < nrow(w_x) || length(common) < nrow(w_y)) {
        message("Aligning on ", length(common), " shared features (",
                nrow(w_x), " in x, ", nrow(w_y), " in y).")
    }

    list(x = w_x[common, , drop = FALSE], y = w_y[common, , drop = FALSE])
}

# Column-wise similarity between two bases, as a named k_x by k_y matrix.
# @noRd
.programSimilarity <- function(w_x, w_y, method) {
    sim <- if (identical(method, "cosine")) {
        RcppML::cosine(w_x, w_y)
    } else {
        stats::cor(w_x, w_y)
    }

    sim <- matrix(as.numeric(sim), nrow = ncol(w_x), ncol = ncol(w_y),
                  dimnames = list(colnames(w_x), colnames(w_y)))

    if (anyNA(sim)) {
        stop("The ", method, " similarity between the two bases contains ",
             "missing values, which usually means a program has zero ",
             "variance across the shared features. Drop the empty programs ",
             "or use more shared features.")
    }

    sim
}
