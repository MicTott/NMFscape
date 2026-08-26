#' Run label-guided (semi-supervised) NMF on a SingleCellExperiment
#'
#' Fits NMF while steering the cell embedding \code{H} with known labels.
#' Set \code{mode = "enrich"} to pull \code{H} toward the label structure
#' (semi-supervised factorization) or \code{mode = "remove"} to suppress
#' structure correlated with the labels (adversarial removal, e.g. batch).
#'
#' Guidance is expressed as a target \code{H} built from class centroids by
#' \code{\link[RcppML]{compute_target}}, which needs an embedding to start
#' from. \code{runGuidedNMF} therefore fits \strong{twice}: an unguided base
#' model, then \code{compute_target()} on its \code{H}, then the guided model.
#' Expect roughly double the runtime of \code{\link{runNMFscape}}. Set
#' \code{init_name} to also keep the base fit, which makes the guided/unguided
#' comparison free rather than costing a third fit.
#'
#' Cells whose label is \code{NA} are kept in the factorization but receive a
#' zero target column, so they are left unguided. This is what makes the
#' method semi-supervised: label the cells you are confident about, pass
#' \code{NA} for the rest, and let those cells float.
#'
#' Guidance is only available through the single-layer
#' \code{\link[RcppML]{nmf}} path. The per-factor \code{W()}/\code{H()}
#' configuration objects used by \code{\link{runDeepNMF}} and friends accept
#' \code{target} arguments but silently ignore them, so there is no guided
#' equivalent for the multi-layer recipes and no \code{W}-side (marker gene)
#' anchoring.
#'
#' @param x A SingleCellExperiment or SpatialExperiment object
#' @param k Integer, number of factors for NMF (rank)
#' @param label_col Character, column name in \code{colData(x)} holding the
#'   guiding labels. May be character or factor. \code{NA} entries are
#'   retained and left unguided.
#' @param mode Character, either "enrich" (attract \code{H} toward the label
#'   structure) or "remove" (adversarially suppress label-correlated
#'   structure). Guidance strength is always given as a positive
#'   \code{strength}; \code{mode} chooses the sign internally.
#' @param strength Numeric scalar > 0, guidance strength. Passed to RcppML as
#'   \code{target_lambda = strength} for "enrich" and
#'   \code{target_lambda = -strength} for "remove" (default 0.5).
#' @param assay Character or integer, which assay to use (default "logcounts")
#' @param name Character, name for the reducedDim slot (default "GuidedNMF")
#' @param subset_row Vector specifying which features to use
#' @param whiten Logical, whether \code{\link[RcppML]{compute_target}} applies
#'   OAS-shrinkage ZCA whitening when building class centroids (default TRUE)
#' @param init_name Character, optional name under which to also store the
#'   unguided base fit that the target was built from (default NULL, discard it)
#' @param absorb_d Logical, whether to absorb the diagonal scaling factor into
#'   W and H symmetrically (default TRUE)
#' @param seed Integer, random seed for reproducibility. Both the base and the
#'   guided fit use it, so the two differ only by the guidance.
#' @param verbose Logical, whether to print progress (default TRUE)
#' @param ... Additional arguments passed to both calls to
#'   \code{\link[RcppML]{nmf}}, including \code{tol}, \code{maxit}, \code{L1},
#'   \code{L2}, \code{loss} and \code{mask}
#'
#' @return The input object with guided NMF results stored in:
#'   \itemize{
#'     \item \code{reducedDim(x, name)}: cells x k coefficient matrix
#'     \item \code{metadata(x)[[paste0(name, "_basis")]]}: genes x k basis matrix
#'     \item \code{metadata(x)[[paste0(name, "_model")]]}: raw RcppML S4 nmf object
#'     \item \code{metadata(x)[[paste0(name, "_guidance")]]}: list recording
#'       \code{label_col}, \code{mode}, \code{target_lambda}, \code{whiten}
#'       and the number of guided cells
#'   }
#'   When \code{init_name} is supplied the unguided base fit is stored under
#'   the same three slot names built from \code{init_name}.
#'
#' @seealso \code{\link{refineNMF}} for post hoc correction of an already
#'   fitted model, and \code{\link{runConditionedNMF}} for factoring a
#'   covariate out by conditioning rather than adversarially.
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 100)
#' sce <- logNormCounts(sce)
#' sce$celltype <- rep(c("TypeA", "TypeB"), length.out = ncol(sce))
#'
#' # Semi-supervised: pull programs toward the annotation
#' sce <- runGuidedNMF(sce, k = 5, label_col = "celltype", verbose = FALSE)
#' head(reducedDim(sce, "GuidedNMF"))
#'
#' # Adversarial: suppress a nuisance covariate
#' sce$batch <- rep(c("b1", "b2"), length.out = ncol(sce))
#' sce <- runGuidedNMF(sce, k = 5, label_col = "batch", mode = "remove",
#'                     name = "DebatchedNMF", verbose = FALSE)
#'
#' @export
#' @importFrom RcppML nmf compute_target
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom S4Vectors metadata<-
runGuidedNMF <- function(x, k, label_col, mode = c("enrich", "remove"),
                         strength = 0.5, assay = "logcounts",
                         name = "GuidedNMF", subset_row = NULL,
                         whiten = TRUE, init_name = NULL, absorb_d = TRUE,
                         seed = NULL, verbose = TRUE, ...) {

    .validateSCE(x, assay)

    mode <- match.arg(mode)

    if (!is.numeric(strength) || length(strength) != 1 || is.na(strength) ||
        strength <= 0) {
        stop("strength must be a single positive number. Use ",
             "mode = \"remove\" for adversarial guidance rather than a ",
             "negative strength.")
    }

    if (!label_col %in% colnames(colData(x))) {
        stop("label_col '", label_col, "' not found in colData(x)")
    }

    labels <- colData(x)[[label_col]]

    if (!is.character(labels) && !is.factor(labels)) {
        stop("label_col '", label_col, "' must be a character or factor ",
             "column; guided NMF groups cells into classes")
    }

    labels <- factor(as.character(labels))
    n_guided <- sum(!is.na(labels))

    if (nlevels(labels) < 2) {
        stop("label_col '", label_col, "' must have at least 2 non-NA ",
             "levels; found ", nlevels(labels))
    }

    if (n_guided < 2) {
        stop("label_col '", label_col, "' has fewer than 2 non-NA labels; ",
             "there is nothing to guide toward")
    }

    target_lambda <- if (mode == "enrich") strength else -strength

    feature_names <- .subsetRowNames(x, subset_row)
    mat <- .extractAssayMatrix(x, assay, subset_row)

    if (verbose) {
        message("Fitting unguided base NMF with k=", k,
                " to initialise the guidance target...")
    }

    base_fit <- RcppML::nmf(data = mat, k = k, seed = seed,
                            verbose = FALSE, ...)

    if (verbose) {
        message("Building guidance target from '", label_col, "' (",
                n_guided, "/", ncol(x), " cells labelled, ",
                nlevels(labels), " classes)...")
    }

    target_h <- RcppML::compute_target(base_fit@h, labels, whiten = whiten)

    if (verbose) {
        message("Running guided NMF (mode=", mode,
                ", target_lambda=", target_lambda, ")...")
    }

    guided_fit <- RcppML::nmf(data = mat, k = k, seed = seed, verbose = FALSE,
                              target_H = target_h,
                              target_lambda = target_lambda, ...)

    x <- .storeNMFFit(x, guided_fit, name, feature_names, k, absorb_d)

    metadata(x)[[paste0(name, "_guidance")]] <- list(
        label_col = label_col,
        mode = mode,
        target_lambda = target_lambda,
        whiten = whiten,
        n_guided = n_guided
    )

    if (!is.null(init_name)) {
        x <- .storeNMFFit(x, base_fit, init_name, feature_names, k, absorb_d)
        if (verbose) {
            message("Unguided base fit stored in reducedDim(x, '",
                    init_name, "')")
        }
    }

    if (verbose) {
        message("Guided NMF completed. Results stored in reducedDim(x, '",
                name, "')")
    }

    return(x)
}

# Write an RcppML S4 nmf fit into the standard NMFscape SCE slots.
# @keywords internal
# @noRd
.storeNMFFit <- function(x, fit, name, feature_names, k, absorb_d) {
    absorbed <- .absorbDiagonal(fit@w, fit@h, fit@d, absorb_d)
    named <- .setFactorNames(absorbed$basis, absorbed$coeff,
                             feature_names, colnames(x), k)

    reducedDim(x, name) <- named$coeff
    metadata(x)[[paste0(name, "_basis")]] <- named$basis
    metadata(x)[[paste0(name, "_model")]] <- fit
    x
}
