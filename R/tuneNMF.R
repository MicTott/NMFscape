#' Tune FactorNet hyperparameters by held-out cross-validation
#'
#' Generalises \code{\link{selectRank}} from "one rank, one layer" to any
#' NMFscape recipe and any set of \code{\link[RcppML]{nmf_layer}} parameters.
#' Wraps \code{\link[RcppML]{cross_validate_graph}}, which masks a random
#' fraction of matrix entries, fits the graph on the remainder and scores the
#' reconstruction of the masked entries -- the same held-out statistic
#' \code{selectRank()} already uses, applied to a whole FactorNet graph.
#'
#' The \code{recipe} presets build the FactorNet graph for you, using the same
#' topology as the corresponding \code{run*()} function, so the selected values
#' transfer directly:
#' \describe{
#'   \item{\code{"standard"}}{one \code{nmf_layer} on the assay matrix, as in
#'     \code{\link{runNMFscape}}.}
#'   \item{\code{"deep"}}{a chain of \code{nmf_layer} calls, as in
#'     \code{\link{runDeepNMF}}. See "Deep recipe" below.}
#'   \item{\code{"multimodal"}}{\code{\link[RcppML]{factor_shared}} across
#'     modalities feeding one \code{nmf_layer}, as in
#'     \code{\link{runMultiModalNMF}}.}
#'   \item{\code{"conditioned"}}{\code{\link[RcppML]{factor_condition}} feeding
#'     one \code{nmf_layer}, as in \code{\link{runConditionedNMF}}.}
#' }
#'
#' @section Parameter grids:
#' \code{params} is a named list of candidate values. Names are passed straight
#' through to \code{\link[RcppML]{nmf_layer}}, so anything that function accepts
#' can be tuned, e.g. \code{list(k = c(5, 10, 20), L1 = c(0, 0.01))}. Note that
#' FactorNet layers take a \emph{scalar} \code{L1}/\code{L2} per layer, unlike
#' \code{\link{runNMFscape}}, which takes a length-2 c(w, h) vector.
#'
#' @section Deep recipe:
#' \code{runDeepNMF()} takes a per-layer rank \emph{vector}, so a joint grid
#' would have to range over candidate vectors. That is not possible with
#' RcppML 1.0.0: multi-layer graphs diverge (training loss goes to
#' \code{NaN}) when \code{test_fraction > 0}, and every layer reports
#' \code{test_loss = 0}, so a joint held-out score cannot be computed. Rather
#' than return meaningless numbers, \code{tuneNMF()} tunes a deep chain
#' \strong{sequentially, one layer at a time}: layer 1 is cross-validated on the
#' assay matrix, the winning rank is fitted without masking, and layer 2 is then
#' cross-validated on the transpose of that layer's H -- exactly the matrix
#' layer 2 factorizes inside \code{runDeepNMF()}. Every fit scored this way is a
#' single layer, which is the case that works upstream. Accordingly,
#' \code{params$k} for the deep recipe is a \emph{list of per-layer candidate
#' grids}, one element per layer, e.g.
#' \code{params = list(k = list(c(10, 20, 30), c(3, 5, 8)))} searches three
#' ranks for layer 1 and three for layer 2. Non-\code{k} parameters are selected
#' at layer 1 and then held fixed for the deeper layers, matching the scalar
#' \code{L1}/\code{L2} that \code{runDeepNMF()} applies to every layer. Because
#' the search is greedy rather than joint, cost is additive in the per-layer
#' grid sizes rather than multiplicative.
#'
#' Penalties do not rescale between layers: layer 1 sees the assay matrix while
#' layer 2 sees the layer above's H, which is typically orders of magnitude
#' smaller. A penalty that suits layer 1 can therefore collapse a deeper layer,
#' in which case \code{tuneNMF()} stops with an error naming the carried-over
#' settings rather than reporting a degenerate fit. Tuning \code{k} on its own
#' avoids this entirely.
#'
#' @param x A SingleCellExperiment or SpatialExperiment object
#' @param params Named list of candidate parameter values passed to
#'   \code{\link[RcppML]{nmf_layer}}, e.g. \code{list(k = c(5, 10, 20))}. For
#'   \code{recipe = "deep"}, \code{params$k} must instead be a list of per-layer
#'   candidate grids (see "Deep recipe").
#' @param recipe Character, which graph topology to tune. One of "standard",
#'   "deep", "multimodal" or "conditioned". Ignored when \code{layer_fn} is
#'   supplied.
#' @param assay Character, which assay to use (default "logcounts"). For
#'   \code{recipe = "multimodal"} a character vector of two or more assay names
#'   selects same-feature-space modalities.
#' @param subset_row Vector specifying which features to use
#' @param reps Integer, number of replicate fits per parameter combination
#'   (default 3). Replicates use different held-out masks and give the standard
#'   errors reported in \code{$summary}.
#' @param strategy Character, "grid" to evaluate every combination or "random"
#'   to sample \code{n_random} of them (default "grid")
#' @param n_random Integer, number of combinations to sample when
#'   \code{strategy = "random"} (default 20)
#' @param test_fraction Numeric, fraction of matrix entries held out
#'   (default 0.2)
#' @param condition_col Character, column of \code{colData(x)} to condition on.
#'   Required for \code{recipe = "conditioned"}.
#' @param alt_exps Character vector of altExp names, optionally including
#'   "main" for the main experiment. Used by \code{recipe = "multimodal"}.
#' @param modality_names Character vector of short modality names. Defaults to
#'   the assay or altExp names.
#' @param layer_fn Optional escape hatch. A function taking \code{(params)}, or
#'   \code{(params, inputs)}, and returning the output
#'   \code{\link[RcppML]{nmf_layer}} node of a custom graph. When supplied it
#'   replaces the \code{recipe} topology; \code{recipe} then only decides which
#'   matrices become inputs. \code{inputs} is the list of
#'   \code{\link[RcppML]{factor_input}} nodes \code{tuneNMF()} built from
#'   \code{x}, and a custom graph must be assembled from those nodes.
#' @param seed Integer, random seed for the held-out masks (default 42)
#' @param verbose Logical, whether to report progress (default TRUE)
#' @param ... Additional arguments passed to
#'   \code{\link[RcppML]{factor_config}}, e.g. \code{tol}, \code{maxit},
#'   \code{loss}, \code{mask_zeros}
#'
#' @return An object of class \code{nmfscape_tuning} (which inherits from
#'   \code{factor_net_cv}) with elements:
#'   \itemize{
#'     \item \code{summary}: data.frame of parameter combinations ranked by
#'       \code{mean_test_loss}, with \code{se_test_loss}, \code{mean_train_loss}
#'       and \code{n_valid}. For \code{recipe = "deep"} it also has a
#'       \code{layer} column.
#'     \item \code{results}: per-replicate data.frame with \code{test_loss},
#'       \code{train_loss}, \code{iterations} and \code{converged}.
#'     \item \code{best_params}: named list of the winning values, shaped so it
#'       can be passed straight to the matching \code{run*()} function, e.g.
#'       \code{do.call(runDeepNMF, c(list(x = sce), fit$best_params))}.
#'     \item \code{recipe}, \code{params}, \code{reps}, \code{strategy},
#'       \code{n_fits}, \code{config}.
#'   }
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 80, ncells = 60)
#' sce <- logNormCounts(sce)
#'
#' # Standard single-layer recipe: tune k
#' tuned <- tuneNMF(sce, params = list(k = c(2, 4)), reps = 2, verbose = FALSE)
#' tuned$summary
#' sce <- do.call(runNMFscape, c(list(x = sce), tuned$best_params))
#'
#' # Deep recipe: one candidate grid per layer
#' deep_tuned <- tuneNMF(sce, recipe = "deep",
#'                       params = list(k = list(c(6, 8), c(2, 3))),
#'                       reps = 1, verbose = FALSE)
#' deep_tuned$best_params$k
#'
#' @seealso \code{\link{plotTuning}} to visualise the result,
#'   \code{\link{selectRank}} for the single-layer special case,
#'   \code{\link[RcppML]{cross_validate_graph}} for the underlying engine.
#'
#' @export
#' @importFrom RcppML cross_validate_graph factor_input factor_shared
#'   factor_condition nmf_layer factor_net factor_config fit
#' @importFrom utils modifyList
tuneNMF <- function(x, params,
                    recipe = c("standard", "deep", "multimodal", "conditioned"),
                    assay = "logcounts", subset_row = NULL, reps = 3,
                    strategy = c("grid", "random"), n_random = 20,
                    test_fraction = 0.2, condition_col = NULL,
                    alt_exps = NULL, modality_names = NULL,
                    layer_fn = NULL, seed = 42, verbose = TRUE, ...) {

    recipe <- match.arg(recipe)
    strategy <- match.arg(strategy)

    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment or SpatialExperiment object")
    }

    reps <- as.integer(reps)
    if (length(reps) != 1 || is.na(reps) || reps < 1) {
        stop("reps must be a single positive integer")
    }
    if (length(test_fraction) != 1 || !is.numeric(test_fraction) ||
        is.na(test_fraction) || test_fraction <= 0 || test_fraction >= 1) {
        stop("test_fraction must be a single number strictly between 0 and 1")
    }

    use_custom <- !is.null(layer_fn)
    if (use_custom && !is.function(layer_fn)) {
        stop("layer_fn must be a function of (params) or (params, inputs)")
    }

    .checkTuneParams(params, recipe = if (use_custom) "standard" else recipe)

    config <- .tuneConfig(test_fraction = test_fraction, seed = seed, ...)

    if (recipe == "deep" && !use_custom) {
        return(.tuneDeep(x = x, params = params, assay = assay,
                         subset_row = subset_row, config = config, reps = reps,
                         strategy = strategy, n_random = n_random,
                         seed = seed, verbose = verbose))
    }

    inputs <- .tuneInputs(x, recipe = recipe, assay = assay,
                          subset_row = subset_row, alt_exps = alt_exps,
                          modality_names = modality_names)

    if (use_custom) {
        fn <- if (length(formals(layer_fn)) >= 2) {
            function(p) layer_fn(p, inputs)
        } else {
            layer_fn
        }
    } else {
        node <- switch(recipe,
            standard = inputs[[1]],
            multimodal = do.call(RcppML::factor_shared, inputs),
            conditioned = RcppML::factor_condition(
                inputs[[1]], .tuneConditionMatrix(x, condition_col))
        )
        layer_name <- switch(recipe, standard = "nmf",
                             multimodal = "shared_nmf", conditioned = "nmf")
        fn <- function(p) {
            do.call(RcppML::nmf_layer,
                    c(list(node, name = layer_name), as.list(p)))
        }
    }

    n_fits <- .tuneFitCount(params, reps, strategy, n_random)
    .reportCost(n_fits, verbose)

    cv <- RcppML::cross_validate_graph(
        inputs = inputs, layer_fn = fn, params = params, config = config,
        reps = reps, strategy = strategy, n_random = n_random,
        seed = as.integer(seed), verbose = FALSE
    )

    if (all(is.na(cv$summary$mean_test_loss)) ||
        isTRUE(all(cv$summary$mean_test_loss == 0))) {
        stop("Cross-validation produced no usable held-out loss (every fit ",
             "returned test_loss = 0 or NA). This graph cannot be ",
             "cross-validated with the installed version of RcppML.")
    }

    best <- .tuneBestParams(cv$best_params, recipe = recipe,
                            use_custom = use_custom)

    if (verbose) {
        message("Best: ", .describeParams(best), " (mean test loss = ",
                signif(min(cv$summary$mean_test_loss, na.rm = TRUE), 5), ")")
    }

    .tuningObject(summary = cv$summary, results = cv$results,
                  best_params = best, recipe = recipe, params = params,
                  reps = reps, strategy = strategy, n_fits = n_fits,
                  config = config,
                  param_names = names(params))
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Validate the parameter grid. Deep grids carry a list of per-layer vectors
# under `k`; every other recipe wants plain atomic vectors throughout.
# @noRd
.checkTuneParams <- function(params, recipe) {
    if (!is.list(params) || length(params) == 0) {
        stop("params must be a non-empty named list of candidate values, ",
             "e.g. list(k = c(5, 10, 20))")
    }
    if (is.null(names(params)) || any(!nzchar(names(params))) ||
        anyDuplicated(names(params)) > 0) {
        stop("All elements of params must have unique, non-empty names")
    }

    if (recipe == "deep") {
        if (!"k" %in% names(params)) {
            stop("params must include 'k' for recipe = \"deep\"")
        }
        k_grids <- params$k
        if (!is.list(k_grids) || length(k_grids) < 2) {
            stop("For recipe = \"deep\", params$k must be a list with one ",
                 "candidate grid per layer and at least 2 layers, e.g. ",
                 "list(k = list(c(10, 20, 30), c(3, 5, 8))). Use ",
                 "recipe = \"standard\" for single-layer NMF.")
        }
        for (i in seq_along(k_grids)) {
            g <- k_grids[[i]]
            if (!is.numeric(g) || length(g) == 0 || anyNA(g) || any(g < 1)) {
                stop("params$k[[", i, "]] must be a non-empty numeric vector ",
                     "of positive ranks")
            }
        }
        others <- params[setdiff(names(params), "k")]
    } else {
        others <- params
    }

    for (nm in names(others)) {
        v <- others[[nm]]
        if (!is.atomic(v) || length(v) == 0 || all(is.na(v))) {
            stop("params$", nm, " must be a non-empty atomic vector of ",
                 "candidate values")
        }
    }
    invisible(TRUE)
}

# factor_config() with NMFscape's defaults, overridable through `...`.
# @noRd
.tuneConfig <- function(test_fraction, seed, ...) {
    defaults <- list(tol = 1e-5, maxit = 100, loss = "mse", verbose = FALSE,
                     seed = seed, test_fraction = test_fraction)
    do.call(RcppML::factor_config, modifyList(defaults, list(...)))
}

# Build the factor_input node(s) a recipe factorizes.
# @noRd
.tuneInputs <- function(x, recipe, assay, subset_row, alt_exps,
                        modality_names) {
    if (recipe != "multimodal") {
        .validateSCE(x, assay)
        mat <- .extractAssayMatrix(x, assay, subset_row)
        return(list(RcppML::factor_input(mat, name = "input")))
    }

    if (!is.null(alt_exps)) {
        if (length(alt_exps) < 2) {
            stop("recipe = \"multimodal\" requires at least 2 modalities in ",
                 "alt_exps")
        }
        if (is.null(modality_names)) modality_names <- alt_exps
        mats <- lapply(alt_exps, function(ae) {
            if (ae == "main") {
                .validateSCE(x, assay[1])
                .extractAssayMatrix(x, assay[1], subset_row)
            } else {
                if (!ae %in% altExpNames(x)) {
                    stop("altExp '", ae, "' not found in x")
                }
                ae_sce <- altExp(x, ae)
                which_assay <- if ("counts" %in% assayNames(ae_sce)) {
                    "counts"
                } else {
                    1L
                }
                mat <- assay(ae_sce, which_assay)
                if (any(mat < 0)) {
                    warning("Negative values in altExp '", ae,
                            "'. Setting to 0.")
                    mat[mat < 0] <- 0
                }
                mat
            }
        })
    } else {
        if (length(assay) < 2) {
            stop("recipe = \"multimodal\" needs either alt_exps, or two or ",
                 "more assay names in `assay`")
        }
        if (is.null(modality_names)) modality_names <- assay
        mats <- lapply(assay, function(a) {
            .validateSCE(x, a)
            .extractAssayMatrix(x, a, subset_row)
        })
    }

    if (length(modality_names) != length(mats)) {
        stop("modality_names must have one entry per modality")
    }
    n_cells <- vapply(mats, ncol, integer(1))
    if (length(unique(n_cells)) != 1) {
        stop("All modalities must have the same number of cells. Found: ",
             paste(paste0(modality_names, "=", n_cells), collapse = ", "))
    }

    lapply(seq_along(mats), function(i) {
        RcppML::factor_input(mats[[i]], name = modality_names[i])
    })
}

# One-hot / numeric design matrix for the conditioning variable. Mirrors the
# encoding runConditionedNMF() uses so the tuned k transfers.
# @noRd
.tuneConditionMatrix <- function(x, condition_col) {
    if (is.null(condition_col)) {
        stop("condition_col is required for recipe = \"conditioned\"")
    }
    if (!condition_col %in% colnames(colData(x))) {
        stop("condition_col '", condition_col, "' not found in colData(x)")
    }
    condition_var <- colData(x)[[condition_col]]
    if (any(is.na(condition_var))) {
        stop("condition_col '", condition_col, "' contains NA values")
    }
    if (is.character(condition_var) || is.factor(condition_var)) {
        condition_var <- as.factor(condition_var)
        Z <- stats::model.matrix(~ condition_var - 1)
        colnames(Z) <- levels(condition_var)
    } else if (is.numeric(condition_var)) {
        Z <- matrix(condition_var, ncol = 1)
        colnames(Z) <- condition_col
    } else {
        stop("condition_col must be character, factor, or numeric")
    }
    Z
}

# @noRd
.tuneFitCount <- function(params, reps, strategy, n_random) {
    n_combos <- prod(vapply(params, length, integer(1)))
    if (strategy == "random") n_combos <- min(n_combos, n_random)
    as.integer(n_combos) * reps
}

# @noRd
.reportCost <- function(n_fits, verbose, prefix = "Cross-validating ") {
    if (verbose) {
        message(prefix, n_fits, " model fits...")
    }
    if (n_fits > 100) {
        warning(n_fits, " model fits requested. This may take a long time; ",
                "consider a coarser grid, fewer reps, or ",
                "strategy = \"random\".", call. = FALSE)
    }
    invisible(n_fits)
}

# @noRd
.describeParams <- function(p) {
    paste(vapply(names(p), function(nm) {
        paste0(nm, " = ", paste(p[[nm]], collapse = ", "))
    }, character(1)), collapse = "; ")
}

# runNMFscape() takes length-2 c(w, h) penalties while nmf_layer() takes one
# scalar, so widen them for the standard recipe to keep best_params callable.
# @noRd
.tuneBestParams <- function(best, recipe, use_custom) {
    if (use_custom || recipe != "standard") return(best)
    for (nm in intersect(c("L1", "L2"), names(best))) {
        if (length(best[[nm]]) == 1) best[[nm]] <- rep(best[[nm]], 2)
    }
    best
}

# @noRd
.tuningObject <- function(summary, results, best_params, recipe, params, reps,
                          strategy, n_fits, config, param_names) {
    rownames(summary) <- NULL
    rownames(results) <- NULL
    structure(
        list(summary = summary, results = results, best_params = best_params,
             recipe = recipe, params = params, reps = reps,
             strategy = strategy, n_fits = n_fits, config = config,
             param_names = param_names),
        class = c("nmfscape_tuning", "factor_net_cv")
    )
}

# Sequential, layer-by-layer cross-validation of a deep chain.
#
# RcppML 1.0.0 cannot score a multi-layer graph: with test_fraction > 0 the
# chain diverges (train loss -> NaN) and every layer reports test_loss = 0.
# Single layers are fine, so each layer is cross-validated on its own input:
# layer 1 on the assay matrix, layer i > 1 on t(H) of the layer above, which is
# exactly the matrix runDeepNMF() feeds it.
# @noRd
.tuneDeep <- function(x, params, assay, subset_row, config, reps, strategy,
                      n_random, seed, verbose) {

    .validateSCE(x, assay)
    mat <- .extractAssayMatrix(x, assay, subset_row)

    k_grids <- lapply(params$k, function(g) sort(unique(as.integer(g))))
    n_layers <- length(k_grids)
    others <- params[setdiff(names(params), "k")]

    n_fits <- .tuneFitCount(c(list(k = k_grids[[1]]), others), reps,
                            strategy, n_random)
    for (i in seq_len(n_layers)[-1]) {
        n_fits <- n_fits + length(k_grids[[i]]) * reps
    }
    .reportCost(n_fits, verbose,
                prefix = paste0("Cross-validating ", n_layers,
                                " layers sequentially, "))

    full_config <- config
    full_config$test_fraction <- 0

    summaries <- vector("list", n_layers)
    all_results <- vector("list", n_layers)
    chosen_k <- integer(n_layers)
    fixed <- others[FALSE]
    current <- mat

    for (i in seq_len(n_layers)) {
        grid <- k_grids[[i]]
        max_k <- min(dim(current)) - 1L
        keep <- grid <= max_k
        if (!any(keep)) {
            stop("No candidate rank for layer ", i, " fits its input matrix ",
                 "(", nrow(current), " x ", ncol(current), "). Supply ranks ",
                 "below ", max_k + 1L, ".")
        }
        if (!all(keep) && verbose) {
            message("Layer ", i, ": dropping rank(s) ",
                    paste(grid[!keep], collapse = ", "),
                    " that exceed the input dimension")
        }
        grid <- grid[keep]

        layer_params <- c(list(k = grid),
                          if (i == 1L) others else fixed)
        node <- RcppML::factor_input(current, name = "input")
        layer_name <- paste0("layer_", i)
        fn <- function(p) {
            do.call(RcppML::nmf_layer,
                    c(list(node, name = layer_name), as.list(p)))
        }

        if (verbose) {
            message("Layer ", i, ": cross-validating k = ",
                    paste(grid, collapse = ", "), " on a ", nrow(current),
                    " x ", ncol(current), " matrix")
        }

        cv <- RcppML::cross_validate_graph(
            inputs = list(node), layer_fn = fn, params = layer_params,
            config = config, reps = reps, strategy = strategy,
            n_random = n_random, seed = as.integer(seed + 1000L * i),
            verbose = FALSE
        )

        if (all(is.na(cv$summary$mean_test_loss)) ||
            isTRUE(all(cv$summary$mean_test_loss == 0))) {
            reason <- if (i > 1L && length(fixed) > 0) {
                paste0("The layer-1 settings carried over here (",
                       .describeParams(fixed), ") are the likely cause: this ",
                       "layer factorizes the layer above's H, which is on a ",
                       "much smaller scale than the assay, so a penalty ",
                       "tuned on the assay can collapse it. Try smaller ",
                       "penalties, or tune k on its own.")
            } else {
                paste0("This layer cannot be cross-validated with the ",
                       "installed version of RcppML.")
            }
            stop("Layer ", i, " produced no usable held-out loss (every fit ",
                 "returned NA or zero). ", reason)
        }

        summaries[[i]] <- cbind(layer = i, cv$summary)
        all_results[[i]] <- cbind(layer = i, cv$results)

        k_i <- as.integer(cv$best_params$k)
        chosen_k[i] <- k_i
        if (i == 1L) {
            fixed <- cv$best_params[setdiff(names(cv$best_params), "k")]
        }

        if (i < n_layers) {
            out <- do.call(RcppML::nmf_layer,
                           c(list(node, k = k_i, name = layer_name), fixed))
            res <- RcppML::fit(RcppML::factor_net(inputs = list(node),
                                                  output = out,
                                                  config = full_config))
            current <- t(res$layers[[layer_name]]$H)
        }
    }

    if (any(diff(chosen_k) >= 0)) {
        warning("Selected ranks are not strictly decreasing (k = ",
                paste(chosen_k, collapse = ", "),
                "); hierarchical factorization normally needs k to shrink ",
                "with depth.", call. = FALSE)
    }

    summary_df <- do.call(rbind, summaries)
    results_df <- do.call(rbind, all_results)
    best <- c(list(k = chosen_k), fixed)

    if (verbose) {
        message("Best: ", .describeParams(best))
    }

    .tuningObject(summary = summary_df, results = results_df,
                  best_params = best, recipe = "deep", params = params,
                  reps = reps, strategy = strategy, n_fits = n_fits,
                  config = config, param_names = names(params))
}
