#' Run batch-conditioned NMF on a SingleCellExperiment
#'
#' Performs NMF while factoring out a conditioning variable (e.g., batch,
#' donor, sample) so that learned programs capture biological variation
#' rather than technical effects. Uses
#' \code{\link[RcppML]{factor_condition}} internally.
#'
#' @param x A SingleCellExperiment or SpatialExperiment object
#' @param k Integer, number of factors
#' @param condition_col Character, column name in \code{colData(x)} for the
#'   conditioning variable. Factor/character values are one-hot encoded;
#'   numeric values are used directly.
#' @param assay Character, which assay to use (default "logcounts")
#' @param name Character, name for the reducedDim slot (default "CondNMF")
#' @param subset_row Vector specifying which features to use
#' @param tol Numeric, tolerance for convergence (default 1e-5)
#' @param maxit Integer, maximum iterations (default 100)
#' @param L1 Numeric vector of length 2, L1 regularization (default c(0,0))
#' @param L2 Numeric vector of length 2, L2 regularization (default c(0,0))
#' @param distribution Character, loss function (default "mse")
#' @param absorb_d Logical, whether to absorb diagonal scaling (default TRUE)
#' @param seed Integer, random seed for reproducibility
#' @param verbose Logical, whether to print progress (default TRUE)
#' @param ... Additional arguments passed to \code{\link[RcppML]{factor_config}}
#'
#' @return The input object with conditioned NMF results stored in:
#'   \itemize{
#'     \item \code{reducedDim(x, name)}: cells x k coefficient matrix
#'     \item \code{metadata(x)[[paste0(name, "_basis")]]}: genes x k basis
#'     \item \code{metadata(x)[[paste0(name, "_model")]]}: raw FactorNet result
#'     \item \code{metadata(x)[[paste0(name, "_condition")]]}: list with the
#'       condition column name and model matrix used
#'   }
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#' sce$batch <- sample(c("A", "B"), ncol(sce), replace = TRUE)
#' sce <- runConditionedNMF(sce, k = 5, condition_col = "batch")
#'
#' @export
#' @importFrom RcppML factor_input factor_condition nmf_layer factor_net
#'   factor_config fit
runConditionedNMF <- function(x, k, condition_col, assay = "logcounts",
                              name = "CondNMF", subset_row = NULL,
                              tol = 1e-5, maxit = 100,
                              L1 = c(0, 0), L2 = c(0, 0),
                              distribution = "mse", absorb_d = TRUE,
                              seed = NULL, verbose = TRUE, ...) {

    .validateSCE(x, assay)

    if (!condition_col %in% colnames(colData(x))) {
        stop("condition_col '", condition_col, "' not found in colData(x)")
    }

    condition_var <- colData(x)[[condition_col]]

    if (any(is.na(condition_var))) {
        stop("condition_col '", condition_col, "' contains NA values; ",
             "please remove or impute before running conditioned NMF")
    }

    # Build condition model matrix
    if (is.character(condition_var) || is.factor(condition_var)) {
        condition_var <- as.factor(condition_var)
        if (nlevels(condition_var) < 2) {
            warning("All cells have the same condition value; ",
                    "conditioning will have no effect")
        }
        Z <- stats::model.matrix(~ condition_var - 1)
        colnames(Z) <- levels(condition_var)
    } else if (is.numeric(condition_var)) {
        if (length(unique(condition_var)) < 2) {
            warning("All cells have the same condition value; ",
                    "conditioning will have no effect")
        }
        Z <- matrix(condition_var, ncol = 1)
        colnames(Z) <- condition_col
    } else {
        stop("condition_col must be character, factor, or numeric")
    }

    feature_names <- .subsetRowNames(x, subset_row)
    mat <- .extractAssayMatrix(x, assay, subset_row)

    if (verbose) {
        message("Running conditioned NMF with k=", k,
                " (conditioning on '", condition_col, "')...")
    }

    # Build FactorNet graph
    input <- RcppML::factor_input(mat, name = "input")
    conditioned <- RcppML::factor_condition(input, Z)
    output <- RcppML::nmf_layer(conditioned, k = k, name = "nmf")

    config <- RcppML::factor_config(
        tol = tol, maxit = maxit, loss = distribution,
        seed = seed, verbose = FALSE, ...
    )
    net <- RcppML::factor_net(inputs = list(input), output = output,
                              config = config)
    result <- RcppML::fit(net)

    # Extract results
    layer <- result$layers$nmf
    absorbed <- .absorbDiagonal(layer$W, layer$H, layer$d, absorb_d)

    named <- .setFactorNames(absorbed$basis, absorbed$coeff,
                             feature_names, colnames(x), k)

    # Store
    reducedDim(x, name) <- named$coeff
    metadata(x)[[paste0(name, "_basis")]] <- named$basis
    metadata(x)[[paste0(name, "_model")]] <- result
    metadata(x)[[paste0(name, "_condition")]] <- list(
        col = condition_col, Z = Z
    )

    if (verbose) {
        message("Conditioned NMF completed. Results stored in ",
                "reducedDim(x, '", name, "')")
    }

    return(x)
}
