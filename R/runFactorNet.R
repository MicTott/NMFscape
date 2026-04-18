#' Store a fitted FactorNet result in a SingleCellExperiment
#'
#' Power-user function for storing results from a custom RcppML FactorNet
#' graph in standard NMFscape slots. Users build and fit their own graph,
#' then use this function for SCE integration.
#'
#' @param x A SingleCellExperiment or SpatialExperiment object
#' @param result A fitted FactorNet result (output of \code{fit(factor_net(...))})
#' @param layer_name Character, name of the layer in the result to use as
#'   the primary output
#' @param name Character, name for the reducedDim slot (default "FactorNet")
#' @param absorb_d Logical, whether to absorb diagonal scaling (default TRUE)
#' @param verbose Logical, whether to print progress (default TRUE)
#'
#' @return The input object with results stored in:
#'   \itemize{
#'     \item \code{reducedDim(x, name)}: cells x k coefficient matrix
#'     \item \code{metadata(x)[[paste0(name, "_basis")]]}: features x k basis
#'       (or per-modality if multi-modal)
#'     \item \code{metadata(x)[[paste0(name, "_model")]]}: raw FactorNet result
#'   }
#'
#' @examples
#' \dontrun{
#' library(RcppML)
#' library(scuttle)
#' sce <- mockSCE(ngenes = 100, ncells = 50)
#' sce <- logNormCounts(sce)
#'
#' # Build custom graph
#' mat <- assay(sce, "logcounts")
#' inp <- factor_input(mat)
#' layer <- nmf_layer(inp, k = 5, name = "my_nmf")
#' net <- factor_net(inputs = list(inp), output = layer)
#' result <- fit(net)
#'
#' # Store in SCE
#' sce <- runFactorNet(sce, result, layer_name = "my_nmf")
#' }
#'
#' @export
runFactorNet <- function(x, result, layer_name, name = "FactorNet",
                         absorb_d = TRUE, verbose = TRUE) {

    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment or SpatialExperiment object")
    }

    available <- names(result$layers)
    if (!layer_name %in% available) {
        stop("Layer '", layer_name, "' not found in result. ",
             "Available layers: ", paste(available, collapse = ", "))
    }

    layer <- result$layers[[layer_name]]
    w_data <- layer$W
    h_mat <- layer$H
    d_vec <- layer$d
    k <- length(d_vec)

    # Handle multi-modal case (W is a list of per-modality matrices)
    if (is.list(w_data) && !is.matrix(w_data)) {
        if (absorb_d) {
            sqrt_d <- sqrt(d_vec)
            coeff_matrix <- t(diag(sqrt_d, nrow = k) %*% h_mat)
        } else {
            coeff_matrix <- t(h_mat)
        }

        factor_names <- paste0("NMF_", seq_len(k))
        rownames(coeff_matrix) <- colnames(x)
        colnames(coeff_matrix) <- factor_names
        reducedDim(x, name) <- coeff_matrix

        modality_names <- names(w_data)
        for (mod in modality_names) {
            w_mat <- w_data[[mod]]
            if (absorb_d) {
                basis <- w_mat %*% diag(sqrt_d, nrow = k)
            } else {
                basis <- w_mat
            }
            colnames(basis) <- factor_names
            metadata(x)[[paste0(name, "_basis_", mod)]] <- basis
        }
        metadata(x)[[paste0(name, "_modalities")]] <- modality_names
    } else {
        # Standard single-input case
        absorbed <- .absorbDiagonal(w_data, h_mat, d_vec, absorb_d)
        named <- .setFactorNames(absorbed$basis, absorbed$coeff,
                                 rownames(w_data), colnames(x), k)
        reducedDim(x, name) <- named$coeff
        metadata(x)[[paste0(name, "_basis")]] <- named$basis
    }

    metadata(x)[[paste0(name, "_model")]] <- result

    if (verbose) {
        message("FactorNet results stored in reducedDim(x, '", name, "')")
    }

    return(x)
}
