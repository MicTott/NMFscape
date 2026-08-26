#' Diagnose the fit of an NMF model
#'
#' Runs the model diagnostics that RcppML computes but that NMFscape has not
#' previously surfaced, in one call: whether the residual dispersion is
#' structured by gene or by cell, whether the data carry more zeros than the
#' fitted distribution predicts, which variance-power family best matches the
#' observed mean-variance relationship, and how sparse each program is.
#'
#' The distribution score test compares the power-variance family
#' V(mu) = mu^p without refitting, so it is a cheap second opinion on the
#' \code{distribution} argument of \code{\link{runNMFscape}}: it is most
#' meaningful for a model fit with \code{distribution = "mse"}, whose
#' residuals carry no distributional assumption of their own.
#'
#' Note that RcppML's \code{variance_explained()} is deliberately not reported.
#' Despite being a generic, RcppML 1.0.0 defines a method only for
#' \code{svd} objects, not for \code{nmf} models, so there is nothing to wrap.
#' Use \code{\link{evaluateNMF}} for reconstruction loss instead.
#'
#' @param x A SingleCellExperiment or SpatialExperiment object with NMF results
#'   from \code{\link{runNMFscape}}, \code{\link{consensusNMF}} or
#'   \code{\link{refineNMF}}, which are the functions that store an S4
#'   \code{nmf} model. The FactorNet recipes store a \code{factor_net_result}
#'   instead and are not supported here.
#' @param name Character, name of the NMF result to diagnose (default "NMF").
#' @param assay Character, which assay the NMF was run on
#'   (default "logcounts").
#' @param cv_threshold Numeric, coefficient-of-variation threshold above which
#'   dispersion is declared structured by row or by column (default 0.5).
#' @param zi_threshold Numeric, minimum excess zero fraction needed to declare
#'   zero inflation (default 0.05).
#' @param powers Numeric vector of variance powers to score
#'   (default c(0, 1, 2, 3), i.e. Gaussian, Poisson-like, Gamma and inverse
#'   Gaussian).
#' @param test_nb Logical, whether to also run the negative-binomial
#'   overdispersion diagnostic (default TRUE).
#' @param min_mu Numeric, floor applied to predicted values to avoid dividing
#'   by zero (default 1e-6).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{dispersion}}{From \code{\link[RcppML]{diagnose_dispersion}}:
#'       the recommended dispersion \code{mode}, the \code{global_phi}
#'       estimate, and the per-row and per-column CVs.}
#'     \item{\code{zero_inflation}}{From
#'       \code{\link[RcppML]{diagnose_zero_inflation}}: the
#'       \code{excess_zero_rate}, a \code{has_zi} flag, the recommended
#'       \code{zi_mode}, and per-row and per-column excess rates.}
#'     \item{\code{distribution}}{From
#'       \code{\link[RcppML]{score_test_distribution}}: a \code{scores}
#'       data.frame, the \code{best_power} and \code{best_distribution}, and
#'       the negative-binomial diagnostic when requested.}
#'     \item{\code{sparsity}}{From \code{\link[RcppML]{sparsity}}: a data.frame
#'       with the sparsity of each factor in the \code{w} and \code{h}
#'       matrices.}
#'     \item{\code{k}}{The rank of the diagnosed model.}
#'   }
#'
#' @seealso \code{\link{evaluateNMF}} for reconstruction loss,
#'   \code{\link{selectDistribution}} for BIC-based distribution selection.
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ngenes = 200, ncells = 100)
#' sce <- logNormCounts(sce)
#' sce <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)
#'
#' diag <- diagnoseNMF(sce)
#' diag$dispersion$mode
#' diag$distribution$best_distribution
#' head(diag$sparsity)
#'
#' @export
#' @importFrom RcppML diagnose_dispersion diagnose_zero_inflation
#'   score_test_distribution
#' @importMethodsFrom RcppML sparsity
diagnoseNMF <- function(x, name = "NMF", assay = "logcounts",
                        cv_threshold = 0.5, zi_threshold = 0.05,
                        powers = c(0, 1, 2, 3), test_nb = TRUE,
                        min_mu = 1e-6) {

    if (!is(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment or SpatialExperiment object")
    }

    if (!assay %in% assayNames(x)) {
        stop("assay '", assay, "' not found in x")
    }

    model <- getModel(x, name)

    if (.isFactorNet(model)) {
        stop("diagnoseNMF() requires an S4 nmf model, but '", name,
             "' was produced by a FactorNet recipe (runDeepNMF, ",
             "runMultiModalNMF, runConditionedNMF or runFactorNet). ",
             "RcppML's diagnostics take a single W and H, which a ",
             "multi-layer graph does not have. Per-layer W, d and H are ",
             "available as metadata(x)[['", name, "_model']]$layers, and the ",
             "fitted loss as $total_loss.")
    }

    mat <- assay(x, assay)

    if (any(mat < 0)) {
        mat[mat < 0] <- 0
    }

    if (nrow(model@w) != nrow(mat)) {
        stop("The model for '", name, "' was fit on ", nrow(model@w),
             " features but assay '", assay, "' has ", nrow(mat),
             " rows. If the fit used subset_row, diagnose it on the same ",
             "features.")
    }

    list(
        dispersion = RcppML::diagnose_dispersion(
            data = mat, model = model,
            cv_threshold = cv_threshold, min_mu = min_mu
        ),
        zero_inflation = RcppML::diagnose_zero_inflation(
            data = mat, model = model, threshold = zi_threshold
        ),
        distribution = RcppML::score_test_distribution(
            data = mat, model = model, powers = powers,
            test_nb = test_nb, min_mu = min_mu
        ),
        sparsity = RcppML::sparsity(model),
        k = ncol(model@w)
    )
}
