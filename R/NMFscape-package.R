#' NMFscape: Non-negative Matrix Factorization for Single Cell and Spatial Data
#'
#' High-performance NMF for \code{SingleCellExperiment} and
#' \code{SpatialExperiment} objects, built on the RcppML backend.
#'
#' @section Known RcppML limitations:
#' NMFscape is developed against RcppML 1.0.0. Several arguments in that
#' version are documented but do not do what they say, and NMFscape works
#' around them rather than passing them through:
#' \itemize{
#'   \item \code{graph_W}/\code{graph_H} do not perform graph-aware smoothing
#'     (a randomly rewired graph changes the fit as much as the true one), so
#'     \code{\link{runSpatialNMF}} implements the penalty itself.
#'   \item The per-factor \code{W()}/\code{H()} \code{target} argument is
#'     ignored by \code{fit()}, so guidance is only available through
#'     \code{\link{runGuidedNMF}}, which uses the single-layer path.
#'   \item Adversarial target removal collapses the factorization, so
#'     \code{runGuidedNMF(mode = "remove")} is disabled.
#'   \item \code{predict()} fails on multi-layer graphs, so
#'     \code{\link{predictNMF}} refuses them with an explanation.
#'   \item Multi-layer graphs report a held-out loss of exactly zero, so
#'     \code{\link{tuneNMF}} refuses to cross-validate them.
#'   \item \code{bipartiteMatch()} errors on negative costs, so
#'     \code{\link{alignPrograms}} floors its cost matrix at zero.
#'   \item \code{variance_explained()} has no method for \code{nmf} objects, so
#'     \code{\link{diagnoseNMF}} omits it.
#' }
#' Loading NMFscape with a newer RcppML emits a startup message, because these
#' workarounds become incorrect if the underlying behaviour changes.
#'
#' @keywords internal
"_PACKAGE"
