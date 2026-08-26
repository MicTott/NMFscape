#' Report GPU availability for the RcppML backend
#'
#' Reports whether the installed RcppML build can offload factorization to a
#' GPU, and what device it would use. GPU support depends on how RcppML was
#' compiled, so it is easy to miss that it is either available or absent; this
#' is the one call that answers the question.
#'
#' @param force_recheck Logical, whether to re-probe the device rather than
#'   reuse RcppML's cached answer (default FALSE).
#'
#' @return A list with \code{available}, a logical flag, and \code{device},
#'   the device information reported by \code{\link[RcppML]{gpu_info}} or
#'   \code{NULL} when no GPU is usable.
#'
#' @seealso \code{\link{runNMFscape}}, which notes GPU availability when
#'   \code{verbose = TRUE}.
#'
#' @examples
#' gpuInfo()
#'
#' @export
#' @importFrom RcppML gpu_available gpu_info
gpuInfo <- function(force_recheck = FALSE) {
    available <- isTRUE(RcppML::gpu_available(force_recheck = force_recheck))
    list(
        available = available,
        device = if (available) RcppML::gpu_info() else NULL
    )
}

# TRUE when the RcppML build can use a GPU; never errors, so verbose
# reporting cannot break a factorization.
# @noRd
.gpuAvailable <- function() {
    isTRUE(tryCatch(RcppML::gpu_available(), error = function(e) FALSE))
}
