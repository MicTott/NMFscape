# The version of RcppML this package has been tested against. Several RcppML
# arguments are documented but non-functional in 1.0.0, and NMFscape works
# around them: runSpatialNMF() implements its own graph penalty because
# graph_W/graph_H do not smooth, alignPrograms() floors the bipartiteMatch()
# cost because negative entries error, tuneNMF() refuses multi-layer graphs
# because their held-out loss is identically zero, and runGuidedNMF() disables
# mode = "remove". If those are fixed upstream the workarounds become wrong,
# so a newer RcppML should be checked rather than silently trusted.
.RCPPML_TESTED <- "1.0.0"

.onLoad <- function(libname, pkgname) {
    installed <- tryCatch(utils::packageVersion("RcppML"),
                          error = function(e) NULL)
    if (!is.null(installed) && installed > package_version(.RCPPML_TESTED)) {
        packageStartupMessage(
            "NMFscape was tested against RcppML ", .RCPPML_TESTED,
            " but ", as.character(installed), " is installed. NMFscape works ",
            "around several arguments that are non-functional in ",
            .RCPPML_TESTED, "; see ?NMFscape-package. Please verify results ",
            "if you rely on runSpatialNMF(), alignPrograms(), tuneNMF() or ",
            "runGuidedNMF().")
    }
    invisible()
}
