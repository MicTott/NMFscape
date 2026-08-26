setup_viz_sce <- function(ngenes = 200, ncells = 120) {
    sce <- scuttle::mockSCE(ngenes = ngenes, ncells = ncells)
    sce <- scuttle::logNormCounts(sce)
    sce <- runNMFscape(sce, k = 4, seed = 1, verbose = FALSE)
    sce$celltype <- rep(c("A", "B", "C"), length.out = ncol(sce))
    sce
}

test_that("plotProgramDots summarizes usage by mean, not sum", {
    skip_if_not_installed("scuttle")
    sce <- setup_viz_sce()
    # deliberately unbalanced groups: a sum would track group size
    sce$grp <- c(rep("big", 100), rep("small", 20))
    p <- plotProgramDots(sce, group = "grp")

    expect_true("meanWeight" %in% colnames(p$data))
    coeffs <- reducedDim(sce, "NMF")
    expect_equal(p$data$meanWeight[p$data$Program == "NMF_1" &
                                       p$data$Group == "small"],
                 mean(coeffs[sce$grp == "small", "NMF_1"]))
})

test_that("the plotting and DEP layer builds", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scran")
    skip_if_not_installed("pheatmap")
    sce <- setup_viz_sce()

    expect_s3_class(plotProgramHeatmap(sce, n_genes = 5), "pheatmap")
    deps <- FindAllDEPs(sce, cell_type_col = "celltype")
    expect_named(deps, sort(unique(sce$celltype)), ignore.order = TRUE)
    expect_s3_class(plotDEPsVolcano(deps), "ggplot")
    expect_no_error(plotDEPsHeatmap(deps))

    set.seed(1)
    reducedDim(sce, "UMAP") <- matrix(rnorm(ncol(sce) * 2), ncol = 2,
                                      dimnames = list(colnames(sce),
                                                      c("UMAP1", "UMAP2")))
    expect_s3_class(vizDimRed(sce, program = 1), "ggplot")
    expect_s3_class(vizUMAP(sce, program = 1), "ggplot")
})
