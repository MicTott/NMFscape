test_that("plotProgramDots summarizes usage by mean, not sum", {
    sce <- runNMFscape(makeSCE(ncells = 120), k = 4, seed = 1, verbose = FALSE)
    # deliberately unbalanced: summing weights would track group size instead
    sce$grp <- c(rep("big", 100), rep("small", 20))
    p <- plotProgramDots(sce, group = "grp")

    expect_true("meanWeight" %in% colnames(p$data))
    coeffs <- reducedDim(sce, "NMF")
    expect_equal(p$data$meanWeight[p$data$Program == "NMF_1" &
                                       p$data$Group == "small"],
                 mean(coeffs[sce$grp == "small", "NMF_1"]))
})

test_that("the plotting and differential-program layer builds", {
    skip_if_not_installed("pheatmap")
    sce <- runNMFscape(makeSCE(ngenes = 300), k = 4, seed = 1, verbose = FALSE)

    expect_s3_class(plotProgramHeatmap(sce, n_genes = 5), "pheatmap")

    deps <- FindAllDEPs(sce, cell_type_col = "celltype")
    expect_named(deps, c("A", "B", "C"), ignore.order = TRUE)
    expect_s3_class(plotDEPsVolcano(deps), "ggplot")
    expect_no_error(plotDEPsHeatmap(deps))

    set.seed(1)
    reducedDim(sce, "UMAP") <- matrix(rnorm(ncol(sce) * 2), ncol = 2,
                                      dimnames = list(colnames(sce),
                                                      c("UMAP1", "UMAP2")))
    expect_s3_class(vizDimRed(sce, program = 1), "ggplot")
    expect_s3_class(vizUMAP(sce, program = 1), "ggplot")
})

test_that("plotProgramEnrichment shows the top gene sets per program", {
    res <- data.frame(program = rep(c("NMF_1", "NMF_2"), each = 3),
                      pathway = paste0("SET_", seq_len(6)),
                      pval = 1e-3,
                      padj = c(0.01, 0.2, 0.3, 0.02, 0.4, 0.5),
                      size = 20L, ES = 0.5, NES = c(2, 1, -1, 2, 1, -1),
                      leading_edge = "a,b")
    p <- plotProgramEnrichment(res, top_n = 2)

    expect_s3_class(p, "ggplot")
    # padj_cutoff drops the non-significant sets rather than filling to top_n
    expect_equal(nrow(p$data), 2L)
})
