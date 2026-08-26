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
    expect_s3_class(p, "ggplot")
    expect_true("meanWeight" %in% colnames(p$data))

    coeffs <- reducedDim(sce, "NMF")
    for (prog in colnames(coeffs)) {
        for (g in c("big", "small")) {
            expected <- mean(coeffs[sce$grp == g, prog])
            actual <- p$data$meanWeight[p$data$Program == prog &
                                            p$data$Group == g]
            expect_equal(actual, expected)
        }
    }
})

test_that("plotProgramDots honors color_palette and program selection", {
    skip_if_not_installed("scuttle")
    sce <- setup_viz_sce()

    p <- plotProgramDots(sce, group = "celltype",
                         color_palette = c("lightgrey", "blue"))
    expect_s3_class(p, "ggplot")
    expect_silent(invisible(ggplot2::ggplot_build(p)))

    p2 <- plotProgramDots(sce, group = "celltype", programs = c(1, 3))
    expect_equal(length(unique(as.character(p2$data$Program))), 2)

    expect_error(plotProgramDots(sce, group = "nope"), "not found")
})

test_that("plotProgramHeatmap builds for both gene-selection modes", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("pheatmap")
    sce <- setup_viz_sce()

    expect_s3_class(plotProgramHeatmap(sce, n_genes = 5), "pheatmap")
    expect_s3_class(plotProgramHeatmap(sce, n_genes = 5, make_unique = FALSE),
                    "pheatmap")
    expect_s3_class(plotProgramHeatmap(sce, programs = c(1, 2), n_genes = 4),
                    "pheatmap")
    expect_error(plotProgramHeatmap(sce, nmf_name = "nope"), "not found")
})

test_that("FindAllDEPs returns per-group statistics", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scran")
    sce <- setup_viz_sce()

    deps <- FindAllDEPs(sce, cell_type_col = "celltype")
    expect_named(deps, sort(unique(sce$celltype)), ignore.order = TRUE)
    for (nm in names(deps)) {
        df <- deps[[nm]]
        expect_true(all(c("fold_enrichment", "mean_usage_in",
                          "mean_usage_out") %in% colnames(df)))
        expect_equal(nrow(df), ncol(reducedDim(sce, "NMF")))
    }
    expect_error(FindAllDEPs(sce, cell_type_col = "nope"), "not found")
})

test_that("DEP plots build from FindAllDEPs output", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scran")
    sce <- setup_viz_sce()
    deps <- FindAllDEPs(sce, cell_type_col = "celltype")

    expect_s3_class(plotDEPsVolcano(deps), "ggplot")
    expect_s3_class(plotDEPsVolcano(deps, cell_types = "A"), "ggplot")
    expect_no_error(plotDEPsHeatmap(deps))
})

test_that("vizDimRed and vizUMAP validate their inputs", {
    skip_if_not_installed("scuttle")
    sce <- setup_viz_sce()

    expect_error(vizDimRed(sce, program = 1), "not found in reducedDims")
    expect_error(vizUMAP(sce, program = 1), "not found in reducedDims")

    # supply coordinates directly rather than depending on uwot
    set.seed(1)
    reducedDim(sce, "UMAP") <- matrix(rnorm(ncol(sce) * 2), ncol = 2,
                                      dimnames = list(colnames(sce),
                                                      c("UMAP1", "UMAP2")))
    expect_s3_class(vizDimRed(sce, program = 1), "ggplot")
    expect_s3_class(vizDimRed(sce, program = "NMF_2"), "ggplot")
    expect_s3_class(vizUMAP(sce, program = 1), "ggplot")
    expect_error(vizDimRed(sce, program = 1, dims = c(1, 5)), "not available")
})
