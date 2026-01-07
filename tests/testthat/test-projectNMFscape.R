test_that("projectNMFscape works with basic projection", {
    skip_if_not_installed("scuttle")
    
    library(scuttle)
    
    # Create source data with NMF results
    sce_source <- mockSCE(ngenes = 200, ncells = 100)
    sce_source <- logNormCounts(sce_source)
    sce_source <- runNMFscape(sce_source, k = 5, verbose = FALSE)
    
    # Create target data with overlapping features
    sce_target <- mockSCE(ngenes = 180, ncells = 80)
    rownames(sce_target) <- rownames(sce_source)[1:180]  # Ensure overlap
    sce_target <- logNormCounts(sce_target)
    
    # Test basic projection
    sce_target <- projectNMFscape(target = sce_target, source = sce_source, verbose = FALSE)
    
    # Check that results are stored correctly
    expect_true("projected_NMF" %in% reducedDimNames(sce_target))
    expect_true("projected_NMF_basis" %in% names(metadata(sce_target)))
    expect_true("projected_NMF_info" %in% names(metadata(sce_target)))
    
    # Check dimensions
    projected_coeffs <- reducedDim(sce_target, "projected_NMF")
    expect_equal(nrow(projected_coeffs), ncol(sce_target))
    expect_equal(ncol(projected_coeffs), 5)
    
    # Check projection info
    proj_info <- metadata(sce_target)$projected_NMF_info
    expect_equal(proj_info$n_factors, 5)
    expect_equal(proj_info$source_name, "NMF")
    expect_equal(proj_info$target_assay, "logcounts")
})

test_that("projectNMFscape handles feature matching correctly", {
    skip_if_not_installed("scuttle")
    
    library(scuttle)
    
    # Create source and target with different feature sets
    sce_source <- mockSCE(ngenes = 150, ncells = 50)
    rownames(sce_source) <- paste0("Gene_", 1:150)
    sce_source <- logNormCounts(sce_source)
    sce_source <- runNMFscape(sce_source, k = 3, verbose = FALSE)
    
    sce_target <- mockSCE(ngenes = 120, ncells = 60)
    rownames(sce_target) <- paste0("Gene_", 51:170)  # 50% overlap
    sce_target <- logNormCounts(sce_target)
    
    # Test projection with feature matching
    sce_target <- projectNMFscape(target = sce_target, source = sce_source, 
                                 match_features = TRUE, verbose = FALSE)
    
    # Should succeed with partial overlap
    expect_true("projected_NMF" %in% reducedDimNames(sce_target))
    
    # Check that common features were used
    proj_info <- metadata(sce_target)$projected_NMF_info
    expect_true(proj_info$n_common_features > 0)
    expect_true(proj_info$n_common_features < 150)
})

test_that("projectNMFscape fails appropriately with no common features", {
    skip_if_not_installed("scuttle")
    
    library(scuttle)
    
    # Create source and target with no overlapping features
    sce_source <- mockSCE(ngenes = 100, ncells = 50)
    rownames(sce_source) <- paste0("SourceGene_", 1:100)
    sce_source <- logNormCounts(sce_source)
    sce_source <- runNMFscape(sce_source, k = 3, verbose = FALSE)
    
    sce_target <- mockSCE(ngenes = 100, ncells = 60)
    rownames(sce_target) <- paste0("TargetGene_", 1:100)
    sce_target <- logNormCounts(sce_target)
    
    # Should fail with no common features
    expect_error(
        projectNMFscape(target = sce_target, source = sce_source, 
                       match_features = TRUE, verbose = FALSE),
        "No common features found"
    )
})

test_that("projectNMFscape validates inputs correctly", {
    skip_if_not_installed("scuttle")
    
    library(scuttle)
    
    sce <- mockSCE(ngenes = 100, ncells = 50)
    sce <- logNormCounts(sce)
    
    # Test with object without NMF results
    sce_no_nmf <- mockSCE(ngenes = 100, ncells = 50)
    sce_no_nmf <- logNormCounts(sce_no_nmf)
    
    expect_error(
        projectNMFscape(target = sce, source = sce_no_nmf, verbose = FALSE),
        "NMF results .* not found in source object"
    )
    
    # Test with non-SCE objects
    expect_error(
        projectNMFscape(target = "not_sce", source = sce, verbose = FALSE),
        "target must be a SingleCellExperiment"
    )
    
    expect_error(
        projectNMFscape(target = sce, source = "not_sce", verbose = FALSE),
        "source must be a SingleCellExperiment"
    )
})

test_that("projectNMFscape works with custom names", {
    skip_if_not_installed("scuttle")
    
    library(scuttle)
    
    # Create source with custom NMF name
    sce_source <- mockSCE(ngenes = 150, ncells = 50)
    sce_source <- logNormCounts(sce_source)
    sce_source <- runNMFscape(sce_source, k = 4, name = "custom_NMF", verbose = FALSE)
    
    # Create target
    sce_target <- mockSCE(ngenes = 150, ncells = 60)
    rownames(sce_target) <- rownames(sce_source)
    sce_target <- logNormCounts(sce_target)
    
    # Project with custom names
    sce_target <- projectNMFscape(target = sce_target, source = sce_source,
                                 source_name = "custom_NMF", 
                                 result_name = "my_projection",
                                 verbose = FALSE)
    
    # Check custom names are used
    expect_true("my_projection" %in% reducedDimNames(sce_target))
    expect_true("my_projection_basis" %in% names(metadata(sce_target)))
    expect_true("my_projection_info" %in% names(metadata(sce_target)))
})