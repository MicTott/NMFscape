test_that("gpuInfo reports availability without erroring", {
    info <- gpuInfo()
    expect_type(info, "list")
    expect_equal(names(info), c("available", "device"))
    expect_type(info$available, "logical")
    expect_length(info$available, 1)
    expect_false(is.na(info$available))
    if (!info$available) {
        expect_null(info$device)
    }
    expect_type(gpuInfo(force_recheck = TRUE)$available, "logical")
})
