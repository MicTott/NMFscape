# Changes in version 0.99.3 (2026-08-25)

## Bug Fixes

* `runNMFscape()`, `consensusNMF()`, `runDeepNMF()` and `runConditionedNMF()`
  now honor a character `subset_row`. Previously `rownames(x)[subset_row]`
  assumed integer indexing, so passing gene names (the usual HVG workflow)
  silently produced `NA` rownames in the stored basis, breaking `getBasis()`,
  `getTopFeatures()` and `plotProgramHeatmap()`. Unknown feature names now
  error instead of being subset out.
* `runNMFscape(distribution = "auto")` now actually selects a distribution.
  It read `$best` from `auto_nmf_distribution()`, which returns `$loss`, so
  the chosen loss was `NULL` and silently fell back to `mse`.
* `runDeepNMF()`, `runMultiModalNMF()` and `runConditionedNMF()` now apply
  `L1`/`L2`. They accepted and documented the arguments but never passed them
  to `nmf_layer()`, so regularization was ignored.
* `runDeepNMF()` supports three or more layers. The `length(k) > 2` branch was
  unfinished and errored with non-conformable arguments. Results for two
  layers are unchanged.
* `getDiagonal()`, `predictNMF()` and `reconstructNMF()` now work with results
  from the FactorNet recipes, which store a `factor_net_result` rather than an
  S4 `nmf`. `evaluateNMF()` reports why it cannot and where the fitted loss
  lives.
* `plotProgramDots()` colors by mean rather than summed program weight, so the
  scale no longer tracks group size, and honors `color_palette`.

## Breaking Changes

* `L1`/`L2` in `runDeepNMF()`, `runMultiModalNMF()` and `runConditionedNMF()`
  are now scalars rather than length-2 vectors, matching
  `RcppML::nmf_layer()`, which applies one penalty per layer. Passing a
  length-2 value errors with an explanation.
* The data frame behind `plotProgramDots()` renames `sumWeight` to
  `meanWeight`.

## Vignettes

* All seven vignettes now execute at build time. Six previously set
  `eval = FALSE` globally, so none of the documented workflows were verified
  by the build -- which is how the `L1`/`L2`, `distribution = "auto"` and
  deep-NMF bugs above went unnoticed.
* The mechanics vignettes (*Choosing k*, *Deep NMF*, *Multi-Modal NMF*,
  *Transfer Learning*) now use simulated data with known structure instead of
  ExperimentHub downloads. They build in seconds, need no network, and can
  show recovery against ground truth: `selectRank()` finds the true rank, deep
  NMF recovers the simulated hierarchy one-to-one, and projection maps each
  query cell type onto the matching reference program.
* The interpretation vignettes (*Getting Started*, *Discovering Cell Type
  Programs*) keep the real Zeisel and Tasic datasets.
* Fixed duplicate chunk labels in the cell-type vignette, which `eval = FALSE`
  had been masking.
* Figures render at `fig.retina = 1`, `dpi = 96`. BiocStyle's retina default
  quadrupled every raster for no benefit in an HTML vignette. Building the
  vignettes still raises the installed size to ~11Mb, most of which is the
  ~0.9Mb of self-contained BiocStyle assets each of the seven vignettes
  carries.

## Tests

* Test suite expanded from 115 to 198. New regression coverage for character
  `subset_row`, `distribution = "auto"`, FactorNet regularization, deep NMF at
  three and four layers, and the FactorNet accessor paths.
* First tests for the visualization and differential-program layer, which had
  none: `plotProgramDots()`, `plotProgramHeatmap()`, `FindAllDEPs()`,
  `plotDEPsVolcano()`, `plotDEPsHeatmap()`, `vizDimRed()` and `vizUMAP()`.

## Internal

* Dropped unused `Matrix` and `parallel` imports, declared `grDevices`, and
  required `ggplot2 (>= 3.4.0)` for `linewidth`.
* `NEWS` converted to `NEWS.md` so R can parse it.

# Changes in version 0.99.1 (2026-04-08)

## Major Changes

* Upgraded to RcppML (>= 1.0.0), which introduces an S4 nmf class with
  slots @w, @d, @h, @misc and a three-factor decomposition A = W * diag(d) * H.
  The diagonal scaling factor is absorbed symmetrically into W and H by default
  (via the new absorb_d parameter, default TRUE).
* Raw RcppML nmf models are now stored in metadata(x)[[paste0(name, "_model")]]
  so downstream functions like predictNMF() and refineNMF() can use them.

## New Features

Rank selection and diagnostics:
* selectRank() - cross-validation for optimal k via held-out reconstruction
  error. Supports replicates for error-bar estimates.
* plotRankSelection() - elbow plot visualization of CV results.
* selectDistribution() - compare loss functions (MSE, Poisson, NB, etc.) and
  recommend the best via BIC.
* evaluateNMF() - compute reconstruction loss for a fitted model.
* assessNMF() - evaluate embedding quality against known labels (ARI, NMI,
  silhouette, classifier accuracy).

Clustering and stability:
* consensusNMF() - consensus NMF across replicates with cophenetic stability
  score, cluster assignments, and consensus matrix.

Transfer learning:
* predictNMF() - project new SCE data onto an existing NMF model, with
  automatic feature alignment.
* refineNMF() - label-guided refinement of an NMF model using cell type
  annotations, with optional batch correction.

FactorNet recipe wrappers (multi-layer and multi-modal NMF):
* runDeepNMF() - hierarchical NMF with multiple layers. Stores per-layer W/H/d
  plus an effective basis mapping genes directly to meta-programs.
* runMultiModalNMF() - joint NMF across modalities via shared cell embedding.
  Supports SCE assays (same feature space) and altExps (different feature
  spaces, e.g., RNA + ATAC).
* runConditionedNMF() - batch-conditioned NMF that factors out a categorical
  or numeric covariate during factorization.
* runFactorNet() - power-user escape hatch for storing custom FactorNet graph
  results in standard NMFscape SCE slots.

Accessor enhancements:
* getModel() - extract the raw RcppML S4 nmf model from metadata.
* getDiagonal() - extract the diagonal scaling vector d.
* getBasis() and getTopFeatures() gained a modality parameter for multi-modal
  results, with helpful error messages when a modality is required but omitted.

New parameters on existing functions:
* runNMFscape() gained L2, distribution (with "auto" option that calls
  auto_nmf_distribution()), and absorb_d parameters. Supports all 6 RcppML
  loss distributions (mse, gp, nb, gamma, inverse_gaussian, tweedie).

## Vignettes

New vignette suite focused on snRNA-seq brain data:
* "Getting Started with NMFscape" - basic workflow (updated).
* "NMFscape: Package Overview" - full function reference and vignette map.
* "Choosing k: Rank Selection and Diagnostics" - CV, distribution selection,
  consensus stability, label-based assessment.
* "Discovering Cell Type Programs" - biological interpretation using Tasic
  visual cortex data, including label-guided refinement.
* "Deep NMF: Hierarchical Gene Programs" - multi-layer factorization.
* "Multi-Modal NMF Analysis" - joint RNA + chromatin factorization.
* "Transfer Learning and Batch Correction" - reference mapping with predictNMF()
  and batch-corrected joint embedding with runConditionedNMF().

## Bug Fixes

* Fixed R CMD check errors in vizDimRed() and vizUMAP() examples by wrapping
  runUMAP-dependent code in \dontrun{}.
* Fixed duplicate vignette title NOTE by repurposing NMFscape.Rmd as a package
  overview (distinct from getting-started.Rmd).

## Internal

* Added R/factornet-utils.R with shared internal helpers (.validateSCE,
  .extractAssayMatrix, .absorbDiagonal, .setFactorNames) to reduce
  duplication across FactorNet wrappers.
* Test suite expanded from 19 to 115 tests covering all new functionality.

## Deprecations / Breaking Changes

* Direct access to nmf_result$w and nmf_result$h from RcppML::nmf() output no
  longer works. Use getBasis() and getCoefficients() accessors, or access S4
  slots via getModel(sce)@w and getModel(sce)@h.

# Changes in version 0.99.0 (2025-09-09)

## Initial Release

* First Bioconductor submission of NMFscape
* Fast NMF implementation using RcppML backend for SingleCellExperiment objects
* Core functions for running NMF analysis on single-cell data
* Accessor functions: getBasis(), getCoefficients(), getTopFeatures(), reconstructNMF()
* Visualization functions: vizUMAP(), vizDimRed(), plotProgramHeatmap(), plotProgramDots()
* Differential program analysis: FindAllDEPs(), plotDEPsVolcano(), plotDEPsHeatmap()
* Support for identifying cell type-specific gene programs
* Integration with Bioconductor SingleCellExperiment framework
