# MitoThreshold reproducible analysis pipeline
#
# PURPOSE
#   1. Load one unintegrated Seurat object containing raw RNA counts.
#   2. Apply an initial nFeature_RNA filter.
#   3. Split the object by a user-selected metadata column (typically by run, dataset or batch).
#   4. Run SCTransform, clustering and DoubletFinder independently per split.
#   5. Retain singlets, rerun SCTransform and merge the processed splits.
#   6. Evaluate mitochondrial-threshold response trajectories for the complete
#      dataset and, optionally, for levels of a biological grouping variable.
#
# HOW TO USE
#   A. Install packages once, if needed:
#      install.packages(c(
#        "Seurat", "DoubletFinder", "dplyr", "purrr", "tidyr", "stringr",
#        "tibble", "ggplot2", "Matrix", "slider", "readr", "msigdbr", "future"
#      ))
#      install.packages("BiocManager")
#      BiocManager::install("AUCell")
#   B. Edit only the USER CONFIGURATION section below.
#   C. Start a fresh R session.
#   D. Run: Rscript MitoThreshold_Reproducible_Pipeline.R
#      Or, in RStudio: source("MitoThreshold_Reproducible_Pipeline.R")
#
# INPUT
#   Supported formats are: (1) a 10x feature-barcode .h5 file, (2) an .rds
#   file containing one Seurat object, or (3) an .RData/.rda file containing
#   a Seurat object. For .RData/.rda, set input_object_name if the file
#   contains more than one Seurat object.
#
#   A raw 10x .h5 file does not usually contain sample-level metadata. Supply
#   metadata either through h5_cell_metadata (one value applied to every cell)
#   or h5_metadata_file (a CSV/TSV with one row per barcode). The split_by and
#   sample_id columns must exist after import; group_by may be NULL.
#
# OUTPUT
#   01_preprocessing/   DoubletFinder diagnostics, singlet objects and summary
#   02_mito_analysis/   threshold trajectories, transition calls and plots
#   processed_singlets.rds
#   run_configuration.rds and sessionInfo.txt
#
# IMPORTANT
#   - Mitochondrial values are converted to fractions internally (0-1).
#   - No integrated assay is required or used.
#   - PCA threshold trajectories use cells within each mitochondrial bin,
#     matching the original analysis. Other trajectories are cumulative and
#     use cells below each cutoff; bins are left-closed and right-open.
#   - Only gene-variance and PCA-variance analyses require at least
#     min_cells_per_bin cells in an individual bin. Cell recovery, library
#     complexity, gene expression, lineage markers, cell-level AUCell and
#     pseudobulk AUCell retain every observed bin through the maximum bin.
#   - Existing non-empty output directories are never overwritten.

set.seed(2)

# USER CONFIGURATION ---------------------------------------------------------

cfg <- list(
  input_file ="/path/to/input/file",
  input_object_name = NULL,
  output_dir = "/path/to/save/output",
  
  # 10x .h5 import options (ignored for .rds/.RData/.rda inputs)
  h5_feature_type = NULL, # used when the H5 contains multiple modalities
  h5_project = NULL,
  h5_cell_metadata = list(
    PID = NULL, # optional metadata to add, if desired
    Tissue = NULL # optional metadata to add, if desired
  ),
  # Optional CSV/TSV with a barcode column and metadata columns. This is useful
  # when cells in one H5 require different sample/batch annotations.
  h5_metadata_file = NULL,
  h5_metadata_barcode_col = "barcode",
  
  # To resume only the mitochondrial analysis after preprocessing has already
  # completed, supply the earlier processed_singlets.rds path here and choose
  # a new, empty output_dir. Leave NULL for a complete run.
  resume_processed_file = NULL,
  
  # Metadata columns
  split_by = "Dataset",      # DoubletFinder is run independently per level (typically Run or Dataset)
  sample_id = "PID",         # required for pseudobulk pathway summaries to pseudobulk by each participant/donor
  group_by = "Tissue",       # e.g. "Tissue"; NULL runs total-dataset analysis
  group_levels = c("TLN", "NTLN", "PBMC"),       # optional plotting order, e.g. c("TLN", "NTLN")
  
  # Initial cell filter
  nfeature_min = 100L,
  nfeature_max = 7000L,
  
  # Mitochondrial metadata
  mito_col = "percent.mito", # calculated from mito_pattern if absent
  mito_pattern = "^MT-",     # use "^mt-" for many mouse annotations
  mito_scale = "auto",       # "auto", "percent" (0-100), or "fraction" (0-1)
  mito_bin_width = 0.05,
  mito_cutoff_min = 0.05,
  mito_cutoff_max = NULL,    # NULL uses the next observed bin, capped at 1
  
  # SCTransform and DoubletFinder
  seed = 2L,
  sct_vst_flavor = "v2",
  dims_use = 1:30,
  clustering_resolution = 0.8,
  doublet_rate = 0.06,        # scalar or named vector keyed by split level
  doubletfinder_pN = 0.25,
  save_umap = TRUE,
  
  # Threshold-response metrics
  # Only gene-variance and PCA-variance analyses require at least this many
  # cells in the individual mitochondrial bin. Every other trajectory uses all
  # observed bins through the maximum bin.
  min_cells_per_bin = 20L,
  pca_npcs = 15L,
  pca_summary_npcs = 10L,
  pca_variable_features = 2000L,
  
  # Optional biological metrics
  run_gene_expression = TRUE,
  run_aucell = TRUE,
  run_pseudobulk_aucell = TRUE,
  species = "Homo sapiens",
  msigdb_collection = "H",
  auc_fraction = 0.05,
  auc_batch_size = 1000L,
  
  housekeeping_genes = c(
    "ACTB", "GAPDH", "RPLP0", "RPS18", "B2M", "EEF1A1", "TUBA1B"
  ),
  apoptosis_genes = c(
    "BAX", "BAK1", "CASP3", "CASP8", "TP53", "FAS", "DDIT3",
    "BBC3", "PMAIP1", "GADD45A", "CYCS"
  ),
  mitochondrial_genes = c(
    "MT-CO1", "MT-CO2", "MT-CO3", "MT-ND1", "MT-ND2", "MT-CYB",
    "MT-ATP6"
  ),
  lineage_markers = c(
    "CD3D", "CD3E", "CD4", "CD8A",
    "MS4A1", "CD79A", "CD79B", "IGHM",
    "NKG7", "GZMA", "GZMB", "KLRD1",
    "LYZ", "CD14", "S100A8", "S100A9", "FCER1A", "CLEC9A", "CD1C",
    "IGKC", "MZB1", "XBP1", "JCHAIN",
    "MKI67", "TOP2A", "STMN1"
  ),
  lineage_groups = c(
    CD3D = "T cell", CD3E = "T cell", CD4 = "T cell", CD8A = "T cell",
    MS4A1 = "B cell", CD79A = "B cell", CD79B = "B cell", IGHM = "B cell",
    NKG7 = "NK cell", GZMA = "NK cell", GZMB = "NK cell", KLRD1 = "NK cell",
    LYZ = "Myeloid", CD14 = "Myeloid", S100A8 = "Myeloid",
    S100A9 = "Myeloid", FCER1A = "Myeloid", CLEC9A = "Myeloid",
    CD1C = "Myeloid", IGKC = "Plasma cell", MZB1 = "Plasma cell",
    XBP1 = "Plasma cell", JCHAIN = "Plasma cell",
    MKI67 = "Proliferating", TOP2A = "Proliferating", STMN1 = "Proliferating"
  ),
  hallmark_sets = c(
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_HYPOXIA",
    "HALLMARK_APOPTOSIS",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_PEROXISOME",
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
  )
)

# DEPENDENCIES ---------------------------------------------------------------

required_packages <- c(
  "Seurat", "DoubletFinder", "dplyr", "purrr", "tidyr", "stringr",
  "tibble", "ggplot2", "Matrix", "slider", "readr", "future"
)

if (tolower(tools::file_ext(cfg$input_file)) %in% c("h5", "hdf5")) {
  required_packages <- c(required_packages, "hdf5r")
}

if (isTRUE(cfg$run_aucell)) {
  required_packages <- c(required_packages, "AUCell", "msigdbr")
}

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running the pipeline: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(Matrix)
  library(slider)
})

# GENERAL HELPERS ------------------------------------------------------------

safe_name <- function(x) {
  x <- stringr::str_replace_all(as.character(x), "[/\\\\:*?\"<>|]", "_")
  stringr::str_replace_all(x, "\\s+", "_")
}

create_empty_output_dir <- function(path) {
  if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop(
      "Output directory is not empty: ", normalizePath(path),
      "\nChoose a new cfg$output_dir. Existing results will not be overwritten."
    )
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

make_dir <- function(...) {
  path <- file.path(...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

assert_metadata <- function(obj, columns) {
  columns <- columns[!is.na(columns) & nzchar(columns)]
  missing <- setdiff(columns, colnames(obj[[]]))
  if (length(missing) > 0L) {
    stop("Missing metadata column(s): ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

validate_config <- function(cfg) {
  # Backward-compatible alias used in the first polished draft.
  if (is.null(cfg$min_cells_per_bin) && !is.null(cfg$min_cells_per_cutoff)) {
    warning(
      "cfg$min_cells_per_cutoff has been renamed to cfg$min_cells_per_bin; ",
      "using the supplied value. Please update the configuration name."
    )
    cfg$min_cells_per_bin <- cfg$min_cells_per_cutoff
  }
  if (
    is.null(cfg$min_cells_per_bin) || length(cfg$min_cells_per_bin) != 1L ||
    !is.finite(cfg$min_cells_per_bin) || cfg$min_cells_per_bin < 1
  ) {
    stop("cfg$min_cells_per_bin must be one finite number greater than or equal to 1.")
  }
  cfg$min_cells_per_bin <- as.integer(cfg$min_cells_per_bin)
  cfg
}

read_barcode_metadata <- function(path, barcode_col) {
  if (!file.exists(path)) stop("H5 metadata file does not exist: ", path)
  extension <- tolower(tools::file_ext(path))
  if (extension == "csv") {
    metadata <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (extension %in% c("tsv", "txt")) {
    metadata <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    stop("cfg$h5_metadata_file must be a .csv, .tsv or .txt file.")
  }
  if (!barcode_col %in% colnames(metadata)) {
    stop("H5 metadata table has no barcode column named: ", barcode_col)
  }
  barcodes <- as.character(metadata[[barcode_col]])
  if (anyNA(barcodes) || any(!nzchar(barcodes)) || anyDuplicated(barcodes)) {
    stop("The H5 metadata barcode column must be complete and unique.")
  }
  metadata
}

add_h5_metadata <- function(obj, constant_metadata = NULL,
                            metadata_file = NULL, barcode_col = "barcode") {
  if (!is.null(constant_metadata)) {
    if (!is.list(constant_metadata) || is.null(names(constant_metadata)) ||
        any(!nzchar(names(constant_metadata)))) {
      stop("cfg$h5_cell_metadata must be NULL or a named list.")
    }
    for (column in names(constant_metadata)) {
      value <- constant_metadata[[column]]
      if (!length(value) %in% c(1L, ncol(obj))) {
        stop("H5 constant metadata '", column,
             "' must have length 1 or equal the number of cells.")
      }
      obj[[column]] <- value
    }
  }
  
  if (!is.null(metadata_file)) {
    metadata <- read_barcode_metadata(metadata_file, barcode_col)
    matched_rows <- match(colnames(obj), as.character(metadata[[barcode_col]]))
    if (anyNA(matched_rows)) {
      missing_barcodes <- colnames(obj)[is.na(matched_rows)]
      stop(
        "The H5 metadata table is missing ", length(missing_barcodes),
        " cell barcode(s). First missing barcode(s): ",
        paste(utils::head(missing_barcodes, 5L), collapse = ", ")
      )
    }
    metadata_columns <- setdiff(colnames(metadata), barcode_col)
    collisions <- intersect(metadata_columns, colnames(obj[[]]))
    if (length(collisions) > 0L) {
      stop(
        "H5 metadata columns would overwrite existing Seurat metadata: ",
        paste(collisions, collapse = ", ")
      )
    }
    for (column in metadata_columns) {
      obj[[column]] <- metadata[[column]][matched_rows]
    }
  }
  obj
}

load_seurat_object <- function(path, object_name = NULL,
                               h5_feature_type = "Gene Expression",
                               h5_project = "SeuratProject",
                               h5_cell_metadata = NULL,
                               h5_metadata_file = NULL,
                               h5_metadata_barcode_col = "barcode") {
  if (!file.exists(path)) stop("Input file does not exist: ", path)
  extension <- tolower(tools::file_ext(path))
  
  if (extension %in% c("h5", "hdf5")) {
    h5_data <- Seurat::Read10X_h5(
      path, use.names = TRUE, unique.features = TRUE
    )
    if (inherits(h5_data, "Matrix") || is.matrix(h5_data)) {
      counts <- h5_data
    } else if (is.list(h5_data)) {
      if (length(h5_data) == 1L &&
          (is.null(h5_feature_type) || !h5_feature_type %in% names(h5_data))) {
        counts <- h5_data[[1L]]
      } else {
        if (is.null(h5_feature_type) || !h5_feature_type %in% names(h5_data)) {
          stop(
            "The H5 file contains these feature types: ",
            paste(names(h5_data), collapse = ", "),
            ". Set cfg$h5_feature_type to the RNA feature type."
          )
        }
        counts <- h5_data[[h5_feature_type]]
      }
    } else {
      stop("Read10X_h5 returned an unsupported object type: ", class(h5_data)[[1L]])
    }
    if (nrow(counts) == 0L || ncol(counts) == 0L) {
      stop("The selected H5 count matrix is empty.")
    }
    obj <- Seurat::CreateSeuratObject(
      counts = counts, project = h5_project, assay = "RNA",
      min.cells = 0L, min.features = 0L
    )
    obj <- add_h5_metadata(
      obj,
      constant_metadata = h5_cell_metadata,
      metadata_file = h5_metadata_file,
      barcode_col = h5_metadata_barcode_col
    )
  } else if (extension == "rds") {
    obj <- readRDS(path)
  } else if (extension %in% c("rdata", "rda")) {
    input_env <- new.env(parent = emptyenv())
    loaded <- load(path, envir = input_env)
    if (!is.null(object_name)) {
      if (!object_name %in% loaded) {
        stop("cfg$input_object_name was not found in the input file.")
      }
      obj <- input_env[[object_name]]
    } else {
      candidates <- loaded[vapply(loaded, function(nm) {
        inherits(input_env[[nm]], "Seurat")
      }, logical(1))]
      if (length(candidates) != 1L) {
        stop(
          "The .RData file must contain exactly one Seurat object, or set ",
          "cfg$input_object_name. Candidates: ", paste(candidates, collapse = ", ")
        )
      }
      obj <- input_env[[candidates]]
    }
  } else {
    stop("Input must be a 10x .h5/.hdf5, .rds, .RData or .rda file.")
  }
  
  if (!inherits(obj, "Seurat")) stop("The loaded object is not a Seurat object.")
  if (!"RNA" %in% Seurat::Assays(obj)) stop("The object has no RNA assay.")
  obj
}

join_assay_layers <- function(obj, assay = "RNA") {
  if (inherits(obj[[assay]], "Assay5") && length(SeuratObject::Layers(obj[[assay]])) > 1L) {
    obj[[assay]] <- SeuratObject::JoinLayers(obj[[assay]])
  }
  obj
}

get_layer <- function(obj, assay, layer) {
  SeuratObject::LayerData(obj, assay = assay, layer = layer)
}

write_csv_safe <- function(x, path) {
  readr::write_csv(as.data.frame(x), path, na = "")
}

mito_mid_from_bin <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  values <- stringr::str_extract_all(as.character(x), "[0-9.]+")
  purrr::map_dbl(values, function(parts) {
    numbers <- suppressWarnings(as.numeric(parts))
    if (length(numbers) >= 2L) return(numbers[[2]])
    if (length(numbers) == 1L) return(numbers[[1]])
    NA_real_
  })
}

resolve_doublet_rate <- function(rate, split_name) {
  if (length(rate) == 1L && is.null(names(rate))) return(as.numeric(rate))
  if (!is.null(names(rate)) && split_name %in% names(rate)) {
    return(as.numeric(rate[[split_name]]))
  }
  stop("No cfg$doublet_rate was supplied for split level: ", split_name)
}

# PREPROCESSING --------------------------------------------------------------

run_sct <- function(obj, cfg) {
  DefaultAssay(obj) <- "RNA"
  if ("SCT" %in% Seurat::Assays(obj)) obj[["SCT"]] <- NULL
  Seurat::SCTransform(
    obj,
    assay = "RNA",
    new.assay.name = "SCT",
    vst.flavor = cfg$sct_vst_flavor,
    vars.to.regress = NULL,
    verbose = FALSE
  )
}

valid_dimensions <- function(obj, requested) {
  max_dim <- min(max(requested), ncol(obj) - 1L, nrow(obj[["SCT"]]) - 1L)
  dims <- requested[requested <= max_dim]
  if (length(dims) < 2L) stop("Too few cells/features for the requested PCA dimensions.")
  dims
}

run_doubletfinder_split <- function(obj, split_name, cfg, plot_dir, object_dir) {
  message("\n--- SCTransform and DoubletFinder: ", split_name, " ---")
  
  obj <- run_sct(obj, cfg)
  dims <- valid_dimensions(obj, cfg$dims_use)
  # The original workflow used RunPCA's 50-PC default and then supplied PCs
  # 1:30 to clustering/DoubletFinder. Retaining that rank is important because
  # randomized PCA can change slightly when a different number of PCs is fitted.
  pca_npcs <- min(max(50L, max(dims)), ncol(obj) - 1L, nrow(obj[["SCT"]]) - 1L)
  obj <- Seurat::RunPCA(obj, npcs = pca_npcs, verbose = FALSE)
  obj <- Seurat::FindNeighbors(obj, dims = dims, verbose = FALSE)
  obj <- Seurat::FindClusters(
    obj, resolution = cfg$clustering_resolution, verbose = FALSE
  )
  obj$mitothresh_df_cluster <- as.character(Seurat::Idents(obj))
  
  if (isTRUE(cfg$save_umap)) {
    obj <- Seurat::RunUMAP(obj, dims = dims, verbose = FALSE)
    qc_plot <- Seurat::DimPlot(obj, reduction = "umap", label = TRUE) +
      ggplot2::ggtitle(paste(split_name, "pre-DoubletFinder clusters"))
    ggplot2::ggsave(
      file.path(plot_dir, paste0(safe_name(split_name), "_clusters.pdf")),
      qc_plot, width = 8, height = 6
    )
  }
  
  sweep_results <- DoubletFinder::paramSweep(obj, PCs = dims, sct = TRUE)
  sweep_stats <- DoubletFinder::summarizeSweep(sweep_results, GT = FALSE)
  pk_table <- DoubletFinder::find.pK(sweep_stats)
  pk_table$BCmetric <- as.numeric(as.character(pk_table$BCmetric))
  pk_table$pK <- as.numeric(as.character(pk_table$pK))
  valid_pk <- is.finite(pk_table$BCmetric) & is.finite(pk_table$pK)
  if (!any(valid_pk)) stop("DoubletFinder could not estimate pK for ", split_name)
  best_pk <- pk_table$pK[valid_pk][which.max(pk_table$BCmetric[valid_pk])]
  
  rate <- resolve_doublet_rate(cfg$doublet_rate, split_name)
  expected_raw <- round(rate * ncol(obj))
  homotypic <- DoubletFinder::modelHomotypic(obj$mitothresh_df_cluster)
  expected_adjusted <- max(1L, round(expected_raw * (1 - homotypic)))
  
  before_cols <- colnames(obj[[]])
  obj <- DoubletFinder::doubletFinder(
    obj,
    PCs = dims,
    pN = cfg$doubletfinder_pN,
    pK = best_pk,
    nExp = expected_adjusted,
    reuse.pANN = NULL,
    sct = TRUE
  )
  new_cols <- setdiff(colnames(obj[[]]), before_cols)
  class_col <- grep("^DF.classifications_", new_cols, value = TRUE)
  pann_col <- grep("^pANN_", new_cols, value = TRUE)
  if (length(class_col) != 1L || length(pann_col) != 1L) {
    stop("Could not identify new DoubletFinder output columns for ", split_name)
  }
  obj$DF.class <- obj[[class_col]][, 1]
  obj$pANN <- obj[[pann_col]][, 1]
  
  if (isTRUE(cfg$save_umap)) {
    doublet_plot <- Seurat::DimPlot(
      obj,
      reduction = "umap",
      group.by = "DF.class",
      cols = c(Singlet = "grey75", Doublet = "firebrick")
    ) + ggplot2::ggtitle(paste(split_name, "DoubletFinder"))
    ggplot2::ggsave(
      file.path(plot_dir, paste0(safe_name(split_name), "_doublets.pdf")),
      doublet_plot, width = 8, height = 6
    )
  }
  
  n_total <- ncol(obj)
  n_doublets <- sum(obj$DF.class == "Doublet", na.rm = TRUE)
  singlets <- subset(obj, subset = DF.class == "Singlet")
  if (ncol(singlets) < 3L) stop("Fewer than three singlets remain in ", split_name)
  singlets <- run_sct(singlets, cfg)
  saveRDS(
    singlets,
    file.path(object_dir, paste0(safe_name(split_name), "_singlets.rds")),
    compress = FALSE
  )
  
  summary <- tibble::tibble(
    split = as.character(split_name),
    optimal_pK = best_pk,
    n_cells = n_total,
    expected_doublets_raw = expected_raw,
    homotypic_proportion = homotypic,
    expected_doublets_adjusted = expected_adjusted,
    called_doublets = n_doublets,
    called_doublet_proportion = n_doublets / n_total,
    retained_singlets = ncol(singlets)
  )
  list(object = singlets, summary = summary)
}

merge_singlets <- function(objects) {
  if (length(objects) == 1L) return(objects[[1]])
  merge(x = objects[[1]], y = objects[-1], merge.data = TRUE)
}

# MITOCHONDRIAL VALUES AND STRATA -------------------------------------------

prepare_mito_fraction <- function(obj, cfg) {
  if (!cfg$mito_col %in% colnames(obj[[]])) {
    message("Calculating ", cfg$mito_col, " from feature pattern ", cfg$mito_pattern)
    obj[[cfg$mito_col]] <- Seurat::PercentageFeatureSet(
      obj, assay = "RNA", pattern = cfg$mito_pattern
    )
  }
  values <- as.numeric(obj[[cfg$mito_col]][, 1])
  if (any(!is.finite(values))) stop("The mitochondrial metadata contains missing/non-finite values.")
  if (any(values < 0)) stop("Mitochondrial values cannot be negative.")
  
  scale <- cfg$mito_scale
  if (identical(scale, "auto")) scale <- if (max(values) > 1) "percent" else "fraction"
  if (!scale %in% c("percent", "fraction")) {
    stop("cfg$mito_scale must be 'auto', 'percent' or 'fraction'.")
  }
  fraction <- if (scale == "percent") values / 100 else values
  if (any(fraction > 1)) stop("Converted mitochondrial fractions exceed 1; check cfg$mito_scale.")
  obj$mito_fraction <- fraction
  # Match cut(..., right = FALSE): 0.05 belongs to [0.05, 0.10), not
  # [0, 0.05). Keeping the numeric upper bound avoids fragile factor labels.
  obj$mito_bin_upper <- round(
    (floor(fraction / cfg$mito_bin_width) + 1) * cfg$mito_bin_width,
    10
  )
  obj
}

make_cutoffs <- function(values, cfg) {
  observed_upper <- round(
    (floor(values / cfg$mito_bin_width) + 1) * cfg$mito_bin_width,
    10
  )
  maximum <- if (is.null(cfg$mito_cutoff_max)) {
    max(observed_upper)
  } else {
    cfg$mito_cutoff_max
  }
  cutoffs <- sort(unique(observed_upper[
    observed_upper >= cfg$mito_cutoff_min & observed_upper <= maximum
  ]))
  if (length(cutoffs) < 3L) {
    stop("At least three mitochondrial cutoffs are required; adjust the cutoff range.")
  }
  cutoffs
}

valid_cutoffs_for_cells <- function(meta, cells, cutoffs, cfg) {
  if (is.null(cfg$min_cells_per_bin)) {
    stop(
      "Missing cfg$min_cells_per_bin. If your configuration still uses ",
      "min_cells_per_cutoff, rename it to min_cells_per_bin."
    )
  }
  bin_values <- meta$mito_bin_upper[match(cells, rownames(meta))]
  counts <- table(factor(
    bin_values, levels = cutoffs
  ))
  cutoffs[as.integer(counts) >= cfg$min_cells_per_bin]
}

observed_cutoffs_for_cells <- function(meta, cells, cutoffs) {
  bin_values <- meta$mito_bin_upper[match(cells, rownames(meta))]
  cutoffs[cutoffs %in% unique(bin_values[is.finite(bin_values)])]
}

make_strata <- function(obj, cfg) {
  meta <- obj[[]]
  strata <- list(Total = rownames(meta))
  if (!is.null(cfg$group_by) && nzchar(cfg$group_by)) {
    assert_metadata(obj, cfg$group_by)
    observed <- unique(as.character(meta[[cfg$group_by]]))
    observed <- observed[!is.na(observed)]
    if (!is.null(cfg$group_levels)) {
      absent <- setdiff(cfg$group_levels, observed)
      if (length(absent) > 0L) {
        warning("Configured group levels absent from data: ", paste(absent, collapse = ", "))
      }
      observed <- c(intersect(cfg$group_levels, observed), setdiff(observed, cfg$group_levels))
    }
    grouped <- lapply(observed, function(level) rownames(meta)[as.character(meta[[cfg$group_by]]) == level])
    names(grouped) <- observed
    strata <- c(strata, grouped)
  }
  strata
}

mean_gene_variance <- function(matrix, cells) {
  n <- length(cells)
  if (n < 2L) return(NA_real_)
  x <- matrix[, cells, drop = FALSE]
  sums <- Matrix::rowSums(x)
  sums_sq <- Matrix::rowSums(x ^ 2)
  variances <- (sums_sq - sums ^ 2 / n) / (n - 1)
  mean(variances[is.finite(variances)], na.rm = TRUE)
}

metric_trajectories <- function(obj, strata, cutoffs, cfg) {
  meta <- obj[[]]
  rna <- get_layer(obj, "RNA", "data")
  sct <- get_layer(obj, "SCT", "data")
  
  purrr::imap_dfr(strata, function(stratum_cells, stratum_name) {
    observed_cutoffs <- observed_cutoffs_for_cells(
      meta, stratum_cells, cutoffs
    )
    valid_cutoffs <- valid_cutoffs_for_cells(meta, stratum_cells, cutoffs, cfg)
    # The original pipeline used every observed bin for cumulative cell
    # recovery and complexity, but removed sparse bins before accumulating
    # gene-variance sufficient statistics. Keeping those rules separate is
    # essential because the rising-elbow call depends on the full curve tail.
    purrr::map_dfr(observed_cutoffs, function(cutoff) {
      cells <- intersect(
        stratum_cells,
        rownames(meta)[meta$mito_fraction < cutoff]
      )
      n <- length(cells)
      bin_cells <- intersect(
        stratum_cells,
        rownames(meta)[meta$mito_bin_upper == cutoff]
      )
      # The legacy variance analysis removed sparse bins before accumulating
      # sufficient statistics. Preserve that behavior for reproducibility.
      variance_is_valid <- cutoff %in% valid_cutoffs
      variance_cells <- if (variance_is_valid) {
        intersect(
          stratum_cells,
          rownames(meta)[
            meta$mito_bin_upper %in% valid_cutoffs[valid_cutoffs <= cutoff]
          ]
        )
      } else {
        character()
      }
      cell_meta <- meta[cells, , drop = FALSE]
      tibble::tibble(
        stratum = stratum_name,
        cutoff = cutoff,
        n_cells = n,
        n_cells_bin = length(bin_cells),
        retained_fraction = n / length(stratum_cells),
        mean_nFeature = mean(cell_meta$nFeature_RNA, na.rm = TRUE),
        mean_nCount = mean(cell_meta$nCount_RNA, na.rm = TRUE),
        mean_complexity = mean(
          cell_meta$nFeature_RNA / pmax(cell_meta$nCount_RNA, 1), na.rm = TRUE
        ),
        mean_variance_RNA = if (variance_is_valid) {
          mean_gene_variance(rna, variance_cells)
        } else {
          NA_real_
        },
        mean_variance_SCT = if (variance_is_valid) {
          mean_gene_variance(sct, variance_cells)
        } else {
          NA_real_
        }
      )
    })
  })
}

# PCA is evaluated within non-overlapping mitochondrial bins, matching the
# original scripts. It is deliberately kept separate from cumulative metrics.
pca_bin_trajectories <- function(obj, strata, cutoffs, cfg) {
  meta <- obj[[]]
  obj <- Seurat::FindVariableFeatures(
    obj, assay = "SCT", selection.method = "vst",
    nfeatures = cfg$pca_variable_features, verbose = FALSE
  )
  features <- intersect(
    Seurat::VariableFeatures(obj[["SCT"]]),
    rownames(get_layer(obj, "SCT", "scale.data"))
  )
  if (length(features) < 2L) {
    stop("Fewer than two SCT variable features have computed residuals for PCA.")
  }
  
  purrr::imap_dfr(strata, function(stratum_cells, stratum_name) {
    valid_cutoffs <- valid_cutoffs_for_cells(meta, stratum_cells, cutoffs, cfg)
    purrr::map_dfr(valid_cutoffs, function(cutoff) {
      lower <- max(0, cutoff - cfg$mito_bin_width)
      eligible <- rownames(meta)[meta$mito_bin_upper == cutoff]
      cells <- intersect(stratum_cells, eligible)
      npcs <- min(cfg$pca_npcs, length(cells) - 1L, length(features) - 1L)
      if (npcs < 1L) return(NULL)
      obj_bin <- subset(obj, cells = cells)
      Seurat::VariableFeatures(obj_bin[["SCT"]]) <- features
      obj_bin <- Seurat::RunPCA(
        obj_bin, assay = "SCT", features = features,
        npcs = npcs, verbose = FALSE
      )
      eigenvalues <- obj_bin[["pca"]]@stdev ^ 2
      percent <- eigenvalues / sum(eigenvalues) * 100
      summary_n <- min(cfg$pca_summary_npcs, length(percent))
      tibble::tibble(
        stratum = stratum_name,
        cutoff = cutoff,
        bin_lower = lower,
        n_cells_bin = length(cells),
        variance_PC1 = percent[[1]],
        cumulative_variance_PC1_to_N = sum(percent[seq_len(summary_n)]),
        N = summary_n
      )
    })
  })
}

# GENE EXPRESSION ------------------------------------------------------------

gene_category_table <- function(cfg) {
  lineage_labels <- unname(cfg$lineage_groups[cfg$lineage_markers])
  lineage_labels[is.na(lineage_labels)] <- "Other"
  dplyr::bind_rows(
    tibble::tibble(
      gene = cfg$housekeeping_genes, category = "Housekeeping",
      lineage = NA_character_
    ),
    tibble::tibble(
      gene = cfg$apoptosis_genes, category = "Apoptosis",
      lineage = NA_character_
    ),
    tibble::tibble(
      gene = cfg$mitochondrial_genes, category = "Mitochondrial",
      lineage = NA_character_
    ),
    tibble::tibble(
      gene = cfg$lineage_markers, category = "Lineage marker",
      lineage = lineage_labels
    )
  ) %>% dplyr::distinct(gene, .keep_all = TRUE)
}

gene_expression_trajectories <- function(obj, strata, cutoffs, cfg) {
  gene_table <- gene_category_table(cfg)
  rna <- get_layer(obj, "RNA", "data")
  gene_table <- dplyr::filter(gene_table, gene %in% rownames(rna))
  if (nrow(gene_table) == 0L) {
    warning("None of the configured QC/marker genes are present; gene trajectories skipped.")
    return(tibble::tibble(
      stratum = character(), cutoff = double(), gene = character(),
      mean_expression = double(), n_cells = integer(),
      n_cells_bin = integer(), category = character(), lineage = character()
    ))
  }
  rna <- rna[gene_table$gene, , drop = FALSE]
  meta <- obj[[]]
  
  purrr::imap_dfr(strata, function(stratum_cells, stratum_name) {
    observed_cutoffs <- observed_cutoffs_for_cells(
      meta, stratum_cells, cutoffs
    )
    purrr::map_dfr(observed_cutoffs, function(cutoff) {
      cells <- intersect(stratum_cells, rownames(meta)[meta$mito_fraction < cutoff])
      tibble::tibble(
        stratum = stratum_name,
        cutoff = cutoff,
        gene = rownames(rna),
        mean_expression = Matrix::rowMeans(rna[, cells, drop = FALSE]),
        n_cells = length(cells),
        n_cells_bin = sum(
          meta$mito_bin_upper[match(stratum_cells, rownames(meta))] == cutoff
        )
      )
    })
  }) %>% dplyr::left_join(gene_table, by = "gene")
}

# AUCELL --------------------------------------------------------------------

run_aucell_batched <- function(expression, gene_sets, batch_size, auc_fraction, seed) {
  gene_sets <- lapply(gene_sets, intersect, y = rownames(expression))
  gene_sets <- gene_sets[lengths(gene_sets) > 0L]
  if (length(gene_sets) == 0L) stop("No AUCell gene-set genes are present.")
  cells <- colnames(expression)
  batches <- split(seq_along(cells), ceiling(seq_along(cells) / batch_size))
  output <- matrix(
    NA_real_, nrow = length(gene_sets), ncol = length(cells),
    dimnames = list(names(gene_sets), cells)
  )
  set.seed(seed)
  for (i in seq_along(batches)) {
    batch_cells <- cells[batches[[i]]]
    message("AUCell batch ", i, "/", length(batches))
    rankings <- AUCell::AUCell_buildRankings(
      expression[, batch_cells, drop = FALSE],
      plotStats = FALSE, splitByBlocks = TRUE, nCores = 1
    )
    auc <- AUCell::AUCell_calcAUC(
      gene_sets,
      rankings,
      aucMaxRank = ceiling(auc_fraction * nrow(expression)),
      nCores = 1
    )
    values <- as.matrix(AUCell::getAUC(auc))
    output[rownames(values), colnames(values)] <- values
  }
  output
}

make_aucell_gene_sets <- function(expression, cfg) {
  hallmark <- msigdbr::msigdbr(
    species = cfg$species,
    collection = cfg$msigdb_collection
  ) %>%
    dplyr::filter(.data$gs_name %in% cfg$hallmark_sets) %>%
    dplyr::select("gs_name", "gene_symbol") %>%
    dplyr::group_by(.data$gs_name) %>%
    dplyr::summarise(genes = list(unique(.data$gene_symbol)), .groups = "drop") %>%
    tibble::deframe()
  hallmark$MT_GENES <- cfg$mitochondrial_genes
  hallmark <- lapply(hallmark, intersect, y = rownames(expression))
  hallmark[lengths(hallmark) > 0L]
}

aucell_cell_trajectories <- function(obj, strata, cutoffs, cfg) {
  expression <- get_layer(obj, "RNA", "data")
  sets <- make_aucell_gene_sets(expression, cfg)
  auc <- run_aucell_batched(
    expression, sets, cfg$auc_batch_size, cfg$auc_fraction, cfg$seed
  )
  meta <- obj[[]]
  values <- purrr::imap_dfr(strata, function(stratum_cells, stratum_name) {
    observed_cutoffs <- observed_cutoffs_for_cells(
      meta, stratum_cells, cutoffs
    )
    purrr::map_dfr(observed_cutoffs, function(cutoff) {
      cells <- intersect(stratum_cells, rownames(meta)[meta$mito_fraction < cutoff])
      cells <- intersect(cells, colnames(auc))
      tibble::tibble(
        stratum = stratum_name,
        cutoff = cutoff,
        gene_set = rownames(auc),
        mean_auc = Matrix::rowMeans(auc[, cells, drop = FALSE]),
        n_cells = length(cells),
        n_cells_bin = sum(
          meta$mito_bin_upper[match(stratum_cells, rownames(meta))] == cutoff
        )
      )
    })
  })
  list(trajectories = values, gene_sets = sets)
}

aucell_pseudobulk_trajectories <- function(obj, cutoffs, gene_sets, cfg) {
  assert_metadata(obj, c(cfg$sample_id, cfg$group_by))
  meta <- obj[[]]
  
  score_one_mode <- function(stratum_values, id_prefix) {
    keep <- !is.na(meta[[cfg$sample_id]]) & !is.na(stratum_values)
    mode_meta <- meta[keep, , drop = FALSE]
    stratum_values <- as.character(stratum_values[keep])
    composite_key <- paste(
      as.character(mode_meta[[cfg$sample_id]]), stratum_values,
      sprintf("%.10f", mode_meta$mito_bin_upper), sep = "\r"
    )
    unique_keys <- unique(composite_key)
    key_to_id <- stats::setNames(
      sprintf("%s%08d", id_prefix, seq_along(unique_keys)), unique_keys
    )
    mode_meta$mitothresh_pb_id <- unname(key_to_id[composite_key])
    first_cell <- match(unique_keys, composite_key)
    pb_meta <- tibble::tibble(
      pseudobulk = unname(key_to_id[unique_keys]),
      sample = as.character(mode_meta[[cfg$sample_id]][first_cell]),
      stratum = stratum_values[first_cell],
      bin_upper = mode_meta$mito_bin_upper[first_cell]
    )
    
    mode_obj <- subset(obj, cells = rownames(mode_meta))
    mode_obj$mitothresh_pb_id <- mode_meta$mitothresh_pb_id
    aggregate <- Seurat::AggregateExpression(
      mode_obj,
      group.by = "mitothresh_pb_id",
      assays = "RNA",
      slot = "data",
      return.seurat = FALSE,
      verbose = FALSE
    )$RNA
    rankings <- AUCell::AUCell_buildRankings(
      aggregate, plotStats = FALSE, nCores = 1
    )
    auc <- as.matrix(AUCell::getAUC(AUCell::AUCell_calcAUC(
      gene_sets,
      rankings,
      aucMaxRank = ceiling(cfg$auc_fraction * nrow(aggregate)),
      nCores = 1
    )))
    
    auc_long <- as.data.frame(t(auc)) %>%
      tibble::rownames_to_column("pseudobulk") %>%
      tidyr::pivot_longer(
        cols = -"pseudobulk", names_to = "gene_set", values_to = "auc"
      ) %>%
      dplyr::left_join(pb_meta, by = "pseudobulk")
    
    bin_summary <- auc_long %>%
      dplyr::group_by(.data$stratum, .data$bin_upper, .data$gene_set) %>%
      dplyr::summarise(
        bin_mean_auc = mean(.data$auc, na.rm = TRUE),
        n_pseudobulk_profiles = dplyr::n(),
        .groups = "drop"
      )
    bin_cell_counts <- tibble::tibble(
      stratum = stratum_values,
      bin_upper = mode_meta$mito_bin_upper
    ) %>%
      dplyr::count(.data$stratum, .data$bin_upper, name = "n_cells_bin")
    
    # Match the original pseudobulk calculation: first average sample-level
    # pseudobulks within each bin, then form a cell-count-weighted cumulative
    # mean across bins. Equal weighting of every sample-bin profile gives a
    # materially different curve when sample contributions are unbalanced.
    bin_summary %>%
      dplyr::left_join(bin_cell_counts, by = c("stratum", "bin_upper")) %>%
      dplyr::arrange(.data$stratum, .data$gene_set, .data$bin_upper) %>%
      dplyr::group_by(.data$stratum, .data$gene_set) %>%
      dplyr::mutate(
        weighted_auc_sum = cumsum(.data$bin_mean_auc * .data$n_cells_bin),
        cumulative_cells = cumsum(.data$n_cells_bin),
        mean_auc = .data$weighted_auc_sum / .data$cumulative_cells,
        n_pseudobulk_profiles_cumulative = cumsum(.data$n_pseudobulk_profiles)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(
        .data$bin_upper %in% cutoffs
      ) %>%
      dplyr::transmute(
        stratum = .data$stratum,
        cutoff = .data$bin_upper,
        gene_set = .data$gene_set,
        mean_auc = .data$mean_auc,
        n_cells = .data$cumulative_cells,
        n_cells_bin = .data$n_cells_bin,
        n_pseudobulk_profiles = .data$n_pseudobulk_profiles,
        n_pseudobulk_profiles_cumulative = .data$n_pseudobulk_profiles_cumulative
      )
  }
  
  output <- list(score_one_mode(rep("Total", nrow(meta)), "PBT"))
  if (!is.null(cfg$group_by) && nzchar(cfg$group_by)) {
    output[[2]] <- score_one_mode(meta[[cfg$group_by]], "PBG")
  }
  dplyr::bind_rows(output)
}

# TRANSITION DETECTION -------------------------------------------------------

# The original transition-selection function is inserted below without
# dataset-specific code. It accepts a long table and returns one transition
# call per stratum plus the smoothed composite trajectory.

compute_elbow <- function(df_long,
                          group_col,
                          metric_col,
                          value_col,
                          mito_bin_col = "mito_bin") {
  
  # HARDEN TYPES
  df_long <- df_long %>%
    mutate(
      "{group_col}"  := as.character(.data[[group_col]]),
      "{metric_col}" := as.character(.data[[metric_col]])
    )
  
  val_raw <- df_long[[value_col]]
  
  # 1) If already numeric, keep it
  if (is.numeric(val_raw)) {
    df_long[[value_col]] <- val_raw
  } else {
    # 2) Try as.numeric first (handles scientific notation like 1e-5)
    val_num <- suppressWarnings(as.numeric(as.character(val_raw)))
    
    # 3) If that completely fails, fall back to parse_number (handles "12%", "1,234")
    if (all(is.na(val_num))) {
      val_num <- readr::parse_number(as.character(val_raw))
    }
    
    df_long[[value_col]] <- val_num
  }
  
  # fail fast with helpful context
  if (all(is.na(df_long[[value_col]]))) {
    bad <- unique(head(as.character(val_raw), 20))
    stop(
      "compute_elbow(): '", value_col, "' became all NA after coercion. ",
      "Example raw values: ", paste(bad, collapse = " | ")
    )
  }
  
  df_long <- df_long %>%
    mutate(mito_mid = mito_mid_from_bin(.data[[mito_bin_col]]))
  
  # normalise within group + metric (safe when range == 0) 
  df_norm <- df_long %>%
    group_by(.data[[group_col]], .data[[metric_col]]) %>%
    mutate(
      vmin = min(.data[[value_col]], na.rm = TRUE),
      vmax = max(.data[[value_col]], na.rm = TRUE),
      vrng = vmax - vmin,
      value_norm = ifelse(is.finite(vrng) & vrng > 0,
                          (.data[[value_col]] - vmin) / vrng,
                          0)  # flat line -> 0 everywhere
    ) %>%
    ungroup() %>%
    select(-vmin, -vmax, -vrng)
  
  # SNR weights
  snr_df <- df_norm %>%
    group_by(.data[[group_col]], .data[[metric_col]]) %>%
    arrange(mito_mid) %>%
    summarise(
      signal = IQR(.data[[value_col]], na.rm = TRUE),
      noise  = mad(diff(.data[[value_col]]), constant = 1, na.rm = TRUE),
      snr = signal / (noise + 1e-6),
      .groups = "drop"
    )
  
  df_weighted <- df_norm %>%
    left_join(snr_df, by = c(group_col, metric_col)) %>%
    group_by(.data[[group_col]]) %>%
    mutate(weight = snr / max(snr, na.rm = TRUE)) %>%
    ungroup()
  
  df_composite <- df_weighted %>%
    group_by(.data[[group_col]], mito_mid) %>%
    summarise(
      composite = weighted.mean(value_norm, w = weight, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_composite_smooth <- df_composite %>%
    group_by(.data[[group_col]]) %>%
    arrange(mito_mid) %>%
    mutate(
      composite_smooth = slide_dbl(
        composite, ~ mean(.x, na.rm = TRUE),
        .before = 1, .after = 1, .complete = TRUE
      ),
      composite_smooth = ifelse(is.na(composite_smooth), composite, composite_smooth)
    ) %>%
    ungroup()
  
  # Features to drive method choice
  curve_features <- df_composite_smooth %>%
    group_by(.data[[group_col]]) %>%
    arrange(mito_mid) %>%
    mutate(
      d1 = composite_smooth - lag(composite_smooth),
      idx = row_number(),
      n = n()
    ) %>%
    summarise(
      has_neg = any(d1 < 0, na.rm = TRUE),
      has_pos = any(d1 > 0, na.rm = TRUE),
      
      max_drop = abs(min(d1, na.rm = TRUE)),
      max_rise = max(d1, na.rm = TRUE),
      med_abs_step = mean(abs(d1), na.rm = TRUE),
      mad_abs_step = mad(abs(d1), constant = 1, na.rm = TRUE),
      
      start_value = first(composite_smooth),
      peak_value  = max(composite_smooth, na.rm = TRUE),
      end_value   = dplyr::last(composite_smooth),
      
      total_gain = peak_value - start_value,
      plateau_drop_frac = (peak_value - end_value) / (total_gain + 1e-6),
      
      early_mean_abs_slope = mean(abs(d1[idx <= ceiling(n / 3)]), na.rm = TRUE),
      late_mean_abs_slope  = mean(abs(d1[idx > floor(2 * n / 3)]), na.rm = TRUE),
      late_flat_ratio = late_mean_abs_slope / (early_mean_abs_slope + 1e-6),
      
      u_shaped = has_neg & has_pos &
        (max_drop > 0.15) &
        (max_rise > 0.15),
      
      drop_z = max_drop / (mad_abs_step + 1e-6),
      
      plateau_like = has_pos &
        (total_gain > 0.05) &
        (late_flat_ratio < 0.5) &
        (plateau_drop_frac < 0.15),
      
      .groups = "drop"
    )
  
  # add flatline-after-knee score by joining knee index (recover the knee position)
  knee_pos <- df_composite_smooth %>%
    group_by(.data[[group_col]]) %>%
    arrange(mito_mid) %>%
    mutate(d1 = composite_smooth - lag(composite_smooth), idx = row_number()) %>%
    filter(!is.na(d1)) %>%
    summarise(knee_idx = idx[which.min(d1)], .groups = "drop")
  
  # add flatline-after-knee score (knee_idx must be present BEFORE summarise)
  df_steps <- df_composite_smooth %>%
    group_by(.data[[group_col]]) %>%
    arrange(mito_mid) %>%
    mutate(
      d1  = composite_smooth - lag(composite_smooth),
      idx = row_number()
    ) %>%
    filter(!is.na(d1)) %>%
    left_join(knee_pos, by = group_col)  # <-- now knee_idx exists in each group
  
  flatline_df <- df_steps %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      knee_idx = first(knee_idx),
      pre_mad  = mad(d1[idx <= knee_idx], constant = 1, na.rm = TRUE),
      post_mad = mad(d1[idx >  knee_idx], constant = 1, na.rm = TRUE),
      flatline_ratio = post_mad / (pre_mad + 1e-6),
      .groups = "drop"
    )
  
  curve_features <- curve_features %>%
    left_join(knee_pos, by = group_col) %>%
    left_join(flatline_df %>% select(-knee_idx), by = group_col) %>%
    mutate(
      knee_conf = drop_z / (flatline_ratio + 1e-6)
    )
  
  # Knee (max drop) threshold ----
  boundary_k <- 2  # exclude first/last 1 bins if knee lands there
  
  knee_threshold <- df_composite_smooth %>%
    group_by(.data[[group_col]]) %>%
    arrange(mito_mid) %>%
    mutate(
      d1 = composite_smooth - lag(composite_smooth),
      idx = row_number(),
      n   = n(),
      is_boundary = idx <= boundary_k | idx > (n - boundary_k)
    ) %>%
    filter(!is.na(d1)) %>%
    summarise(
      raw_idx = idx[which.min(d1)],
      raw_mid = mito_mid[which.min(d1)],
      interior_mid = mito_mid[which.min(ifelse(is_boundary, Inf, d1))],
      n = first(n),
      .groups = "drop"
    ) %>%
    transmute(
      !!group_col := .data[[group_col]],
      at_boundary = raw_idx <= boundary_k | raw_idx > (n - boundary_k),
      mito_threshold = if_else(at_boundary, interior_mid, raw_mid),
      method = if_else(at_boundary, "max_drop_boundary_fallback", "max_drop")
    )
  
  # Valley threshold (interior minimum) ----
  valley_threshold <- df_composite_smooth %>%
    group_by(.data[[group_col]]) %>%
    arrange(mito_mid) %>%
    mutate(is_boundary = mito_mid == min(mito_mid) | mito_mid == max(mito_mid)) %>%
    filter(!is_boundary) %>%
    filter(composite_smooth == min(composite_smooth, na.rm = TRUE)) %>%
    summarise(
      mito_threshold = max(mito_mid),   # (existing choice)
      method = "valley",
      .groups = "drop"
    )
  
  # Transition threshold (curvature-based) ----
  transition_threshold <- {
    df_curvature <- df_composite_smooth %>%
      group_by(.data[[group_col]]) %>%
      arrange(mito_mid) %>%
      mutate(
        delta = composite_smooth - lag(composite_smooth),
        curvature = abs((lead(composite_smooth) - composite_smooth) -
                          (composite_smooth - lag(composite_smooth)))
      ) %>%
      ungroup()
    
    transition_bins <- df_curvature %>%
      group_by(.data[[group_col]]) %>%
      mutate(
        curvature_z = (curvature - mean(curvature, na.rm = TRUE)) /
          mad(curvature, na.rm = TRUE)
      ) %>%
      filter(curvature_z > 1, delta < 0) %>%
      ungroup()
    
    if (nrow(transition_bins) == 0) {
      df_curvature %>%
        distinct(.data[[group_col]]) %>%
        transmute(
          !!group_col := .data[[group_col]],
          mito_threshold = NA_real_,
          method = NA_character_
        )
    } else {
      transition_bins %>%
        group_by(.data[[group_col]]) %>%
        summarise(
          mito_threshold = max(mito_mid, na.rm = TRUE),
          method = "transition_end",
          .groups = "drop"
        )
    }
  }
  
  # Rising-elbow threshold (for monotone increasing curves that flatten) ----
  rising_elbow_threshold <- {
    df_rise <- df_composite_smooth %>%
      group_by(.data[[group_col]]) %>%
      arrange(mito_mid) %>%
      mutate(idx = row_number()) %>%
      group_modify(~{
        d <- .x
        
        if (nrow(d) < 3) {
          return(tibble(
            mito_threshold = NA_real_,
            method = NA_character_
          ))
        }
        
        x1 <- d$mito_mid[1]
        y1 <- d$composite_smooth[1]
        x2 <- d$mito_mid[nrow(d)]
        y2 <- d$composite_smooth[nrow(d)]
        
        # line between first and last points
        denom <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
        
        if (!is.finite(denom) || denom == 0) {
          return(tibble(
            mito_threshold = NA_real_,
            method = NA_character_
          ))
        }
        
        # perpendicular distance from each point to the chord
        d$dist_to_chord <- abs(
          (y2 - y1) * d$mito_mid -
            (x2 - x1) * d$composite_smooth +
            x2 * y1 - y2 * x1
        ) / denom
        
        # avoid choosing first/last points
        d_inner <- d[2:(nrow(d) - 1), , drop = FALSE]
        
        tibble(
          mito_threshold = d_inner$mito_mid[which.max(d_inner$dist_to_chord)],
          method = "rising_elbow"
        )
      }) %>%
      ungroup()
    
    df_rise
  }
  
  # Plateau-start threshold (better for gradual rise-to-plateau curves) ----
  plateau_threshold <- {
    PLATEAU_PROP <- 0.98   # first point reaching 98% of plateau height
    FUTURE_GAIN_TOL <- 0.005  # allow at most 0.5% more gain after threshold
    
    df_plateau <- df_composite_smooth %>%
      group_by(.data[[group_col]]) %>%
      arrange(mito_mid) %>%
      mutate(
        peak_value = max(composite_smooth, na.rm = TRUE),
        plateau_cutoff = peak_value * PLATEAU_PROP,
        reaches_plateau = composite_smooth >= plateau_cutoff
      ) %>%
      ungroup()
    
    plateau_hits <- df_plateau %>%
      group_by(.data[[group_col]]) %>%
      mutate(
        future_max = rev(cummax(rev(composite_smooth))),
        future_gain = future_max - composite_smooth
      ) %>%
      filter(reaches_plateau, future_gain <= FUTURE_GAIN_TOL) %>%
      summarise(
        mito_threshold = min(mito_mid, na.rm = TRUE),
        method = "plateau_start",
        .groups = "drop"
      )
    
    if (nrow(plateau_hits) == 0) {
      df_composite_smooth %>%
        distinct(.data[[group_col]]) %>%
        transmute(
          !!group_col := .data[[group_col]],
          mito_threshold = NA_real_,
          method = NA_character_
        )
    } else {
      plateau_hits
    }
  }
  
  # Plateau-elbow threshold ----
  plateau_elbow_threshold <- {
    GAIN_FRAC_REQ   <- 0.75   # need to have achieved at least 75% of total rise
    FUTURE_GAIN_FRAC <- 0.12  # allow at most 10% of total rise to remain
    SLOPE_FRAC      <- 0.50   # current slope <= 40% of early slope
    USE_SMOOTHED_D1 <- FALSE  # use raw local slope from composite, not over-smoothed slope
    
    df_plateau_elbow <- df_composite_smooth %>%
      group_by(.data[[group_col]]) %>%
      arrange(mito_mid) %>%
      group_modify(~{
        d <- .x
        
        if (nrow(d) < 5) {
          return(tibble(
            mito_threshold = NA_real_,
            method = NA_character_
          ))
        }
        
        # Use local slope from raw composite or smoothed composite
        if (USE_SMOOTHED_D1) {
          d <- d %>%
            mutate(
              d1 = composite_smooth - lag(composite_smooth)
            )
        } else {
          d <- d %>%
            mutate(
              d1 = composite - lag(composite)
            )
        }
        
        # Early reference slope from first few positive steps
        early_steps <- d$d1[2:min(4, nrow(d))]
        early_steps <- early_steps[is.finite(early_steps) & early_steps > 0]
        
        if (length(early_steps) == 0) {
          return(tibble(
            mito_threshold = NA_real_,
            method = NA_character_
          ))
        }
        
        early_ref <- max(early_steps, na.rm = TRUE)
        
        start_val <- d$composite_smooth[1]
        peak_val  <- max(d$composite_smooth, na.rm = TRUE)
        total_gain <- peak_val - start_val
        
        if (!is.finite(total_gain) || total_gain <= 0) {
          return(tibble(
            mito_threshold = NA_real_,
            method = NA_character_
          ))
        }
        
        d <- d %>%
          mutate(
            gain_frac = (composite_smooth - start_val) / total_gain,
            future_max = rev(cummax(rev(composite_smooth))),
            future_gain = future_max - composite_smooth,
            future_gain_frac = future_gain / total_gain,
            slope_frac = d1 / (early_ref + 1e-6),
            flat_enough = !is.na(slope_frac) & slope_frac <= SLOPE_FRAC
          )
        
        hits <- d %>%
          filter(
            gain_frac >= GAIN_FRAC_REQ,
            future_gain_frac <= FUTURE_GAIN_FRAC,
            flat_enough
          )
        
        if (nrow(hits) == 0) {
          tibble(
            mito_threshold = NA_real_,
            method = NA_character_
          )
        } else {
          tibble(
            mito_threshold = min(hits$mito_mid, na.rm = TRUE),
            method = "plateau_elbow"
          )
        }
      }) %>%
      ungroup()
    
    df_plateau_elbow
  }
  
  # Final choice logic
  # Tune these once on your plots
  KNEE_CONF_CUTOFF <- 6     # "in your face" threshold
  DROPZ_CUTOFF     <- 4     # backup if flatline_ratio is weird
  
  final_threshold <- curve_features %>%
    left_join(
      plateau_elbow_threshold %>%
        dplyr::rename(
          mito_threshold_plateau_elbow = mito_threshold,
          method_plateau_elbow = method
        ),
      by = group_col
    ) %>%
    left_join(
      plateau_threshold %>%
        dplyr::rename(
          mito_threshold_plateau = mito_threshold,
          method_plateau = method
        ),
      by = group_col
    ) %>%
    left_join(
      rising_elbow_threshold %>%
        dplyr::rename(
          mito_threshold_rising = mito_threshold,
          method_rising = method
        ),
      by = group_col
    ) %>%
    left_join(knee_threshold, by = group_col) %>%
    left_join(valley_threshold, by = group_col, suffix = c("_knee", "_valley")) %>%
    left_join(
      transition_threshold %>%
        dplyr::rename(
          mito_threshold_transition = mito_threshold,
          method_transition = method
        ),
      by = group_col
    ) %>%
    transmute(
      !!group_col := .data[[group_col]],
      mito_threshold = case_when(
        plateau_like & !is.na(mito_threshold_rising) ~ mito_threshold_rising,
        plateau_like & !is.na(mito_threshold_plateau_elbow) ~ mito_threshold_plateau_elbow,
        plateau_like & !is.na(mito_threshold_plateau) ~ mito_threshold_plateau,
        
        !plateau_like & !is.na(mito_threshold_knee) &
          (knee_conf >= KNEE_CONF_CUTOFF | drop_z >= DROPZ_CUTOFF) ~ mito_threshold_knee,
        
        u_shaped & !is.na(mito_threshold_valley) ~ mito_threshold_valley,
        !is.na(mito_threshold_transition) ~ mito_threshold_transition,
        !is.na(mito_threshold_valley) ~ mito_threshold_valley,
        !is.na(mito_threshold_knee) ~ mito_threshold_knee,
        TRUE ~ NA_real_
      ),
      method = case_when(
        plateau_like & !is.na(mito_threshold_rising) ~ "rising_elbow",
        plateau_like & !is.na(mito_threshold_plateau_elbow) ~ "plateau_elbow",
        plateau_like & !is.na(mito_threshold_plateau) ~ "plateau_start",
        
        !plateau_like & !is.na(mito_threshold_knee) &
          (knee_conf >= KNEE_CONF_CUTOFF | drop_z >= DROPZ_CUTOFF) ~ "max_drop_obvious",
        
        u_shaped & !is.na(mito_threshold_valley) ~ "valley_override",
        !is.na(mito_threshold_transition) ~ method_transition,
        !is.na(mito_threshold_valley) ~ "valley",
        !is.na(mito_threshold_knee) ~ "max_drop_fallback",
        TRUE ~ NA_character_
      )
    )
  
  # bin lookup for overlay (use upper bound)
  bin_lookup <- df_long %>%
    distinct(.data[[mito_bin_col]]) %>%
    mutate(mito_mid = mito_mid_from_bin(.data[[mito_bin_col]]))
  
  df_composite_binned <- df_composite_smooth %>%
    left_join(bin_lookup, by = "mito_mid") %>%
    mutate(
      !!mito_bin_col := factor(.data[[mito_bin_col]], levels = unique(bin_lookup[[mito_bin_col]]))
    )
  
  list(
    final_threshold = final_threshold,
    composite_smooth = df_composite_binned
  )
  
}

detect_all_transitions <- function(metrics, pca_metrics) {
  metric_long <- metrics %>%
    dplyr::select(
      "stratum", "cutoff", "n_cells", "mean_complexity",
      "mean_variance_RNA", "mean_variance_SCT"
    ) %>%
    tidyr::pivot_longer(
      cols = -c("stratum", "cutoff"),
      names_to = "metric", values_to = "value"
    )
  
  structural <- dplyr::filter(
    metric_long,
    grepl("mean_variance", .data$metric),
    is.finite(.data$value)
  )
  complexity <- dplyr::filter(metric_long, .data$metric == "mean_complexity")
  retained <- dplyr::filter(metric_long, .data$metric == "n_cells")
  pca_long <- pca_metrics %>%
    dplyr::select(
      "stratum", "cutoff", "variance_PC1",
      "cumulative_variance_PC1_to_N"
    ) %>%
    tidyr::pivot_longer(
      cols = -c("stratum", "cutoff"),
      names_to = "metric", values_to = "value"
    )
  
  inputs <- list(
    gene_variance = structural,
    pca_variance = pca_long,
    complexity = complexity,
    retained_cells = retained
  )
  
  results <- purrr::imap(inputs, function(data, analysis) {
    if (nrow(data) == 0L) return(NULL)
    result <- compute_elbow(
      df_long = data,
      group_col = "stratum",
      metric_col = "metric",
      value_col = "value",
      mito_bin_col = "cutoff"
    )
    result$final_threshold$analysis <- analysis
    result$composite_smooth$analysis <- analysis
    result
  })
  
  calls <- purrr::map_dfr(results, "final_threshold")
  composites <- purrr::map_dfr(results, "composite_smooth")
  consensus <- calls %>%
    dplyr::group_by(.data$stratum) %>%
    dplyr::summarise(
      consensus_median = stats::median(.data$mito_threshold, na.rm = TRUE),
      consensus_mean = mean(.data$mito_threshold, na.rm = TRUE),
      consensus_min = min(.data$mito_threshold, na.rm = TRUE),
      consensus_max = max(.data$mito_threshold, na.rm = TRUE),
      n_metrics = sum(is.finite(.data$mito_threshold)),
      .groups = "drop"
    )
  list(calls = calls, composites = composites, consensus = consensus)
}

# PLOTS ---------------------------------------------------------------------

save_standard_plots <- function(metrics, pca_metrics, gene_expression,
                                auc_cell, auc_pseudobulk, transitions, plot_dir) {
  metric_plot_data <- metrics %>%
    dplyr::select(
      "stratum", "cutoff", "retained_fraction",
      "mean_complexity", "mean_variance_RNA", "mean_variance_SCT"
    ) %>%
    tidyr::pivot_longer(
      cols = -c("stratum", "cutoff"),
      names_to = "metric", values_to = "value"
    )
  p <- ggplot2::ggplot(
    metric_plot_data,
    ggplot2::aes(x = .data$cutoff, y = .data$value, colour = .data$metric)
  ) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$metric),
                        cols = ggplot2::vars(.data$stratum), scales = "free") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::labs(x = "Cumulative mitochondrial cutoff", y = NULL, colour = NULL)
  ggplot2::ggsave(file.path(plot_dir, "cumulative_metric_trajectories.pdf"),
                  p, width = 12, height = 10)
  
  if (nrow(pca_metrics) > 0L) {
    pca_long <- pca_metrics %>%
      dplyr::select("stratum", "cutoff", "variance_PC1",
                    "cumulative_variance_PC1_to_N") %>%
      tidyr::pivot_longer(-c("stratum", "cutoff"),
                          names_to = "metric", values_to = "value")
    p <- ggplot2::ggplot(
      pca_long,
      ggplot2::aes(.data$cutoff, .data$value, colour = .data$metric)
    ) + ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::facet_wrap(ggplot2::vars(.data$stratum), scales = "free") +
      ggplot2::theme_classic() +
      ggplot2::labs(x = "Mitochondrial bin upper bound", y = "% variance", colour = NULL)
    ggplot2::ggsave(file.path(plot_dir, "pca_bin_trajectories.pdf"), p, width = 12, height = 6)
  }
  
  if (!is.null(gene_expression) && nrow(gene_expression) > 0L) {
    qc_expression <- dplyr::filter(
      gene_expression, .data$category != "Lineage marker"
    )
    p <- ggplot2::ggplot(
      qc_expression,
      ggplot2::aes(.data$cutoff, .data$mean_expression,
                   colour = .data$gene, group = .data$gene)
    ) + ggplot2::geom_line() +
      ggplot2::facet_grid(rows = ggplot2::vars(.data$category),
                          cols = ggplot2::vars(.data$stratum), scales = "free") +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::labs(
        title = "Housekeeping, apoptosis and mitochondrial genes",
        x = "Cumulative mitochondrial cutoff", y = "Mean expression"
      )
    ggplot2::ggsave(file.path(plot_dir, "qc_gene_expression_trajectories.pdf"),
                    p, width = 14, height = 11)
    
    lineage_expression <- dplyr::filter(
      gene_expression, .data$category == "Lineage marker"
    )
    if (nrow(lineage_expression) > 0L) {
      lineage_expression$lineage <- factor(
        lineage_expression$lineage,
        levels = c(
          "T cell", "B cell", "NK cell", "Myeloid", "Plasma cell",
          "Proliferating", "Other"
        )
      )
      p_lineage <- ggplot2::ggplot(
        lineage_expression,
        ggplot2::aes(
          .data$cutoff, .data$mean_expression,
          colour = .data$gene, group = .data$gene
        )
      ) + ggplot2::geom_line() + ggplot2::geom_point(size = 1.5) +
        ggplot2::facet_grid(
          rows = ggplot2::vars(.data$lineage),
          cols = ggplot2::vars(.data$stratum), scales = "free"
        ) +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::labs(
          title = "Lineage-marker expression",
          x = "Cumulative mitochondrial cutoff", y = "Mean expression"
        )
      ggplot2::ggsave(
        file.path(plot_dir, "lineage_marker_trajectories.pdf"),
        p_lineage, width = 16, height = 14
      )
    }
  }
  
  if (!is.null(auc_cell) && nrow(auc_cell) > 0L) {
    p <- ggplot2::ggplot(
      auc_cell,
      ggplot2::aes(.data$cutoff, .data$mean_auc,
                   colour = .data$gene_set, group = .data$gene_set)
    ) + ggplot2::geom_line() +
      ggplot2::facet_wrap(ggplot2::vars(.data$stratum), scales = "free") +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::labs(x = "Cumulative mitochondrial cutoff", y = "Mean AUCell AUC")
    ggplot2::ggsave(file.path(plot_dir, "aucell_cell_trajectories.pdf"),
                    p, width = 13, height = 7)
  }
  
  if (!is.null(auc_pseudobulk) && nrow(auc_pseudobulk) > 0L) {
    p <- ggplot2::ggplot(
      auc_pseudobulk,
      ggplot2::aes(.data$cutoff, .data$mean_auc,
                   colour = .data$gene_set, group = .data$gene_set)
    ) + ggplot2::geom_line() +
      ggplot2::facet_wrap(ggplot2::vars(.data$stratum), scales = "free") +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::labs(x = "Cumulative mitochondrial cutoff", y = "Mean pseudobulk AUC")
    ggplot2::ggsave(file.path(plot_dir, "aucell_pseudobulk_trajectories.pdf"),
                    p, width = 13, height = 7)
  }
  
  if (nrow(transitions$calls) > 0L) {
    p <- ggplot2::ggplot(
      transitions$calls,
      ggplot2::aes(.data$mito_threshold, .data$analysis)
    ) + ggplot2::geom_point(size = 3) +
      ggplot2::geom_vline(
        data = transitions$consensus,
        ggplot2::aes(xintercept = .data$consensus_mean),
        colour = "firebrick", linetype = "dashed"
      ) +
      ggplot2::facet_wrap(ggplot2::vars(.data$stratum)) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = "Inferred mitochondrial transition", y = NULL)
    ggplot2::ggsave(file.path(plot_dir, "transition_summary.pdf"), p, width = 10, height = 6)
  }
}

# MAIN WORKFLOW --------------------------------------------------------------
# Allow large objects during SCTransform/DoubletFinder.
# This raises the permitted object size; it does not allocate 200 GB immediately.
options(future.globals.maxSize = 200 * 1024^3)

# Avoid copying the large Seurat object to multiple parallel workers.
future::plan(future::sequential)

run_pipeline <- function(cfg) {
  cfg <- validate_config(cfg)
  set.seed(cfg$seed)
  output_dir <- create_empty_output_dir(cfg$output_dir)
  preprocessing_dir <- make_dir(output_dir, "01_preprocessing")
  split_object_dir <- make_dir(preprocessing_dir, "split_objects")
  singlet_dir <- make_dir(preprocessing_dir, "singlet_objects")
  df_plot_dir <- make_dir(preprocessing_dir, "doubletfinder_plots")
  analysis_dir <- make_dir(output_dir, "02_mito_analysis")
  table_dir <- make_dir(analysis_dir, "tables")
  plot_dir <- make_dir(analysis_dir, "plots")
  
  saveRDS(cfg, file.path(output_dir, "run_configuration.rds"))
  capture.output(utils::sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  
  if (!is.null(cfg$resume_processed_file)) {
    message("Resuming from processed singlets; SCTransform/DoubletFinder will be skipped...")
    obj <- load_seurat_object(cfg$resume_processed_file)
    if (!"SCT" %in% Seurat::Assays(obj)) {
      stop("The resume object has no SCT assay and cannot be used.")
    }
    assert_metadata(obj, c(cfg$sample_id, cfg$group_by))
    obj <- join_assay_layers(obj, "RNA")
    obj <- prepare_mito_fraction(obj, cfg)
  } else {
    message("Loading input Seurat object...")
    obj <- load_seurat_object(
      path = cfg$input_file,
      object_name = cfg$input_object_name,
      h5_feature_type = cfg$h5_feature_type,
      h5_project = cfg$h5_project,
      h5_cell_metadata = cfg$h5_cell_metadata,
      h5_metadata_file = cfg$h5_metadata_file,
      h5_metadata_barcode_col = cfg$h5_metadata_barcode_col
    )
    DefaultAssay(obj) <- "RNA"
    assert_metadata(obj, c(cfg$split_by, cfg$sample_id, cfg$group_by))
    if (!all(c("nFeature_RNA", "nCount_RNA") %in% colnames(obj[[]]))) {
      stop("Input metadata must contain nFeature_RNA and nCount_RNA.")
    }
    raw_counts <- get_layer(obj, "RNA", "counts")
    if (nrow(raw_counts) == 0L || ncol(raw_counts) == 0L) stop("RNA raw counts are empty.")
    
    obj <- prepare_mito_fraction(obj, cfg)
    keep <- rownames(obj[[]])[
      obj$nFeature_RNA > cfg$nfeature_min & obj$nFeature_RNA < cfg$nfeature_max
    ]
    if (length(keep) == 0L) stop("No cells remain after the nFeature_RNA filter.")
    obj <- subset(obj, cells = keep)
    
    split_values <- unique(as.character(obj[[cfg$split_by]][, 1]))
    split_values <- split_values[!is.na(split_values)]
    if (length(split_values) == 0L) stop("No valid levels found in cfg$split_by.")
    message("Processing ", length(split_values), " split level(s): ",
            paste(split_values, collapse = ", "))
    
    # Save then process one split at a time. This is slower than keeping every
    # object in memory but is safer for large datasets and portable to HPC use.
    split_files <- setNames(character(length(split_values)), split_values)
    for (split_name in split_values) {
      cells <- rownames(obj[[]])[as.character(obj[[cfg$split_by]][, 1]) == split_name]
      split_obj <- subset(obj, cells = cells)
      path <- file.path(split_object_dir, paste0(safe_name(split_name), "_input.rds"))
      saveRDS(split_obj, path, compress = FALSE)
      split_files[[split_name]] <- path
    }
    rm(obj)
    invisible(gc())
    
    processed <- vector("list", length(split_values))
    names(processed) <- split_values
    summaries <- vector("list", length(split_values))
    names(summaries) <- split_values
    for (split_name in split_values) {
      split_obj <- readRDS(split_files[[split_name]])
      result <- run_doubletfinder_split(
        split_obj, split_name, cfg, df_plot_dir, singlet_dir
      )
      processed[[split_name]] <- result$object
      summaries[[split_name]] <- result$summary
      rm(split_obj, result)
      invisible(gc())
    }
    doublet_summary <- dplyr::bind_rows(summaries)
    write_csv_safe(doublet_summary, file.path(preprocessing_dir, "doublet_summary.csv"))
    
    message("Merging singlet objects...")
    obj <- merge_singlets(processed)
    rm(processed, summaries)
    invisible(gc())
    obj <- join_assay_layers(obj, "RNA")
    obj <- prepare_mito_fraction(obj, cfg)
    saveRDS(obj, file.path(output_dir, "processed_singlets.rds"), compress = FALSE)
  }
  
  message("Preparing log-normalized RNA data...")
  DefaultAssay(obj) <- "RNA"
  obj <- Seurat::NormalizeData(
    obj, normalization.method = "LogNormalize", scale.factor = 10000,
    verbose = FALSE
  )
  obj <- join_assay_layers(obj, "RNA")
  
  cutoffs <- make_cutoffs(obj$mito_fraction, cfg)
  strata <- make_strata(obj, cfg)
  message("Calculating cumulative threshold-response metrics...")
  metrics <- metric_trajectories(obj, strata, cutoffs, cfg)
  pca_metrics <- pca_bin_trajectories(obj, strata, cutoffs, cfg)
  variance_available <- any(is.finite(metrics$mean_variance_RNA)) ||
    any(is.finite(metrics$mean_variance_SCT)) || nrow(pca_metrics) > 0L
  if (!variance_available) {
    warning(
      "No mitochondrial bins passed cfg$min_cells_per_bin = ",
      cfg$min_cells_per_bin, " for the variance analyses. Non-variance ",
      "trajectories will still be saved for every observed bin."
    )
  }
  write_csv_safe(metrics, file.path(table_dir, "cumulative_metrics.csv"))
  write_csv_safe(pca_metrics, file.path(table_dir, "pca_bin_metrics.csv"))
  
  gene_expression <- NULL
  if (isTRUE(cfg$run_gene_expression)) {
    message("Calculating cumulative gene-expression trajectories...")
    gene_expression <- gene_expression_trajectories(obj, strata, cutoffs, cfg)
    write_csv_safe(gene_expression, file.path(table_dir, "gene_expression_trajectories.csv"))
    write_csv_safe(
      dplyr::filter(gene_expression, .data$category != "Lineage marker"),
      file.path(table_dir, "qc_gene_expression_trajectories.csv")
    )
    write_csv_safe(
      dplyr::filter(gene_expression, .data$category == "Lineage marker"),
      file.path(table_dir, "lineage_marker_trajectories.csv")
    )
  }
  
  auc_cell <- NULL
  auc_pseudobulk <- NULL
  if (isTRUE(cfg$run_aucell)) {
    message("Calculating AUCell trajectories...")
    auc_result <- aucell_cell_trajectories(obj, strata, cutoffs, cfg)
    auc_cell <- auc_result$trajectories
    write_csv_safe(auc_cell, file.path(table_dir, "aucell_cell_trajectories.csv"))
    if (isTRUE(cfg$run_pseudobulk_aucell)) {
      auc_pseudobulk <- aucell_pseudobulk_trajectories(
        obj, cutoffs, auc_result$gene_sets, cfg
      )
      write_csv_safe(
        auc_pseudobulk,
        file.path(table_dir, "aucell_pseudobulk_trajectories.csv")
      )
    }
  }
  
  message("Detecting transition points...")
  transitions <- detect_all_transitions(metrics, pca_metrics)
  write_csv_safe(transitions$calls, file.path(table_dir, "transition_calls.csv"))
  write_csv_safe(transitions$composites, file.path(table_dir, "transition_composites.csv"))
  write_csv_safe(transitions$consensus, file.path(table_dir, "consensus_transitions.csv"))
  
  save_standard_plots(
    metrics, pca_metrics, gene_expression, auc_cell, auc_pseudobulk,
    transitions, plot_dir
  )
  saveRDS(
    list(
      metrics = metrics,
      pca_metrics = pca_metrics,
      gene_expression = gene_expression,
      aucell_cell = auc_cell,
      aucell_pseudobulk = auc_pseudobulk,
      transitions = transitions
    ),
    file.path(analysis_dir, "mitothreshold_results.rds")
  )
  message("\nPipeline completed successfully. Results: ", output_dir)
  invisible(list(object = obj, transitions = transitions))
}

if (sys.nframe() == 0L && Sys.getenv("MITOTHRESH_NO_RUN", unset = "0") != "1") {
  pipeline_result <- run_pipeline(cfg)
}
