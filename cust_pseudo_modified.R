#!/usr/bin/env Rscript

# --- Global error handler ---
options(error = function() {
  msg <- geterrmessage()
  cat("✖ GLOBAL ERROR —", msg, "\n", file = stderr())
  q(save = "no", status = 1)
})

# --- Helper for safe steps ---
safe_step <- function(expr, step_name) {
  tryCatch(
    expr,
    error = function(e) {
      cat("✖ ERROR in", step_name, ":", e$message, "\n", file = stderr())
      q(save = "no", status = 1)
    }
  )
}

# --- Memory reporting ---
report_memory <- function(label) {
  m <- gc(reset = TRUE)[2, 2]
  cat(sprintf("📊 Memory after %s: %.1f MB\n", label, m))
}

# --- Load libraries ---
safe_step(library(Seurat),        "load Seurat")
safe_step(library(monocle3),      "load monocle3")
safe_step(library(ggplot2),       "load ggplot2")
safe_step(library(dplyr),         "load dplyr")
safe_step(library(SeuratWrappers),"load SeuratWrappers")
report_memory("libs loaded")

# --- Paths & cell types ---
input_file    <- "/mnt/NFS/pisces/chlu0571/10XSC/cellranger/cust_pseudo/rds/human_ncku_pdac_scRNA_merge_seurat4_20241121.rds"
output_folder <- "/mnt/NFS/pisces/chlu0571/10XSC/cellranger/cust_pseudo/trajectory_results_merged"

# Make sure this directory exists and is writable!
safe_step({
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
}, "create output directory")

cell_types <- c(
  "Acinar", "Ductal", "Epithelial_cells",
  "Tumor", "Tumor_EMT", "Tumor_liver", "Tumor_untreated"
)

# --- Read & subset using a different approach ---
safe_step({
  # Load the RDS file
  cat("📂 Loading RDS file:", input_file, "\n")
  scrna <- readRDS(input_file)
  
  # Print top metadata columns
  cat("\nMetadata columns:\n")
  print(head(colnames(scrna@meta.data)))
  
  # Check if orig.ident exists (common merge identifier)
  if ("orig.ident" %in% colnames(scrna@meta.data)) {
    cat("\nSamples in merged data (orig.ident):\n")
    print(table(scrna$orig.ident))
  }

  # Check if it's a valid Seurat object
  if (!inherits(scrna, "Seurat")) {
    stop("The loaded object is not a Seurat object")
  }
  
  # Print information about the object
  cat("📊 Loaded Seurat object with", nrow(scrna@meta.data), "cells and", nrow(scrna), "features\n")
  
  # Check if cell_type column exists
  if (!"cell_type" %in% colnames(scrna@meta.data)) {
    cat("⚠️ 'cell_type' column not found. Available columns are:", 
        paste(colnames(scrna@meta.data), collapse=", "), "\n")
    stop("cell_type column not found in metadata")
  }
  
  # Create a subset expression to use with the Seurat subset function
  cat("📊 Creating a new Seurat object with selected cell types\n")
  
  # APPROACH 1: Create a new Seurat object from scratch
  # This avoids using subset() which seems to cause errors
  
  # Create a logical vector for selecting cells
  selected_cells <- scrna@meta.data$cell_type %in% cell_types & !is.na(scrna@meta.data$cell_type)
  
  # Get count of selected cells
  num_selected <- sum(selected_cells)
  cat("📊 Selected", num_selected, "cells with specified cell types\n")
  
  if (num_selected == 0) {
    stop("No cells found with the specified cell types")
  }
  
  # Get cell barcodes
  cell_barcodes <- rownames(scrna@meta.data)[selected_cells]
  
  # Print first few barcodes for debugging
  cat("📊 First few barcodes:", paste(head(cell_barcodes, 3), collapse=", "), "...\n")
  
  # Instead of using subset(), manually extract the data
  # Create a new Seurat object with only the selected cells
  cat("📊 Extracting counts for selected cells...\n")
  counts <- scrna@assays$RNA@counts[, cell_barcodes]
  
  # Create a new Seurat object
  cat("📊 Creating new Seurat object...\n")
  new_scrna <- CreateSeuratObject(counts = counts)
  
  # Copy over necessary metadata
  for (col in colnames(scrna@meta.data)) {
    if (col %in% c("cell_type")) {  # You can add more columns as needed
      new_scrna@meta.data[[col]] <- scrna@meta.data[cell_barcodes, col]
    }
  }
  
  # Transfer reduction embeddings if they exist
  if ("umap" %in% names(scrna@reductions)) {
    cat("📊 Transferring UMAP embeddings...\n")
    new_scrna[["umap"]] <- CreateDimReducObject(
      embeddings = scrna@reductions$umap@cell.embeddings[cell_barcodes, ],
      key = "UMAP_",
      assay = DefaultAssay(scrna)
    )
  }
  
  # Set the active assay
  DefaultAssay(new_scrna) <- "RNA"
  
  # Replace the original object
  scrna <- new_scrna
  rm(new_scrna, counts)
  gc()
  
  cat("📊 Created new Seurat object with", ncol(scrna), "cells and", nrow(scrna), "features\n")
  
  report_memory("after subsetting")
}, "load + subset")

# --- Convert to Monocle3 CDS ---
safe_step({
  cat("📊 Converting to Monocle3 cell_data_set...\n")
  cds <- as.cell_data_set(scrna)
  cds$cell_type <- factor(cds@colData$cell_type)
  report_memory("after conversion")
}, "to CDS")

# --- Down-sample to ≤10k cells to avoid segfault ---
safe_step({
  max_cells <- 10000
  total     <- ncol(cds)
  if (total > max_cells) {
    set.seed(42)
    
    # Check if cell_type is available for stratification
    if ("cell_type" %in% names(colData(cds))) {
      # Stratified sampling by cell type to preserve proportions
      keep <- c()
      cell_types_present <- unique(cds$cell_type)
      
      cat("📊 Performing stratified down-sampling...\n")
      
      for (ct in cell_types_present) {
        ct_cells <- colnames(cds)[cds$cell_type == ct]
        ct_prop <- min(1, max_cells * length(ct_cells) / total)
        
        # Make sure we don't try to sample more cells than available
        sample_size <- min(ceiling(length(ct_cells) * ct_prop), length(ct_cells))
        ct_sample <- sample(ct_cells, sample_size)
        
        cat(sprintf("  - %s: %d → %d cells (%.1f%%)\n", 
            ct, length(ct_cells), length(ct_sample), 100 * length(ct_sample) / length(ct_cells)))
        keep <- c(keep, ct_sample)
      }
      
      # If we still have too many cells, do a final random sample
      if (length(keep) > max_cells) {
        keep <- sample(keep, max_cells)
        cat("  - Final random sample to", length(keep), "cells\n")
      }
    } else {
      # Simple random sampling if cell_type is not available
      cat("⚠️ cell_type not found in colData, performing random sampling\n")
      keep <- sample(colnames(cds), max_cells)
    }
    
    cds  <- cds[, keep]
    cat("⚡ Down-sampled from", total, "→", length(keep), "cells\n")
    report_memory("after down-sampling")
  }
}, "down-sample")

# --- Learn trajectory ---
safe_step({
  cat("📊 Running cluster_cells...\n")
  cds <- cluster_cells(cds, resolution = 0.5)
  
  cat("📊 Running learn_graph...\n")
  cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)
  
  # Check if we have Ductal cells
  ductal_cells <- colnames(cds)[cds$cell_type == "Ductal"]
  
  # If we have Ductal cells, use them as root
  if (length(ductal_cells) > 0) {
    cat("🔍 Using", length(ductal_cells), "Ductal cells as root\n")
    cds <- order_cells(cds, root_cells = ductal_cells)
  } else {
    # If no Ductal cells, use the first available cell type as root
    first_type <- levels(cds$cell_type)[1]
    root_cells <- colnames(cds)[cds$cell_type == first_type]
    cat("⚠️ No Ductal cells found. Using", length(root_cells), first_type, "cells as root instead\n")
    cds <- order_cells(cds, root_cells = root_cells)
  }
  
  report_memory("trajectory learned")
}, "learn_graph + order_cells")

# --- Plot all cells in one UMAP ---
safe_step({
  # by Cell Type
  cat("📊 Plotting cells by cell type...\n")
  p1 <- plot_cells(
    cds,
    color_cells_by = "cell_type",
    show_trajectory_graph   = TRUE,
    label_cell_groups       = FALSE,
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_roots = FALSE,
    label_branch_points = FALSE,
    cell_size = 1
  )
  ggsave(
    file.path(output_folder, "all_cell_types_trajectory.png"),
    p1, width = 10, height = 8
  )
  
  # by Pseudotime
  cat("📊 Plotting cells by pseudotime...\n")
  p2 <- plot_cells(
    cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph   = TRUE,
    label_cell_groups       = FALSE,
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_roots = FALSE,
    label_branch_points = FALSE,
    cell_size = 1
  )
  ggsave(
    file.path(output_folder, "all_pseudotime_trajectory.png"),
    p2, width = 10, height = 8
  )
}, "plot")

cat("✓ All trajectory analyses complete!\n")