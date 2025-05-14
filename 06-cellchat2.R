# Cellchat pipeline with debug
# Updated 06-cellchat.R for Seurat 5.0.2 and CellChat 2.1.2

Sys.setenv(R_MAX_VSIZE = "100GB")
options(bitmapType = "cairo")
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("CellChat"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("Cairo"))

option_list <- list(
    make_option(c("-i", "--input"),
                type = "character",
                help = "Input Seurat RDS file"),
    make_option(c("-s", "--species"),
                type = "character",
                help = "Species of cells (human or mouse)"),
    make_option(c("-t", "--tissue"),
                type = "character",
                help = "Tissue type for the sc-type analysis"),
    make_option(c("-o", "--output"),
                type = "character",
                help = "Output file path for cell type scoring results"),
    make_option(c("-d", "--outputDir"),
                type = "character",
                help = "Output directory for plots")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$species)) {
    print_help(opt_parser)
    stop("All arguments must be supplied (input file, species, tissue, output file, and output directory).\n", call. = FALSE)
}

message("Starting CellChat analysis")

# Load the Seurat object
scrna_seurat <- readRDS(opt$input)

# Debug the matrix
print("Assays in the Seurat object:")
print(Assays(scrna_seurat))
print("Metadata columns:")
print(colnames(scrna_seurat@meta.data))
print("Head of RNA assay data:")
print(head(scrna_seurat@assays$RNA@data[, 1:5]))

# Create the CellChat object.
data.input <- scrna_seurat@assays$RNA@data
cellchat <- createCellChat(object = data.input, meta = scrna_seurat@meta.data, group.by = "customclassif")

# Set the CellChat database based on species
if (opt$species == "human") {
    CellChatDB.use <- CellChatDB.human
} else if (opt$species == "mouse") {
    CellChatDB.use <- CellChatDB.mouse
} else {
    stop("Invalid species. Please specify 'human' or 'mouse'.")
}
cellchat@DB <- CellChatDB.use

# Debug: Check signaling genes in the database
print("Signaling genes in the database:")
print(head(CellChatDB.use$gene))

# Check if signaling genes are present in the data
signaling_genes <- intersect(rownames(cellchat@data), CellChatDB.use$gene$Symbol)
print("Number of signaling genes found in the data:")
print(length(signaling_genes))

# Subset the data to include only signaling genes
cellchat <- subsetData(cellchat)

# Debug: Check dimensions and summary of the signaling data
print("Dimensions of data.signaling matrix:")
print(dim(cellchat@data.signaling))
print("Summary of data.signaling matrix:")
print(summary(cellchat@data.signaling))

# Ensure metadata is assigned correctly
cellchat@meta <- scrna_seurat@meta.data

# Identify overexpressed genes and interactions, compute communication probability, and filter communications
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, 
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, 
                 title.name = "Interaction weights/strength")

# Create output directory for CellChat results
output_folder <- "./"
cellchat_folder <- file.path(output_folder, "Cellchat")
if (!dir.exists(cellchat_folder)) {
  dir.create(cellchat_folder, recursive = TRUE)
}

ptm <- Sys.time()
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
execution.time <- Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd = TRUE)

# Use CairoPNG for interaction network (Number)
CairoPNG(filename = file.path(cellchat_folder, 'interaction_network_Number.png'),
         width = 1600, height = 1600, res = 150)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge = FALSE, vertex.label.cex = 1.5, title.name = "Number of interactions")
dev.off()

# Use CairoPNG for interaction network (Weights)
CairoPNG(filename = file.path(cellchat_folder, 'interaction_network_weights.png'),
         width = 1600, height = 1600, res = 150)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, 
                 vertex.label.cex = 1.5, title.name = "Interaction weights/strength")
dev.off()

# Generate PDF version of the interaction networks
pdf(file.path(cellchat_folder, "interaction_network.pdf"), width = 50, height = 25)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, 
                 vertex.label.cex = 1.5, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, 
                 vertex.label.cex = 1.5, title.name = "Interaction weights/strength")
dev.off()

# Use CairoPNG for the detailed interaction network plots
CairoPNG(filename = file.path(cellchat_folder, 'interaction_network1.png'),
         width = 3000, height = 2000, res = 150)
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), 
                   vertex.label.cex = 1.5, title.name = rownames(mat)[i])
}
dev.off()

# Generate PDF for the detailed interaction network plots
pdf(file.path(cellchat_folder, "interaction_network1.pdf"), width = 50, height = 25)
par(mfrow = c(3,4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), 
                   vertex.label.cex = 1.5, title.name = rownames(mat)[i])
}
dev.off()

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver <- seq(1, 4)

# Loop over each pathway and use tryCatch to safely generate the hierarchy plots
for (i in 1:length(pathways.show.all)) {
  tryCatch({
    CairoPNG(filename = file.path(cellchat_folder, paste0(pathways.show.all[i], "_hierarchy_individual.png")),
             width = 3000, height = 2000, res = 100)
    netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
    dev.off()
    message("Successfully generated hierarchy plot for ", pathways.show.all[i])
  }, error = function(e) {
    message("Error generating hierarchy plot for ", pathways.show.all[i], ": ", e)
    # Optionally, you can write the error to a log file or take additional actions.
  })
  
  # Compute and save the contribution plot regardless of hierarchy plot errors
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename = file.path(cellchat_folder, paste0(pathways.show.all[i], "_L-R_contribution.pdf")),
         plot = gg, width = 6, height = 2, units = 'in', dpi = 150)
}

# Generate bubble plot for all significant interactions (after removing isolates and flipping axes)
num_cell_types <- length(levels(cellchat@idents))
pdf(file.path(cellchat_folder, "All_the_significant_interactions_Bubble_plot_filtered.pdf"), width = 10, height = 20)
netVisual_bubble(cellchat, sources.use = 1:num_cell_types, targets.use = 1:num_cell_types, remove.isolate = TRUE) +
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
        axis.text.y = element_text(size = 5))
dev.off()

message("CellChat analysis completed.")

