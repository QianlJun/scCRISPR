plot_ldr <- function(filepath, nontarget = NULL, output, sgRNAcut = 30, sel.sgRNA = NULL, sel.target = NULL){
  if (is.null(output)){
    output = paste(getwd(), "/output/plot_ldr", sep = "")
    dir.create(output)
  } else {
    output = paste(output, "/output/plot_ldr", sep = "")
    dir.create(output)
  }

  if (is.null(nontarget)){
    print("Please allocate the name of non-targert gene use parameter:nontarget")
  }

  seurat.obj <- readRDS(filepath)

  # filter sgRNA  < 30 cells support
  tmp <- table(seurat.obj$sgRNA_name)
  sgRNA.keep <- names(tmp)[tmp >= sgRNAcut]
  cells.keep <- Seurat::WhichCells(seurat.obj, expression = sgRNA_name %in% sgRNA.keep)
  seurat.obj <- subset(seurat.obj, cells = cells.keep)
  cellInfo <- seurat.obj@meta.data

  scenicOptions=readRDS(file="int/scenicOptions.Rds")
  aucell_regulonAUC <- SCENIC::loadInt(scenicOptions, "aucell_regulonAUC")
  aucell.mt <- aucell_regulonAUC@assays@data$AUC
  aucell.mt <- aucell.mt[,cells.keep]

  aucell.obj <- Seurat::CreateAssayObject(aucell.mt)
  seurat.obj[["regulonAUC"]] <- aucell.obj
  DefaultAssay(seurat.obj) <- "regulonAUC"

  # Feature selection
  seurat.obj <- Seurat::FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 50, assay = "regulonAUC", verbose = FALSE)
  # Scaling the data
  seurat.obj <- Seurat::ScaleData(seurat.obj, assay = "regulonAUC")
  seurat.obj <- Seurat::RunPCA(seurat.obj, npcs = 30, reduction.name = "PCA_on_regulonAUC", assay = "regulonAUC", verbose = F)
  seurat.obj <- Seurat::RunUMAP(seurat.obj, reduction = "PCA_on_regulonAUC", dims = 1:20, reduction.name = "UMAP_on_regulonsAUC")
  seurat.obj <- Seurat::RunTSNE(seurat.obj, reduction = "PCA_on_regulonAUC", dims = 1:20, reduction.name = "TSNE_on_regulonsAUC")

  # figure umap
  pdf(paste(output, "/umap.pdf"),width = 8,height = 5)
  Seurat::DimPlot(seurat.obj, reduction = "UMAP_on_regulonsAUC", group.by = "target_gene", shuffle = TRUE)
  dev.off()

  # figure tsne
  pdf(paste(output, "/tsne.pdf"),width = 8,height = 5)
  Seurat::DimPlot(seurat.obj, reduction = "TSNE_on_regulonsAUC", group.by = "target_gene", shuffle = TRUE)
  dev.off()

  if (!is.null(sel.sgRNA)){
    seurat.obj <- Seurat::SetIdent(seurat.obj, value = "sgRNA_name")
    subdata <- subset(seurat.obj, idents = sel.sgRNA)

    # figure umap
    pdf(paste(output, "/select sgRNA umap.pdf"),width = 8,height = 5)
    DimPlot(subdata, reduction = "UMAP_on_regulonsAUC", group.by = "sgRNA_name", shuffle = TRUE)
    dev.off()

    # figure tsne
    pdf(paste(output, "/select sgRNA tsne.pdf"),width = 8,height = 5)
    DimPlot(subdata, reduction = "TSNE_on_regulonsAUC", group.by = "sgRNA_name", shuffle = TRUE)
    dev.off()
  }

  if (!is.null(sel.target)){
    seurat.obj <- Seurat::SetIdent(seurat.obj, value = "target_gene")
    subdata <- subset(seurat.obj, idents = sel.target)

    # figure umap
    pdf(paste(output, "/select target gene umap.pdf"),width = 8,height = 5)
    Seurat::DimPlot(subdata, reduction = "UMAP_on_regulonsAUC", group.by = "target_gene", shuffle = TRUE)
    dev.off()

    # figure tsne
    pdf(paste(output, "/select target gene tsne.pdf"),width = 8,height = 5)
    Seurat::DimPlot(subdata, reduction = "TSNE_on_regulonsAUC", group.by = "target_gene", shuffle = TRUE)
    dev.off()
  }

  saveRDS(seurat.obj, paste(output, "seuratobj.rds", sep = ""))
}
