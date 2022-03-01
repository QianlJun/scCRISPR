plot_geneffect <- function(filepath, nontarget = NULL, output, sgRNAcut = 30, sgRNA = NULL,
                           targetgene = NULL, gene = NULL, test.use = "wilcox"){
  if (is.null(output)){
    output = paste(getwd(), "/output/plot_geneffect", sep = "")
    dir.create(output)
  } else {
    output = paste(output, "/output/plot_geneffect", sep = "")
    dir.create(output)
  }

  if (is.null(nontarget)){
    print("Please allocate the name of non-targert gene use parameter:nontarget")
  }

  seurat.obj <- readRDS(filepath)
  DefaultAssay(seurat.obj) <- "RNA"

  # filter sgRNA with < 30 cells support
  tmp <- table(seurat.obj$sgRNA_name)
  sgRNA.keep <- names(tmp)[tmp >= sgRNAcut]
  cells.keep <- Seurat::WhichCells(seurat.obj, expression = sgRNA_name %in% sgRNA.keep)
  seurat.obj <- subset(seurat.obj, cells = cells.keep)
  cellInfo <- seurat.obj@meta.data

  seurat.obj <- subset(seurat.obj, cells = cells.keep)
  cellInfo <- seurat.obj@meta.data

  nontargetsgRNA <- unique(seurat.obj$sgRNA_name[seurat.obj$target_gene == nontarget])

  genes <- unique(seurat.obj$target_gene)
  genes <- genes[-which(genes == nontarget)]

  rna.data <- GetAssayData(seurat.obj,assay = "RNA",slot = "data")
  rna.data <- as.data.frame(rna.data)


  rna.mt <- as.data.frame(rowMeans(rna.data[,rownames(cellInfo)[cellInfo$sgRNA_name == unqiue(cellInfo$sgRNA_name)[1]]]))
  for (i in unique(cellInfo$sgRNA_name)){
    a <- as.data.frame(rowMeans(rna.data[,rownames(cellInfo)[cellInfo$sgRNA_name == i]]))
    colnames(a) <- i
    rna.mt <- cbind(rna.mt,a)
  }
  rna.mt <- rna.mt[,-1]

  non.target <- as.data.frame(rowMeans(rna.mt[genes,nontargetsgRNA]))
  colnames(non.target) <- c("Non_target")

  bargene <- data.frame(barcode=cellInfo$sgRNA_name,gene=cellInfo$gene)
  bargene <- bargene[!duplicated(bargene$sgRNA_name),]
  bargene.filt <- bargene[bargene$target_gene %in% genes,]
  rna.barcode <- rna.mt[genes,bargene.filt$sgRNA_name[bargene.filt$target_gene %in% genes]]

  # removed sgRNAs
  sg.rm <- c()
  bargene.filt <- bargene.filt[!bargene.filt$barcode %in% sg.rm,]
  rna.barcode <- rna.mt[genes,bargene.filt$sgRNA_name[bargene.filt$target_gene %in% genes]]
  # sgRNA affect target gene expression
  for (i in genes){
    bc <- seurat.obj$barcode[seurat.obj$target_gene == i]
    df <- as.data.frame(rna.barcode[i,bc])
    rownames(df) <- i
    colnames(df) <- bc
    non.tar <- as.data.frame(non.target[i,])
    colnames(non.tar) <- "Non_target"
    rownames(non.tar) <- i
    tmp <- (cbind(non.tar,df) + 0.0001) / (non.target[i,] + 0.0001)
    tmp <- as.data.frame(t(tmp))
    colnames(tmp) <- "value"
    tmp$sgRNA <- factor(rownames(tmp),levels = rownames(tmp))
    p1 <- ggplot(tmp,aes(x=sgRNA,y=value)) + geom_point() + theme_bw() + geom_hline(aes(yintercept=1),colour="Red",linetype="dashed") + xlab("")
    ggsave(paste(output,"/",i,".png",sep = ""),p1,width = 4,height = 4)
    }

  # violin plot
  seurat.obj$tmp_label <- "Non-Targeting"
  seurat.obj$tmp_label[seurat.obj$target_gene == targetgene] <- paste(targetgene," KO",sep = "")
  seurat.obj$tmp_label <- factor(seurat.obj$tmp_label, levels = c("Non-Targeting", paste(targetgene," KO",sep = "")))
  pdf(paste(output,"/", targetgene,"_violinplot.pdf",sep = ""), width = 6, height = 4)
  VlnPlot(subset(seurat.obj, cells = WhichCells(seurat.obj, expression = target_gene %in% c(targetgene,"Non-Targeting"))), features = c(targetgene),group.by = "orig.ident", adjust = 2, assay = "RNA", pt.size = 0.1, split.by = "tmp_label", split.plot = TRUE)
  dev.off()

  # ko vs nontarget
  markers <- FindMarkers(seurat.obj, ident.1=sgRNA, ident.2=nontargetsgRNA ,test.use = test.use, min.pct = 0, logfc.threshold = 0, only.pos = FALSE, assay = "RNA")
  markers$gene <- rownames(markers)
  markers <- markers[,c(6,1:5)]
  markers$effect <- "none"
  markers$effect[which(markers$avg_log2FC > 0)] <- "up"
  markers$effect[which(markers$avg_log2FC < 0)] <- "down"
  markers$reliability <- 1
  markers$reliability[which(markers$p_val < 0.05)] <- 2
  markers$reliability[which(markers$p_val < 0.05 & markers$p_val_adj < 0.05)] <- 4
  markers$reliability[which(markers$p_val < 0.05 & markers$p_val_adj < 0.05 & markers$avg_log2FC > 0 & markers$`pct.1` > 0.2)] <- 5
  markers$reliability[which(markers$p_val < 0.05 & markers$p_val_adj < 0.05 & markers$avg_log2FC < 0 & markers$`pct.2` > 0.2)] <- 5
  markers$reliability[which(markers$p_val < 0.05 & markers$p_val_adj >= 0.05 & markers$avg_log2FC > 0 & markers$`pct.1` > 0.2)] <- 3
  markers$reliability[which(markers$p_val < 0.05 & markers$p_val_adj >= 0.05 & markers$avg_log2FC < 0 & markers$`pct.2` > 0.2)] <- 3
  write.table(markers, paste(output, "/markers.csv", sep = ""), sep=",", col.names = T, row.names = F, quote = F)
  # if (markers$p_val < 0.05){
  #     if (markers$p_val_adj < 0.05){
  #         if (markers$avg_log2FC > 0 && markers$`pct.1` > 0.2){
  #             markers$reliability <- 5
  #         } else if (markers$avg_log2FC < 0 && markers$`pct.2` > 0.2){
  #             markers$reliability <- 5
  #         } else {
  #             markers$reliability <- 4
  #         }
  #     } else {
  #         if (markers$avg_log2FC > 0 && markers$`pct.1` > 0.2){
  #             markers$reliability <- 3
  #         } else if (markers$avg_log2FC < 0 && markers$`pct.2` > 0.2){
  #             markers$reliability <- 3
  #         } else {
  #             markers$reliability <- 2
  #         }
  #     }
  # }

}
