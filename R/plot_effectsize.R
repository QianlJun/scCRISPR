plot_effectsize <- function(filepath, nontarget = NULL, output, sgRNAcut = 30){
  if (is.null(output)){
    output = paste(getwd(), "/output/plot_effectsize", sep = "")
    dir.create(output)
  } else {
    output = paste(output, "/output/plot_effectsize", sep = "")
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

  sg.mt <- as.data.frame(rowMeans(aucell.mt[,rownames(cellInfo)[cellInfo$sgRNA_name == unique(cellInfo$sgRNA_name)[1]]]))
  for (i in unique(cellInfo$sgRNA_name)){
    if (length(rownames(cellInfo)[cellInfo$sgRNA_name == i]) > 1){
      a <- as.data.frame(rowMeans(aucell.mt[,rownames(cellInfo)[cellInfo$sgRNA_name == i]]))
      colnames(a) <- i
      sg.mt <- cbind(sg.mt,a)
    }
  }
  sg.mt <- sg.mt[,-1]

  annotation_col <- seurat.obj@meta.data[,c(4,6)]
  anno.col <- annotation_col[!duplicated(annotation_col$sgRNA_name),]
  rownames(anno.col) <- anno.col$sgRNA_name

  sg.mt.nodup <- sg.mt
  regu.tmp <- strsplit(rownames(sg.mt.nodup),split = "_")
  regu.ve <- c()
  for (i in 1:length(regu.tmp)){
    regu.ve <- c(regu.ve, regu.tmp[[i]][1])
  }
  regu.tmp2 <- strsplit(regu.ve,split = " ")
  regu.ve2 <- c()
  for (i in 1:length(regu.tmp2)){
    regu.ve2 <- c(regu.ve2, regu.tmp2[[i]][1])
  }
  sg.mt.nodup$regulon <- regu.ve2
  sg.mt.nodup <- sg.mt.nodup[!duplicated(sg.mt.nodup$regulon),]
  rownames(sg.mt.nodup) <- sg.mt.nodup$regulon
  sg.mt.nodup <- sg.mt.nodup[,-ncol(sg.mt.nodup)]

  sg.mt.nodup.gene <- as.data.frame(rowMeans(sg.mt.nodup[,anno.col$target_gene == unique(anno.col$target_gene)[1]]))
  for (i in unique(anno.col$target_gene)){
    if (length(anno.col$barcode[anno.col$target_gene == i]) > 1){
      a <- as.data.frame(rowMeans(sg.mt.nodup[,anno.col$target_gene == i]))
      colnames(a) <- i
      sg.mt.nodup.gene <- cbind(sg.mt.nodup.gene,a)
    } else {
      a <- as.data.frame(sg.mt.nodup[,anno.col$target_gene == i])
      colnames(a) <- i
      sg.mt.nodup.gene <- cbind(sg.mt.nodup.gene,a)
    }
  }
  sg.mt.nodup.gene <- sg.mt.nodup.gene[,-1]

  regulon.sel <- unique(anno.col$target_gene)

  for (i in regulon.sel){
    tmp <- sg.mt.nodup.gene[i,]
    tmp <- tmp - tmp[,which(colnames(tmp) == nontarget)]
    tmp <- as.data.frame(t(tmp))
    tmp$gene <- rownames(tmp)
    tmp$gene <- factor(tmp$gene, levels = tmp$gene[order(tmp[,1])])
    p4 <- ggpubr::ggscatter(tmp, x = "gene", y = i,
                    color = i,
                    size = 1,
                    label = tmp$gene,
                    font.label = 9,
                    repel = T,
                    xlab = "Perturbations",
                    ylab = "Effect size") + ggplot2::theme_base() + ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0, size = 0)) + ggplot2::geom_hline(yintercept = 0, linetype="dashed") +
      ggplot2::scale_color_gradientn(values = c(0,0.2,0.4,0.6,0.8,1), colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25")) + ggplot2::ggtitle(i)
    pdf(paste(output,"/",i," effect_size.pdf",sep = ""), width = 10, height = 6)
    print(p4)
    dev.off()
  }

}
