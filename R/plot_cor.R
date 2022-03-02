# compute correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
}

plot_cor <- function(filepath, nontarget = NULL, output, sgRNAcut = 30){
  if (is.null(output)){
    output = paste(getwd(), "/output/plot_cor", sep = "")
    dir.create(output)
  } else {
    output = paste(output, "/output/plot_cor", sep = "")
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

  cormat.p <- round(cor(sg.mt.nodup, method = "pearson"),2)
  upper_tri.p <- get_upper_tri(cormat.p)
  melted_cormat.p <- reshape2::melt(upper_tri.p,na.rm = TRUE)
  cormat.s <- round(cor(sg.mt.nodup, method = "spearman"),2)
  upper_tri.s <- get_upper_tri(cormat.s)
  melted_cormat.s <- reshape2::melt(upper_tri.s,na.rm = TRUE)

  # Pearson
  ggheatmap.p <- ggplot2::ggplot(melted_cormat.p, aes(Var2, Var1, fill = value))+
    ggplot2::geom_tile(color = "white")+
    ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0.75, limit = c(0.5,1), space = "Lab",
                         name="Pearson\nCorrelation") +
    ggplot2::theme_minimal()+ # minimal theme
    ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1))+
    ggplot2::coord_fixed()
  p1<-ggheatmap.p +
    ggplot2::geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    ggplot2::theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          legend.justification = c(1,0))+ggplot2::scale_y_discrete(position = "right")+ggplot2::theme(legend.direction="horizontal",legend.position=c(0.4,0.9))
  # Spearman
  ggheatmap.s <- ggplot2::ggplot(melted_cormat.s, aes(Var2, Var1, fill = value))+
    ggplot2::geom_tile(color = "white")+
    ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0.75, limit = c(0.5,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1))+
    ggplot2::coord_fixed()
  p2<-ggheatmap.s +
    ggplot2::geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    ggplot2::theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          legend.justification = c(1,0)) + ggplot2::scale_y_discrete(position = "right")+ggplot2::theme(legend.direction="horizontal",legend.position=c(0.4,0.9))

  pdf(paste(output, "/sgRNA_correlation.pdf"),width = 12,height = 6)
  cowplot::plot_grid(ncol = 2,p1,p2)
  dev.off()
}
