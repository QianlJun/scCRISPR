cells_need <- function(filepath, output = NULL, genecut = 30){
  if (is.null(output)){
    output = paste(getwd(), "/output/cells_need", sep = "")
    dir.create(output)
  } else {
    output = paste(output, "/output/cells_need", sep = "")
    dir.create(output)
  }

  if (is.null(nontarget)){
    print("Please tell the name of non-targert gene use parameter:nontarget")
  }

  scenicOptions=readRDS(file="int/scenicOptions.Rds")
  aucell_regulonAUC <- SCENIC::loadInt(scenicOptions, "aucell_regulonAUC")
  aucell.mt <- aucell_regulonAUC@assays@data$AUC

  # default keep target genes(not sgRNA) with >= 30 cells support
  seurat.obj <- readRDS(filepath)
  DefaultAssay(seurat.obj) <- "RNA"
  cellInfo <- seurat.obj@meta.data
  tmp <- table(seurat.obj$gene)
  genes.keep <- names(tmp)[tmp >= genecut]
  cells.keep <- WhichCells(seurat.obj, expression = gene %in% genes.keep)

  aucell.mt <- aucell.mt[,cells.keep]

  sgRNA <- unique(as.vector(seurat.obj$sgRNA_name))

  for (x in 1:length(sgRNA)){
    i = sgRNA[x]
    aucell.mt.tmp <- as.data.frame(aucell.mt[,sample(WhichCells(seurat.obj, expression = sgRNA_name == i), 1)])
    colnames(aucell.mt.tmp) <- 1
    for (j in c(3,5,10,15,20,30,40,50,60,70,80,90,100,120,130,150,170,185,200,220)){
      if (j <= length(WhichCells(seurat.obj, expression = sgRNA_name == i))){
        mt.tmp <- as.data.frame(rowMeans(as.data.frame(aucell.mt[,sample(WhichCells(seurat.obj, expression = sgRNA_name == i), j)])))
        colnames(mt.tmp) <- j
        aucell.mt.tmp <- cbind(aucell.mt.tmp, mt.tmp)
      }
    }
    diff.abs.sum <- as.data.frame(sum(abs(aucell.mt.tmp[,2] - aucell.mt.tmp[,1])))
    for (m in 2:ncol(aucell.mt.tmp)){
      tmp <- as.data.frame(sum(abs(aucell.mt.tmp[,m] - aucell.mt.tmp[,m-1])))
      colnames(tmp) <- m
      diff.abs.sum <- cbind(diff.abs.sum,tmp)
    }
    diff.abs.sum <- diff.abs.sum[,-1]
    diff.abs.sum <- as.data.frame(t(diff.abs.sum))
    colnames(diff.abs.sum) <- "Difference"
    diff.abs.sum$cells <- rownames(diff.abs.sum)
    diff.abs.sum$cells <- factor(diff.abs.sum$cells, levels = diff.abs.sum$cells)
    p1 <- ggplot(diff.abs.sum, aes(x=cells,y=Difference)) + geom_point() + theme_bw()
    ggsave(paste(output, "/",i,".png",sep = ""),p1,width = 6,height = 4)
  }

}
