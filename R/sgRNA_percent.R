sgRNA_percent <- function(sgcount, barcode, output = NULL){
  if (is.null(output)){
    output = paste(getwd(), "/output/sgRNA_percent", sep = "")
    dir.create(output)
  } else {
    output = paste(output, "/output/sgRNA_percent", sep = "")
    dir.create(output)
  }

  `%notin%` <- Negate(`%in%`)
  sgcount <- read.table(sgcount, header = T, sep = "\t")
  barcodes <- read.table(barcode)
  sgcount <- sgcount[which(sgcount$cell_barcode %in% barcodes),]
  counts <- as.data.frame(table(table(sgcount$cell)))
  counts <- rbind(data.frame(Var1=0, Freq=length(which(barcodes %notin% sgcount$cell_barcode))), counts)
  counts$Var1 <- paste("cells get ", counts$Var1, " sgRNAs", sep = "")
  Labels <- as.vector(counts$Var1)
  Labels <- paste(Labels, "(", round(counts$Freq/sum(counts$Freq)*100, 2), "%)   ", sep = "")
  write.table(counts, paste(output, "/cells get sgRNA percentage.txt", sep = ""), col.names = T, row.names = F, sep = "\t")
  p1 <- ggplot2::ggplot(counts, ggplot2::aes(x = "", y=Freq, fill=Var1)) + ggplot2::geom_bar(stat="identity", width = 0.5, position = "stack") + ggplot2::coord_polar("y",start=0) + ggplot2::theme_minimal()
  p2 <- p1 + ggplot2::theme(axis.ticks = ggplot2::element_blank()) + ggplot2::theme(legend.title=ggplot2::element_blank(), legend.position = "right") + ggplot2::scale_fill_discrete(labels = Labels)
  p3 <- p2 + ggplot2::theme(axis.text.x = ggplot2::element_blank()) + ggplot2::theme(panel.grid=ggplot2::element_blank()) + ggplot2::theme(panel.border=ggplot2::element_blank()) + ggplot2::xlab("") + ggplot2::ylab(paste("Total cells: ", sum(counts$Freq), sep = ""))
  pdf(paste(output, "/cells get sgRNA percentage.pdf", sep = ""), width = 5, height = 5)
  p3
  dev.off()
}
