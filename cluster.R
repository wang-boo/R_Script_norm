x = get0(x = "go_ppi", ifnotfound = "all")
require(ComplexHeatmap,lib.loc = lib.loc)
sig_mt <- sig[[x]]
# names(pro_quant_df)
Group.sub <- names(sig_mt)[which(names(sig_mt) %in% Group.all)]
hm <- as.matrix(sig_mt[c(Group.sub)])
row.names(hm) <- sig_mt$Symbol
scl <- t(scale(t(log(hm))))
anno <- melt(scl, varnames = c("name", "group"))
df.. <- data.frame(
  'Bio.Group' = stri_extract(
    str = colnames(hm),
    regex =  sample.list$group.pattern,
    simplify = T,
    mode = sample.list$group.mode
  )
)
anno.col <-
  rainbow(
    n = length(sample.list$group.names.order),
    s = 0.7,
    v = 0.6,
    start = 0.2,
    end = 0.9
  )
names(anno.col) <- sample.list$group.names.order
ha <-
  HeatmapAnnotation(
    df = df..,
    show_legend = T,
    col = list(Bio.Group = anno.col),
    annotation_legend_param = list(
      grid_height = unit(6,"mm"),
      labels_gp = gpar(fontsize = 12),
      title_gp = gpar(fontsize = 12)
    )
  )

if(nrow(scl) > 70){
  rname <- 6
  height <- 9
  show.rnames <- F
}else if(nrow(scl) > 50){
  rname <- 7
  show.rnames <- T
  height <- 9
}else{
  rname <- 8
  height <- max(6,nrow(scl)/7)
  show.rnames <- T
}

if(ncol(scl) > 40){
  show.cnames <- F
}else if(ncol(scl) > 20){
  show.cnames <- T
  cname <- 11
}else{
  show.cnames <- T
  cname <- 12
}

cluster.heat <- Heatmap(
  heatmap_legend_param = list(
    legend_direction = c("vertical"),
    title_position = c("topcenter"),
    grid_height = unit(4,"mm"),
    labels_gp = gpar(fontsize = 12),
    title_gp = gpar(fontsize = 12)
  ),
  scl,
  width = unit(10,"cm"),
  # rect_gp = gpar(col = "white", lty = 1, lwd = ifelse(nrow(scl)>=50,NA,0.01)),
  name = "Z-score",
  row_names_gp = gpar(fontsize = rname),
  show_row_names = show.rnames,
  show_column_names = show.cnames,
  column_names_gp = gpar(fontsize = cname),
  top_annotation = ha
)

write.table(scl,"cluster.txt",sep="\t",row.names=FALSE)

pdf(
  file = paste(dir_2, "cluster.pdf", sep = ""),
  height = unit(height,"cm")
)
cluster.heat
dev.off()


tiff(
  res = 144,
  file = paste(dir_2, "cluster.tiff", sep = ""),
  units = "px",
  width = 1000,
  height = height*150
)
cluster.heat
dev.off()

