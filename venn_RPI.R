require("VennDiagram", lib.loc = lib.loc)
  sig.. <- sig$class
if(is.null(up.vs.dn)) venn_ <- cpr$names else venn_ <- up.vs.dn
for (v in venn_) {
  subsig <- sig..[[v]]
  
  category <- 
    as.data.frame(t(stri_split(str = v, fixed = "-")[[1]]), stringsAsFactors = F)
  
  venn.plot <- draw.pairwise.venn(
    area1 = nrow(subsig[[1]]) + nrow(subsig[[3]]),
    area2 = nrow(subsig[[2]]) + nrow(subsig[[3]]),
    cross.area = nrow(subsig[[3]]),
    category = category,
    fill = c("lightblue", "green3"),
    lty = "blank",
    euler.d = T,
    # cex = 1,
    # cat.cex = 1,
    #   cat.dist = 0.09,
    # cat.just = list(c(-1, -1), c(1, 1)),
    # # ext.pos = 30,
    # ext.dist = -0.05,
    # ext.length = 0.85,
    # ext.line.lwd = 2,
    ext.line.lty = "dashed"
  )
  
  pdf(file = paste(dir_1,paste("venn", v, "pdf", sep = "."), sep = ""))
  grid.draw(venn.plot)
  dev.off()
  tiff(file = paste(dir_1,paste("venn", v, fig.type, sep = "."), sep = ""))
  grid.draw(venn.plot)
  dev.off()
}
