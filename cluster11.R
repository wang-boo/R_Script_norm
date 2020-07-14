library(ComplexHeatmap)
library(circlize)   # for colorRamp2

setwd("D:/WD_1711125/combined/txt/")

cluster1 <- function(){
    mat <- read.table("mat.txt", sep="\t", header=T)
    mat <- t(scale(t(mat)))
    samples <- c("blank", "model", "positive", "treat")
    print(samples)
    write.table(mat, "mat.txt",sep="\t")
    anno.col = rainbow(
        n = length(samples),
        s = 0.7,
        v = 0.6,
        start = 0.2,
        end = 0.9
        )
    names(anno.col)=samples
    ha_column = HeatmapAnnotation(
        df = data.frame("Bio.Group" = samples),
        col = list(Bio.Group = anno.col),
        annotation_legend_param = list(
        grid_height = unit(6,"mm"),
        labels_gp = gpar(fontsize = 12),
        title_gp = gpar(fontsize = 12))
        )
    print(anno.col)
    print(ha_column)
    cluster.heat <- Heatmap(
                            heatmap_legend_param = list(
                            legend_direction = c("vertical"),
                            title_position = c("topcenter"),
                            grid_height = unit(4, "mm"),
                            labels_gp = gpar(fontsize = 6),
                            title_gp = gpar(fontsize = 6)
                            ),
                            mat,
                            col = colorRamp2(c(min(mat),max(mat)),c("green","red")),
                            name = "Z-score",
                            row_names_gp = gpar(fontsize = 6),
                            show_column_names = TRUE,
                            show_row_names = FALSE,
                            column_names_gp = gpar(fontsize = 8),
                            top_annotation = ha_column)
    height = 12
    savePlot(
    filename = "cluster",
    type = "pdf",
    device=dev.cur(),
    restoreConsole = TRUE
    )

}

cluster1()
