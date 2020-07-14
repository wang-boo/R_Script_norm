cor_mat = cor(mar)
od = hclust(dist(cor_mat))$order
cor_mat = cor_mat[od,od]
nm = rownames(cor_mat)
col_fun = circlize::colorRamp2(c(-1,0,1),c("green","white","red"))

# `col = col_fun` here is used to generate the legend
Heatmap(cor_mat,
        name="correlation",
        col=col_fun, 
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width = width, height= height, fill){
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col="grey", fill = NA))
            if(i == j){
                grid.text(nm[i], x = x, y = y)
            } else if (i > j){
                grid.circle(x = x, y = y, r = abs(cor_mat[i,j])/2 * min(unit.c(width,height)),
                gp = gpar(fill = col_fun(cor_mat[i,j]), col = NA))
            }else{
                grid.text(sprintf("%.1f", cor_mat[i,j]), x, y, gp = gpar(fontsize = 8))
            }
},
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names=FALSE
)