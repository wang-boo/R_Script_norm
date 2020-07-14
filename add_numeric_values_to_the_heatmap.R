Heatmap(mat, name="foo", cell_fun = function(j,i,x,y,width,height,fill){
    grid.text(sprintf("%.1f",mat[i,j]), x, y, gp=gpar(fontsize=10))
})
