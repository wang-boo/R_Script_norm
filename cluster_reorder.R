pushViewport(viewport(layout=grid.layout(nr=1,nc=3)))

pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
draw(Heatmap(mat,name="foo", row_dend_reorder=FALSE, column_title="no reordering"), newpage=FALSE)
upViewport()

pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
draw(Heatmap(mat,name="foo", row_dend_reorder=TRUE, column_title="applied reordering"), newpage=FALSE)
upViewport()

library(dendsort)
dend = dendsort(hcluster(dist(mat)))
pushViewport(viewport(layout.pos.row=1, layout.pos.col=3))
draw(Heatmap(mat,name="foo", cluster_rows=dend, row_dend_reorder=FALSE,
    column_title = "reordering by dendsort"), newpage=FALSE)
upViewport(2)