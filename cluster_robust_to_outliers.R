robust_dist=function(x,y){
 qx = quantile(x,c(0.1,0.9))
 qy=quantile(y,c(0.1,0.9))
 l=x>qx[1] & x < qx[2] & y>qy[1] & y < qy[2]
 x =x[l]
 y=y[l]
 sqrt(sum((x-y)^2))
 }

Heatmap(mat_with_outliers,name="foo",
col=colorRamp2(c(-3,0,3),
c("green","white","red")),
clustering_distance_rows=robust_dist,
clustering_distance_columns=robust_dist)
