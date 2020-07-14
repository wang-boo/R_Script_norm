library(circlize)
library(ComplexHeatmap)

set.seed(123)
a1 <- matrix(rnorm(16,-1),4)

a2 <- matrix(rnorm(32,1),8)
a <- rbind(a1,a2)
b1 <- matrix(rnorm(24,1),4)
b2 <- matrix(rnorm(48,-1),8)
b <- rbind(b1,b2)
mat <- cbind(a,b)
sample(nrow(mat),nrow(mat))
mat <- mat[sample(nrow(mat),nrow(mat)),sample(ncol(mat),ncol(mat))]
rownames(mat) = paste0("R",1:12)
colnames(mat) = paste0("R",1:10)
Heatmap(mat, col=colorRamp2(c(-3,0,3), c("green","white","red")), cluster_rows=FALSE, cluster_columns=FALSE)

Heatmap(mat,col=rev(rainbow(10))

colors = structure(circlize::rand_color(4),names=c("1","2","3","4"))
Heatmap(discrete_mat,col=colors)
