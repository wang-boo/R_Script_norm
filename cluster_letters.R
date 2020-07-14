mat_letters = matrix(sample(letters[1:4],100,replace=TRUE),10)

# distance in the ASCII table
dist_letters=function(x,y){
    x = strtoi(charToRaw(paste(x,collapse="")),base=16)
    y = strtoi(charToRaw(paste(y,collapse="")),base=16)
    sqrt(sum((x-y)^2))
}
Heatmap(mat_letters,name="foo",col=structure(2:5,names=letters[1:4]),
        clustering_distance_rows=dist_letters,clustering_distance_columns=dist_letters,
        cell_fun=function(j,i,x,y,w,h,col){
            grid.text(mat_letters[i,j],x,y)
        })
