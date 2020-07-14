#### correlation #######################################################
library("GGally", lib.loc = lib.loc)
aft <- pro_quant_df[, Group.all]
bf <- pro_quant_df[, sample.list$raw]
names(bf) <- names(aft)

my_aft <- function(data, mapping, ...) {
  ggplot(data = log2(aft), mapping = mapping) +
    geom_point(
      color = "tomato",
      alpha = 0.5,
      size = 0.1,
      shape = 15
    ) +
    geom_density2d(
      alpha = 0.3,
      color = "grey10",
      bins = 5,
      size = 0.1
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = "tomato4",
      size = 0.25
    )
}
my_bf <- function(data, mapping, ...) {
  ggplot(data = log2(bf), mapping = mapping) +
    geom_point(
      ...,
      color = "olivedrab3",
      alpha = 0.5,
      size = 0.1,
      shape = 15
    ) +
    geom_density2d(
      alpha = 0.3,
      color = "grey10",
      bins = 5,
      size = 0.1
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = "olivedrab4",
      size = 0.25
    )
}

 correlation <-  ggpairs(
  data = log(bf),
  lower = list(continuous = my_bf),
  upper = list(continuous = my_aft),
  diag = NULL
  
) + theme_bw()
  
# ggsave("correlation.svg",cor,width = 6,height = 6,device = "svg")
 ggsave(
   paste(dir_1,"correlation.pdf",sep = ""),
   correlation,
   width = length(sample.list$bioID),
   height = length(sample.list$bioID),
   device = "pdf"
 )
 # ggsave(
 #   paste(dir_1,paste("correlation",fig.type,sep = "."),sep = ""),
 #   correlation,
 #   width = length(sample.list$bioID),
 #   # dpi = 144,
 #   height = length(sample.list$bioID),
 #   device = "tiff"
 # )
 

