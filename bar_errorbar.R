require(ggplot2, lib.loc = lib.loc)
require(Cairo, lib.loc = lib.loc)
# options()
# bar.plot.symbol <- sig$`Case-Ctrl`$Symbol[8:16]
bar_plot <-
  subset.data.frame(
    pro_quant_df,
    subset = Symbol %in% bar.plot.symbol,
    select = c("Symbol", sample.list$bioID)
  )
rownames(bar_plot) <- bar_plot$Symbol


bar_m <- melt(bar_plot)
bar_m$group <-
  stri_extract(
    str = bar_m$variable,
    regex = sample.list$group.pattern,
    simplify = T,
    mode = sample.list$group.mode
  )


mean_mt <-
  acast(
    data = na.omit(bar_m),
    formula = Symbol ~ group,
    fun.aggregate = mean
  )
bar_plot_ <- bar_plot[rownames(mean_mt), ]
bar_plot_[sample.list$bioID] <- bar_plot_[-1] / mean_mt[, group.1]
###

bar_m_ <- melt(bar_plot_)
bar_m_$group <-
  stri_extract(
    str = bar_m_$variable,
    regex = sample.list$group.pattern,
    simplify = T,
    mode = sample.list$group.mode
  )


sd_mt_ <-
  acast(
    data = na.omit(bar_m_),
    formula = Symbol ~ group,
    fun.aggregate = sd
  )

mean_mt_ <-
  acast(
    data = na.omit(bar_m_),
    formula = Symbol ~ group,
    fun.aggregate = mean
  )
bar.plot.df <-
  data.frame(avg = melt(data.frame(gene = rownames(mean_mt_), mean_mt_)),
             sd = melt(data.frame(gene = rownames(sd_mt_), sd_mt_)))
bar.plot.df$avg.variable <-
  factor(
    bar.plot.df$avg.variable,
    levels = sample.list$group.names.order,
    ordered = T
  )
bar.plot.df$sd.variable <-
  factor(bar.plot.df$sd.variable,
         levels = sample.list$group.names.order,
         ordered = T)


b <-
  ggplot(bar.plot.df, aes(x = avg.gene, y = avg.value, fill = avg.variable))
b_plot <- b +
  geom_errorbar(
    aes(
      ymax = avg.value + sd.value,
      ymin = avg.value ,
      group = avg.variable
    ),
    position = position_dodge(0.6),
    color = "black",
    alpha = 0.75,
    size = 0.2,
    width = 0.4
  ) +
  scale_fill_brewer(type = "qual",
                    direction = 1,
                    palette = 2) +
  geom_bar(
    stat = "identity",
    color = "white",
    alpha = 1,
    size = 0.25,
    position = position_dodge(0.6),
    width = 0.6
  ) +
  coord_trans(y = "sqrt", limy = c(0, max(
    bar.plot.df$avg.value + bar.plot.df$sd.value
  ))) +
  xlab("Gene symbol") +
  ylab("Fold change") +
  guides(alpha = "none",
         size = "none",
         fill = guide_legend(title = ""))+
  theme(
    legend.text = element_text(size = 12,family = "Times New Roman"),
    axis.title = element_text(size = 14,family = "Times New Roman"),
    axis.text.x = element_text(size = 12,angle = 45,hjust = 1)
    
  )


ggsave(
  device = cairo_pdf,
  filename = "error.bar.plot.pdf",
  plot = b_plot,
  width = max(8,length(bar.plot.symbol)),
  height = 5,
  dpi = dpi.set
)
ggsave(
  filename = paste("error.bar.plot",fig.type,sep = "."),
  plot = b_plot,
  width = max(8,length(bar.plot.symbol)),
  height = 5,
  dpi = dpi.set
)

# sig$`Case-Ctrl`$Symbol
