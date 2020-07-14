# x = get0(x = "go_ppi", ifnotfound = 2)

pro_quant_df_fc <- pro_quant_df
pro_quant_df_fc$log2.ratio <-
  log2(pro_quant_df_fc[paste("fc..", x, sep = "")])[, 1]
pro_quant_df_fc$regulation <- 0
pro_quant_df_fc$regulation[pro_quant_df_fc$log2.ratio >  lr_co &
                             pro_quant_df_fc$Unique.peptides >= unique_co] <- 1
pro_quant_df_fc$regulation[pro_quant_df_fc$log2.ratio < -lr_co &
                             pro_quant_df_fc$Unique.peptides >= unique_co] <- -1
fc.col.list <- list(col.low = "forestgreen",
                    col.mid = "black",
                    col.high = "tomato")

fc_plot <-
  ggplot(pro_quant_df_fc,
         aes(
           Rank,
           log2.ratio,
           fill = regulation,
           alpha = abs(regulation)
         )) +
  geom_point(
    shape = 21,
    color = "white",
    size = 3,
    stroke = 0.5
  ) +
  scale_fill_gradient2(
    low = fc.col.list$col.low,
    mid = fc.col.list$col.mid,
    high = fc.col.list$col.high
  ) +
  scale_alpha_continuous(range = c(0.1, 0.8)) +
  guides(
    alpha = "none",
    size = "none",
    fill = "none",
    color = "none"
  ) +
  geom_hline(
    yintercept = lr_co,
    color = fc.col.list$col.high,
    size = 0.5,
    alpha = 0.8
  ) +
  geom_hline(
    yintercept = -lr_co,
    color = fc.col.list$col.low,
    size = 0.5,
    alpha = 0.8
  ) +
  ylim(c(-max(abs(
    pro_quant_df_fc$log2.ratio
  )), max(abs(
    pro_quant_df_fc$log2.ratio
  )))) +
  ylab(bquote(italic(log[2] ~ ratio))) +
  xlab(bquote(italic(iBAQ ~ Rank)))+
  theme(axis.title = element_text(size = 14))



ggsave(
  dpi = dpi.set,
  units = "cm",
  scale = 0.5,
  paste(dir_2,paste("FC",fig.type,sep = "."),sep = ""),
  plot = fc_plot,
  width = 40,
  height = 30
)

ggsave(
  dpi = dpi.set,
  units = "cm",
  scale = 0.5,
  paste(dir_2,"FC.pdf",sep = ""),
  plot = fc_plot,
  width = 40,
  height = 30
)
