#### rank ###########################################################################################
x = get0(x = "go_ppi", ifnotfound = "all")

rank_plot <- x
rank <-
  cbind.data.frame(pro_quant_df, as.data.frame(sig.cutoff.list), label = sig.cutoff.list[[rank_plot]])
ylab <- bquote(log[10] ~ (iBAQ))
legend_title <-
  bquote(atop(Significant ~ 'in', italic(N) ~ Compare))
r <-
  ggplot(rank,
         aes(
           Rank,
           log10(iBAQ),
           color = factor(x = rank$label,ordered = T),
           size = label,
           alpha = label
         )) +
  ylab(ylab) +
  xlab("Protein rank") +
  ggtitle("Protein Abundance Range") +
  geom_rug() +
  geom_point(color = "black",
             alpha = 0.5,
             size = 0.5) +
  # facet_grid(facets = label ~ .)+
  # scale_shape_identity() +
  scale_color_discrete(direction = -1,
                       h = c(0, 360),
                       # c = 10,
                       # l = 60,
                       h.start = 110) +
  scale_alpha(range = c(0.0001, 1)) +
  scale_size(breaks = c(0, 1), range = c(0, 0.25)) +
  theme_light() +
  theme(
    aspect.ratio = 1,
    legend.background = element_rect(),
    legend.text = element_text(size = 5)
  ) +
  guides(
    alpha = "none",
    size = "none",
    color = guide_legend(title = legend_title)
  )

ggsave(
  device = cairo_pdf,
  filename = paste(dir_2,"rank.pdf",sep = ""),
  r,
  height = 4,
  scale = 1.5
)
ggsave(
  filename = paste(dir_2,"rank.", fig.type, sep = ""),
  dpi = dpi.set,
  r,
  height = 4,
  scale = 1.5
)
# setwd(dir_1)
