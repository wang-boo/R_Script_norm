# x = get0(x = "go_ppi", ifnotfound = 2)
require(ggrepel, lib.loc = lib.loc)
require(dplyr, lib.loc = lib.loc)
select_var <- function(compare.vc, data, rm.idx = F) {
  ii = compare.vc #<- x
  {
    sig_i <- data #<- pro_quant_df
    # names(sig_i)
    rm.. <-
      unique(sample.list$group[-which(sample.list$group %in% stri_split_fixed(ii, pattern = "-", simplify = T)[1,])])
    # rm.idx = rm.idx
    
    
    # sig_i <- sig$all[sig.cutoff.df[ii]>0,]
    # rm.. <- unique(sample.list$group[-which(sample.list$group %in% stri_split_fixed(ii,pattern = "-",simplify = T)[1,])])
    
    if (length(rm..) > 0){
      
      rs0 <- stri_split_fixed(ii, pattern = "-", simplify = T)[1, ]
      xx..rs <-
        !stri_detect_regex(str = names(sig_i), pattern = "^[:alnum:]{2}[.]{2}") ### pattern is 'fc../cv../pv../qv../fd..'
      rs <-
        gsub(rs0, pattern = "([+]|[~]|[&])", replacement = "\\\\\\1")
      rs.pattern <-
        paste(
          "^[:alnum:]{2}[.]{2}(",
          rs[1],
          "-",
          rs[2],
          "$)",
          "|",
          "^[:alnum:]{2}[.]{2}(",
          rs[2],
          "-",
          rs[1],
          "$)",
          sep = ""
        )
      # names(sig_i)[rs.idx]
      rs.idx.1 <-
        stri_detect_regex(str = names(sig_i), pattern = rs.pattern)
      rs.idx.2 <- xx..rs | rs.idx.1
      cpr.idx.rm <- which(sample.list$group %in% rm..)
      rm.idx.1 <-
        names(sig_i) %in% c(
          sample.list$raw[cpr.idx.rm],
          sample.list$bioID[cpr.idx.rm],
          paste("im..", sample.list$bioID[cpr.idx.rm], sep = ""),
          rm..
        )
      cpr.idx.imp <- which(sample.list$group %in% rs0)
      ### single tag variable names
      rs.idx.imp <-
        names(sig_i) %in% c(paste("im..", sample.list$bioID[cpr.idx.imp], sep = ""))
      rs.idx.cv <-
        names(sig_i) %in% c(paste("cv..", sample.list$group[cpr.idx.imp], sep = ""))
      
      if (!quant.type %in% c("LFQ", "iBAQ", "TMT")) {
        rs.idx.imp[] <- F
      }
      
      rs.idx <- (!rm.idx.1) & (rs.idx.2 | rs.idx.imp | rs.idx.cv)
      sig_rsv <- sig_i[rs.idx]
    } else{
      sig_rsv <- sig_i
    }
  }
  return(sig_rsv)
}

if(nrow(sig[[x]]) > 0){
vp.cpr <- x
vp.sig <- sig[[x]]
if (is.null(get0("vp.text.label", ifnotfound = NULL))) {
  text.label <- vp.sig$Symbol
} else{
  text.label <- vp.text.label
}

vp.df <- select_var(compare.vc = vp.cpr, data = pro_quant_df)
vp.df$Symbol[-which(vp.df$Symbol %in% text.label)] <- NA
{
  fc.idx <- names(vp.df)[grep(x = names(vp.df), pattern = "fc..")]
  pv.idx <- names(vp.df)[grep(x = names(vp.df), pattern = "pv..")]
  color. <- (log2(vp.df[, fc.idx])) * (-log10(vp.df[, pv.idx]))
  
  
  valcano <- data.frame(
    symbol = vp.df$Symbol,
    log.fc = log2(vp.df[, fc.idx]) ,
    color = color.,
    p.value = -log10(vp.df[, pv.idx]),
    label = as.numeric(sig.cutoff.list[[x]])
  )
  if (nrow(vp.sig) > label.num & is.null(get0("vp.text.label", ifnotfound = NULL))) {
    top20 <- subset.data.frame(valcano, label == 1, select = color)
    the20 <- rank(-abs(top20$color)) == label.num
    min.col <- abs(top20$color[the20])
    valcano$symbol[abs(valcano$color) < min.col] <- NA
  }
  
  valcano$symbolL <- valcano$symbol
  valcano$symbolR <- valcano$symbol
  valcano$symbolL[valcano$log.fc>0] <- NA
  valcano$symbolR[valcano$log.fc<0] <- NA
  
  vp.col.list <- list(col.low = "limegreen",
                      col.mid = "grey30",
                      col.high = "red")
  
}



v.label <-
  data.frame(
    x.label = c(-1, ceiling(-max(
      abs(valcano$log.fc[valcano$label > 0])
    )) * 0.9, 1),
    y.label = c(-0.2, -log10(pv_co) * 0.9, -0.2),
    col.label = c(
      vp.col.list$col.low,
      vp.col.list$col.mid,
      vp.col.list$col.high
    ),
    label = paste(
      "italic(",
      c("fc-dn==", "p-value==", "fc-up=="),
      c(2 ^ lr_co, pv_co, 2 ^ lr_co),
      ")",
      sep = ''
    ),
    stringsAsFactors = F
  )

{
  ylab <- bquote(-log[10] ~ (p.value))
  xlab <- bquote(log[2] ~ (ratio))
  legend_title_vp <-  bquote(Significant ~ legend)
  xmin <- min((valcano$color[valcano$label > 0]))
  xmax <- max((valcano$color[valcano$label > 0]))
  valcano$color[which(valcano$color > xmax)] <- xmax
  valcano$color[which(valcano$color < xmin)] <- xmin
  }
attach(vp.col.list,warn.conflicts = F)
valcano.plot <-
  ggplot(valcano,
         aes(
           log.fc,
           p.value,
           fill = color,
           size = p.value,
           alpha = label
         )) +
  ylab(ylab) +  xlab(xlab) + ggtitle(vp.cpr) +
  geom_point(shape = 21,
             stroke = 0.1,
             color = "white"
             ) +
  geom_hline(
    yintercept = -log10(pv_co),
    linetype = 5,
    size = 0.3,
    color = alpha(colour = col.mid, 0.5)
  ) +
  geom_vline(
    xintercept = lr_co,
    linetype = 5,
    size = 0.3,
    color = alpha(colour = col.high, 0.5)
  ) +
  geom_vline(
    xintercept = -lr_co,
    linetype = 5,
    size = 0.3,
    color = alpha(colour = col.low, 0.75)
  ) +
  geom_label(
    nudge_x = c(-1, 0, 1),
    parse = T,
    inherit.aes = F,
    data = v.label ,
    aes(x = x.label, y = y.label, label = label),
    color = "gold4",
    alpha = 0.75,
    size = 2
  ) +
  geom_text_repel(
    box.padding = 0.1,
    xlim = c(NA,-lr_co),
    # ylim = c(-log10(pv_co),NA),
    point.padding = 0.2,
    segment.size = 0.1,
    segment.alpha = 0.2,
    aes(label = symbolL, family = "Times New Roman"),
    color = "darkgreen",
    alpha = 0.8,
    size = 1.7
  ) +
  geom_text_repel(
    box.padding = 0.1,
    xlim = c(lr_co,NA),
    # ylim = c(-log10(pv_co),NA),
    point.padding = 0.2,
    segment.size = 0.1,
    segment.alpha = 0.2,
    aes(label = symbolR, family = "Times New Roman"),
    color = "firebrick4",
    alpha = 0.8,
    size = 1.7
  ) +
  xlim(-max(abs(valcano$log.fc[valcano$label > 0])), max(abs(valcano$log.fc[valcano$label > 0]))) +
  scale_alpha(range = c(0.1, 1)) +
  scale_size(range = c(0.1, 2.5)) +
  scale_fill_gradient2(
    na.value = "black",
    # limits=c(-10,10),
    midpoint = 0,
    mid = NULL,
    low = vp.col.list$col.low,
    high = vp.col.list$col.high,
    breaks = c(-min(abs(xmin),abs(xmax)),0,min(abs(xmin),abs(xmax))),
    label = c("Down-regulated", "Not-significant", "Up-regulated")
  ) +
  theme( 
    aspect.ratio = 1,
    panel.grid.major = element_line(size=0.2),
    panel.grid.minor = element_line(size=0.2),
    text = element_text(face = "italic",family = "Times New Roman"),
    # legend.background = element_rect(),
    # legend.key.size = unit(.8,"cm"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  guides(
    alpha = "none",
    size = "none",
    color = "none",
    fill = guide_legend(title = legend_title_vp)
  )



ggsave(
  device = cairo_pdf,
  paste(dir_2,"v-plot.pdf",sep = ""),
  dpi = dpi.set,
  valcano.plot,
  height = 6
  # scale = 2
)
ggsave(
  filename = paste(dir_2,paste("v-plot", fig.type, sep = "."),sep = ""),
  dpi = dpi.set,
  valcano.plot,
  height = 6
  )
# setwd(dir_1)
}

