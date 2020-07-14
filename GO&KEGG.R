#### GO & KEGG ############

x = get0(x = "go_ppi", ifnotfound = names(sig)[1])
{
  library("clusterProfiler", lib.loc = lib.loc)
  library("STRINGdb", lib.loc = lib.loc)
  # library("pathview", lib.loc = lib.loc)
  library("stringi", lib.loc = lib.loc)
  require(clusterProfiler, lib.loc = lib.loc)
  require(ggplot2, lib.loc = lib.loc)
  require(Cairo, lib.loc = lib.loc)
  require(RColorBrewer, lib.loc = lib.loc)
}
if(custom){
  write.csv(x = pro_quant_df,file = paste(dir_2,"annotated list.csv",sep = ""),row.names = F)
  x <- "custom"
  sig <- list()
  gene_id <- sig[[x]] <- pro_quant_df[grep(pro_quant_df$Symbol,pattern = "^[WXYZ]"),]
}else if (fca)  {
  gene_id <- TC_df
} else{
  gene_id <- sig[[x]]
}
if (length(gene_id$ENTREZID) <= 1) {
  print("NO significant gene")
} else{
  cat(paste("\n\n\n==============",  x, "==============\n\n\n"))
  print(gene_id$ENTREZID)
  
###############################################################################################################
  
  enGO <- function(gene, ...) {
    cc <-
      enrichGO(
        gene = gene,
        OrgDb = OrgDb,
        ont = "CC",
        pvalueCutoff = 0.05,
        minGSSize = 2,
        readable = T,
        pAdjustMethod = "fdr"
      )
    bp <-
      enrichGO(
        gene = gene,
        OrgDb = OrgDb,
        ont = "BP",
        pvalueCutoff = 0.05,
        minGSSize = 2,
        readable = T,
        pAdjustMethod = "fdr"
      )
    mf <-
      enrichGO(
        gene = gene,
        OrgDb = OrgDb,
        ont = "MF",
        pvalueCutoff = 0.05,
        minGSSize = 2,
        readable = T,
        pAdjustMethod = "fdr"
      )
    list(cc = cc, bp = bp, mf = mf)
  }
  listGO <- enGO(gene = gene_id$ENTREZID)   ### CONFIG ###
  
  
  #####################################################################################################
  {
    if (all(is.na.data.frame(x = listGO$cc[1]))) {
      cc <- NULL
    } else{
      cc <- listGO$cc
    }
    
    if (all(is.na.data.frame(x = listGO$bp[1]))) {
      bp <- NULL
    } else{
      bp <- listGO$bp
    }
    
    if (all(is.na.data.frame(x = listGO$mf[1]))) {
      mf <- NULL
    } else{
      mf <- listGO$mf
    }
    
  }
  
  
  num = GO.max     ### < CONFIG > ###
  
  if (!is.null(cc)) {
    ccL <- nchar(cc@result$Description)[1:num]
  } else{
    ccL <- 0
  }
  if (!is.null(bp)) {
    bpL <- nchar(bp@result$Description)[1:num]
  } else{
    bpL <- 0
  }
  if (!is.null(mf)) {
    mfL <- nchar(mf@result$Description)[1:num]
  } else{
    mfL <- 0
  }
  maxL <- max(na.omit(c(ccL, bpL, mfL)))
  
 
  ### wrap string to multi row in axis label
  dsct2 <- function(x, width = GO.dscpt.length , size = 10, reduce = 0.75, ...){
    size <- ifelse(nchar(x,type = "width")>=width,size*reduce,size)
    x_2 <- paste(strwrap(x = x, width = width),collapse = "\n")
    return(c(x_2,size))
  }
  
  if (!is.null(cc)) {
    cc2 <-
      sapply(
        cc@result$Description,
        FUN = dsct2,
        size = 10,
        USE.NAMES = F
      )
    cc@result$Description <- cc2[1, ]
    cc@result$size <- cc2[2, ]
    cc@result$type <- "Cellular component"

  }
  
  
  
  if (!is.null(bp)) {
    bp2 <-
      sapply(
        bp@result$Description,
        FUN = dsct2,
        size = 10,
        USE.NAMES = F
      )
    bp@result$Description <- bp2[1, ]
    bp@result$size <- bp2[2, ]
    bp@result$type <- "Biological process"

  }
  
  if (!is.null(mf)) {
    mf2 <-
      sapply(
        mf@result$Description,
        FUN = dsct2,
        size = 10,
        USE.NAMES = F
      )
    mf@result$Description <- mf2[1, ]
    mf@result$size <- mf2[2, ]
    mf@result$type <- "Molecular function"

  }
  
  
  
  {
    eg_cc <- as.data.frame(cc)
    eg_bp <- as.data.frame(bp)
    eg_mf <- as.data.frame(mf)
    eg_go <- rbind.data.frame(eg_bp, eg_cc, eg_mf)
  }
  
  eg_go.list <- list(
    'Biological process' = as.data.frame(bp),
    'Cellular component' = as.data.frame(cc),
    'Molecular function' = as.data.frame(mf)
  )
  
  ######################################################################################################
  {
    write.csv(
      x = eg_mf[-10],
      file = paste(dir_2,"enrichGO.",  x, ".MF.csv", sep = ""),
      row.names = F
    )
    write.csv(
      x = eg_cc[-10],
      file = paste(dir_2,"enrichGO.",  x, ".CC.csv", sep = ""),
      row.names = F
    )
    write.csv(
      x = eg_bp[-10],
      file = paste(dir_2,"enrichGO.",  x, ".BP.csv", sep = ""),
      row.names = F
    )
    write.csv(
      x = eg_go[-10],
      file = paste(dir_2,"enrichGO.",  x, ".all", ".csv", sep = ""),
      row.names = F
    )
    
  }
  ###################################
  {
    ### < CONFIG > ###
    iterm.N <- function(df, max = num) {
      if (nrow(df) == 0) {
        df.l <- NULL
      } else{
        df.l <- 1:ifelse(max(nrow(df) > max), max, nrow(df))
      }
      return(df.l)
    }
    # iterm.N(df = eg_mf)
    col <- data.frame(
      bar =   brewer.pal(n = 8 , name = "Paired")[1:4 * 2 - 1],
      point = brewer.pal(n = 8 , name = "Paired")[1:4 * 2],
      line =  brewer.pal(n = 8 , name = "Paired")[1:4 * 2]
    )
  }
  
  {
    eg_cc_ <- eg_cc[iterm.N(df = eg_cc), ]
    eg_bp_ <- eg_bp[iterm.N(df = eg_bp), ]
    eg_mf_ <- eg_mf[iterm.N(df = eg_mf), ]
    eg_go_ <- rbind.data.frame(eg_bp_, eg_cc_,  eg_mf_)
    eg_go_$col.bar <- col$bar[as.factor(eg_go_$type)]
    eg_go_$col.point <- col$point[as.factor(eg_go_$type)]
  }
  
  if (nrow(eg_go) > 0) {
    eg_go_$Description <-
      factor(eg_go_$Description, levels = rev(x = eg_go_$Description))
    
    sec_scale <-
      round(max(-log10(eg_go_$p.adjust))) / max(eg_go_$Count)
    
    gg_go <-
      ggplot(eg_go_, aes(x = Description, y = log10(p.adjust)))
    
    {
      fig.go.p... <-  gg_go +
        
        geom_bar(
          color = "white",
          stat = "identity",
          size = 0.25,
          width = 0.75,
          fill = eg_go_$col.point,
          alpha = 0.6
        ) +
        # theme_bw() +
        
        xlab("") +
        geom_path(
          aes(x = c(
            rev(iterm.N(df = eg_bp_)),
            rev(iterm.N(df = eg_cc_)),
            rev(iterm.N(df = eg_mf_))
          ),
          y = -Count * sec_scale),
          size = 0.25,
          stat = "identity",
          linetype = 1,
          color = alpha(eg_go_$col.point, alpha = 0.75)
        ) +
        
        facet_grid(type ~ .,switch = "y",
                   scales = "free_y",
                   space = "free_y" 
                   ) +
        
        geom_point(
          aes(y = -Count * sec_scale),
          stat = "identity",
          color = "white",
          fill = alpha(eg_go_$col.point, alpha = 1),
          shape = 21,
          # stroke = 0.5,
          size = 2
        ) +
        
        coord_flip(expand = T) + 
        scale_x_discrete(position = "top")+
        
        scale_y_continuous(
          sec.axis = sec_axis(
            trans = ~ . / sec_scale,
            name = "Count"
            # breaks = seq(0,ceiling(max(eg_go_$Count)),4)
          ),
          position = "right",
          name = bquote(log[10]~(p.adjust))
          )
      
      
      
      fig.go.p.0 <-  fig.go.p... + 
        theme(
          panel.grid.major = element_line(size = 0.1),
          panel.grid.minor = element_line(size = 0.1),
          axis.title.x = element_text(
            size = 14,
            family = "Times New Roman",
            face = "italic"
          ),
          axis.text.y = element_text(
            family = "Segoe UI",
            lineheight = 0.8,
            size = 11,
            angle = 0,
            hjust = 1,
            vjust = 0.4
          ),
          strip.text.y = element_text(
            family = "Times New Roman",
            size = 14,
            angle = 180
          )
        )
      
      fig.go.p.90 <-  fig.go.p... + 
        theme(
          panel.grid.major = element_line(size = 0.1),
          panel.grid.minor = element_line(size = 0.1),
          axis.title.x = element_text(
            size = 14,
            family = "Times New Roman",
            face = "italic"
          ),
          axis.text.y = element_text(
            family = "Segoe UI",
            lineheight = 0.8,
            size = 11,
            angle = 0,
            hjust = 1,
            vjust = 0.4
          ),
          strip.text.y = element_text(
            size = 14,
            family = "Times New Roman",
            angle = 270
          )
        )
      
      }
    
    {
      fig.go.asp.0 <- (nrow(eg_go_) * 1 + 4)/(max(GO.dscpt.length,50) * 0.1 + 25)
      ggsave(
        device = CairoPDF,
        dpi = dpi.set,
        units = "cm",
        scale = 0.8,
        paste(dir_2,"GO-all.h.pdf",sep = ""),
        plot = fig.go.p.0,
        width = max(GO.dscpt.length,50) * 0.1 + 28,
        height = nrow(eg_go_) * 1 + 6
      )
      ggsave(
        dpi = dpi.set,
        units = "cm",
        scale = 0.8,
        paste(dir_2,paste("GO-all.h",fig.type,sep = "."),sep = ""),
        plot = fig.go.p.0,
        width = max(GO.dscpt.length,50) * 0.1 + 28,
        height = nrow(eg_go_) * 1 + 6
      )
    }
    {
      fig.go.asp.90 <- (nrow(eg_go_) * 1 + 4)/(max(GO.dscpt.length,50) * 0.1 + 20)
      ggsave(
        device = CairoPDF,
        dpi = dpi.set,
        units = "cm",
        scale = 0.8,
        paste(dir_2,"GO-all.v.pdf",sep = ""),
        plot = fig.go.p.90,
        width = max(GO.dscpt.length,50) * 0.1 + 24,
        height = nrow(eg_go_) * 1 + 6
      )
      
      ggsave(
        dpi = dpi.set,
        units = "cm",
        scale = 0.8,
        paste(dir_2,paste("GO-all.v",fig.type,sep = "."),sep = ""),
        plot = fig.go.p.90,
        width = max(GO.dscpt.length,50) * 0.1 + 24,
        height = nrow(eg_go_) * 1 + 6
      )
    }
  } else{
    cat("============ No enriched iterm ================")
  }
  
  
  
  ######################################################################################################
  # log_ratio <- log2(as.vector(gene_id$[, 1]))
  # names(log_ratio) <- gene_id$ENTREZID
  
  ek <-
    enrichKEGG(
      gene = na.omit(gene_id$ENTREZID),
      organism = organism,
      pvalueCutoff = 0.05,
      pAdjustMethod = "fdr",
      minGSSize = 2
    )
  if (!is.null(ek)) {
    if (nrow(ek@result) != 0) {
      ek2 <-
        sapply(
          ek@result$Description,
          FUN = dsct2,
          size = 10,
          USE.NAMES = F
        )
      ek@result$Description <- ek2[1,]
      ek@result$size <- ek2[2,]
      
    }
  }
  
  {
    ek_summary <- as.data.frame(ek)
    write.csv(
      x = ek_summary,
      file = paste(dir_2,"enrichKEGG.",  x, ".csv", sep = ""),
      row.names = F
    )
    }
  
  if (sum(ek_summary$p.adjust < 0.05) & F) {
    ### < CONFIG >  ###
    pv.out <- pathview(
      gene.data = na.omit(gene_id$ENTREZID),
      pathway.id = ek_summary$ID,
      species = organism,
      out.suffix = "enrich_gene",
      kegg.native = T,
      sign.pos = "topright",
      same.layer = T,
      new.signature = T
    )
  }
  
  
  if (nrow(ek_summary) >= 5) {
    sec_scale.ek <-
      round(max(-log10(ek_summary$p.adjust))) / max(ek_summary$Count)
    ek_summary$Description <-
      factor(ek_summary$Description,
             levels = rev(x = ek_summary$Description))
    
    ek_summary_ <- ek_summary[iterm.N(ek_summary), ]
    
    gg_ek <-
      ggplot(ek_summary_, aes(x = Description, y = log10(p.adjust)))
    
    
    
    ek... <- gg_ek + coord_flip() + 
      geom_bar(
        color = "white",
        stat = "identity",
        width = 0.6,
        size = 0.25,
        fill = col$bar[4],
        alpha = 0.8
      ) +
      geom_path(
        aes(x = length(Description):1, y = -Count * sec_scale.ek),
        stat = "identity",
        size = 0.25,
        color = alpha(col$line[4], 0.25)
      ) +
      geom_point(
        aes(x = length(Description):1, y = -Count * sec_scale.ek),
        stat = "identity",
        fill = alpha(col$point[4], alpha = 0.9),
        shape = 21,
        stroke = 0.5,
        size = 2,
        color = "white"
      ) +
      ggtitle("Enrich KEGG") + xlab("") +
      scale_x_discrete(position = "top")+
      scale_y_continuous(
        sec.axis = sec_axis(
          trans = ~ . / sec_scale.ek,
          name = "Count",
          breaks = (1:ceiling(max(
            ek_summary_$Count
          ) / 2)) * 2
        ),
        position = "right",
        name = expression(log[10] ~ (p.adjust))
      )
      
    fig.ek <-  ek... + 
      theme(
        panel.grid.major = element_line(size = 0.1),
        panel.grid.minor = element_line(size = 0.1),
        axis.title.x = element_text(
          size = 14,
          family = "Times New Roman",
          face = "italic"
        ),
        axis.text.y = element_text(
          family = "Arial Narrow",
          lineheight = 0.8,
          size = 11,
          angle = 0,
          hjust = 1,
          vjust = 0.4
        ),
        strip.text.y = element_text(
          size = 14,
          family = "DengXian",
          angle = 0
        )
      )
}

  if (nrow(ek_summary) >= 5) {
    fig.ek.asp <- (nrow(eg_go_) * 0.5 + 2)/(max(GO.dscpt.length,50) * 0.1 + 20)
    
    ggsave(
      dpi = dpi.set,
      units = "cm",
      scale = 0.7,
      paste(dir_2,paste("KEGG",fig.type,sep = "."),sep = ""),
      plot = fig.ek,
      width = 21 + GO.dscpt.length * 0.1,
      height = nrow(eg_go_) * 0.6 + 2
    )
    ggsave(
      device = CairoPDF,
      dpi = dpi.set,
      units = "cm",
      scale = 0.7,
      paste(dir_2,"KEGG.pdf",sep=""),
      plot = fig.ek,
      width = 21 + GO.dscpt.length * 0.1,
      height = nrow(eg_go_) * 0.6 + 2
    )
  }
}
#### GO pie chart ####
if(F){
  require("plotly", lib.loc = lib.loc)
  require("webshot", lib.loc = lib.loc)
  if (!require("webshot"))   {
    install.packages("webshot")
    webshot::install_phantomjs()
  }
  
  plot.ly.eg <- function(eg.df, sub.title = eg.df[, "type"]) {
    if (nrow(eg.df) == 0){ return(NULL)}
    p <-
      plot_ly(
        width = 500,
        height = 500,
        data = eg.df[1:10,],
        labels = ~ Description,
        values = ~ Count,
        type = 'pie',
        # sizes = c(50, 50),
        textposition = 'inside',
        textinfo = 'label+percent',
        insidetextfont = list(color = 'white'),
        # hoverinfo = 'text',
        # text = ~
        # marker = list(colors = colors, line = list(color = '#FFFFFF', width = 1)),
        pull = 0.005,
        showlegend = F
      ) %>% plotly::layout(
        margin = list(
          l = 0,
          r = 0,
          b = 0,
          t = 50,
          pad = 0
        ),
        
        title = paste("Gene ontology (GO)", sub.title),
        xaxis = list(
          showgrid = T,
          zeroline = F,
          showticklabels = F
        ),
        yaxis = list(
          showgrid = T,
          zeroline = F,
          showticklabels = F
        )
      ) #%>% offline(out_dir = getwd())
    return(p)
  }
  
  
  
  
  
  fig.go.pie <- list(
    cc = plot.ly.eg(
      eg.df = eg_go.list$`Cellular component`,
      sub.title = "-- Cellular component"
    ),
    bp = plot.ly.eg(
      eg.df = eg_go.list$`Biological process`,
      sub.title = "-- Biological process"
    ),
    mf = plot.ly.eg(
      eg.df = eg_go.list$`Molecular function`,
      sub.title = "-- Molecular function"
    )
  )
  
  for (f in names(fig.go.pie)) {
    f. <- fig.go.pie[[f]]
    if (!is.null(f.))
      export(f., file = paste(dir_2,"GO pie chart.", f, ".pdf", sep = ""))
  }
}


# setwd(dir_1)
