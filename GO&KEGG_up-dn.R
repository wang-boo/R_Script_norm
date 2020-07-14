#### GO & KEGG ############

{
  library("clusterProfiler", lib.loc = lib.loc)
  library("STRINGdb", lib.loc = lib.loc)
  # library("pathview", lib.loc = lib.loc)
  library("stringi", lib.loc = lib.loc)
  require(clusterProfiler,lib.loc=lib.loc)
  require(ggplot2,lib.loc=lib.loc)
  require(RColorBrewer,lib.loc=lib.loc)
}

########
# x = get0(x = "go_ppi", ifnotfound = 2)
go_ppi = get0(x = "go_ppi", ifnotfound = 2)
x <- go_ppi
vs <- paste(cpr$vs.1[x-1],cpr$vs.2[x-1],sep = "-")
if (!dir.exists(paste(dir_1, vs, sep = "/"))) {
  dir.create(paste(dir_1,vs, sep = "/"))
}
if (!dir.exists(paste(dir_1,vs, "up",sep = "/"))) {
  dir.create(paste(dir_1,vs, "up",sep = "/"))
  dir.create(paste(dir_1,vs, "dn",sep = "/"))
}


setwd(paste(dir_1, vs, sep = "/"))
sub_sig <- sig[[x]]
gene_id0 <- sub_sig

if(up.dn) {
  gene_id_up <-
    subset.data.frame(x = gene_id0, subset = gene_id0[paste("avg.ratio", names(sig)[x], sep = ".")] > 1)
  gene_id_dn <-
    subset.data.frame(x = gene_id0, subset = gene_id0[paste("avg.ratio", names(sig)[x], sep = ".")] < 1)
  regulation <- list(up = gene_id_up, dn = gene_id_dn)
}
for(r in names(regulation)) {
  gene_id <- regulation[[r]]
  setwd(paste(dir_1, vs, r, sep = "/"))
  
  if (length(gene_id$ENTREZID) <= 1) {
    print("NO significant gene")
    
  } else{
    cat(
      paste(
        "\n\n\n==============",
        names(sig)[go_ppi],
        "==============\n\n\n",
        "\n\n\n==============",
        r,
        "running ==============\n\n\n"
      )
    )
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
    
    dsct <- function (x,
                      CUT = 60,
                      size = 13,
                      reduce = 0.75) {
      ### < CONFIG > ###
      x <- stri_trim_both(str = x)
      if (maxL > CUT & nchar(x, type = "width") > CUT) {
        sps <- gregexpr(
          pattern = " ",
          text = paste(x, " ", sep = ""),
          fixed = F
        )[[1]]
        if (length(which(sps - CUT >= 0)) == 1) {
          x_2 <- x
          size <- size
        } else{
          index <- sps[which(sps - CUT >= 0)[1]]
          x_2 <-
            paste(substr(x, 1, index - 1),
                  "\n",
                  substr(x, index + 1, nchar(x)),
                  sep = "")
          size <- size * reduce
        }
      } else {
        x_2 <- x
        size <- size
      }
      c(x_2, size)
    }
    
    if (!is.null(cc)) {
      cc2 <-
        sapply(
          cc@result$Description,
          FUN = dsct,
          size = 10,
          USE.NAMES = F
        )
      cc@result$Description <- cc2[1, ]
      cc@result$size <- cc2[2, ]
      cc@result$type <- "CC"
    }
    
    if (!is.null(bp)) {
      bp2 <-
        sapply(
          bp@result$Description,
          FUN = dsct,
          size = 10,
          USE.NAMES = F
        )
      bp@result$Description <- bp2[1, ]
      bp@result$size <- bp2[2, ]
      bp@result$type <- "BP"
    }
    
    if (!is.null(mf)) {
      mf2 <-
        sapply(
          mf@result$Description,
          FUN = dsct,
          size = 10,
          USE.NAMES = F
        )
      mf@result$Description <- mf2[1, ]
      mf@result$size <- mf2[2, ]
      mf@result$type <- "MF"
    }
    
    
    
    {
      eg_cc <- as.data.frame(cc)
      eg_bp <- as.data.frame(bp)
      eg_mf <- as.data.frame(mf)
      eg_go <- rbind.data.frame(eg_bp, eg_cc, eg_mf)
    }
    ######################################################################################################
    {
      write.csv(
        x = eg_mf[-10],
        file = paste("enrichGO.",  names(sig)[x], ".MF.csv", sep = ""),
        row.names = F
      )
      write.csv(
        x = eg_cc[-10],
        file = paste("enrichGO.",  names(sig)[x], ".CC.csv", sep = ""),
        row.names = F
      )
      write.csv(
        x = eg_bp[-10],
        file = paste("enrichGO.",  names(sig)[x], ".BP.csv", sep = ""),
        row.names = F
      )
      write.csv(
        x = eg_go[-10],
        file = paste("enrichGO.",  names(sig)[x], ".all", ".csv", sep = ""),
        row.names = F
      )
      
    }
    ###################################
    {
      ### < CONFIG > ###
      
      col <- data.frame(
        bar =   brewer.pal(n = 8 , name = "Paired")[1:4 * 2 - 1],
        point = brewer.pal(n = 8 , name = "Paired")[1:4 * 2],
        line =  brewer.pal(n = 8 , name = "Paired")[1:4 * 2]
      )
    }
    
    {
      eg_cc_ <- eg_cc[1:ifelse(nrow(eg_cc) > num, num, nrow(eg_cc)), ]
      eg_bp_ <- eg_bp[1:ifelse(nrow(eg_bp) > num, num, nrow(eg_bp)), ]
      eg_mf_ <- eg_mf[1:ifelse(nrow(eg_mf) > num, num, nrow(eg_mf)), ]
      eg_go_ <- rbind.data.frame(eg_bp_, eg_cc_,  eg_mf_)
      eg_go_$col.bar <- col$bar[as.factor(eg_go_$type)]
      eg_go_$col.point <- col$point[as.factor(eg_go_$type)]
      
      eg_go_$Description <-
        factor(eg_go_$Description, levels = rev(x = eg_go_$Description))
    }
    if (nrow(eg_go) > 0) {
      gg_go <-
        ggplot(eg_go_, aes(x = Description, y = -log10(p.adjust)))
      p <- gg_go +
        geom_bar(
          color = "white",
          #alpha(eg_go_$col.point,alpha = 0.75),
          stat = "identity",
          size = 0.25,
          width = 0.6,
          fill = eg_go_$col.point,
          alpha = 0.6
        ) +
        # theme_bw() +
        theme(axis.text.y = element_text(
          # size = as.numeric(eg_go_$size) + 1,
          size = 11,
          angle = 0,
          hjust = 1,
          vjust = 0.3
        )) + xlab("") +
        facet_grid(type ~ ., scales = "free_y", space = "free_y") +
        geom_point(
          aes(y = Count),
          stat = "identity",
          color = "white",
          fill = alpha(eg_go_$col.point, alpha = 1),
          shape = 21,
          stroke = 0.25,
          size = 2
        ) + coord_flip(expand = T) +
        scale_y_continuous(
          sec.axis = sec_axis(trans = ~ ., name = "Count"),
          position = "right",
          name = expression(-log[10] ~ (p.adjust))
        )
      
      
      ggsave(
        dpi = 600,
        units = "cm",
        scale = 0.75,
        "GO-all.pdf",
        plot = p,
        width = 25 + maxL * 0.1,
        height = nrow(eg_go_) * 1.2 + 2
      )
      ggsave(
        dpi = 600,
        units = "cm",
        scale = 0.75,
        "GO-all.png",
        plot = p,
        width = 25 + maxL * 0.1,
        height = nrow(eg_go_) * 1.2 + 2
      )
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
    
    
    {
      ek_summary <- as.data.frame(ek)
      write.csv(
        x = ek_summary,
        file = paste("enrichKEGG.",  names(sig[x]), ".csv", sep = ""),
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
      ek_summary_ <-
        ek_summary[1:ifelse(nrow(ek_summary) > num, num, nrow(ek_summary)), ]
      
      gg_ek <-
        ggplot(ek_summary_, aes(x = reorder(Description, X = -p.adjust), -log10(p.adjust)))
      
      
      
      e <- gg_ek + coord_flip() +
        geom_bar(
          color = "white",
          stat = "identity",
          width = 0.6,
          size = 0.25,
          fill = col$bar[4],
          alpha = 0.8
        ) +
        # theme_bw() +
        theme(axis.text.y = element_text(
          size = 11,
          angle = 0,
          hjust = 1,
          vjust = 0.3
        )) +
        geom_path(
          aes(x = length(Description):1, y = Count),
          stat = "identity",
          color = alpha(col$line[4], 0.25)
        ) +
        geom_point(
          aes(x = length(Description):1, y = Count),
          stat = "identity",
          fill = alpha(col$point[4], alpha = 0.9),
          shape = 21,
          stroke = 0.5,
          size = 2,
          color = "white"
        ) +
        ggtitle("Enrich KEGG") + xlab("") +
        scale_y_continuous(
          sec.axis = sec_axis(trans = ~ ., name = "Count"),
          position = "right",
          name = expression(-log[10] ~ (p.adjust))
        )
      
      
    }
    if (nrow(ek_summary) >= 5) {
      ggsave(
        dpi = 600,
        units = "cm",
        scale = 0.7,
        "KEGG.png",
        plot = e,
        width = 20 + maxL * 0.1,
        height = nrow(eg_go_) * 0.5 + 2
      )
      ggsave(
        dpi = 600,
        units = "cm",
        scale = 0.7,
        "KEGG.pdf",
        plot = e,
        width = 20 + maxL * 0.1,
        height = nrow(eg_go_) * 0.5 + 2
      )
    }
  }
}


setwd(dir_1)
