{
  to.win.dir <- function(dir) {
    return(stringi::stri_replace_all(dir, fixed = "\\", replacement = "/"))
  }
  if (exists("config.list"))
    attach(config.list, warn.conflicts = F)
  
  dir_0 <-
    to.win.dir(paste(choose.dir(default = "Z:/"), "/combined/txt/", sep = ""))
  # config.dir <- paste(dir_0,"config.list",sep = "")
  # setwd(dir_0)
  if (!exists("config.list")) {
    config.file <-
      list.files(path = dir_0,
                 recursive = T,
                 full.names = T)[grep(
                   list.files(
                     path = dir_0,
                     recursive = T,
                     full.names = T
                   ),
                   pattern = "config.list",
                   fixed = T
                 )]
    if (length(config.file) != 0) {
      message("file < config.list > exists, choose file ?\n\n")
      for (f in 1:length(config.file)) {
        message(paste("  [", f, "] ", config.file[f], "\n\n", sep = ""))
      }
      n <- as.numeric(readline(prompt = "  choose = "))
      if (!is.na(n)) {
        load(file = config.file[n], verbose = T)
      } else{
        message("\nOR source new configuration\n")
      }
    } else{
      message("\nOR source new configuration\n")
    }
  }
}
if (exists("config.list"))
  attach(config.list, warn.conflicts = F)
{
  species_info <- list(
    HUMAN = list(
      OrgDb = "org.Hs.eg.db",
      organism = "hsa",
      species = 9606
    ),
    MOUSE = list(
      OrgDb = "org.Mm.eg.db",
      organism = "mmu",
      species = 10090
    ),
    RAT  = list(
      OrgDb = "org.Rn.eg.db",
      organism = "rno",
      species = 10116
    ),
    DOG  = list(
      OrgDb = "org.Cf.eg.db",
      organism = "cfa",
      species = 9615
    )
    
  )
  OrgDb = species_info[[SPECIES_INFO]][[1]]
  organism = species_info[[SPECIES_INFO]][[2]]
  species = species_info[[SPECIES_INFO]][[3]]
  fca <- quant.type == "count"
  dir_1 <- paste(dir_0, run, "/", sep = "")
  if (!dir.exists(dir_1))
  {
    dir.create(dir_1)
  }
  }
{
  lib.loc <- "C:/Program Files/R/R-3.4.3/library/"
  library("Biobase", lib.loc = lib.loc)
  library("reshape2", lib.loc = lib.loc)
  library("ggplot2", lib.loc = lib.loc)
  library("gplots", lib.loc = lib.loc)
  library("stringi", lib.loc = lib.loc)
  library("RColorBrewer", lib.loc = lib.loc)
  library("extrafont", lib.loc  =  lib.loc)
  suppressMessages(loadfonts(device = "win"))
  library("siggenes", lib.loc  =  lib.loc)
  library("dplyr", lib.loc  =  lib.loc)
  
}
{
  # setwd(dir_0 )
  protein <-
    read.table(
      blank.lines.skip = T,
      file = paste(dir_0, "proteinGroups.txt", sep  =  ""),
      header = T,
      sep = "\t",
      quote = "",
      stringsAsFactors = F
    )
  major.n <-
    stri_split_fixed(protein$Majority.protein.IDs,
                     pattern = ";",
                     simplify = T)
  major.500 <- major.n[, 1:min(ncol(major.n), 500)]
  protein$Majority.protein.IDs <-
    apply(major.500, 1, function(x)
      paste(x[x  !=  ""], collapse = ";"))
  
  peptide <-
    read.table(
      blank.lines.skip = T,
      file = paste(dir_0, "peptides.txt", sep  =  ""),
      header = T,
      sep = "\t",
      quote = "",
      stringsAsFactors = F
    )
  # names(protein)
}
{
  pt <-
    c("Reporter.intensity.corrected.",
      "iBAQ.",
      "LFQ.intensity.",
      "MS.MS.count.")
  names(pt) <- c("TMT", "iBAQ", "LFQ", "count")
  reporter <-
    names(protein)[grep(pattern = pt[quant.type], x = names(x = protein))]
}
#### BIO GROUP ############################################################################################

# setwd(dir_1)
sample.list.dir <- paste(dir_1, "sample.list", sep = "")
if (file.exists(sample.list.dir))
{
  message("\n< sample.list exsits >, loading it? [Y] / [N]\n")
  load.file <- readline(prompt = ">>\n")
  
  if (load.file %in% c("y", "Y", ""))
  {
    load(file = sample.list.dir, verbose = T)
    Group.all <- sample.list$bioID
    reporter <- sample.list$raw
  }
} else{
  load.file <- "N"
}

if (!file.exists(sample.list.dir) | load.file %in% c("N", "n"))
{
  check <- T
  while (check)
  {
    samples.info <- data.frame(raw = reporter, stringsAsFactors = F)
    message("======= sample ID structure =======\n")
    print(samples.info)
    print(data.frame(up = config.list$up, down = config.list$dn))
    message("===================================\n")
    message(
      "\n please input bioGroup pattern regular expression ### '[:alnum:]+[.]{1}' for 'BioGroup.' ###\n"
    )
    pattern.read <-
      readline(prompt = "new regular expression pattern = \ndefault pattern = [:alnum:]+")
    if (pattern.read == "")
    {
      group.pattern <- "[:alnum:]+"
      message(paste("\n=====", group.pattern, "=====\n", sep = " "))
    } else
    {
      group.pattern <- pattern.read
      message(paste("\n=====", group.pattern, "=====\n", sep = " "))
    }
    
    # message(
    #   "\n please input bioGroup pattern regular expression position ### '[Enter] for prefix, [0] for postfix' ###\n"
    # )
    # position.read <-  readline(prompt = "pattern position = \ndefault position = prefix")
    if (T) {
      group.mode <- "first" ### first/last
      # message("\n ===== prefix =====")
    } else{
      group.mode <- "last"
      message("\n ===== postfix =====")
    }
    
    sample.number <- length(reporter)
    message(
      "\nplease input [y] to assign new bioID for each reporter, OR press [Enter] directly if use raw ID"
    )
    if (readline() == "")
    {
      samples.info$bioID <- reporter
      samples.info$group <-
        stri_extract(str = samples.info$bioID,
                     regex = group.pattern,
                     mode = group.mode)
      
      message("========== use raw bioID ==========\n")
      print(samples.info)
      message("===================================\n")
    } else
    {
      if (is.null(samples.info$bioID))
      {
        samples.info$bioID <- ""
      }
      
      info.file <- paste(dir_1, ".info.file", sep  =  "")
      if (!file.exists(info.file))
      {
        write.table(
          x = samples.info,
          file = info.file,
          quote = F,
          sep = "\t",
          row.names = F
        )
      }
      shell(cmd = paste('start excel ', info.file, sep  =  ""))
      .readline <- readline(prompt = "Enter to continue ...")
      samples.info <-
        read.table(
          info.file,
          sep = "\t" ,
          header = T,
          stringsAsFactors = F
        )
      
      
      # samples.info <- na.omit(edit(name = samples.info))
      samples.info.na <-
        stri_detect(samples.info$bioID, regex = "[:graph:]")
      # samples.info.na <- !stri_detect(samples.info$bioID,regex = "[:space:]")
      reporter <- reporter[samples.info.na]
      
      samples.info.1 <- samples.info[samples.info.na,]
      
      samples.info.1$group <-
        stri_extract(str = samples.info.1$bioID,
                     regex = group.pattern,
                     mode = group.mode)
      message("========== use new bioID ==========\n")
      print(samples.info.1)
      message("===================================\n")
    }
    
    message(
      "\ncheck bioID assignment, press [Enter] if done, else press [N] to reassign bioID\n"
    )
    check <- ifelse(readline() == "", F, T)
    if (!check)
    {
      message("======== assign bioID done ========\n")
    }
  }
  sample.list <- as.list(samples.info.1)
  sample.list[["group.mode"]] <- group.mode
  sample.list[["group.pattern"]] <- group.pattern
  Group.all <- sample.list$bioID
}

{
  sample.list[["group.names.order"]] <- unique(c(up, dn))
  sample.list[["replicates"]] <- table(sample.list$group)
}
#### configure pv.mt structure ####
group.sep="-"
if (length(sample.list$group.names.order) >= 2)
{
  up.vs.dn <- paste(up, dn, sep = group.sep) ### compare selected
  sample.list[["up.vs.dn"]] <-  up.vs.dn
  cpr <-
    cbind.data.frame(
      vs = stri_split_fixed(
        str = up.vs.dn,
        pattern = group.sep,
        simplify = T
      ),
      stringsAsFactors = F
    )
  cpr$fc <- paste(cpr$vs.1, group.sep, cpr$vs.2, sep = "")
  cpr$imputation <-
    paste("im..", cpr$vs.1, group.sep, cpr$vs.2, sep = "")
  cpr$names <- sample.list$up.vs.dn
  
  cpr$pv <- paste("pv..", sample.list$up.vs.dn, sep = "")
  cpr$class <- paste("class..", cpr$fc, sep = "")
  # cpr$ratio <-
  sample.list[["contrast"]] <- cpr
}


################################################################################################
################################################################################################

{
  # names(protein)
  if (all(c("Gene.names") %in% names(protein)))
  {
    column <-
      c(
        "Majority.protein.IDs",
        "Gene.names",
        # "Protein.names",
        "Number.of.proteins",
        "Unique.peptides",
        "Sequence.coverage....",
        "MS.MS.count",
        "Q.value",
        "Score",
        reporter,
        "iBAQ"
      )
    prot <- protein[, column]
    nrow(prot)
    gene.names.na.bl <- prot$Gene.names == "" |
      !isUnique(
        stri_split_fixed(
          str = prot$Gene.names,
          pattern = ";",
          n = 2,
          simplify = T
        )[, 1])
    
    prot$Symbol <-
      stri_split_fixed(
        str = prot$Gene.names,
        pattern = ";",
        n = 2,
        simplify = T
      )[, 1]
    prot$Symbol[gene.names.na.bl] <-
      stri_split_fixed(
        str = prot$Majority.protein.IDs,
        pattern = ";",
        n = 2,
        simplify = T
      )[, 1][gene.names.na.bl]
    
    prot$Rank <-  rank(x = -prot$iBAQ)
    
}else
  {
    column <-
      c(
        "Majority.protein.IDs",
        # "Gene.names",
        # "Protein.names",
        "Number.of.proteins",
        "Unique.peptides",
        "Sequence.coverage....",
        "MS.MS.count",
        "Q.value",
        "Score",
        reporter,
        "iBAQ"
      )
    prot <- protein[, column]
    prot$Symbol <-
      stri_extract_first(prot$Majority.protein.IDs, regex = "[:alnum:]+[-]{0,1}")
    prot$Rank <-  rank(x = -prot$iBAQ)
    
}
  
  if (!fca)
  {
    if (length(grep("CON_", x = prot$Majority.protein.IDs)) > 0)
    {
      prot <- prot[-grep("CON_", x = prot$Majority.protein.IDs), ]
    } else
    {
      prot <- prot
    }
  }
  
  if (length(grep("REV_", x = prot$Majority.protein.IDs)) > 0)
  {
    prot <- prot[-grep("REV_", x = prot$Majority.protein.IDs), ]
  } else
  {
    prot <- prot
  }
    
  prot_bak <- prot
  if (quant.type %in% c("count")) {
    prot <-  prot[rowSums(prot[reporter] > 0) >= 1,] ## 20180116 hq
  } else{
    prot <-
      prot[rowSums(prot[reporter] > 0) >= length(reporter) * missing.ratio,] ## 20171228 hq
  }
  prot[Group.all] <- prot[c(reporter)]
  # names(prot)
  
  {
    if (quant.type %in% c("TMT", "iBAQ")) {
      ###  CAL FACTORS: sum median custom ####
      prot_rpt_sub <- prot[reporter]
      median_factor <-
        apply(X = prot_rpt_sub,
              MARGIN = 2,
              FUN = median) /
        (apply(
          X = prot_rpt_sub,
          MARGIN = 2,
          FUN = median
        )[1])
        
      nor_mt <-
        t(t(prot_rpt_sub) / median_factor)          ### < CONFIG > ###
      colMeans(nor_mt, na.rm = T)
      
      prot[sample.list$bioID] <- nor_mt
    }
  }

nrow(prot)
#### get ratio ####
  

#### imputation ####
  if (quant.type %in% c("LFQ", "iBAQ", "TMT"))
  {
    prot_imp <- log2(prot[Group.all])
    prot_imp.bl <-
      as.data.frame(prot[Group.all] == 0, optional = T)
    prot_imp.df <- as.data.frame(melt(prot_imp), optional = T)
    
    prot_imp.df$value[is.infinite(prot_imp.df$value)] <- NA
    prot_imp.df$imputation <- F
    prot_imp.df$imputation[is.na(prot_imp.df$value)] <- T
    
    if (imp.mode == "s")
    {
      imp.avg <-
        acast(
          prot_imp.df,
          formula = . ~ variable,
          fun.aggregate = mean,
          na.rm = T
        )
      imp.sd <-
        acast(
          prot_imp.df,
          formula = . ~ variable,
          fun.aggregate = sd,
          na.rm = T
        )
      
      m <- imp.avg - imp.downshift * imp.sd
      s <- imp.width * imp.sd
      
      for (im in sample.list$bioID)
      {
        prot_imp[prot_imp.bl[, im], im] <-
          rnorm(n = sum(prot_imp.bl[, im]),
                sd = s[, im],
                mean = m[, im])
      }
      
    } else
      if (imp.mode == "a")
      {
        imp.avg <- mean(prot_imp.df$value, na.rm = T)
        imp.sd  <-   sd(prot_imp.df$value, na.rm = T)
        
        m <- imp.avg - imp.downshift * imp.sd
        s <- imp.width * imp.sd
        
        for (im in sample.list$bioID)
        {
          prot_imp[prot_imp.bl[, im], im] <-
            rnorm(n = sum(prot_imp.bl[, im]),
                  sd = s,
                  mean = m)
        }
      }
    row.names(prot_imp.bl) <- prot$Majority.protein.IDs
    prot_imp.bl$Majority.protein.IDs <- row.names(prot_imp.bl)
    
    prot_in <- prot[c("Majority.protein.IDs", sample.list$bioID)]
    prot_in[Group.all] <- 2 ^ (prot_imp[Group.all])
    
  } else
    if (quant.type %in% c("count"))
    {
      {
        prot_imp <- log2(prot[Group.all])
        prot_imp.bl <-
          as.data.frame(prot[Group.all] == 0, optional = T)
        prot_imp.df <- as.data.frame(melt(prot_imp), optional = T)
        prot_imp.df$imputation <- F
        row.names(prot_imp.bl) <- prot$Majority.protein.IDs
        prot_imp.bl$Majority.protein.IDs <- row.names(prot_imp.bl)
        
      }
      prot_in <- prot[c("Majority.protein.IDs", sample.list$bioID)]
      names(prot_in)
    }
  rownames(prot) <- prot$Majority.protein.IDs
  rownames(prot_in) <- prot_in$Majority.protein.IDs
}
nrow(prot_in)

nrow(prot_imp.bl)
#### FC / imputation AND/OR CV ####

{
  {
    pro_imp_m <- melt(
      data = prot_imp.bl,
      id.vars = "Majority.protein.IDs",
      variable.name = "replicates",
      as.is = T
    )
    
    pro_imp_m$group <-
      stri_extract(
        str = pro_imp_m$replicates,
        regex = sample.list$group.pattern,
        simplify = T,
        mode = sample.list$group.mode
      )
    pro_imp_m$group <-
      factor(pro_imp_m$group,
             levels = sample.list$group.names.order,
             ordered = T)
    
    imp_mt <-
      acast(
        data = pro_imp_m,
        formula = Majority.protein.IDs ~ group,
        fun.aggregate = mean
      )
    imp_df <- as.data.frame(x = imp_mt, optional = T)
  }
  
  {
    pro_m_ <- melt(data = prot_in,
                   variable.name = "replicates",
                   as.is = T)
    pro_m <- pro_m_
    
    pro_m$group <-
      stri_extract(
        str = pro_m$replicates,
        regex = sample.list$group.pattern,
        simplify = T,
        mode = sample.list$group.mode
      )
    pro_m$group <-
      factor(pro_m$group,
             levels = sample.list$group.names.order,
             ordered = T)
    
    avg_mt <-
      acast(
        data = pro_m,
        formula = Majority.protein.IDs ~ group,
        fun.aggregate = mean
      )
    avg_df <- as.data.frame(x = avg_mt, optional = T)
    # names(avg_df)
    test.num <- length(unique(pro_m$group))
  }
  
  {
    if (all(sample.list$replicates > 2))
    {
      cv_mt <-
        acast(
          data = na.omit(pro_m),
          formula = Majority.protein.IDs ~ group,
          fun.aggregate = function(x)
            100 * sd(x) / mean(x)
        )
      cv_df <- as.data.frame(cv_mt, optional = T)
      
      sd_mt <-
        acast(
          data = na.omit(pro_m),
          formula = Majority.protein.IDs ~ group,
          fun.aggregate = sd
        )
      sd_df <- as.data.frame(sd_mt, optional = T)
    }
  }
}

#### T.TEST OR ANOVA #########################################################
if (all(sample.list$replicates > 2))
{
  pv_list <- list()
  qv_list <- list()
  fdr_list <- list()
  sam.out_list <- list()
  
  fdr.cl.all <-
    factor(
      stri_extract(
        str = names(prot_in)[-1],
        regex = sample.list$group.pattern,
        simplify = T,
        mode = sample.list$group.mode
      ),
      levels = sample.list$group.names.order,
      ordered = T
    )
  # test.num <- length(unique(pro_m$group))
  ## tow case SAM
  if (test.num == 2 | use.t == use.t)
  {
    for (cpr.name in cpr$names)
    {
      fdr.cl.idx.x <- fdr.cl.all %in% cpr[match(cpr.name, cpr$names), 1:2]
      fdr.cl.x <- fdr.cl.all[fdr.cl.idx.x]
      fdr.dt.x <- prot_in[-1][fdr.cl.idx.x]
      sam.out_list[[cpr.name]] <-
        sam(
          data = scale(log(fdr.dt.x)),
          var.equal =  F,
          cl = fdr.cl.x,
          rand = 123
        )
      pv_list[[cpr.name]] <- sam.out_list[[cpr.name]]@p.value
      qv_list[[cpr.name]] <- sam.out_list[[cpr.name]]@q.value
      if (sum(summary(sam.out_list[[cpr.name]])@mat.fdr[, "Called"])  ==  0) {
        fdr_co_idx <- rep(T, nrow(fdr.dt.x))
      } else{
        delta.fdr.0 <- findDelta(sam.out_list[[cpr.name]], fdr = fdr_co)
        delta.fdr <-
          rbind(delta.fdr.0, delta.fdr.0) ### bug fixed for if delta.fdr.0 has only one row
        delta_min <- min(as.data.frame(delta.fdr)$Delta)
        row.sig_vt <-
          summary(sam.out_list[[cpr.name]], delta_min)@row.sig.genes
        fdr_co_idx <- rep(F, nrow(prot_in))
        fdr_co_idx[row.sig_vt] <- T
      }
      fdr_list[[cpr.name]] <- fdr_co_idx
    }
    pv_df <- as.data.frame(pv_list, optional = T)
    rownames(pv_df) <- prot_in[, "Majority.protein.IDs"]
    qv_df <- as.data.frame(qv_list, optional = T)
    rownames(qv_df) <- prot_in[, "Majority.protein.IDs"]
    fdr_df <- as.data.frame(fdr_list, optional = T)
    rownames(fdr_df) <- prot_in[, "Majority.protein.IDs"]
  } else  if (test.num > 2) {
    ## three or more case SAM & anova-TukeyHSD
    cpr.name <- paste(levels(x = fdr.cl.all), collapse = "-")
    
    fdr.dt.all <- prot_in[-1][names(prot_in)[-1]]
    sam.out_list[[cpr.name]] <-
      sam(data = log(fdr.dt.all),
          cl = fdr.cl.all,
          rand = 123)
    plot(sam.out_list[[cpr.name]], 5)
    pv_list[[cpr.name]] <- sam.out_list[[cpr.name]]@p.value
    # pv_df0 <- as.data.frame(pv_list, optional = T)
    # rownames(pv_df0) <- prot_in["Majority.protein.IDs"]
    qv_list[[cpr.name]] <- sam.out_list[[cpr.name]]@q.value
    
    if (sum(summary(sam.out_list[[cpr.name]])@mat.fdr[, "Called"]) == 0) {
      fdr_co_idx <- rep(T, nrow(fdr.dt.all))
    } else{
      delta.fdr <- findDelta(sam.out_list[[cpr.name]], fdr = fdr_co)
      delta_min <- min(as.data.frame(delta.fdr)$Delta)
      row.sig_vt <-
        summary(sam.out_list[[cpr.name]], delta_min)@row.sig.genes
      fdr_co_idx <- rep(F, nrow(prot_in))
      fdr_co_idx[row.sig_vt] <- T
    }
    fdr_list[[cpr.name]] <- fdr_co_idx
    
    qv_df <- as.data.frame(qv_list, optional = T)
    rownames(qv_df) <- prot_in[, "Majority.protein.IDs"]
    fdr_df <- as.data.frame(fdr_list, optional = T)
    rownames(fdr_df) <- prot_in[, "Majority.protein.IDs"]
    
    anova.df <- function(data.df, cl.vt, ...) {
      Tukey.pv <- as.data.frame(optional = T, t(apply(
        X = data.df,
        MARGIN = 1,
        FUN = function(.x = X,
                       .cl.vt = cl.vt,
                       ...) {
          aov.lm <-
            aov(
              formula = value ~ variable,
              data = data.frame(value = .x, variable = .cl.vt)
            )
          TukeyHSD(aov.lm)$variable[, "p adj"]
        }
      )))
      return(Tukey.pv)
    }
    pv_df <- anova.df(fdr.dt.all, fdr.cl.all)
    # anova.df(data.df = fdr.dt.all[1:14, ], factor(c("a", "b", "c", "a", "b", "c")))
    
  }
}

#### average fc ####
fc_list <- list()
im_list <- list()

if (fca) {
  FCa <- function(idx = 1, df = df) {
    #### Nj and Nx with coexist protiens ####
    coexist <-
      rowSums(df[c(cpr$vs.1[idx], cpr$vs.2[idx])] >= fca.spc) == 2
    
    Nj <- sum(df[coexist, cpr$vs.1[idx]])
    Nx <- sum(df[coexist, cpr$vs.2[idx]])
    Cix <- df[cpr$vs.2[idx]] / Nx
    n <- 1
    b <- 1
    a <- b / ave(Nx)
    Ci <- (1 / n) * Cix
    Tij <- df[cpr$vs.1[idx]] / Nj
    FCij <-
      ((Tij + a) / (Ci + a)) * (rowSums((df[c(cpr$vs.1[idx], cpr$vs.2[idx])])) > 0)
    
    return(log2(FCij))
  }
  
  for (j in 1:nrow(cpr)) {
    fc_list[cpr$fc[j]] <-  FCa(idx = j, df = avg_df)
  }
  # names(fc_list)
  fc_df <-
    as.data.frame(x = fc_list,
                  optional = T,
                  row.names = row.names(avg_df))
  # names(fc_df)
  class_df <- fc_df
  prot.o <- prot[rownames(fc_df),]
  for (.f in 1:nrow(cpr)) {
    stopifnot(all.equal.character(rownames(fc_df), rownames(prot.o), rowname.vt))
    
    f <- cpr$fc[.f]
    vs.1 <- cpr$vs.1[.f]
    vs.2 <- cpr$vs.2[.f]
    
    idxT   <-
      which(
        fc_df[, f] >= lr_co  &
          prot.o$Unique.peptides >= unique_co &
          (avg_df[, vs.1] > spc_co | avg_df[, vs.2] > spc_co)
      )
    idxC   <-
      which((fc_df[, f]) <= -lr_co &
              prot.o$Unique.peptides >= unique_co &
              (avg_df[, vs.1] > spc_co | avg_df[, vs.2] > spc_co)
      )
    idxNA  <- which((fc_df[, f]) == -Inf)
    
    class_df[, f] <- "T&C"
    class_df[idxT , f] <- "T"
    class_df[idxC , f] <- "C"
    class_df[idxNA, f] <-  NA
    
  }
  # names(class_df)
} else{
  for (j in 1:nrow(cpr)) {
    fc_list[[cpr$fc[j]]] <-
      avg_df[, cpr$vs.1[j]] / avg_df[, cpr$vs.2[j]]
    im_list[[cpr$im[j]]] <-
      imp_df[, cpr$vs.1[j]] < imp_co |
      imp_df[, cpr$vs.2[j]] < imp_co
  }
  # names(fc_list)
  fc_df <-
    as.data.frame(x = fc_list,
                  optional = T,
                  row.names = row.names(avg_df))
  im_df.rs <-
    as.data.frame(x = im_list,
                  optional = T,
                  row.names = row.names(imp_df))
}

####### merge data #######
rowname.vt <- rownames(avg_df)
prot.od <- prot[rowname.vt,]
prot_in.od <- prot_in[rowname.vt,]

prot_imp.od <- prot_imp.bl[rowname.vt, -length(prot_imp.bl)]
names(prot_imp.od) <- paste("im.", names(prot_imp), sep = ".")

if (all(sample.list$replicates > 2)) {
  # pv_df.od <- pv_df[rowname.vt, ]
  fc_df.od <- fc_df
  cv_df.od <- cv_df
  names(fc_df.od) <- paste("fc.", names(fc_df), sep = ".")
  names(cv_df.od) <- paste("cv.", names(cv_df), sep = ".")
  
  pv_df.. <- pv_df
  qv_df.. <- qv_df
  fdr_df.. <- fdr_df
  names(pv_df..) <- paste("pv.", names(pv_df), sep = ".")
  names(qv_df..) <- paste("qv.", names(qv_df), sep = ".")
  names(fdr_df..) <- paste("fd.", names(fdr_df), sep = ".")
  pv_qv_fdr.od <- cbind.data.frame(pv_df.., qv_df..)[rowname.vt, ]
  
  
  if (all(
    all.equal.character(rowname.vt, rownames(cv_df.od)),
    all.equal.character(rowname.vt, rownames(fc_df.od)),
    # all.equal.character(rowname.vt, rownames(pv_df.od)),
    all.equal.character(rowname.vt, rownames(pv_qv_fdr.od)),
    all.equal.character(rowname.vt, rownames(prot.od)),
    all.equal.character(rowname.vt, rownames(prot_imp.od)),
    all.equal.character(rowname.vt, rownames(prot_in.od))
  )) {
    message("\n===== all rownames equal =====\n")
  }
  
  pro_quant_df <-
    cbind.data.frame(prot.od[!names(prot.od) %in% sample.list$bioID]
                     ,
                     prot_in.od[Group.all]
                     ,
                     prot_imp.od
                     ,
                     avg_df
                     ,
                     cv_df.od
                     ,
                     fc_df.od
                     # , pv_df.od
                     ,
                     pv_qv_fdr.od)
} else if (fca) {
  fc_df.od <- fc_df
  class_df.od <- class_df
  names(fc_df.od) <- paste("fc.", names(fc_df), sep = ".")
  names(class_df.od) <- paste("class.", names(class_df), sep = ".")
  if (all(
    all.equal.character(rowname.vt, rownames(prot_imp.od)),
    all.equal.character(rowname.vt, rownames(fc_df.od)),
    all.equal.character(rowname.vt, rownames(class_df.od)),
    all.equal.character(rowname.vt, rownames(prot.od)),
    all.equal.character(rowname.vt, rownames(prot_in.od))
  )) {
    message("\n===== all rownames equal =====\n")
  }
  
  pro_quant_df <-
    cbind.data.frame(prot.od[!names(prot.od) %in% sample.list$bioID]
                     ,
                     prot_in.od[Group.all]
                     ,
                     prot_imp.od
                     ,
                     avg_df
                     ,
                     fc_df.od
                     ,
                     class_df.od)
} else
{
  fc_df.od <- fc_df
  names(fc_df.od) <- paste("fc.", names(fc_df), sep = ".")
  
  if (all(
    all.equal.character(rowname.vt, rownames(prot_imp.od)),
    all.equal.character(rowname.vt, rownames(fc_df.od)),
    all.equal.character(rowname.vt, rownames(prot.od)),
    all.equal.character(rowname.vt, rownames(prot_in.od)))
  ) {
    message("\n===== all rownames equal =====\n")
  }
  
  pro_quant_df <-
    cbind.data.frame(prot.od[!names(prot.od) %in% sample.list$bioID]
                     , prot_in.od[Group.all]
                     , prot_imp.od
                     , avg_df
                     , fc_df.od)
}

#### UNIPROT MAPPING ####
{
  if (!exists("uniprot.tab"))
  {
    .tab <-
      paste("Y:/anno/", dir("Y:/anno/")[grep(dir("Y:/anno/"),
                                             pattern = SPECIES_INFO,
                                             ignore.case = T)], sep = "")
    if (length(.tab) == 1) {
      uniprot.tab <-
        read.table(
          file = .tab,
          header = T,
          sep = "\t",
          stringsAsFactors = F)
    } else if (length(.tab) >  1) {
      uniprot.tab <-
        read.table(
          file = choose.files(
            default = "Y:/anno/*.tab",
            caption = paste("Choose uniprot annotation file", SPECIES_INFO, sep = "--"))
          ,
          header = T,
          sep = "\t",
          stringsAsFactors = F
        )
    } else{
      uniprot.tab <-
        read.table(
          file = choose.files(
            default = "*.tab",
            caption = paste("Choose uniprot annotation file", SPECIES_INFO, sep = "--")
          ),
          header = T,
          sep = "\t",
          stringsAsFactors = F
        )
    }
  }
  # my_gene2(str = pro_quant_df$Majority.protein.IDs[4])
  my_gene2 <- function(str, ...)
  {
    vt <-
      stri_extract_all_regex(
        str = str,
        pattern = "[:alnum:]+(-*)[:alnum:]*",
        simplify = T,
        omit_no_match = F
      )[1, ]
    # vt <- vt0[idx,]
    gi <-
      subset.data.frame(x = uniprot.tab,
                        subset = Entry %in% vt &
                          Cross.reference..GeneID. != "")
    acc.order <- vt[which(vt %in% gi$Entry)]
    choose <- subset.data.frame(gi, subset = Entry == acc.order[1])
    gn <-
      subset.data.frame(x = uniprot.tab,
                        subset = Entry %in% vt &
                          (Gene.names...primary.. != ""))
    GS_ <-
      paste(names(table(gn$Gene.names...primary..))[order(table(gn$Gene.names...primary..), decreasing = T)], collapse = ";")
    EG_ <-
      stri_extract_first(choose$Cross.reference..GeneID., regex = "[:alnum:]+")
    if (length(GS_) == 0)
      GS_ <- NA
    if (length(EG_) == 0)
      EG_ <- NA
    eg.gs <- data.frame(EG = EG_,
                        GS = GS_,
                        stringsAsFactors = F)
    return(eg.gs)
  }
  {
    time0 <- Sys.time()
    str <- pro_quant_df$Majority.protein.IDs
    sp.anno <-
      as.data.frame(t(
        sapply(
          X = str,
          FUN = my_gene2,
          simplify = T
        )))
    
    pro_quant_df$ENTREZID <- unlist(sp.anno$EG)
    pro_quant_df$Alias <- unlist(sp.anno$GS)
    Sys.time() - time0
  }
}


data.df <- pro_quant_df

#### significant function define ####
sig.cutoff <-
  function(data.df,
           imp.rs = im_df.rs,
           quant.type = quant.type,
           ...) {
    ratio.cpr <-
      names(data.df)[grep(x = names(data.df), pattern = "fc..")]
    pv.cpr <-
      names(data.df)[grep(x = names(data.df), pattern = "pv..")]
    qv.cpr <-
      names(data.df)[grep(x = names(data.df), pattern = "qv..")]
    fdr.cpr <-
      names(data.df)[grep(x = names(data.df), pattern = "fd..")]
    unique.cpr <-
      names(data.df)[grep(x = names(data.df), pattern = "^Unique.peptides")]
    imputation.cpr <- names(imp.rs)
    
    sig.cutoff.list <- list()
    
    if (any(sample.list$replicates <  3)) {
      sig.cutoff.df <- cbind.data.frame(imp.rs,
                                        abs(log2(data.df[ratio.cpr])) > lr_co,
                                        data.df[unique.cpr] >= unique_co)
    } else{
      sig.cutoff.df <- cbind.data.frame(
        imp.rs,
        abs(log2(data.df[ratio.cpr])) > lr_co,
        data.df[pv.cpr] < pv_co,
        data.df[qv.cpr] < qv_co,
        data.df[fdr.cpr],
        data.df[unique.cpr] >= unique_co
      )
    }
    
    for (i in 1:nrow(cpr)) {
      sig.cutoff.list[[cpr$names[i]]] <-
        apply(X = sig.cutoff.df[na.omit(c(
          imputation.cpr[i],
          ratio.cpr[i],
          unique.cpr,
          qv.cpr[ifelse(test.num == 2 |
                          use.t, i, 1)],
          fdr.cpr[ifelse(test.num == 2 |
                           use.t, i, 1)],
          pv.cpr[i]
        ))],
        MARGIN = 1,
        FUN = prod)
    }
    sig.cutoff.list[["all"]] <-
      apply(X = as.data.frame(sig.cutoff.list),
            MARGIN = 1,
            FUN = sum)
    return(sig.cutoff.list)
}
#### significant mapping ####
if (fca)
{
  sig.cutoff.list <- list()
  for (i in 1:nrow(cpr))
  {
    sig.cutoff.list[[cpr$fc[i]]] <- !is.na(pro_quant_df[, cpr$class[i]])
  }
  sig.cutoff.list[["all"]] <-
    apply(X = as.data.frame(sig.cutoff.list),
          MARGIN = 1,
          FUN = sum)
  sig <- list()
  sig[["all"]] <- pro_quant_df[sig.cutoff.list$all > 0,]
} else
{
  sig <- list()
  sig.cutoff.list <- sig.cutoff(data.df = pro_quant_df)
  sig[["all"]] <- pro_quant_df[sig.cutoff.list$all > 0,]
}



#### extract significant protein groups ####

if (nrow(sig$all) > 0)
{
  sig.cutoff.df <-
    as.data.frame(sig.cutoff.list, optional = T)[sig.cutoff.list$all > 0,]
  
  colNameFilter <-
    function(x,
             group.vs,
             .group.sep = group.sep,
             .rep.sep = rep.sep,
             .column = column,
             .sample.list = sample.list)
    {
      require(dplyr, lib.loc = lib.loc)
      rs0 <-
        gsub(group.vs, pattern = "([+]|[~]|[&]|[ ]|[-])", replacement = "\\\\\\1")
      rs.sep <-
        gsub(.group.sep,
             pattern = "([+]|[~]|[&]|[ ]|[-])",
             replacement = "\\\\\\1")
      rs.rep.sep <-
        gsub(.rep.sep, pattern = "([+]|[~]|[&]|[ ]|[-])", replacement = "\\\\\\1")
      rs <-
        stri_split_fixed(rs0, pattern = rs.sep, simplify = T)[1, ]
      
      pattern.list <- list(
        group = .sample.list$group.pattern,
        a_b = paste(paste("^(", c(
          rs0, rs[1], rs[2]
        ), ")$", sep = ""), collapse = "|"),
        a.x = paste(
          paste("^(", c(rs[1], rs[2]), rs.rep.sep, ")[:digit:]+$", sep = ""),
          collapse = "|"
        ),
        xx..a_b = paste(
          paste(
            "^[:alnum:]{2,5}[.]{2}(",
            c(rs0, rs[1], rs[2]),
            ")[",
            rs.rep.sep,
            "]*",
            "[:digit:]*$",
            sep = ""
          )
          ,
          collapse = "|"
        ),
        a_b..x = paste(
          paste(
            "^",
            c(rs0, rs[1], rs[2]),
            "[",
            rs.rep.sep,
            "]{2}[:digit:]+$",
            sep = ""
          ),
          collapse = "|"
        )
      )
      column.rs <- c(column[!column %in% reporter],
                     .sample.list$raw[stri_detect_regex(.sample.list$group, pattern = pattern.list$a_b)],
                     x[stri_detect_regex(x, pattern = pattern.list$a_b)],
                     x[stri_detect_regex(x, pattern = pattern.list$a.x)],
                     x[stri_detect_regex(x, pattern = pattern.list$a_b..x)],
                     x[stri_detect_regex(x, pattern = pattern.list$xx..a_b)],
                     c("Symbol", "ENTREZID", "Alias"))
      return(x[x %in% column.rs])
    }
  
  for (ii in cpr$fc)
  {
    sig_i <- sig$all[sig.cutoff.df[ii] > 0,]
    colName.ii <- names(sig_i) %>% colNameFilter(., group.vs = ii)
    sig[[ii]] <- sig_i[colName.ii]
    sig[[paste(ii, "up", sep  =  "..")]] <-
      sig_i[colName.ii] %>% filter(., .[paste("fc..", ii, sep  =  "")]  >  1)
    sig[[paste(ii, "dn", sep  =  "..")]] <-
      sig_i[colName.ii] %>% filter(., .[paste("fc..", ii, sep  =  "")]  <  1)
  }
  if (up.dn) {
    up.vs.dn.0 <- sample.list[["up.vs.dn"]]
    up.vs.dn <- unique(c(
      up.vs.dn.0,
      paste(up.vs.dn.0, "up", sep = ".."),
      paste(up.vs.dn.0, "dn", sep = "..")
    ))
  }
  
}

#### for ChIRP-MS ####
if (fca) {
  for (i in 1:nrow(cpr))
  {
    class_list <- list()
    for (j in c('T', 'C', 'T&C')) {
      k <- sig[[cpr$fc[i]]][, cpr$class[i]]
      class_list[[j]] <-
        subset.data.frame(x = sig[[cpr$fc[i]]], subset =  k %in% j)
    }
    sig$class[[cpr$fc[i]]] <- class_list
  }
}
#### write protein tables ####
# setwd(dir_1)
nrow(pro_quant_df)
write.csv(
  x = pro_quant_df,
  file = paste(dir_1, "all proteins.csv", sep = ""),
  row.names = F
)
write.csv(
  x = sig$all,
  file = paste(dir_1, "significant proteins.csv", sep = ""),
  row.names = F)
#### write data summary ####
source('Y:/R_Script_norm/v201801/sink.R', echo = F)
save(sample.list, file = paste(dir_1, "sample.list", sep = ""))
save(config.list, file = paste(dir_1, "config.list", sep = ""))

#### quantification CV ####
if (all(sample.list$replicates  >  2)) {
  c_v_ <-
    melt(data = cv_df,
         variable.name = "Bio.Group",
         value.name = "CV")
  
  cv_plot <-
    ggplot(data = c_v_, mapping = aes(x = Bio.Group, y = CV)) +
    geom_boxplot() + ylab("CV %") +
    theme(aspect.ratio = 0.6)
  ggsave(paste(dir_1, "CV.tiff", sep = ""),
         cv_plot,
         width = 6,
         dpi = 144)
  ggsave(paste(dir_1, "CV.pdf", sep = ""), cv_plot, width = 6)
}

# setwd(dir_1)

