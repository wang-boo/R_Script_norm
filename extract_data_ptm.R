
{
  to.win.dir <- function(dir){
    return(stringi::stri_replace_all(dir,fixed = "\\",replacement = "/"))
  }
  if (exists("config.list")) attach(config.list,warn.conflicts = F)
  
  dir_0 <- to.win.dir(paste(choose.dir(default = "Z:/"), "/combined/txt/", sep = ""))
  # config.dir <- paste(dir_0,"config.list",sep = "")
  # setwd(dir_0)
  if (!exists("config.list")) {
    config.file <-
      list.files(path = dir_0,recursive = T,full.names = T)[grep(list.files(path = dir_0,recursive = T,full.names = T),
                                     pattern = "config.list",
                                     fixed = T)]
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
if (exists("config.list")) attach(config.list,warn.conflicts = F)



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
  
  dir_1 <- paste(dir_0, run,"/", sep = "")
  if (!dir.exists(dir_1)) {
    dir.create(dir_1)
  }
}
{
  library("Biobase", lib.loc = lib.loc)
  library("reshape2", lib.loc = lib.loc)
  library("ggplot2", lib.loc = lib.loc)
  library("gplots", lib.loc = lib.loc)
  library("stringi", lib.loc = lib.loc)
  library("RColorBrewer", lib.loc = lib.loc)
  library("extrafont", lib.loc=lib.loc); loadfonts(device = "win")
  library("siggenes", lib.loc=lib.loc)
}

{
  # setwd(dir_0 )
  protein <-
    read.table(
      blank.lines.skip = T,
      file = paste(dir_0,"proteinGroups.txt",sep=""),
      header = T,
      sep = "\t",
      quote = "",
      stringsAsFactors = F
    )
  major.n <- stri_split_fixed(protein$Majority.protein.IDs,pattern = ";",simplify = T)
  major.500 <- major.n[,1:min(ncol(major.n),500)]
  protein$Majority.protein.IDs <- apply(major.500,1,function(x) paste(x[x!=""],collapse = ";"))
  
  ptm.txt <- dir(dir_0)[grep(x = dir(dir_0),pattern = paste("(",ptm.site,")",sep = ""))]
  ptm <-
    read.table(
      blank.lines.skip = T,
      file = paste(dir_0,ptm.txt,sep=""),
      header = T,
      sep = "\t",
      quote = "",
      stringsAsFactors = F
    )
  names(ptm)
}
{
  pt <- c("(Reporter.intensity.corrected.)[0-9][_]{3}[0-9]","(Intensity)[.]{1}[0-9][_]{3}[0-9]")
  names(pt) <- c("TMT","LFQ")
  reporter <- names(ptm)[ grep(pattern = pt[quant.type],x = names(ptm))]
}




#### BIO GROUP ############################################################################################

# setwd(dir_1)
sample.list.dir <- paste(dir_1,"sample.list",sep = "")
if (file.exists(sample.list.dir)) {
  message("\n< sample.list exsits >, loading it? [Y] / [N]\n")
  load.file <- readline(prompt = ">>\n")
  
  if (load.file %in% c("y", "Y","")) {
    load(file = sample.list.dir,verbose = T)
    Group.all <- sample.list$bioID
  }
}else{
  load.file <- "N"
}

if(!file.exists(sample.list.dir) | load.file %in% c("N","n")){
  check <- T
  while (check) {
    samples.info <- data.frame(raw = reporter, stringsAsFactors = F)
    samples.info.df <-
      data.frame(
        stri_split_fixed(
          str = samples.info$raw,
          pattern = "__",
          simplify = T
        ),stringsAsFactors = F
      )
    samples.info.df.p <- subset.data.frame(samples.info.df,subset = X2=="_1",select = X1)
    message("======= sample ID structure =======\n")
    print(samples.info.df.p)
    print(data.frame(up = config.list$up,down = config.list$dn))
    message("===================================\n")
    message(
      "\n please input bioGroup pattern regular expression ### '[:alnum:]+[.]{1}' for 'BioGroup.' ###\n"
    )
    pattern.read <-
      readline(prompt = "new regular expression pattern = \ndefault pattern = [:alnum:]+")
    if (pattern.read == "") {
      group.pattern <- "[:alnum:]+"
      message(paste("\n=====", group.pattern, "=====\n", sep = " "))
    } else{
      group.pattern <- pattern.read
      message(paste("\n=====", group.pattern, "=====\n", sep = " "))
    }
    
    message(
      "\n please input bioGroup pattern regular expression position ### '[Enter] for prefix, [0] for postfix' ###\n"
    )
    position.read <-
      readline(prompt = "pattern position = \ndefault position = prefix")
    if (position.read == "") {
      group.mode <- "first" ### first/last
      message("\n ===== prefix =====")
    } else{
      group.mode <- "last"
      message("\n ===== postfix =====")
    }
    
    sample.number <- length(reporter)/3
    message(
      "\nplease input [y] to assign new bioID for each reporter, OR press [Enter] directly if use raw ID"
    )
    if (readline(prompt = ">>") == "") {
      samples.info$bioID <- reporter
      samples.info$group <-
        stri_extract(str = samples.info$bioID,
                     regex = group.pattern,
                     mode = group.mode)
      
      message("========== use raw bioID ==========\n")
      print(samples.info)
      message("===================================\n")
    } else{
      if(is.null(samples.info$bioID))  samples.info$bioID <- rep("",sample.number) 
{
      .samples.info.cast <-
        as.data.frame(acast(samples.info.df, formula = X1 ~ X2))
      samples.info.raw.df <-
        data.frame(raw = rownames(.samples.info.cast),
                   .samples.info.cast,
                   row.names = NULL)
      samples.info.df <-
        na.omit(edit(name = samples.info.raw.df, edit.row.names = F))
      names(samples.info.df)[5] <- "bioID"
      samples.info.df$group <-
        stri_extract(str = samples.info.df$bioID,
                     regex = group.pattern,
                     mode = group.mode)
      suppressWarnings(samples.info.melt <- melt(samples.info.df, id.vars = c("raw", "bioID", "group")))
      samples.info.short <- data.frame(
        raw =   paste(samples.info.df$raw,  samples.info.df$value,sep = ""),
        bioID = paste(samples.info.df$bioID,samples.info.df$value,sep = ""),
        group = paste(samples.info.df$group,samples.info.df$value,sep = ""),
        stringsAsFactors = F
      )
      samples.info.long <- data.frame(
        raw = paste(samples.info.melt$raw,samples.info.melt$value,sep = ""),
        bioID = paste(samples.info.melt$bioID,samples.info.melt$value,sep = ""),
        group = paste(samples.info.melt$group,samples.info.melt$value,sep = ""),
        stringsAsFactors = F
      )
}      
      message("========== use new bioID ==========\n")
      print(samples.info.long)
      print(samples.info.short)
      message("===================================\n")
    }
    message("\ncheck bioID assignment, press [Enter] if done, else press [N] to reassign bioID\n")
    check <- ifelse(readline() == "", F, T)
    if (!check) {
      message("======== assign bioID done ========\n")
    }
  }
  sample.list <- as.list(samples.info.short)
  sample.list[["group.mode"]] <- group.mode
  sample.list[["group.pattern"]] <- group.pattern
  Group.all <- sample.list$bioID
}

#### bio.compare config ####
{
  u <- table(up)
  d <- -table(dn)
  u.d <- sort(c(d, u), decreasing = T)
  u.d.df <- data.frame(order = u.d)
  agg.cpr <-
    aggregate.data.frame(x = u.d.df,
                         by = list(groupName = names(u.d)),
                         FUN = sum)
  sample.list[["group.names.order"]] <-
    agg.cpr$groupName[order(agg.cpr$order, decreasing = F)]
  sample.list[["replicates"]] <- table(sample.list$group)
}
################################################################################################
{
  # names(ptm)
  if (T) {
    col.choose <-
      c(
        "Protein.group.IDs",
        "Protein",
        "Amino.acid",
        "Position",
        "Score.for.localization",
        "Localization.prob",
        "GlyGly..K..Probabilities",
        "GlyGly..K..Score.diffs",
        "Position.in.peptide",
        "Number.of.GlyGly..K.",
        "PEP",
        "Score"
      )
    column <- c(col.choose,reporter)
    
    prot.ptm <- ptm[!(is.na(ptm$Position)|ptm$Number.of.GlyGly..K.==""), column]
    # names(prot.ptm)
    prot.ptm.melt <- melt(data = prot.ptm,id.vars = col.choose)
    prot.ptm.melt$bioID <- stri_split_fixed(str = prot.ptm.melt$variable,pattern = "___",simplify = T)[,1]
    prot.ptm.melt$siteID <- stri_split_fixed(str = prot.ptm.melt$variable,pattern = "___",simplify = T)[,2]
    # names(prot.ptm.melt)
    prot.ptm.dcast <- dcast(prot.ptm.melt[c(col.choose,"bioID","siteID","value")] ,formula = ...~bioID,value.var = "value",fun.aggregate = mean)
    row.names(prot.ptm.dcast) <-
      paste(
        prot.ptm.dcast$Protein,
        paste(prot.ptm.dcast$Amino.acid,
              prot.ptm.dcast$Position, sep = ""),
        prot.ptm.dcast$siteID,
        sep = "-"
      )
  } 
  prot <- prot.ptm.dcast
  prot <- prot[-grep("REV_", x = prot$Protein), ]
  prot <- prot[-grep("CON_", x = prot$Protein), ]
  prot_bak <- prot
  # prot <- prot_bak
  prot <-
    prot[rowSums(prot[sample.list$raw] > 0) >= length(sample.list$raw) * missing.ratio, ]
  prot[Group.all] <- prot[c(sample.list$raw)]
  # names(prot)
{
  if (T) {
    ###  CAL FACTORS: sum median custom ####
    prot_rpt_sub <- prot[sample.list$raw]
    func <- mean
    median_factor <-
      apply(X = prot_rpt_sub,
            MARGIN = 2,
            FUN = func) /
      (apply(X = prot_rpt_sub,
             MARGIN = 2,
             FUN = func)[1])
    nor_mt <-
      t(t(prot_rpt_sub) / median_factor)          ### < CONFIG > ###
    colMeans(nor_mt, na.rm = T)
    
    prot[sample.list$bioID] <- nor_mt
  }
}
  #### imputation ####
  if (T) {
    prot_imp <- log2(prot[Group.all])
    prot_imp.bl <-
      as.data.frame(prot[Group.all] == 0, optional = T)
    prot_imp.df <- as.data.frame(melt(prot_imp), optional = T)
    
    prot_imp.df$value[is.infinite(prot_imp.df$value)] <- NA
    prot_imp.df$imputation <- F
    prot_imp.df$imputation[is.na(prot_imp.df$value)] <- T
    
    if (imp.mode == "s") {
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
      
      for (im in sample.list$bioID) {
        prot_imp[prot_imp.bl[, im], im] <-
          rnorm(n = sum(prot_imp.bl[, im]),
                sd = s[, im],
                mean = m[, im])
      }
      
    } else if (imp.mode == "a") {
      imp.avg <- mean(prot_imp.df$value, na.rm = T)
      imp.sd  <-   sd(prot_imp.df$value, na.rm = T)
      
      m <- imp.avg - imp.downshift * imp.sd
      s <- imp.width * imp.sd
      
      for (im in sample.list$bioID) {
        prot_imp[prot_imp.bl[, im], im] <-
          rnorm(n = sum(prot_imp.bl[, im]),
                sd = s,
                mean = m)
      }
    }
    # prot_in <- list()
    prot_imp.bl$IDs <- row.names(prot_imp.bl)<- rownames(prot)
    
    # prot_in[sample.list$bioID] <- prot[c(sample.list$bioID)]
    prot_in <- cbind.data.frame(IDs=row.names(prot_imp),2 ^ (prot_imp[Group.all]))
    
  } else if (quant.type %in% c("count")) {
    {
      prot_imp <- log2(prot[Group.all])
      prot_imp.bl <-
        as.data.frame(prot[Group.all] == 0, optional = T)
      prot_imp.df <- as.data.frame(melt(prot_imp), optional = T)
      prot_imp.df$imputation <- F
      row.names(prot_imp.bl) <- prot$IDs
      prot_imp.bl$IDs <- row.names(prot_imp.bl)
      
    }
    prot_in <- prot[c("IDs", sample.list$bioID)]
    names(prot_in)
  }
  # rownames(prot) <- prot$IDs
  # rownames(prot_in) <- prot_in$IDs
}
#### FC / imputation AND/OR CV ####
  # sample.list$group.names.order <- c("a","b","c")
{
  {
    pro_imp_m <- melt(
      data = prot_imp.bl,
      id.vars = "IDs",
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
      acast(data = pro_imp_m,
            formula = IDs ~ group,
            fun.aggregate = mean)
    imp_df <- as.data.frame(x = imp_mt, optional = T)
    names(imp_df)
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
      acast(data = pro_m,
            formula = IDs ~ group,
            fun.aggregate = mean)
    avg_df <- as.data.frame(x = avg_mt, optional = T)
    # names(avg_df)
    test.num <- length(unique(pro_m$group))
  }
  
  {
    if (all(sample.list$replicates > 2)) {
      cv_mt <-
        acast(
          data = na.omit(pro_m),
          formula = IDs ~ group,
          fun.aggregate = function(x)
            100 * sd(x) / mean(x)
        )
      cv_df <- as.data.frame(cv_mt, optional = T)
      
      sd_mt <-
        acast(
          data = na.omit(pro_m),
          formula = IDs ~ group,
          fun.aggregate = sd
        )
      sd_df <- as.data.frame(sd_mt, optional = T)
    }
  }
}
#### configure pv.mt structure ####
if (length(sample.list$group.names.order) >= 2) {
  test.data <-
    data.frame(
      group = factor(
        x = rep(x = sample.list$group.names.order, each = 3),
        levels = sample.list$group.names.order
      ),
      value = rnorm(length(sample.list$group.names.order) * 3)
    )
  a0 <- aov(formula = value ~ group, data = test.data)
  aov0 <- TukeyHSD(a0, conf.level = 0.95)$group
  pv.mt <- matrix(nrow = nrow(avg_df), ncol = nrow(aov0))
  colnames(pv.mt) <- rownames(aov0)
  row.names(pv.mt) <- row.names(avg_df)
  cpr <-
    cbind.data.frame(
      vs = stri_split_fixed(
        str = rownames(aov0),
        pattern = "-",
        simplify = T
      ),
      stringsAsFactors = F
    )
  cpr$names <- rownames(aov0)
  cpr$fc <- paste(cpr$vs.1, "-", cpr$vs.2, sep = "")
  cpr$imputation <- paste("im..", cpr$vs.1, "-", cpr$vs.2, sep = "")
  
  cpr$pv <- paste("pv..", rownames(aov0), sep = "")
  cpr$class <- paste("class..", cpr$fc, sep = "")
  sample.list[["contrast"]] <- cpr
  up.vs.dn = paste(up, dn, sep = "-") ### compare selected
  sample.list[["up.vs.dn"]] <-  up.vs.dn ### compare selected
}

#### T.TEST OR ANOVA #########################################################
if (all(sample.list$replicates > 2)) {
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
  if (test.num == 2 | use.t == use.t) {
    for (cpr.name in cpr$names) {
      fdr.cl.idx.x <- fdr.cl.all %in% cpr[match(cpr.name, cpr$names), 1:2]
      fdr.cl.x <- fdr.cl.all[fdr.cl.idx.x]
      fdr.dt.x <- prot_in[-1][fdr.cl.idx.x]
      sam.out_list[[cpr.name]] <-
        sam(data = scale(log(fdr.dt.x)),
            cl = fdr.cl.x,
            rand = 123)
      pv_list[[cpr.name]] <- sam.out_list[[cpr.name]]@p.value
      qv_list[[cpr.name]] <- sam.out_list[[cpr.name]]@q.value
      if(sum(summary(sam.out_list[[cpr.name]])@mat.fdr[,"Called"])==0){
        fdr_co_idx <- rep(T,nrow(fdr.dt.all))
      }else{
        delta.fdr.0 <- findDelta(sam.out_list[[cpr.name]], fdr = fdr_co)
        delta.fdr <- rbind(delta.fdr.0,delta.fdr.0) ### bug fixed for if delta.fdr.0 has only one row
        delta_min <- min(as.data.frame(delta.fdr)$Delta)
        row.sig_vt <- summary(sam.out_list[[cpr.name]], delta_min)@row.sig.genes
        fdr_co_idx <- rep(F,nrow(prot_in))
        fdr_co_idx[row.sig_vt] <- T
      }
      fdr_list[[cpr.name]] <- fdr_co_idx
    }
    pv_df <- as.data.frame(pv_list, optional = T)
    rownames(pv_df) <- prot_in[,"IDs"]
    qv_df <- as.data.frame(qv_list, optional = T)
    rownames(qv_df) <- prot_in[,"IDs"]
    fdr_df <- as.data.frame(fdr_list,optional = T)
    rownames(fdr_df) <- prot_in[,"IDs"]
  } else  if (test.num > 2) {
## three or more case SAM & anova-TukeyHSD
    cpr.name <- paste(levels(x = fdr.cl.all), collapse = "-")
    
    fdr.dt.all <- prot_in[-1][names(prot_in)[-1]]
    sam.out_list[[cpr.name]] <-
      sam(data = log(fdr.dt.all),
          cl = fdr.cl.all,
          rand = 123)
    plot(sam.out_list[[cpr.name]],0.3)
    pv_list[[cpr.name]] <- sam.out_list[[cpr.name]]@p.value
    # pv_df0 <- as.data.frame(pv_list, optional = T)
    # rownames(pv_df0) <- prot_in["IDs"]
    qv_list[[cpr.name]] <- sam.out_list[[cpr.name]]@q.value
    
    if(sum(summary(sam.out_list[[cpr.name]])@mat.fdr[,"Called"])==0){
      fdr_co_idx <- rep(T,nrow(fdr.dt.all))
    }else{
      delta.fdr <- findDelta(sam.out_list[[cpr.name]], fdr = fdr_co)
      delta_min <- min(as.data.frame(delta.fdr)$Delta)
      row.sig_vt <- summary(sam.out_list[[cpr.name]], delta_min)@row.sig.genes
      fdr_co_idx <- rep(F,nrow(prot_in))
      fdr_co_idx[row.sig_vt] <- T
    }
    fdr_list[[cpr.name]] <- fdr_co_idx
    
    qv_df <- as.data.frame(qv_list, optional = T)
    rownames(qv_df) <- prot_in[,"IDs"]
    fdr_df <- as.data.frame(fdr_list,optional = T)
    rownames(fdr_df) <- prot_in[,"IDs"]

    anova.df <- function(data.df, cl.vt, ...) {
      Tukey.pv <- as.data.frame(optional = T, t(
        apply(
          X = data.df,
          MARGIN = 1,
          FUN = function(.x = X,
                         .cl.vt = cl.vt,
                         ...) {
            aov.lm <-
              aov(formula = value ~ variable,
                  data = data.frame(value = .x, variable = .cl.vt))
            TukeyHSD(aov.lm)$variable[, "p adj"]
          }
        )
      ))
      return(Tukey.pv)
    }
    pv_df <- anova.df(fdr.dt.all,fdr.cl.all)
      # anova.df(data.df = fdr.dt.all[1:14, ], factor(c("a", "b", "c", "a", "b", "c")))

  }
  
}

#### average fc ####
fc_list <- list()
im_list <- list()
{
  for (j in 1:nrow(cpr)) {
    fc_list[[cpr$fc[j]]] <-  avg_df[, cpr$vs.1[j]] / avg_df[, cpr$vs.2[j]]
    im_list[[cpr$im[j]]] <-  imp_df[, cpr$vs.1[j]] < imp_co | imp_df[, cpr$vs.2[j]] < imp_co
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
#### merge data #######
  rowname.vt <- rownames(avg_df)
  prot.od <- prot[rowname.vt, ]
  prot_in.od <- prot_in[rowname.vt, ]
  
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
  pv_qv_fdr.od <- cbind.data.frame(pv_df..,qv_df..,fdr_df..)[rowname.vt,]
  
  
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
                     , prot_in.od[Group.all]
                     , prot_imp.od
                     , avg_df
                     , cv_df.od
                     , fc_df.od
                     # , pv_df.od
                     , pv_qv_fdr.od
    )
}else {
  fc_df.od <- fc_df
  names(fc_df.od) <- paste("fc.", names(fc_df), sep = ".")
  
  if (all(
    all.equal.character(rowname.vt, rownames(prot_imp.od)),
    all.equal.character(rowname.vt, rownames(fc_df.od)),
    all.equal.character(rowname.vt, rownames(prot.od)),
    all.equal.character(rowname.vt, rownames(prot_in.od))
  )) {
    message("\n===== all rownames equal =====\n")
  }
  
  pro_quant_df <-
    cbind.data.frame(prot.od[!names(prot.od) %in% sample.list$bioID]
                     , prot_in.od[Group.all]
                     , prot_imp.od
                     , avg_df
                     , fc_df.od)
}
#### DEBUG ####
  # names(pro_quant_df)
  # names(prot_imp.od)
#### UNIPROT MAPPING ####
{
  
  if (!exists("uniprot.tab")) {
    .tab <- paste("Y:/anno/",dir("Y:/anno/")[grep(dir("Y:/anno/"),pattern = SPECIES_INFO,ignore.case = T)],sep = "")
    if(length(.tab) == 1){
    uniprot.tab <-
      read.table(
        file = .tab,
        header = T,
        sep = "\t",
        stringsAsFactors = F
      )
    }else if(length(.tab) >1){
    uniprot.tab <-
      read.table(
        file = choose.files(
          default = "Y:/anno/*.tab",
          caption = paste("Choose uniprot annotation file", SPECIES_INFO, sep = "--")
        ),
        header = T,
        sep = "\t",
        stringsAsFactors = F
      )
    }else{
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
  # my_gene2(str = pro_quant_df$IDs[4])
  my_gene2 <- function(str, ...) {
    vt <-
      stri_extract_all_regex(
        str = str,
        pattern = "[:alnum:]+(-*)[:alnum:]*",
        simplify = T,
        omit_no_match = F
      )[1,]
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
    str <- protein$Majority.protein.IDs
    sp.anno <-
      as.data.frame(t(sapply(
        X = str,
        FUN = my_gene2,
        simplify = T
      )))
    
    protein$ENTREZID <- unlist(sp.anno$EG)
    protein$Alias <- unlist(sp.anno$GS)
    Sys.time() - time0
  }
}

  # row.names(protein) <- protein$id
{  assign.map <- function(str, ...) {
    entrezid.fun <- function(str) {
      paste(protein[stri_extract_all(str = str, regex = "[:digit:]+")[[1]], "ENTREZID"], collapse = "|")
    }
    alias.fun <- function(str) {
      paste(protein[stri_extract_all(str = str, regex = "[:digit:]+")[[1]], "Alias"], collapse = "|")
    }
    eid <- sapply(str, entrezid.fun)
    als <- sapply(str, alias.fun)
    data.frame(ENTREZID = eid,Alias = als,id = names(eid))
  }
  group.IDs.map <- assign.map(str = pro_quant_df$Protein.group.IDs)
  .pro_quant_df <- pro_quant_df
  if(all(.pro_quant_df$Protein.group.IDs==group.IDs.map$id)){
    pro_quant_df <- cbind.data.frame(.pro_quant_df,group.IDs.map[-3])
  }
}  
#### DEBUG ####  
  # names(pro_quant_df)
  data.df <- pro_quant_df
#### significant function define ####
sig.cutoff <- function(data.df,imp.rs = im_df.rs,quant.type = quant.type,...){    
  ratio.cpr <- names(data.df)[grep(x = names(data.df),pattern = "fc..")]
  pv.cpr <- names(data.df)[grep(x = names(data.df),pattern = "pv..")]
  qv.cpr <- names(data.df)[grep(x = names(data.df),pattern = "qv..")]
  fdr.cpr <- names(data.df)[grep(x = names(data.df),pattern = "fd..")]
  loc.cpr <- names(data.df)[grep(x = names(data.df),pattern = "^Localization.prob")]
  imputation.cpr <- names(imp.rs)
  
  sig.cutoff.list <- list()
  
  if(any(sample.list$replicates <3) ){
    sig.cutoff.df <- cbind.data.frame(
      imp.rs,
      abs(log2(data.df[ratio.cpr])) > lr_co,
      data.df[unique.cpr] >= unique_co)
  }else{
    sig.cutoff.df <- cbind.data.frame(
      imp.rs,
      abs(log2(data.df[ratio.cpr])) > lr_co,
      data.df[pv.cpr] < pv_co,
      data.df[qv.cpr] < qv_co,
      data.df[fdr.cpr],
      data.df[loc.cpr] >= loc_co*0.99)
  }

  for (i in 1:nrow(cpr)) {
    
    sig.cutoff.list[[cpr$names[i]]] <- apply(X = sig.cutoff.df[na.omit(c(imputation.cpr[i],
                                                                         ratio.cpr[i],
                                                                         loc.cpr,
                                                                         qv.cpr[ifelse(test.num == 2|use.t,i,1)],
                                                                         fdr.cpr[ifelse(test.num == 2|use.t,i,1)],
                                                                         pv.cpr[i]))],
                                             MARGIN = 1, 
                                             FUN = prod)
  }
  sig.cutoff.list[["all"]] <- apply(X = as.data.frame(sig.cutoff.list),MARGIN = 1,FUN = sum)
  return(sig.cutoff.list)
}
#### significant mapping ####
{
  sig <- list()
  sig.cutoff.list <- sig.cutoff(data.df = pro_quant_df)
  sig[["all"]] <- pro_quant_df[sig.cutoff.list$all > 0, ]
}

#### DEBUG ####
  # summary(sig$all)
  # sss <- as.data.frame(sig.cutoff.list,optional = T)
  # names(sss)

#### extract significant protein groups ####

if(nrow(sig$all)>0) {
  sig.cutoff.df <-
    as.data.frame(sig.cutoff.list, optional = T)[sig.cutoff.list$all > 0, ]
  # names(sig.cutoff.df)
  # ii=cpr$names[4]
  for (ii in cpr$names) {
    sig_i <- sig$all[sig.cutoff.df[ii] > 0, ]
    rm.. <-
      unique(sample.list$group[-which(sample.list$group %in% stri_split_fixed(ii, pattern = "-", simplify = T)[1, ])])
    if (length(rm..) > 0) {
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
    sig[[ii]] <- sig_rsv
  }
  # message(ii)
  # print(names(sig_i)[rs.idx])
}



#### DEBUG ####
  # names(sig[[1]])
  # names(pro_quant_df)
#### write protein tables ####
# setwd(dir_1)
write.csv(
  x = pro_quant_df,
  file = paste(dir_1,"all PTM.csv",sep = ""),
  row.names = T
)
write.csv(
  x = sig$all,
  file = paste(dir_1,"significant PTM.csv",sep = ""),
  row.names = T
)
#### write data summary ####
source('Y:/R_Script/sink_ptm.R', echo = F)
save(sample.list,file = paste(dir_1,"sample.list",sep = ""))
save(config.list,file = paste(dir_1,"config.list",sep = ""))

#### quantification CV ####
if(all(sample.list$replicates>2)){
  c_v_ <-
    melt(data = cv_df,
         variable.name = "Bio.Group",
         value.name = "CV")
  
  cv_plot <- ggplot(data = c_v_, mapping = aes(x = Bio.Group, y = CV)) + geom_boxplot() + ylab("CV %")
  ggsave(paste(dir_1,"CV.tiff",sep = ""), cv_plot, width = 4, height = 4,dpi = 144)
  ggsave(paste(dir_1,"CV.pdf",sep = ""), cv_plot, width = 4, height = 4)
}

# setwd(dir_1)

