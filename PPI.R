#### PPI ###########################################################################################
x = get0(x = "go_ppi", ifnotfound = names(sig)[2])
{
  library("STRINGdb")
  options(timeout = 3600)
}
{
  string_db <- STRINGdb$new(version = "10",
                            species = species,
                            # score_threshold = 0,
                            # input_directory = "C:/Users/hq/Documents/R/win-library/3.3/STRINGdb/library"
                            # input_directory = "C:\\Users\\HQ\\Documents\\R\\win-library\\3.3\\STRINGdb\\database"
                            input_directory = "C://Program Files//R//R-3.4.3//library//STRINGdb//R/")
}
if (fca)  {
  gene_id <- TC_df
} else{
  gene_id <- sig[[x]]
}

sub_sig <- gene_id 
if (nrow(sub_sig) > 0  ) {
  if (T) {
    if (!fca)  {
      sub_sig$log.ratio <-
        log(sub_sig[, grep(pattern = "fc..", x = names(sub_sig))])
    }
    # names(sub_sig)
    mapped <-
      string_db$map(
        my_data_frame = sub_sig ,
        my_data_frame_id_col_names = "Symbol",
        removeUnmappedRows = F
      )
    
    sort.id <- names(mapped)[grep(names(mapped),pattern = ifelse(is.null(grep(names(mapped),pattern = "pv..",fixed = T)),"pv..","fc.." ))]
    mapped <- mapped[order(mapped[sort.id]),]
    logCol <-
      stri_replace_all(
        str = ifelse(fca, names(sub_sig)[grep(pattern = "fc..", x = names(sub_sig))], "log.ratio"),
        replacement = ".",
        regex = "[+|-]"
      )
    mapped_col <-
      string_db$add_diff_exp_color(mapped, logFcColStr = logCol)
    # desc <- string_db$add_proteins_description(x)
    payload_id <-
      string_db$post_payload(mapped_col$STRING_id, colors = mapped_col$color)
    tiff(
      filename = paste(dir_2,"PPI.tiff", sep = ""),
      width = 2400,
      height = 2400,
      res = 288
    )
    if (nrow(sub_sig) < 400) {
      string_db$plot_network(
        mapped$STRING_id,
        required_score = 200,
        payload_id = payload_id,
        add_link = F
      )
    } else{
      string_db$plot_ppi_enrichment(mapped$STRING_id, quiet = T)
      
    }
    
  } else {
    mapped <-
      string_db$map(
        my_data_frame = sub_sig ,
        my_data_frame_id_col_names = "Symbol",
        removeUnmappedRows = F
      )
    tiff(
      filename = paste(dir_2, "PPI.tiff", sep = ""),
      width = 2400,
      height = 2400,
      res = 288
    )
    string_db$plot_network(mapped$STRING_id,
                           required_score = 200,
                           add_link = F)
    dev.off()
  }
  
  
  description <- string_db$add_proteins_description(mapped)
############################# write tables #################################################################

  write.csv(
    x = description[c(2:length(description), 1)],
    file = paste(dir_2,paste("significant proteins",  x, "csv", sep = "."),sep = ""),
    row.names = F
  )
}



# setwd(dir_1)

