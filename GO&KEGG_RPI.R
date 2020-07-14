sig_TC <- sig$class
x = get0(x = "go_ppi", ifnotfound = names(sig_TC)[1])
dir_2_bak <- dir_2
for (y in 1:ifelse(xTC, 3, 2)) {
  if (!dir.exists(paste(dir_2_bak, names(sig_TC[[x]])[y], sep = ""))) {
    dir.create(paste(dir_2_bak, names(sig_TC[[x]])[y], sep = ""),recursive = T)
  }
  dir_2 <- paste(dir_2_bak, names(sig_TC[[x]])[y],"/",sep = "")
  TC_df <- sig_TC[[x]][[y]]
  source('Y:/for run/R_Script_norm/v201801/GO&KEGG.R', echo = TRUE)
}
dir_2 <- dir_2_bak
