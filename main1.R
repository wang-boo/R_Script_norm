#### for profiling ####
source('D:/R_Script_norm/v201801/extract_data.R', echo = TRUE)
#### for PTM quant ####
# source('Y:/R_Script_norm/v201801/extract_data_ptm.R', echo = TRUE)


#### main analysis ####

{
  source('D:/R_Script_norm/v201801/correlation.R', echo = TRUE)
  if (sum(sample.list$replicates) >= 3) {
    dir_2 <- dir_1
    if(exists("go_ppi")) remove("go_ppi")
    source('D:/R_Script_norm/v201801/cluster.R', echo = TRUE)
    source('D:/R_Script_norm/v201801/rank.R', echo = TRUE)
  }
}
