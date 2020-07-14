#### for profiling ####
source('Y:/R_Script_norm/v201801/extract_data.R', echo = TRUE)
#### for PTM quant ####
# source('Y:/R_Script_norm/v201801/extract_data_ptm.R', echo = TRUE)


#### main analysis ####

{
  source('Y:/R_Script_norm/v201801/correlation.R', echo = TRUE)
  if (sum(sample.list$replicates) >= 3) {
    dir_2 <- dir_1
    if(exists("go_ppi")) remove("go_ppi")
    source('Y:/R_Script_norm/v201801/cluster.R', echo = TRUE)
    source('Y:/R_Script_norm/v201801/rank.R', echo = TRUE)
  }
}
#### LOCAL ####
for (go_ppi in up.vs.dn) {
  source('D:/R_Script_norm/v201801//setwd.R', echo = TRUE)
  source('D:/R_Script_norm/v201801//rank.R', echo = TRUE)
  if (any(sample.list$replicates < 3)) {
    source('D:/R_Script_norm/v201801//fc-plot.R', echo = TRUE)
  } else{
    source('D:/R_Script_norm/v201801/cluster.R', echo = TRUE)
    source('D:/R_Script_norm/v201801/v-plot.R', echo = TRUE)
  }
}

#### NETWORK 1 ####
for (go_ppi in up.vs.dn) {
  source('D:/R_Script_norm/v201801//setwd.R', echo = TRUE)
  if (up.dn) {
    source('D:/R_Script_norm/v201801//GO&KEGG_up-dn.R', echo = TRUE)
  } else{
    source('D:/R_Script_norm/v201801//GO&KEGG.R', echo = TRUE)
  }
}
#### NETWORK 2 ####
for (go_ppi in up.vs.dn) {
  source('D:/R_Script_norm/v201801//setwd.R', echo = TRUE)
  source('D:/R_Script_norm/v201801//PPI.R', echo = TRUE)
}


#### extra analysis ####
for (go_ppi in up.vs.dn) {
  if (up.dn) {
    source('Y:/R_Script_norm/v201801//GO&KEGG_up-dn.R', echo = TRUE)
  }
}

######
for (go_ppi in up.vs.dn) {
  source('Y:/R_Script_norm/v201801//setwd.R', echo = TRUE)
  source('Y:/R_Script_norm/v201801//v-plot.R', echo = TRUE)
  source('Y:/R_Script_norm/v201801//cluster.R', echo = TRUE)
}
#### bar w/ error bar plot ####
source('Y:/R_Script_norm/v201801//bar_errorbar.R', echo = TRUE)


#### ChIRP-MS ####
# source('Y:/R_Script/run_configuration.R', echo = TRUE)
source('Y:/R_Script_norm/v201801//extract_data1.R', echo = TRUE)
source('Y:/R_Script_norm/v201801//venn_RPI.R', echo = TRUE)
### NETWORK 1 ###
for (go_ppi in up.vs.dn) {
  source('Y:/R_Script_norm/v201801//setwd.R', echo = TRUE)
  source('Y:/R_Script_norm/v201801//GO&KEGG_RPI.R', echo = TRUE)
}

### NETWORK 2 ###
for (go_ppi in up.vs.dn) {
  source('Y:/R_Script_norm/v201801//setwd.R', echo = TRUE)
  source('Y:/R_Script_norm/v201801//PPI_RPI.R', echo = TRUE)
}
