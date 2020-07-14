x = get0(x = "go_ppi", ifnotfound = names(sig)[1])
# setwd(dir_1)
if (x != "all") {
  if (!dir.exists(paste(dir_1, x, sep = ""))) {
    dir.create(paste(dir_1, x, "",sep = ""))
  }
  # setwd(paste(dir_1, x, sep = "/"))
  dir_2 <- paste(dir_1, x,"/",sep = "")
} else{
  dir_2 <- dir_1
  # setwd(dir_1)
}
