config.list <- list(
  lib.loc = "~/R/win-library/3.4",  ### R library
  SPECIES_INFO = "HUMAN", ### MOUSE / HUMAN / RAT / DOG
  run = "RUN-0.1",  ### directory of the current analysis run
  
  quant.type = "LFQ", ### "TMT / iBAQ / LFQ / count
  
  missing.ratio = 1/2, ### ratio of MIN NON-missing value for imputation
  imp.mode = "s", ### s/a for seperate/all column
  imp.downshift = 1.8, ### mean shift of imputation distribution
  imp.width =0.3, ### width of imputation distribution
  
  use.t = T, ### use t.test instead of ANOVA for group > 2
  lr_co = log2(2),   ### log2 ratio cutoff
  pv_co = 0.05,      ### p-value cutoff
  unique_co = 2,     ### unique peptides cutoff
  imp_co = 0.5,      ### imputation ratio less than
  spc_co = 2,        ### for quant.type = "count", min spectra count
  qv_co = 1,      ### local FDR or q-value cutoff
  fdr_co = .9,     ### for FDR control
  delta_co = 0.1,
  up = c("H"),  ### the numerator of compare
  dn = c("NC"),    ### the denominator of compare
  
  up.dn = F, ### separate up/down protein for GO/KEGG
  
  xTC =  F,  ### for IP/ChIRP-MS only, also analyse T&C OR just T / C
  # fca = F,  ### for IP/ChIRP-MS only
  
  vp.text.label =NULL,
  label.num = Inf,  ### max number of labels in valcano plot 
  group.1 = "Ctrl",  ### normalize to 1 in error bar plot
  # vp.text.label = c("Zfr","Zranb2","Zzef1","Zpr1" ), ### manual valcano plot labeling
  # bar.plot.symbol = c("Acadl","Alox15","Lipc","Apoa5","Serpine1","Fasn","Psapl1","Tgfb1"), ### manual bar/error plot
  
  dpi.set = 600,   ### figure resolution
  fig.type = "tiff",   ### figure file type
  
  GO.max = 10,   ###  max  GO iterm for plot
  GO.dscpt.length = 50,   ###  max character number of GO description 
  custom = F,   ### custom analysis
  DATE = Sys.Date()
)
