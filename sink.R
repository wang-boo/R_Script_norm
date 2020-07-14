options(digits = 4)
sink(
  file = paste(dir_1,"REPORT SUMMARY.txt",sep = ""),
  append = F,
  split = FALSE,
  type = "output"
)
print(DATE)
cat(
  "\n==============================================   Service Report Summary   ==============================================\n\n"
)

cat(paste("Protein groups number (with 1% decoy):", nrow(protein)), "\n")
cat(paste("Peptides number (with 1% decoy):", nrow(peptide), "\n"))
cat(paste(
  "Protein groups number (quantifiable, without decoy):",
  nrow(pro_quant_df)
),
"\n")
if (quant.type %in% c("iBAQ", "LFQ")) {
  cat(
    paste(
      "Max NON-missing value allowed for each reserved protein group (imputation):",
      missing.ratio,
      "\n"
    )
  )
}
cat("\n========== Classified protein groups number(compare group) ==========\n")
cat(paste("\n========== Fold change cutoff:", 2 ^ lr_co))
if (all(sample.list$replicates >= 3))
  cat(paste("\n========== p-value cutoff:", pv_co))
cat(paste("\n========== Unique peptide cutoff:", unique_co))

if (fca)
  cat(paste("\n========== spectra count cutoff:", spc_co))
if (is.null(up.vs.dn))
  up.vs.dn <- sample.list$contrast$names
for (s in up.vs.dn) {
  if (fca) {
    cat(paste("\n\n\n\n\n\n", s, "\n"))
    cat(
      "\n========================================================================================================================\n\n"
    )
    s1 <- sig$class[[s]]
    for (ss in names(s1)) {
      cat(paste("\n\n\n\n\n\n", ss, "\n"))
      cat(
        "\n    ================================================================================================================    \n\n"
      )
      print(nrow(s1[[ss]]))
      print(summary(s1[[ss]]))
    }
  } else{
    cat(paste("\n\n\n\n\n\n", s, "\n"))
    s1 <- sig[[s]]
    cat(
      "\n========================================================================================================================\n\n"
    )
    print(nrow(s1))
    print(summary(s1))
  }
}
sink()
