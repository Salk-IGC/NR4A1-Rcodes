# ==========================================================
# 01_RNAseq_preprocess.R
# LKR10 RNA-seq: sample metadata, QC, DESeq2 objects, normalized counts and
# differential expression results for the two sequencing runs used in the paper.
#
#   Run 1: sgNT, sgLkb1, sgCrtc2, sgCreb1; vehicle / IFNg   (36 samples)
#   Run 2: sgNT, sgLkb1, sgNr4a1;          vehicle / IFNg   (24 samples)
#
# Inputs (raw count tables generated from GSE302655):
#   data/run1/raw.txt
#   data/run2/raw.txt
#   Sample metadata is parsed from the count-table column names, which must
#   have the format <sample>.<replicate><treatment>, e.g. "sg1.CREB1.1F".
#   treatment: v = vehicle, F = IFNg; replicate: 1-3
#
# Required R packages: DESeq2, here, qdapRegex, and (used inside myExplore())
#   stringr, gplots, ggplot2, pheatmap, cluster, VennDiagram, openxlsx
#
# Outputs (results/QC/):
#   run1_QC.RDS   list for run 1: [[1]] DESeq2 object, [[2]] transformed (VST) data
#   run2_QC.RDS   list for run 2: [[1]] DESeq2 object, [[2]] transformed (VST) data
#   run1_normalized_counts.csv, run2_normalized_counts.csv
#   run1_QC_plots.pdf, run2_QC_plots.pdf  (written by myExplore())
# Outputs (results/DE/):
#   run1_results.RDS, run2_results.RDS   lists of DE result tables, named by comparison
#   (per-comparison files are in results/DE/run1/ and results/DE/run2/)
# ==========================================================

library(DESeq2)
library(here)
library(qdapRegex)
library(ggplot2)
library(pheatmap)

# Custom helper functions: processRawTxt(), processAnnoTxt(),
# myExplore(), myDiffExp()
source(here("R", "rnaseq_functions.R"))

qc_dir <- here("results", "QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Helper: run QC, save the object, write normalized counts ----
run_qc <- function(raw, design, anno, ref, name, rds_name, norm_csv_name) {
  old_wd <- setwd(qc_dir)
  on.exit(setwd(old_wd), add = TRUE)

  qc <- myExplore(counts     = raw,
                  design     = design,
                  colData1   = "guide_KO",
                  colData2   = "treatment",
                  formulaDE2 = "~ guide_KO + treatment",
                  ref        = ref,
                  name       = name,
                  anno       = anno)

  saveRDS(qc, rds_name)

  norm_counts <- counts(qc[[1]], normalized = TRUE)
  write.csv(norm_counts, norm_csv_name, quote = FALSE)

  invisible(qc)
}

# ==========================================================
# Run 1: 36 samples from GSE302655
# treatment: v = vehicle, F = IFNg; 3 biological replicates per condition
# sg1.CREB1; sg1.CRTC2; sg2.CREB1; sg2.CRTC2; sg1.LKB1; sg1.NT
# ==========================================================
raw1  <- processRawTxt(here("data", "run1", "raw.txt"))
anno1 <- processAnnoTxt(here("data", "run1", "raw.txt"))

# Sample metadata is parsed from the column names: <sample>.<replicate><treatment>
rep_treat1 <- sub(".*\\.", "", colnames(raw1))     # e.g. "1F"
sample1    <- sub("\\.[^.]*$", "", colnames(raw1)) # e.g. "sg1.CREB1", "sgLKB1", "sgNT1"

# KO and guide. LKB1 and NT samples carry no guide number in their names and are guide 1.
KO1 <- ifelse(grepl("^sgNT", sample1),  "NT",
              ifelse(grepl("^sgLKB1", sample1), "LKB1",
                     sub("^sg[12]\\.", "", sample1)))
guide1 <- ifelse(grepl("^sg[12]\\.", sample1),
                 sub("^(sg[12])\\..*", "\\1", sample1),
                 "sg1")

design1 <- data.frame(
  row.names = colnames(raw1),
  guide     = guide1,
  KO        = KO1,
  rep_treat = rep_treat1,
  treatment = substring(rep_treat1, 2, 2),
  Replicate = factor(substring(rep_treat1, 1, 1))
)

design1$guide_KO       <- paste(design1$guide, design1$KO, sep = ".")
design1$guide_KO_treat <- sprintf("%s_%s", design1$guide_KO, design1$treatment)
design1$KO_treat       <- sprintf("%s_%s", design1$KO, design1$treatment)

qc1 <- run_qc(raw1, design1, anno1,
              ref           = "sg1.NT",
              name          = "run1_QC",
              rds_name      = "run1_QC.RDS",
              norm_csv_name = "run1_normalized_counts.csv")

# ==========================================================
# Run 2: 24 samples from GSE302655
# treatment: v = vehicle, F = IFNg; 3 biological replicates per condition
# sg1.Nr4a1; sg2.Nr4a1; sg1.LKB1; sg1.NT1
# ==========================================================
raw2  <- processRawTxt(here("data", "run2", "raw.txt"))
anno2 <- processAnnoTxt(here("data", "run2", "raw.txt"))

# Sample metadata is parsed from the column names: <guide>.<KO>.<replicate><treatment>
parts <- strsplit(colnames(raw2), "[.]+")
field <- function(i) vapply(parts, `[`, i, FUN.VALUE = character(1))

design2 <- data.frame(
  row.names = colnames(raw2),
  guide     = field(1),
  KO        = field(2),
  Replicate = substring(field(3), 1, 1),
  treatment = substring(field(3), 2, 2),
  rep_treat = field(3)
)

design2$guide_KO       <- paste(design2$guide, design2$KO, sep = ".")
design2$guide_KO_treat <- paste(design2$guide, design2$KO, design2$treatment, sep = ".")
design2$KO_treat       <- paste(design2$KO, design2$treatment, sep = ".")

# Sanity check: 4 groups x 2 treatments x 3 replicates
stopifnot(all(design2$treatment %in% c("v", "F")),
          all(table(design2$guide_KO, design2$treatment) == 3))

qc2 <- run_qc(raw2, design2, anno2,
              ref           = "sg1.NT1",
              name          = "run2_QC",
              rds_name      = "run2_QC.RDS",
              norm_csv_name = "run2_normalized_counts.csv")

# ==========================================================
# Run 1: differential expression -> results/DE/run1_results.RDS
# Uses raw1, design1 and anno1 from above.
# The ORDER of this list matters: downstream scripts index the results by
# position (e.g. results1[[1]] = comp1, the IFNg response in sgNT).
# Each entry: c(group 1, group 2, variable tested, reference level,
#               comparison id, design column used to select samples)
# ==========================================================
de_dir1 <- here("results", "DE", "run1")
dir.create(de_dir1, recursive = TRUE, showWarnings = FALSE)

comps1 <- list(
  # sgNT: vehicle vs IFNg
  c("sg1.NT_v", "sg1.NT_F", "treatment", "v", "comp1", "guide_KO_treat"),
  # IFNg-treated: sgNT vs sgLKB1
  c("sg1.NT_F", "sg1.LKB1_F", "guide_KO", "sg1.NT", "comp2", "guide_KO_treat"),
  # IFNg-treated: sgNT vs sg2.CREB1
  c("sg1.NT_F", "sg2.CREB1_F", "guide_KO", "sg1.NT", "comp3b", "guide_KO_treat"),
  # IFNg-treated: sgLKB1 vs sg2.CREB1
  c("sg1.LKB1_F", "sg2.CREB1_F", "guide_KO", "sg1.LKB1", "comp5b", "guide_KO_treat"),
  # IFNg-treated: NT vs CREB1 (both guides)
  c("NT_F", "CREB1_F", "KO_treat", "NT_F", "comp3a", "KO_treat"),
  # IFNg-treated: NT vs CRTC2 (both guides); reference = CRTC2
  c("NT_F", "CRTC2_F", "KO_treat", "CRTC2_F", "comp4", "KO_treat"),
  # IFNg-treated: LKB1 vs CREB1 (both guides)
  c("LKB1_F", "CREB1_F", "KO_treat", "LKB1_F", "comp5a", "KO_treat"),
  # IFNg-treated: LKB1 vs CRTC2 (both guides)
  c("LKB1_F", "CRTC2_F", "KO_treat", "LKB1_F", "comp6", "KO_treat")
)

results1 <- lapply(comps1, function(X) {
  d <- design1[design1[[X[6]]] %in% c(X[1], X[2]), ]
  r <- raw1[, row.names(d)]

  comp_dir <- file.path(de_dir1, paste0(X[5], "_", X[1], "v", X[2]))
  dir.create(comp_dir, showWarnings = FALSE)
  old_wd <- setwd(comp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  myDiffExp(counts      = r,
            design      = d,
            colData1    = X[3],
            colData2    = "Replicate",
            formulaEdgR = NULL,
            formulaDE2  = sprintf("~ Replicate + %s", X[3]),
            ref         = X[4],
            name        = paste0(X[5], "_", X[1], ".v.", X[2], ".deseq2"),
            org         = "mmusculus",
            anno        = anno1,
            db          = "geneontology_Biological_Process",
            runWebG     = FALSE,   # GO enrichment not needed for the results table
            fdr         = 0.05)
})
names(results1) <- vapply(comps1, `[`, 5, FUN.VALUE = character(1))

saveRDS(results1, here("results", "DE", "run1_results.RDS"))

# ==========================================================
# Run 2: differential expression -> results/DE/run2_results.RDS
# Uses raw2, design2 and anno2 from above.
# Each entry: c(group 1, group 2, variable tested, reference level,
#               comparison id, design column used to select samples)
# ==========================================================
de_dir <- here("results", "DE", "run2")
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

comps2 <- list(
  # Ctrl-vehicle vs Ctrl-IFNg
  c("NT1.v", "NT1.F", "treatment", "v", "comp1", "KO_treat"),
  # LKB1-vehicle vs LKB1-IFNg
  c("LKB1.v", "LKB1.F", "treatment", "v", "comp2", "KO_treat"),
  # Nr4a1-vehicle vs Nr4a1-IFNg (both guides, then each guide)
  c("Nr4a1.v", "Nr4a1.F", "treatment", "v", "comp3", "KO_treat"),
  c("sg1.Nr4a1.v", "sg1.Nr4a1.F", "treatment", "v", "comp3a", "guide_KO_treat"),
  c("sg2.Nr4a1.v", "sg2.Nr4a1.F", "treatment", "v", "comp3b", "guide_KO_treat"),
  # Ctrl-IFNg vs LKB1-IFNg
  c("NT1.F", "LKB1.F", "KO_treat", "NT1.F", "comp4", "KO_treat"),
  # Ctrl-IFNg vs Nr4a1-IFNg
  c("NT1.F", "Nr4a1.F", "KO_treat", "NT1.F", "comp5", "KO_treat"),
  c("sg1.NT1.F", "sg1.Nr4a1.F", "KO_treat", "NT1.F", "comp5a", "guide_KO_treat"),
  c("sg1.NT1.F", "sg2.Nr4a1.F", "KO_treat", "NT1.F", "comp5b", "guide_KO_treat"),
  # LKB1-IFNg vs Nr4a1-IFNg
  c("LKB1.F", "Nr4a1.F", "KO_treat", "LKB1.F", "comp6", "KO_treat"),
  c("sg1.LKB1.F", "sg1.Nr4a1.F", "KO_treat", "LKB1.F", "comp6a", "guide_KO_treat"),
  c("sg1.LKB1.F", "sg2.Nr4a1.F", "KO_treat", "LKB1.F", "comp6b", "guide_KO_treat"),
  # Ctrl-vehicle vs LKB1-vehicle
  c("NT1.v", "LKB1.v", "KO_treat", "NT1.v", "comp7", "KO_treat"),
  # Ctrl-vehicle vs Nr4a1-vehicle
  c("NT1.v", "Nr4a1.v", "KO_treat", "NT1.v", "comp8", "KO_treat"),
  c("sg1.NT1.v", "sg1.Nr4a1.v", "KO_treat", "NT1.v", "comp8a", "guide_KO_treat"),
  c("sg1.NT1.v", "sg2.Nr4a1.v", "KO_treat", "NT1.v", "comp8b", "guide_KO_treat"),
  # LKB1-vehicle vs Nr4a1-vehicle
  c("LKB1.v", "Nr4a1.v", "KO_treat", "LKB1.v", "comp9", "KO_treat"),
  c("sg1.LKB1.v", "sg1.Nr4a1.v", "KO_treat", "LKB1.v", "comp9a", "guide_KO_treat"),
  c("sg1.LKB1.v", "sg2.Nr4a1.v", "KO_treat", "LKB1.v", "comp9b", "guide_KO_treat")
)

# myDiffExp() writes its output files to the working directory, so each
# comparison runs inside its own folder and the old directory is restored after.
results2 <- lapply(comps2, function(X) {
  d <- design2[design2[[X[6]]] %in% c(X[1], X[2]), ]
  r <- raw2[, row.names(d)]

  comp_dir <- file.path(de_dir, paste0(X[5], "_", X[1], "v", X[2]))
  dir.create(comp_dir, showWarnings = FALSE)
  old_wd <- setwd(comp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  myDiffExp(counts      = r,
            design      = d,
            colData1    = X[3],
            colData2    = "Replicate",
            formulaEdgR = NULL,
            formulaDE2  = sprintf("~ Replicate + %s", X[3]),
            ref         = X[4],
            name        = paste0(X[5], "_", X[1], ".v.", X[2], ".deseq2"),
            org         = "mmusculus",
            anno        = anno2,
            db          = "geneontology_Biological_Process",
            runWebG     = FALSE,   # GO enrichment not needed for the results table
            fdr         = 0.05)
})
names(results2) <- vapply(comps2, `[`, 5, FUN.VALUE = character(1))

saveRDS(results2, here("results", "DE", "run2_results.RDS"))

# ---- Record package versions ----
writeLines(capture.output(sessionInfo()), here("sessionInfo.txt"))
