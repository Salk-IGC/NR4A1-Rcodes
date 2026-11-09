# NR4A1 directly represses Interferon Response Genes and mediates checkpoint immunotherapy resistance

R code for the RNA-seq and CUT&RUN analyses and figure panels in:

> * NR4A1 directly represses Interferon Response Genes and mediates checkpoint immunotherapy resistance (submitted) *

## Data availability

Raw and processed data are deposited in GEO:

| Data | Accession |
|------|-----------|
| LKR10 RNA-seq | GSE302655 & GSE351227 |
| NR4A1 CUT&RUN | GSE302652 |

## Repository layout

```
R/
  rnaseq_functions.R     RNA-seq helper functions (processRawTxt, processAnnoTxt, myExplore, myDiffExp)
  plot_helpers.R         save_pheatmap_jpeg(), gVolcano4review()
  cutandrun_helpers.R    countmatrix(), plotpairheatmap2()
scripts/
  01_RNAseq_preprocess.R QC, DESeq2 objects, normalized counts, DE results (run first)
  Fig1A.R                Figure 1A
  Fig2D.R                Figure 2D
  Fig3ACG.R              Figure 3A, 3C, 3G
  Fig4ABCD.R             Figure 4A-D (CUT&RUN + RNA-seq)
  FigS1.R                Figure S1A, S1B
  FigS5.R                Figure S5E, S5G
  FigS6.R                Figure S6G
  FigS8.R                Figure S8A (CUT&RUN)
data/                    input files (not tracked here; see below)
results/                 written by 01_RNAseq_preprocess.R
figures/                 written by the figure scripts
```

## How to run

1. Open the repository root as an R project (paths are built with `here::here()`).
2. Put the input files in `data/` (see below).
3. Run `scripts/01_RNAseq_preprocess.R`. This creates `results/QC/run1_QC.RDS`, `results/QC/run2_QC.RDS`,
   normalized count tables, and `results/DE/run1_results.RDS` / `run2_results.RDS`.
4. Run any figure script. Each script header lists its inputs and outputs.

| Script | Needs |
|--------|-------|
| Fig1A | `results/QC/run2_QC.RDS`, `data/tableS1.txt` |
| Fig2D | `results/QC/run1_QC.RDS`, `data/tableS2.txt` |
| Fig3ACG | `run1_QC.RDS`, `run2_QC.RDS`, `data/tableS2.txt` |
| Fig4ABCD | CUT&RUN peak files and bigwigs (GSE302652), `mm10.tss`, `tableS3.txt`, `run2_QC.RDS` |
| FigS1 | `run2_QC.RDS`, `run2_results.RDS`, `data/tableS1.txt` |
| FigS5 | `run1_QC.RDS`, `run2_results.RDS` |
| FigS6 | `run2_QC.RDS`, gene list from Table S2 |
| FigS8 | CUT&RUN peak files and bigwigs (GSE302652) |

Some scripts (Fig4ABCD, FigS1, FigS5, FigS6, FigS8) read files by name from the working directory;
set the working directory to the folder holding those files, or edit the paths at the top of the script.

### Input data

- **RNA-seq raw counts:** `data/run1/raw.txt` and `data/run2/raw.txt`, generated from GSE302655.
  Sample metadata is parsed from the column names, format `<sample>.<replicate><treatment>`
  (e.g. `sg1.CREB1.1F`; `v` = vehicle, `F` = IFNg).
- **Gene lists:** `tableS1.txt` (IFNg-dependent genes, n = 286), `tableS2.txt` (LKB1-dependent genes, n = 166),
  and `tableS3.txt`, from the paper's supplementary tables.
- **Experiments:** run 1 has sgNT, sgLkb1, sgCrtc2 and sgCreb1 cells (36 samples); run 2 has sgNT, sgLkb1 and
  sgNr4a1 cells (24 samples). All were treated with vehicle or IFNg (16 h, 25 ng/mL).

## Requirements

- R 4.3.2 or later
- CRAN: `here`, `ggplot2`, `pheatmap`, `data.table`, `qdapRegex`, `stringr`, `gplots`, `cluster`,
  `VennDiagram`, `openxlsx`, `ggrepel`
- Bioconductor: `DESeq2`, `EnrichedHeatmap`, `ComplexHeatmap`, `rtracklayer`, `GenomicRanges`, `EnhancedVolcano`, `circlize`

`01_RNAseq_preprocess.R` writes the package versions it used to `sessionInfo.txt`.
