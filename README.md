# BAGFOOT2 — differential TF footprinting from ATAC-seq

**BAGFOOT2** is a computational framework to infer **transcription factor (TF) activity changes** between two experimental conditions using **ATAC-seq** data.  
It integrates **sequence-bias–corrected footprinting**, **chromatin accessibility**, **robust outlier statistics**, and **motif enrichment** to provide a robust, interpretable view of TF regulatory dynamics.

BAGFOOT2 is an updated implementation of BaGFoot (Baek et al., *Cell Reports* 2017), designed for **biologists and bioinformaticians** interested in chromatin accessibility, transcriptional regulation, and regulatory genomics.

---

## Abstract
Transcription factors (TFs) are central mediators of gene regulation, interpreting upstream signals into context-specific transcriptional programs.
Chromatin accessibility profiling by ATAC-seq provides a powerful, genome-wide readout of regulatory activity, yet robust inference of TF activity from ATAC-seq remains challenging due to motif redundancy,
limited comparability across conditions, and reliance on a single analytical approach. Here, we present BAGFOOT2,
a computational framework that infers differential TF activity directly from ATAC-seq data in an unbiased and multilayered manner.
BAGFOOT2 quantifies three complementary indicators of differential TF activity for each motif: footprint depth, flanking chromatin accessibility,
and motif enrichment. BAGFOOT2 then compares these indicators between two experimental conditions, identifies TFs whose changes stand out using a robust Mahalanobis distance (Minimum Covariance Determinant), and integrates them into intuitive plots and comprehensive result tables.
To reduce redundancy and improve interpretability, BAGFOOT2 operates on archetypal consensus motif clusters rather than individual, highly similar motifs.
BAGFOOT2 requires only ATAC-seq data, does not depend on prior knowledge of relevant TFs, and does not necessitate focusing exclusively on differentially accessible regions.

---

## Documentation

📘 **HTML vignette (recommended entry point)**  
➡️ https://goldsteinlab-tech.github.io/BAGFOOT2/

The vignette provides:
- a conceptual overview of the BAGFOOT2 workflow,
- detailed explanations of each processing step,
- the command-line usage and all options,
- annotated example R code to run each step,
- guidance on project structure and inputs,
- full description of footprinting, statistics, and output files.

## Motif files (Zenodo)

The genome-wide archetype motif occurrences (**MotifDB**) and hexamer index (**nuccode**) files for **mm10** and **hg38** are available on Zenodo:

**Zenodo record:**  
https://doi.org/10.5281/zenodo.20782782

They are **downloaded automatically** on first use if `--mm_files_dir` is not given. On a shared server, set `BAGFOOT2_MOTIF_DB_DIR` to a folder that already contains them.

---

## Installation

### Conda (recommended)
BAGFOOT2 and its dependencies can be installed via **Conda**:

```bash
conda create -n bagfoot2 -c goldsteinlab-tech -c conda-forge -c bioconda bagfoot2
conda activate bagfoot2
bagfoot2 --help
```

> Make sure `samtools` is available, and `HOMER` if you plan to use motif enrichment (`--homer TRUE`).

### R (from source)

```r
# install devtools if needed
install.packages("devtools")

# install BAGFOOT2 from GitHub
devtools::install_github("goldsteinlab-tech/BAGFOOT2", build_vignettes = TRUE)
```

---

## Quick start

```bash
bagfoot2 \
  --bam_dir    bam_merged/ \
  --peak_dir   macs2/ \
  --output_dir results_fed_fasted/ \
  --control    fed \
  --treatment  fasted \
  --genome     mm10 \
  --cores      4
```

Add `--homer TRUE --homer_bin /path/homer/bin --reps_control 1,2,3 --reps_treatment 1,2,3` to include HOMER motif enrichment.

---

## Input Requirements

BAGFOOT2 operates on ATAC-seq data and requires:

* **Merged BAM files** (one per condition; sorted and indexed): `<condition>_merged.bam`,
* **MACS2 peak files** (BED or narrowPeak): `<condition>_peaks.narrowPeak`,
* Reference genome (`mm10` or `hg38`),
* Optional but recommended: replicate-level BAMs for **HOMER** integration (`<condition>_1.bam`, `<condition>_2.bam`, …; the number of replicates can differ between conditions).

Details on file organization and naming conventions are provided in the vignette.

---

## Key Features

* Bias-corrected ATAC-seq footprinting (hexamer-based correction)
* Explicit modeling of **footprint depth vs. flanking accessibility**
* Archetype motif clusters to reduce motif redundancy
* **Robust Mahalanobis distance (MCD)** with finite-sample p-values and FDR correction
* One-line **command-line** interface, automatic download of motif files
* Batch-wise motif processing with restart capability
* Integration with **HOMER** differential accessibility and motif enrichment
* Publication-ready **PDF / TIFF** figures and a complete **CSV** result table

---

## Output

BAGFOOT2 produces:

* A bagplot of all motifs (Δ flanking accessibility vs. Δ footprint depth) highlighting TFs with differential activity
* Barplots of the top significant TFs ranked by Mahalanobis distance (D2) and by FDR
* A table with all values per motif: ΔFA, ΔFPD, D2, p-value, FDR, direction, HOMER enrichment
* Optional annotation with HOMER motif enrichment results

---

## Citation

If you use BAGFOOT2 in your work, please cite the original BaGFoot method:

> Baek S., Goldstein I., Hager G.L.
> **Bivariate genomic footprinting detects changes in transcription factor activity.**
> *Cell Reports* 19(8):1710–1722 (2017)

A BAGFOOT2 manuscript is in preparation.

---

## License

This project is licensed under the **GNU General Public License (GPL)**.
See the `LICENSE` file for details.

---

## Acknowledgments

This work was supported by the **European Research Council** (ERC-StG #947907, to I.G.).

---

## Contact & Issues

* Maintainer: [ido.goldstein@mail.huji.ac.il](mailto:ido.goldstein@mail.huji.ac.il)
* Bug reports & feature requests:
  [https://github.com/goldsteinlab-tech/BAGFOOT2/issues](https://github.com/goldsteinlab-tech/BAGFOOT2/issues)
