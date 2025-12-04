# OnenessTwoness
## BRCA1/BRCA2-del phenotype random forest classifier



## <font color=black> Introduction </font>
`onenesstwoness` is an R package that implements the **B1+2 classifier** ([Setton, Hadi, Choo et al., Nature 2023](https://www.nature.com/articles/s41586-023-06461-2)), a pan-cancer homologous recombination (HR) deficiency classifier that augments [HRDetect](https://www.nature.com/articles/nm.4292) with structural variant features specific to BRCA1-deficient and BRCA2-deficient tumors. The classifier combines six HRDetect features with five additional SV features (1–100-kb tandem duplications, reciprocal duplications, reciprocal deletion-duplications, reciprocal deletions, and homeologous deletions) to improve HR-deficiency classification, particularly for BRCA2-deficient cancers.


It provides:

- A wrapper to run the HRDetect pipeline from standard VCF/BEDPE/CNV inputs: `run_hrdetect`
- A classifier to compute Oneness/Twoness probabilities: `predict_B1_2`

The typical use case is:

1. Prepare a **complex events gGraph** (e.g. from [JaBbA](https://github.com/mskilab-org/JaBbA) / [gGnome](https://github.com/mskilab-org/gGnome)).
2. Prepare SNV/indel/SV/het pileup data and reference genome.
3. Run `predict_B1_2()` to compute Oneness/Twoness probabilities, optionally letting it call `run_hrdetect()` internally.

---

## <font color=black> Installation </font>

```r
# install.packages("devtools")  # if needed
devtools::install_github("mskilab-org/onenesstwoness-dev")
```

Required R/Bioconductor packages are installed automatically. You must also have the following command‑line tools on your `$PATH` for the HRDetect preprocessing:

`bcftools`, `vcftools`, `bedtools`, `vcf-sort`, `bgzip`, `tabix`

---

## <font color=black> Usage </font>

Minimal example (HRDetect already run):

```r
library(onenesstwoness)

ot <- predict_B1_2(
  complex          = "/path/to/complex_events_ggraph.rds",
  hrdetect_results = "/path/to/hrdetect_results.rds",
  outdir           = "./",
  save             = TRUE
)

```

Full pipeline (run homeology + HRDetect internally):

```r
library(onenesstwoness)

ot <- predict_B1_2(
  complex = "/path/to/complex_events_ggraph.rds",
  snv     = "/path/to/sample.snv.vcf.gz",
  indel   = "/path/to/sample.indel.vcf.gz",   # optional; defaults to snv
  sv      = "/path/to/sample.sv.vcf.gz",
  hets    = "/path/to/sample_hets.tsv",
  mask    = system.file("extdata", "hg19_mask.rds", package = "onenesstwoness"),
  genome  = "/path/to/human_g1k_v37.fasta",
  ref     = "hg19",
  outdir  = "./ot_output",
  cores   = 4, 
  save    = TRUE
)

```

### Key `predict_B1_2` Arguments
| Parameter | Default value | Description/notes |
|-----------|---------------|-------------------|
| `complex` | **required** | Path to complex events `gGraph` RDS (typically from JaBbA/gGnome). |
| `snv` | `NULL` | Somatic SNV VCF (gzipped or plain). Required if `hrdetect_results` is `NULL` and you want to run HRDetect. |
| `indel` | `snv` | Somatic indel VCF. If `NULL`, `snv` is used for both SNVs and indels. |
| `sv` | `NULL` | Structural variant calls (VCF or JaBbA-style RDS). Used for SV signatures in HRDetect. |
| `hets` | `NULL` | Heterozygous SNV pileup/counts table (for allelic CN). Optional but recommended. |
| `homeology` | `NULL` | Precomputed homeology results (list or path to RDS). If `NULL`, computed from `complex`. |
| `homeology_stats` | `NULL` | Precomputed homeology summary statistics. If `NULL`, computed from `homeology`/`complex`. |
| `hrdetect_results` | `NULL` | Path to precomputed HRDetect results RDS. If `NULL`, `run_hrdetect()` is called using `snv`/`indel`/`sv`/`hets`. |
| `genome` | `NULL` | Path to reference genome FASTA (e.g. `human_g1k_v37.fasta`). Required when running HRDetect. |
| `ref` | `"hg19"` | Genome ID string for HRDetect (e.g. `"hg19"`, `"hg38"`). Must be consistent with `genome`. |
| `mask` | `system.file("extdata","hg19_mask.rds", ...)` | RDS or BED of regions to include/exclude in HRDetect preprocessing. |
| `model` | `system.file("model","stash.retrained.model.rds", ...)` | Path to random forest classifier model. Defaults to packaged model. |
| `outdir` | `"./"` | Directory where intermediate files and final results are written. |
| `save` | `TRUE` | If `TRUE`, saves `onenesstwoness_results.rds` in `outdir`. |
| `cores` | `1` | Number of CPU cores for parallelizable steps (e.g. HRDetect preprocessing). |

---

## <font color=black> Usage for HRDetect only </font>

```r
library(onenesstwoness)

res <- run_hrdetect(
  snv    = "/path/to/sample.snv.vcf.gz",
  hets   = "/path/to/sample_hets.tsv",
  jabba  = "/path/to/jabba_simple.rds",   # JaBbA result with CN segments
  indel  = "/path/to/sample.indel.vcf.gz",
  sv     = "/path/to/sample.sv.vcf.gz",
  mask   = system.file("extdata", "hg19_mask.rds", package = "onenesstwoness"),
  genome = "/path/to/human_g1k_v37.fasta",
  ref    = "hg19",
  outdir = "./hrdetect_output",
  save   = TRUE
)

```

## <font color=black> Notes on reproducibility </font>

The pipeline writes multiple temporary files in the working directory (e.g. processed VCFs, BEDPE, CNV tables, region masks).
For reproducibility, use a dedicated directory per sample and set <code>outdir</code> explicitly.
Ensure that <code>genome</code> (FASTA) and <code>ref</code> (e.g. <code>"hg19"</code>) are consistent.

---

## License

[MIT](https://choosealicense.com/licenses/mit/)
````

