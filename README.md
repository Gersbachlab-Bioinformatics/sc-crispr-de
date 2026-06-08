# sc-CRISPR Differential Expression (Negative Binomial)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Perform per-guide differential expression analysis on single-cell CRISPR screen data using a negative binomial GLM. Each sgRNA is tested in a SLURM array job: cells carrying that guide are compared against cells with no detected guide.

## Overview

**Input**
- 10X Genomics output folders for sgRNA counts and gene expression counts
- A gene annotation BED6 file (or explicit target region coordinates) to define which genes to test
- Optionally, a pre-built RData checkpoint to skip data loading

**Output**
One tab-separated file per array job (`<output-basename>.<job>.txt`) with NB model coefficients and p-values for each tested gene.

## Environment setup

### 1. Install Mamba (if not already available)

```bash
conda install -n base -c conda-forge mamba
```

### 2. Create the environment

```bash
mamba create -n sc_crispr_de \
    -c conda-forge -c bioconda \
    r-base \
    r-argparse \
    r-seurat \
    bioconductor-genomicranges \
    r-stringr \
    r-purrr \
    r-mass \
    r-data.table \
    r-matrix
```

Activate it:

```bash
conda activate sc_crispr_de
```

### 3. Optional: bias-reduced NB model (`--apply-bias-reduction`)

```bash
mamba install -n sc_crispr_de -c conda-forge r-brglm2
```

If `r-brglm2` is unavailable via conda, install it from within R:

```r
install.packages("brglm2")
```

## Usage

### Step 1 — (Optional) build an RData checkpoint

Run the script once on a single job to load, filter, and save the processed Seurat object. This avoids re-loading data for every array task:

```bash
Rscript run_negbinom.extract_summary.R \
    --job 1 \
    --sgrnas /path/to/sgrna/data \
    --transcript-data /path/to/gex/data \
    --target-region-chrom chr6 \
    --target-region-start 29700000 \
    --target-region-end 33200000 \
    --gene-annotations-file /path/to/ensembl.hg38.annotations.bed \
    --sctransform \
    --min-nfeature-rna 20 \
    --max-ncount-rna 3000 \
    --max-mt-percent 20 \
    --save-rdata /path/to/output/checkpoint.RData \
    --output-basename /path/to/output/sample.negbinom
```

### Step 2 — Submit the SLURM array

Edit the path variables at the top of `sbatch_run_negbinom.sh`, then:

```bash
bash sbatch_run_negbinom.sh
```

### Step 3 — Run FRACTEL rank aggregation

After all array jobs have completed, run [FRACTEL](https://github.com/Gersbachlab-Bioinformatics/FRACTEL) to aggregate per-guide p-values into element-level significance scores.

**Install FRACTEL** (once):

```bash
pip install --force --no-deps git+https://github.com/Gersbachlab-Bioinformatics/FRACTEL
```

**Configure and run** `run_fractel.sh`:

Edit the variables at the top of `run_fractel.sh` (mirror the paths used in `sbatch_run_negbinom.sh`), then:

```bash
bash run_fractel.sh
```

The script performs three steps automatically:

1. **Merge outputs** — concatenates all `<output-basename>.<job>.txt` files, filters to the `groupGroup2` negbinom component (the perturbation effect row), and renames `Pr(>|z|)` → `pvalue`.
2. **Simulate null** — runs `fractel simulate` to build the RRA null distribution for each unique guide-per-gene count observed in the data.
3. **Run FRACTEL** — runs `fractel run` grouping guides by `gene`, using `grna` as the row identifier and the NB `Estimate` coefficient as the effect size.

**Output** — `<OUTPUT_BASENAME>.tsv.gz`, a compressed TSV with one row per gene and columns:

| Column | Description |
|---|---|
| `gene` | Gene symbol (index) |
| `FRACTEL_pval` | RRA element-level p-value |
| `FRACTEL_pval_fdr_corr` | BH-corrected FDR across all genes |
| `FRACTEL_effect_size` | Weighted-average NB estimate across top guides |

**Optional: calibrate p-values** before running FRACTEL when non-targeting control p-values are not uniformly distributed:

```bash
fractel calibrate \
    --data-frame <OUTPUT_BASENAME>.negbinom_merged.tsv \
    --reference-data-frame <OUTPUT_BASENAME>.negbinom_merged.tsv \
    --reference-df-select-col gene \
    --reference-df-select-value <NTARGETING_KEYWORD> \
    --pval-col pvalue \
    --interpolated-col pvalue_calibrated \
    --output-basename <OUTPUT_BASENAME>.negbinom_calibrated
```

Then re-run Step 3 of `run_fractel.sh` with `--data-frame` pointing to the calibrated file and `--pval-col pvalue_calibrated`.

### Key arguments

| Argument | Description |
|---|---|
| `--job` | Index into the gRNA list (set by `$SLURM_ARRAY_TASK_ID`) |
| `--sgrnas` | 10X folder(s) for gRNA counts |
| `--transcript-data` | 10X folder(s) for gene expression counts |
| `--target-region-bed` | BED file defining genomic region(s) to test |
| `--gene-annotations-file` | BED6 gene annotations to intersect with target region |
| `--gene-symbol-whitelist` | Plain-text list of gene symbols to test (alternative to BED-based selection) |
| `--rdata` | Skip data loading by providing a saved checkpoint |
| `--save-rdata` | Save a processed Seurat object before DE analysis |
| `--sctransform` | Normalize with SCTransform instead of log-normalization |
| `--latentvar` | Covariate(s) to include in the NB model (e.g. `nCount_RNA`) |
| `--min-cells-per-gene-and-group` | Minimum expressing cells per group to test a gene |
| `--sgrnas-binary-mixture-model` | Input gRNA matrix is already a binary cell × guide assignment |
| `--apply-bias-reduction` | Use bias-reduced NB (`brnb`) instead of `glm.nb` |

Run `Rscript run_negbinom.extract_summary.R --help` for the full argument list.

## Dependencies

| Package | Source | Notes |
|---|---|---|
| argparse | CRAN (conda-forge: `r-argparse`) | CLI argument parsing |
| Seurat | CRAN (conda-forge: `r-seurat`) | Single-cell data structures and normalization |
| GenomicRanges | Bioconductor (bioconda: `bioconductor-genomicranges`) | Genomic interval operations |
| stringr | CRAN (conda-forge: `r-stringr`) | String utilities |
| purrr | CRAN (conda-forge: `r-purrr`) | Functional programming helpers |
| MASS | CRAN (conda-forge: `r-mass`) | `glm.nb` for negative binomial regression |
| data.table | CRAN (conda-forge: `r-data.table`) | Fast result aggregation |
| Matrix | CRAN (conda-forge: `r-matrix`) | Sparse matrix operations (Seurat dependency) |
| brglm2 | CRAN (conda-forge: `r-brglm2`) | **Optional** — required for `--apply-bias-reduction` |

## Versioning and releases

The current version is tracked in the [`VERSION`](VERSION) file. To tag and publish a release:

```bash
VERSION=$(cat VERSION)
git tag -a "v${VERSION}" -m "Release v${VERSION}"
git push origin "v${VERSION}"
```

Then create a GitHub release from that tag (via the web UI or `gh release create`).

To bump the version, edit `VERSION` before tagging.

## License

MIT — see [LICENSE](LICENSE).
