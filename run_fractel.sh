#!/usr/bin/env bash
# run_fractel.sh — Merge negbinom DE outputs and run FRACTEL rank aggregation.
#
# Prerequisite: install FRACTEL
#   pip install --force --no-deps git+https://github.com/Gersbachlab-Bioinformatics/FRACTEL
#
# Run after all SLURM array jobs from sbatch_run_negbinom.sh have completed.
# Edit the user-configurable section below, then:
#   bash run_fractel.sh

set -euo pipefail

# --- User-configurable ---
NB_OUTPUT_BASENAME="/path/to/output/sample.negbinom"  # --output-basename used in run_negbinom.extract_summary.R
OUTPUT_DIR="/path/to/fractel/output"
OUTPUT_BASENAME="${OUTPUT_DIR}/sample.fractel"
NUM_SIMULATIONS=100000

# BND: top-k guides used per element in the RRA test.
#   Float in (0,1) → fraction of guides per element (e.g. 0.5 = top 50 %)
#   Integer >= 1   → absolute number of guides
BND=0.5

# Value in the 'gene' column that identifies non-targeting controls.
# Leave empty ("") if no non-targeting controls are present.
NTARGETING_KEYWORD=""
# -------------------------

mkdir -p "${OUTPUT_DIR}"

MERGED_TSV="${OUTPUT_BASENAME}.negbinom_merged.tsv"
SIM_BASENAME="${OUTPUT_BASENAME}.sim"

# ---------------------------------------------------------------------------
# Step 1: merge negbinom per-guide output files and filter to the perturbation
#         effect row (negbinom_component == "groupGroup2").
# ---------------------------------------------------------------------------
echo "[1/3] Merging negbinom outputs and filtering to groupGroup2..."
python3 - <<PYEOF
import glob, sys
import pandas as pd

pattern = "${NB_OUTPUT_BASENAME}.*.txt"
files = sorted(glob.glob(pattern))
if not files:
    sys.exit(f"No files found matching: {pattern}")

df = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
df = df[df["negbinom_component"] == "groupGroup2"].copy()
df = df.rename(columns={"Pr(>|z|)": "pvalue"})
df.to_csv("${MERGED_TSV}", sep="\t", index=False)
print(f"Wrote {len(df):,} rows  ({df['grna'].nunique():,} gRNAs, {df['gene'].nunique():,} genes)  →  ${MERGED_TSV}")
PYEOF

# ---------------------------------------------------------------------------
# Step 2: compute unique per-element guide counts (needed to pre-simulate the
#         null distribution).  Non-targeting controls are excluded because they
#         form the background set and are not aggregated by element.
# ---------------------------------------------------------------------------
echo "[2/3] Running fractel simulate..."
NUM_GUIDES_ARG=$(python3 - <<PYEOF
import pandas as pd
df = pd.read_csv("${MERGED_TSV}", sep="\t")
keyword = "${NTARGETING_KEYWORD}"
if keyword:
    df = df[df["gene"] != keyword]
counts = sorted(df.groupby("gene")["grna"].count().unique().tolist())
print(" ".join(str(c) for c in counts))
PYEOF
)
echo "  Unique guide counts per element: ${NUM_GUIDES_ARG}"

fractel simulate \
    --num-guides ${NUM_GUIDES_ARG} \
    --num-simulations "${NUM_SIMULATIONS}" \
    --bnd "${BND}" \
    --output-basename "${SIM_BASENAME}"

# ---------------------------------------------------------------------------
# Step 3: run FRACTEL element-level rank aggregation.
# ---------------------------------------------------------------------------
echo "[3/3] Running fractel run..."
FRACTEL_ARGS=(
    fractel run
    --data-frame "${MERGED_TSV}"
    --aggregating-cols "gene"
    --pval-col "pvalue"
    --row-id-col "grna"
    --sim-data "${SIM_BASENAME}.npz"
    --bnd "${BND}"
    --effect-size-col "Estimate"
    --pval-type "two-sided"
    --output-basename "${OUTPUT_BASENAME}"
)
if [[ -n "${NTARGETING_KEYWORD}" ]]; then
    FRACTEL_ARGS+=(--keyword-for-background-values "${NTARGETING_KEYWORD}")
fi
"${FRACTEL_ARGS[@]}"

echo "Done.  Results → ${OUTPUT_BASENAME}.tsv.gz"
echo "  Columns: gene | FRACTEL_pval | FRACTEL_pval_fdr_corr | FRACTEL_effect_size"
