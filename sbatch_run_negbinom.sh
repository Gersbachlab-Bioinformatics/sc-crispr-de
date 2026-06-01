#!/bin/bash
# Submit a SLURM array job to run negative binomial DE analysis across all sgRNAs.
# Edit the variables below before running.

# --- User-configurable paths ---
LOGS_DIR="/path/to/logs"
SGRNA_DIR="/path/to/sgrna/data"          # folder(s) with matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz
TRANSCRIPT_DIR="/path/to/transcriptomic/data"
GENE_ANNOTATIONS="/path/to/ensembl.hg38.annotations.bed"
GENE_WHITELIST="/path/to/genes.txt"
RDATA="/path/to/preprocessed.RData"       # set to the checkpoint saved by --save-rdata, or remove flag below
OUTPUT_DIR="/path/to/output"
OUTPUT_BASENAME="${OUTPUT_DIR}/sample.negbinom"
SCRIPT_DIR="$(dirname "$(realpath "$0")")"

# --- SLURM parameters ---
N_SGRNAS=12625   # total number of sgRNAs (sets the array upper bound)

mkdir -p "${LOGS_DIR}"
mkdir -p "${OUTPUT_DIR}"

sbatch \
    -o "${LOGS_DIR}/run_negbinom.%a.out" \
    --array=1-${N_SGRNAS}%500 \
    --mem=8G \
    --cpus-per-task=1 \
    <<EOF
#!/bin/bash
Rscript "${SCRIPT_DIR}/run_negbinom.extract_summary.R" \
    --job \${SLURM_ARRAY_TASK_ID} \
    --sgrnas "${SGRNA_DIR}" \
    --transcript-data "${TRANSCRIPT_DIR}" \
    --de-method negbinom \
    --gene-annotations-file "${GENE_ANNOTATIONS}" \
    --gene-symbol-whitelist "${GENE_WHITELIST}" \
    --sctransform \
    --min-nfeature-rna 20 \
    --max-ncount-rna 3000 \
    --max-mt-percent 20 \
    --rdata "${RDATA}" \
    --min-cells-per-gene-and-group 3 \
    --sgrnas-binary-mixture-model \
    --output-basename "${OUTPUT_BASENAME}"
EOF
