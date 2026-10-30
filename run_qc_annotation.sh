#!/usr/bin/env bash
# Reproducible genome QC + annotation pipeline for 10 Salmonella enterica assemblies.
# Part of the bacterial-genome-qc-annotation portfolio project.
#
# Prerequisites (see README.md):
#   - NCBI Datasets CLI (for data/download_assemblies.sh)
#   - CheckM2 with its database downloaded (run on the 10 assemblies)
#   - Prokka 1.15.6+ with databases set up (prokka --setupdb)
#
# Usage: bash run_qc_annotation.sh
# Outputs: results/quality_report.tsv (from CheckM2), prokka_out/<genome>/ per genome.
set -euo pipefail

ASSEMBLIES="data/assemblies"
RESULTS="results"

# 1. Download assemblies (skip if already present)
if [ ! -f "$ASSEMBLIES/Typhimurium_LT2.fna" ]; then
  echo "Downloading 10 Salmonella assemblies from NCBI..."
  bash data/download_assemblies.sh
else
  echo "Assemblies already present, skipping download."
fi

# 2. CheckM2 quality assessment.
#    The quality_report.tsv bundled in results/ was produced on the author's PC with:
#      checkm2 predict --threads 8 --input "$ASSEMBLIES" --output-directory checkm2_out
#    (CheckM2 version and database date are recorded in README.md.)
#    To reproduce here, uncomment:
# mkdir -p "$RESULTS"
# checkm2 predict --threads 8 --input "$ASSEMBLIES" --output-directory checkm2_out
# cp checkm2_out/quality_report.tsv "$RESULTS/quality_report.tsv"

# 3. Prokka annotation, one run per genome.
mkdir -p prokka_out
for fna in "$ASSEMBLIES"/*.fna; do
  label=$(basename "$fna" .fna)
  if [ -f "prokka_out/$label/$label.txt" ]; then
    echo "Prokka output for $label already present, skipping."
    continue
  fi
  echo "Annotating $label with Prokka..."
  prokka --outdir "prokka_out/$label" --prefix "$label" \
    --genus Salmonella --species enterica \
    --cpus 4 --quiet "$fna"
done

# 4. Summarize annotation stats and build figures.
python3 qc_analysis.py

echo "Done. See results/genome_qc_summary.csv and figures/."
