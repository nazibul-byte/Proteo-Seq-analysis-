#!/usr/bin/env bash
set -euo pipefail

# Proteo-Seq / snATAC-seq alignment wrapper
#
# Usage:
#   ./run_nfcore_atacseq.sh DESIGN.csv OUTDIR
#
# Optional environment variables:
#   NFCORE_ATACSEQ=nf-core/atacseq
#   NFCORE_ATACSEQ_REVISION=2.0
#   GENOME=hg38
#   READ_LENGTH=75
#   MAX_CPUS=1
#   MAX_MEMORY=10.GB
#   PROFILE=singularity

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 DESIGN.csv OUTDIR" >&2
  exit 1
fi

DESIGN_CSV="$1"
OUTDIR="$2"

NFCORE_ATACSEQ="${NFCORE_ATACSEQ:-nf-core/atacseq}"
NFCORE_ATACSEQ_REVISION="${NFCORE_ATACSEQ_REVISION:-2.0}"
GENOME="${GENOME:-hg38}"
READ_LENGTH="${READ_LENGTH:-75}"
MAX_CPUS="${MAX_CPUS:-1}"
MAX_MEMORY="${MAX_MEMORY:-10.GB}"
PROFILE="${PROFILE:-singularity}"

[[ -f "$DESIGN_CSV" ]] || { echo "ERROR: design CSV not found: $DESIGN_CSV" >&2; exit 1; }
mkdir -p "$OUTDIR"

NEXTFLOW_ARGS=(run "$NFCORE_ATACSEQ")
# A revision is used for the canonical remote pipeline. For a local pipeline
# directory, the checked-out local copy is used directly.
if [[ ! -d "$NFCORE_ATACSEQ" ]]; then
  NEXTFLOW_ARGS+=(-r "$NFCORE_ATACSEQ_REVISION")
fi

nextflow "${NEXTFLOW_ARGS[@]}" \
  --input "$DESIGN_CSV" \
  --genome "$GENOME" \
  -profile "$PROFILE" \
  --max_cpus "$MAX_CPUS" \
  --max_memory "$MAX_MEMORY" \
  --outdir "$OUTDIR" \
  --read_length "$READ_LENGTH" \
  --skip_plot_profile \
  --skip_plot_fingerprint \
  --skip_igv \
  --skip_deseq2_qc \
  --skip_peak_annotation \
  --skip_consensus_peaks \
  --skip_peak_qc
