#!/usr/bin/env bash
set -euo pipefail

# Convert per-well paired-end BAM files into one ArchR/10x-style fragments file.
#
# Run from the directory containing the per-well BAM files.
#
# Optional environment variables:
#   THREADS=8
#   JOBS=8
#   MAPQ=30
#   WORKDIR=work_fragments
#   OUT=fragments.tsv.gz

THREADS="${THREADS:-8}"
JOBS="${JOBS:-8}"
MAPQ="${MAPQ:-30}"
WORKDIR="${WORKDIR:-work_fragments}"
OUT="${OUT:-fragments.tsv.gz}"

mkdir -p "$WORKDIR"

have_cmd() { command -v "$1" >/dev/null 2>&1; }
for c in samtools bedtools bgzip tabix sort awk grep zcat uniq; do
  have_cmd "$c" || { echo "ERROR: '$c' not found in PATH" >&2; exit 1; }
done

process_one_bam() {
  local bam="$1"
  local stem
  stem="$(basename "$bam" .bam)"

  # Extract a well barcode from filenames such as d2_r12_c1.
  local bc
  if [[ "$stem" =~ (d[0-9]+_r[0-9]+_c[0-9]+) ]]; then
    bc="${BASH_REMATCH[1]}"
  elif [[ "$stem" =~ (r[0-9]+_c[0-9]+) ]]; then
    bc="${BASH_REMATCH[1]}"
  elif [[ "$stem" =~ (r[0-9]+c[0-9]+) ]]; then
    bc="${BASH_REMATCH[1]}"
  else
    bc="$stem"
  fi

  local tmp="$WORKDIR/$stem"
  echo "[*] $stem (barcode=$bc)"

  # 1. Mark/remove PCR duplicates.
  samtools sort -n -@ "$THREADS" -o "$tmp.namesort.bam" "$bam"
  samtools fixmate -m -@ "$THREADS" "$tmp.namesort.bam" "$tmp.fixmate.bam"
  samtools sort -@ "$THREADS" -o "$tmp.coordsort.bam" "$tmp.fixmate.bam"
  samtools markdup -r -s -@ "$THREADS" "$tmp.coordsort.bam" "$tmp.dedup.bam"

  # 2. Keep properly paired, primary, mapped reads with MAPQ >= threshold.
  #    -f 3: paired + properly paired
  #    -F 3852: remove unmapped, mate-unmapped, secondary, QC-fail,
  #             duplicate, and supplementary alignments.
  samtools view -b -@ "$THREADS" \
    -f 3 -F 3852 -q "$MAPQ" \
    "$tmp.dedup.bam" > "$tmp.filtered.bam"

  # 3. Convert paired alignments to fragment intervals on the same chromosome.
  #    Output columns: chrom, start, end, barcode, count.
  #
  # BEDPE records are restricted to read pairs mapping to the same chromosome.
  bedtools bamtobed -bedpe -i "$tmp.filtered.bam" \
    | awk -v OFS='\t' -v BC="$bc" '$1==$4 { s=$2; if($5<s) s=$5; e=$3; if($6>e) e=$6; if(e>s) print $1, s, e, BC, 1 }' \
    | bgzip -@ "$THREADS" > "$tmp.fragments.tsv.gz"

  tabix -s 1 -b 2 -e 3 "$tmp.fragments.tsv.gz"

  rm -f \
    "$tmp.namesort.bam" \
    "$tmp.fixmate.bam" \
    "$tmp.coordsort.bam" \
    "$tmp.dedup.bam" \
    "$tmp.filtered.bam"
}

export -f process_one_bam
export WORKDIR THREADS MAPQ OUT

shopt -s nullglob
BAMS=( *.bam )
if [[ ${#BAMS[@]} -eq 0 ]]; then
  echo "ERROR: no *.bam files found in $(pwd)" >&2
  exit 1
fi

if command -v parallel >/dev/null 2>&1; then
  parallel -j "$JOBS" process_one_bam ::: "${BAMS[@]}"
else
  for bam in "${BAMS[@]}"; do
    process_one_bam "$bam"
  done
fi

echo "[*] Merging per-well fragments..."
zcat "$WORKDIR"/*.fragments.tsv.gz \
  | sort -T "$WORKDIR" -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 \
  | uniq \
  | bgzip -@ "$THREADS" > "$OUT"

tabix -s 1 -b 2 -e 3 "$OUT"

echo "Done: $OUT and ${OUT}.tbi"
