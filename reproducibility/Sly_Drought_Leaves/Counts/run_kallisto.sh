#!/bin/bash

# kallisto quantification of PRJNA1125861 (Solanum lycopersicum, mature
# leaves, progressive water deficit). 21 paired-end runs, DNBSEQ-G400.
# SRR -> treatment/time mapping in coldata.tsv

set -u

# Reads (Clean_data/) and the SRR list (raw_data/lista) are not in the repository:
# download PRJNA1125861 from SRA and set BASE, INDEX and KALLISTO as needed.
BASE="${BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CLEAN="$BASE/Clean_data"
OUT="$BASE/Counts"
archivo_codigos="$BASE/raw_data/lista"
LOG="$OUT/kallisto.log"

# Ensembl Plants SL3.0 cdna.all, 34658 targets, 1 transcript per gene.
# Index version 10 (kallisto 0.48.x), checked with kallisto inspect.
INDEX="${INDEX:-$BASE/ref_tomato_SL3/index_tomato_09_26}"

# kallisto 0.48.0 from bioconda (with HDF5). The apt build (/usr/bin/kallisto) is
# compiled without HDF5 and does not write abundance.h5.
KALLISTO="${KALLISTO:-kallisto}"

HILOS=16

cd "$OUT" || exit 1

if [ ! -s "$INDEX" ]; then
  echo "ERROR: kallisto index not found at $INDEX" | tee -a "$LOG"
  echo "       set INDEX to the SL3.0 cDNA index" | tee -a "$LOG"
  exit 1
fi

# Read each line of the file
while IFS= read -r nombre
do
  [ -z "$nombre" ] && continue
  case "$nombre" in \#*) continue ;; esac

  R1="$CLEAN/clean_${nombre}_1.fastq.gz"
  R2="$CLEAN/clean_${nombre}_2.fastq.gz"

  # skip if read cleaning has not finished yet
  if [ ! -s "$R1" ] || [ ! -s "$R2" ]; then
    echo "[$(date +%F_%T)] $nombre has no clean data, skipped" | tee -a "$LOG"
    continue
  fi

  # skip samples already quantified
  if [ -s "$nombre/abundance.h5" ]; then
    echo "[$(date +%F_%T)] $nombre already quantified, skipped" | tee -a "$LOG"
    continue
  fi

  echo "[$(date +%F_%T)] Processing $nombre" | tee -a "$LOG"
  mkdir -p "$nombre"

  "$KALLISTO" quant -i "$INDEX" \
                 -o "$nombre" \
                 -t "$HILOS" \
                 "$R1" "$R2" 2>>"$LOG"

  if [ $? -ne 0 ] || [ ! -s "$nombre/abundance.h5" ]; then
    echo "[$(date +%F_%T)] ERROR in $nombre" | tee -a "$LOG"
    rm -rf "$nombre"
  fi
done < "$archivo_codigos"

echo "[$(date +%F_%T)] Aggregating counts" | tee -a "$LOG"
Rscript "$OUT/extraerCounts.r" "$OUT" 2>&1 | tee -a "$LOG"

echo "[$(date +%F_%T)] Done" | tee -a "$LOG"
