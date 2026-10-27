#!/bin/bash

# Cuantificacion con kallisto de PRJNA1125861 (Solanum lycopersicum, hoja
# madura, deficit hidrico progresivo). 21 runs paired-end, DNBSEQ-G400.
# Mapeo SRR -> tratamiento/tiempo en ../raw_data/samples_metadata.tsv

set -u

# Reads (Clean_data/) and the SRR list (raw_data/lista) are not in the repository:
# download PRJNA1125861 from SRA and set BASE, INDEX and KALLISTO as needed.
BASE="${BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CLEAN="$BASE/Clean_data"
OUT="$BASE/Counts"
archivo_codigos="$BASE/raw_data/lista"
LOG="$OUT/kallisto.log"

# Ensembl Plants SL3.0 cdna.all, 34658 targets, 1 transcrito por gen.
# Index version 10 (kallisto 0.48.x), verificado con kallisto inspect.
INDEX="${INDEX:-$BASE/ref_tomato_SL3/index_tomato_09_26}"

# kallisto 0.48.0 de bioconda (con HDF5). El de apt (/usr/bin/kallisto) esta
# compilado sin HDF5 y no escribe abundance.h5.
KALLISTO="${KALLISTO:-kallisto}"

HILOS=16

cd "$OUT" || exit 1

if [ ! -s "$INDEX" ]; then
  echo "ERROR: no se encuentra el indice de kallisto en $INDEX" | tee -a "$LOG"
  echo "       ver README.md en este directorio" | tee -a "$LOG"
  exit 1
fi

# Leer cada linea del archivo
while IFS= read -r nombre
do
  [ -z "$nombre" ] && continue
  case "$nombre" in \#*) continue ;; esac

  R1="$CLEAN/clean_${nombre}_1.fastq.gz"
  R2="$CLEAN/clean_${nombre}_2.fastq.gz"

  # saltar si la limpieza aun no termino
  if [ ! -s "$R1" ] || [ ! -s "$R2" ]; then
    echo "[$(date +%F_%T)] $nombre sin clean data, se omite" | tee -a "$LOG"
    continue
  fi

  # saltar muestras ya cuantificadas
  if [ -s "$nombre/abundance.h5" ]; then
    echo "[$(date +%F_%T)] $nombre ya cuantificado, se omite" | tee -a "$LOG"
    continue
  fi

  echo "[$(date +%F_%T)] Procesando $nombre" | tee -a "$LOG"
  mkdir -p "$nombre"

  "$KALLISTO" quant -i "$INDEX" \
                 -o "$nombre" \
                 -t "$HILOS" \
                 "$R1" "$R2" 2>>"$LOG"

  if [ $? -ne 0 ] || [ ! -s "$nombre/abundance.h5" ]; then
    echo "[$(date +%F_%T)] ERROR en $nombre" | tee -a "$LOG"
    rm -rf "$nombre"
  fi
done < "$archivo_codigos"

echo "[$(date +%F_%T)] Agregando conteos" | tee -a "$LOG"
Rscript "$OUT/extraerCounts.r" "$OUT" 2>&1 | tee -a "$LOG"

echo "[$(date +%F_%T)] Listo" | tee -a "$LOG"
