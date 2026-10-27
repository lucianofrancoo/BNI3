#!/bin/bash
#
# DREAM4 size10 benchmark with BNI3 (paper, Figure 5). Everything runs in series:
#
#   1) Binarization      1.Binarization/BNI3_Binarizator.py (SSD + WCSS + behavior reviewer)
#   2) Method choice     highest Unique_Ratio in the reviewer summary (tie -> WCSS)
#   3) Rule inference    2.Rules_Inference/1.BNI3_Boolean_Rules_Inference.py
#                        (chains the analysis and the rule combination evaluation)
#   4) Gold standard     eval_gold_standard.r (structural metrics vs the DREAM4 gold standard)
#
# Inputs:   ../data/time_series/size10_X/   21-point series and their 5-point subsets (_5times)
#           ../data/gold_standard/
# Outputs:  ${OUT}/binarized/size10_X/<series>/
#           ${OUT}/rules/size10_X/rules_<METHOD>_<series>/
#           ${OUT}/BNI3_size10_complete_evaluation_results.tsv
#
# OUT defaults to ../rerun so the results used in the paper (../results/BNI3) are never
# overwritten. Usage, from any directory:
#   nohup bash run_bni3_dream4.sh > run_bni3_dream4.log 2>&1 &
#
# Resumable: with SKIP_EXISTING=1 finished steps are skipped (a binarization with its
# reviewer summary; a rules folder with rules_by_gene_evaluated.tsv). Delete a rules
# folder left half-way before relaunching.
#
# Runtime on the paper's server (-p 50): ~2.2 h per 21-point series, ~2 min per 5-point one.

# ============================================
# CONFIGURATION
# ============================================

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DREAM4_DIR=$(dirname "$SCRIPT_DIR")
BNI3=$(cd "${DREAM4_DIR}/../.." && pwd)

DATA="${DREAM4_DIR}/data"
OUT="${OUT:-${DREAM4_DIR}/rerun}"

BINARIZATOR="${BNI3}/1.Binarization/BNI3_Binarizator.py"
INFERENCE="${BNI3}/2.Rules_Inference/1.BNI3_Boolean_Rules_Inference.py"
EVAL_R="${SCRIPT_DIR}/eval_gold_standard.r"

PYTHON="python3"

DO_BINARIZATION=1
DO_INFERENCE=1
DO_GOLD_EVAL=1

SKIP_EXISTING=1
NETWORKS=(1 2 3 4 5)

BIN_PROCS=1

# Inference parameters used in the paper (the evaluation is chained by default)
INF_ARGS=(-p 50 -n_rep 10 -pop 50 -gen 50)

# ============================================
# FUNCTIONS
# ============================================

N_OK=0
N_FAIL=0
N_SKIP=0
FAILED_LIST=()

ts() { date '+%Y-%m-%d %H:%M:%S'; }

log() { echo "[$(ts)] $*"; }

fail() {
    N_FAIL=$((N_FAIL + 1))
    FAILED_LIST+=("$1")
    log "  ✗ $1"
}

# Unique_Ratio of one method in the reviewer summary (columns looked up by name)
get_ratio() {
    local summary="$1" method="$2"
    awk -F'\t' -v m="$method" '
        NR==1 { for (i=1; i<=NF; i++) { if ($i=="Method") cm=i; if ($i=="Unique_Ratio") cr=i }; next }
        $cm==m { print $cr; exit }' "$summary"
}

# binarize <raw_file> <out_dir>
binarize() {
    local raw="$1" out="$2"
    local stem
    stem=$(basename "$raw"); stem="${stem%.tsv}"
    local summary="${out}/${stem}_reviewer_summary.tsv"

    if [ "$SKIP_EXISTING" -eq 1 ] && [ -f "$summary" ] \
        && [ -f "${out}/${stem}_binarized_SSD.tsv" ] && [ -f "${out}/${stem}_binarized_WCSS.tsv" ]; then
        log "  ⏭  Binarization already done: ${stem}"
        N_SKIP=$((N_SKIP + 1))
        return 0
    fi

    log "  → Binarizing: ${stem}"
    mkdir -p "$out"
    if "$PYTHON" "$BINARIZATOR" -i "$raw" -O "$out" -p "$BIN_PROCS"; then
        N_OK=$((N_OK + 1))
    else
        fail "binarization ${stem}"
    fi
}

# infer <raw_file> <bin_dir> <rules_dir> <args...>
infer() {
    local raw="$1" bindir="$2" rulesdir="$3"; shift 3
    local stem
    stem=$(basename "$raw"); stem="${stem%.tsv}"
    local summary="${bindir}/${stem}_reviewer_summary.tsv"

    if [ ! -f "$summary" ]; then
        fail "inference ${stem}: ${summary} not found"
        return 1
    fi

    local ssd wcss method
    ssd=$(get_ratio "$summary" SSD)
    wcss=$(get_ratio "$summary" WCSS)
    if [ -z "$ssd" ] || [ -z "$wcss" ]; then
        fail "inference ${stem}: Unique_Ratio not found (SSD='${ssd}', WCSS='${wcss}')"
        return 1
    fi

    # Tie -> WCSS
    if (( $(echo "$wcss >= $ssd" | bc -l) )); then method="WCSS"; else method="SSD"; fi

    local bin_file="${bindir}/${stem}_binarized_${method}.tsv"
    local out="${rulesdir}/rules_${method}_${stem}"

    # With the chained evaluation a run is complete only once the evaluated file exists
    local done_file="${out}/rules_by_gene_evaluated.tsv"
    [[ " $* " == *" --no_evaluation "* ]] && done_file="${out}/rules_by_gene.tsv"

    if [ "$SKIP_EXISTING" -eq 1 ] && [ -f "$done_file" ]; then
        log "  ⏭  Rules already inferred: rules_${method}_${stem}"
        N_SKIP=$((N_SKIP + 1))
        return 0
    fi

    log "  → Inference: ${stem} | SSD=${ssd} WCSS=${wcss} → ${method}"
    mkdir -p "$out"
    if "$PYTHON" "$INFERENCE" -i "$raw" -i_binary "$bin_file" -o "$out" "$@"; then
        N_OK=$((N_OK + 1))
    else
        fail "inference ${stem} (${method})"
    fi
}

# ============================================
# RUN
# ============================================

START=$(date +%s)
log "START DREAM4 size10 with BNI3"
log "BNI3:   ${BNI3}"
log "Output: ${OUT}"

for f in "$BINARIZATOR" "$INFERENCE" "$EVAL_R"; do
    [ -f "$f" ] || { log "ERROR: $f not found"; exit 1; }
done

for net in "${NETWORKS[@]}"; do
    raw_dir="${DATA}/time_series/size10_${net}"
    log "=== size10_${net} ==="
    for raw in "$raw_dir"/*; do
        [ -f "$raw" ] || continue
        stem=$(basename "$raw"); stem="${stem%.tsv}"
        bindir="${OUT}/binarized/size10_${net}/${stem}"

        if [ "$DO_BINARIZATION" -eq 1 ]; then
            binarize "$raw" "$bindir"
        fi
        if [ "$DO_INFERENCE" -eq 1 ]; then
            infer "$raw" "$bindir" "${OUT}/rules/size10_${net}" "${INF_ARGS[@]}"
        fi
    done
done

if [ "$DO_GOLD_EVAL" -eq 1 ]; then
    log "=== Gold standard evaluation ==="
    if Rscript "$EVAL_R" "${OUT}/rules" "${OUT}/binarized" "${DATA}/gold_standard" \
        "${OUT}/BNI3_size10_complete_evaluation_results.tsv" > "${OUT}/eval_gold_standard.log" 2>&1; then
        N_OK=$((N_OK + 1))
        log "  ✓ ${OUT}/BNI3_size10_complete_evaluation_results.tsv"
    else
        fail "gold standard evaluation (see ${OUT}/eval_gold_standard.log)"
    fi
fi

ELAPSED=$(( $(date +%s) - START ))
log "=========================================="
log "END  ($((ELAPSED / 3600))h $(((ELAPSED % 3600) / 60))m)"
log "  OK: ${N_OK}   Skipped: ${N_SKIP}   Failed: ${N_FAIL}"
for f in "${FAILED_LIST[@]}"; do log "    - $f"; done
log "=========================================="

[ "$N_FAIL" -eq 0 ]
