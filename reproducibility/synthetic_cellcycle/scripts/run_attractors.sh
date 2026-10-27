#!/bin/bash
#
# Synthetic cell-cycle network (paper): attractors and figures with BNI3.
# BNI3_Attractorator.py computes the attractors and draws the heatmap and the
# attractor diagram; 4.BNI3_Boolean_network_visualizer.py draws the regulatory network.
#
# Conditions: wild type, CycD knockout (CycD:0) and CycD overexpression (CycD:1).
# The run starts from the rules alone, without a binarized matrix, so no observed
# states are marked and no trajectory is computed (that step needs -b).
#
# Output: ../results/attractors/  and  ../results/network/
#
# Usage, from any directory: bash run_attractors.sh

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BASE=$(dirname "$SCRIPT_DIR")
BNI3=$(cd "${BASE}/../.." && pwd)

ATTRACTORATOR="${BNI3}/3.Attractors/BNI3_Attractorator.py"
NETWORK_VIS="${BNI3}/2.Rules_Inference/4.BNI3_Boolean_network_visualizer.py"
RULES="${BASE}/rules/rules_by_gene.tsv"
OUT="${BASE}/results"

# Conditions: "" = wild type
MUTATIONS=("" "CycD:0" "CycD:1")

# Upstream states drawn in the attractor diagram (0 = none)
PREDECESSORS=2

n_fail=0
for mut in "${MUTATIONS[@]}"; do
    echo ""
    echo "=========================================="
    echo "Condition: ${mut:-wild type}"
    echo "=========================================="

    args=(-i "$RULES" -O "$OUT" --no_path --predecessors "$PREDECESSORS" --svg)
    [ -n "$mut" ] && args+=(-m "$mut")

    if ! python3 "$ATTRACTORATOR" "${args[@]}"; then
        echo "✗ Failed: ${mut:-wild type}"
        n_fail=$((n_fail + 1))
    fi
done

echo ""
echo "=========================================="
echo "Regulatory network (same rules, Position 1)"
echo "=========================================="
if ! python3 "$NETWORK_VIS" -i "$RULES" -o "${OUT}/network/cellcycle_network"; then
    echo "✗ Failed: network drawing"
    n_fail=$((n_fail + 1))
fi

echo ""
echo "Results in: ${OUT}/attractors and ${OUT}/network"
[ "$n_fail" -eq 0 ]
