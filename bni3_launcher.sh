#!/bin/bash
# 
# BNI3 INTERACTIVE LAUNCHER
# ========================================
# Integrated pipeline wrapper for Boolean Network Inference 3
#

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m'

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

print_header() {
    clear
    printf "\n${CYAN}╔═══════════════════════════════════════════════════════════════╗${NC}\n"
    printf "${CYAN}║                                                               ║${NC}\n"
    printf "${CYAN}║                    ${BLUE}🧬 BNI3 PIPELINE LAUNCHER${CYAN}                  ║${NC}\n"
    printf "${CYAN}║                                                               ║${NC}\n"
    printf "${CYAN}╚═══════════════════════════════════════════════════════════════╝${NC}\n\n"
}

# Helper to check file existence
prompt_for_file() {
    local prompt_msg="$1"
    local var_name="$2"
    while true; do
        read -r -p "$prompt_msg" file_path
        if [ -f "$file_path" ]; then
            eval "$var_name=\"$file_path\""
            break
        else
            printf "${RED}ERROR: File '$file_path' does not exist. Please try again.${NC}\n"
        fi
    done
}

run_binarization() {
    printf "\n${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
    printf "${GREEN} MODULE 1: BINARIZATION${NC}\n"
    printf "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n\n"
    printf "This module converts your raw continuous gene expression data into binary states (0/1).\n\n"
    printf "${YELLOW}Available algorithms:${NC}\n"
    printf " - ${CYAN}1) SSD (Short Series Discretization):${NC} Identifies expression step transitions, excellent for time-series.\n"
    printf " - ${CYAN}2) WCSS (Within Cluster Sum of Squares):${NC} Minimizes within-cluster sum of squares for robust clustering.\n\n"
    
    printf "${YELLOW}Input Format Example:${NC}\n"
    printf "The input file must be a tab-separated (.tsv) containing Gene IDs as rows and timepoints as columns.\n"
    printf "ID      T1      T2      T3\n"
    printf "GeneA   0.5     2.1     1.8\n"
    printf "GeneB   1.2     0.3     0.4\n\n"

    prompt_for_file "Enter path to your raw counts TSV file: " input_file

    while true; do
        read -r -p "Choose algorithm (1 for SSD, 2 for WCSS): " algo_choice
        if [[ "$algo_choice" == "1" || "$algo_choice" == "2" ]]; then break; else printf "${RED}Invalid choice.${NC}\n"; fi
    done

    # Generate default output name dynamically based on input_file and algorithm
    input_dir=$(dirname "$input_file")
    input_name=$(basename "$input_file")
    
    if [ "$algo_choice" == "1" ]; then
        suffix="_binarized_SSD"
        algo="SSD"
    else
        suffix="_binarized_WCSS"
        algo="WCSS"
    fi
    
    default_out="${input_dir}/${input_name%.*}${suffix}.${input_name##*.}"

    read -r -p "Enter output FILE path [default: $default_out]: " out_file
    out_file=${out_file:-"$default_out"}

    out_dir=$(dirname "$out_file")
    mkdir -p "$out_dir"
    
    printf "\n${CYAN}>> Running Binarization...${NC}\n"
    if [ "$algo_choice" == "1" ]; then
        python3 "$ROOT_DIR/1.Binarization/BNI3_SSD.py" -i "$input_file" -o "$out_file"
    else
        python3 "$ROOT_DIR/1.Binarization/BNI3_WCSS.py" -i "$input_file" -o "$out_file"
    fi

    if [ $? -eq 0 ]; then
        printf "\n${CYAN}>> Running Binarization Behavior Reviewer...${NC}\n"
        python3 "$ROOT_DIR/1.Binarization/BNI3_behavior_reviewer.py" -i "$out_file" -o "$out_dir/reviewer_summary_${algo}.tsv"
        printf "\n${GREEN}✓ Binarization and review completed successfully! Check outputs in $out_dir${NC}\n"
    else
        printf "\n${RED}✗ Binarization failed! Check the error messages above.${NC}\n"
    fi
}

run_rules_inference() {
    printf "\n${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
    printf "${GREEN} MODULE 2: BOOLEAN RULES INFERENCE & EVALUATION${NC}\n"
    printf "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n\n"
    printf "This module infers logic transition rules among your genes using Gene Expression Programming.\n"
    printf "It automatically runs Inference, Analysis, and Evaluation.\n\n"

    printf "${YELLOW}Input Formats Example:${NC}\n"
    printf "1. Continuous Matrix (Genes as rows, Timepoints as columns):\n"
    printf "ID      T1      T2\nGeneA   0.5     2.1\nGeneB   1.2     0.3\n\n"
    printf "2. Binarized Matrix (Transposed - Genes as columns, Timepoints as rows):\n"
    printf "GeneA   GeneB\n0       1\n1       0\n\n"

    printf "${CYAN}Note on Advanced Control:${NC}\n"
    printf "This launcher uses stable defaults. If you need fine-grained control to change the number of processors (-p), repetitions (-n_rep), population size (-pop), or generations (-gen), please run the script directly from the terminal via:\n"
    printf "  ${YELLOW}python3 2.Rules_Inference/1.BNI3_Boolean_Rules_Inference.py -h${NC}\n\n"

    printf "${YELLOW}Mandatory parameters:${NC}\n"
    prompt_for_file "Enter path to original continuous TSV: " eval_raw
    prompt_for_file "Enter path to binarized TSV: " eval_bin
    
    printf "\n${YELLOW}Advanced parameters (Optional, press Enter to use defaults):${NC}\n"
    read -r -p "Optimal regulators K penalty? [default: 2]: " penalty_k
    penalty_k=${penalty_k:-2}

    # Generate default output directory dynamically based on binarized input file
    bin_dir=$(dirname "$eval_bin")
    bin_name=$(basename "$eval_bin")
    default_out_dir="${bin_dir}/${bin_name%.*}_Inferred_rules/"

    read -r -p "Enter output DIRECTORY [default: $default_out_dir]: " out_dir
    out_dir=${out_dir:-"$default_out_dir"}

    mkdir -p "$out_dir"
    
    printf "\n${CYAN}>> Running Rule Inference...${NC}\n"
    if python3 "$ROOT_DIR/2.Rules_Inference/1.BNI3_Boolean_Rules_Inference.py" -i "$eval_raw" -i_binary "$eval_bin" -o "$out_dir" -reg_optimal "$penalty_k" -p $(nproc); then
        
        printf "\n${CYAN}>> Running Rule Combinations Evaluator...${NC}\n"
        if python3 "$ROOT_DIR/2.Rules_Inference/3.BNI3_Evaluate_rules.py" -i "$out_dir/rules_by_gene.tsv" -m "$eval_bin" -o "$out_dir/evaluation_results.tsv" -n $(nproc) -v; then
            
            printf "\n${CYAN}>> Generating Network Topology Graphic...${NC}\n"
            python3 "$ROOT_DIR/2.Rules_Inference/4.BNI3_Boolean_network_visualizer.py" -i "$out_dir/rules_by_gene_evaluated.tsv" -o "$out_dir/network_topology.png"
            printf "Topology graphic saved at $out_dir/network_topology.png\n"
            printf "\n${GREEN}✓ Rules Inference pipeline completed successfully!${NC}\n"
        else
            printf "\n${RED}✗ Evaluation failed!${NC}\n"
        fi
    else
        printf "\n${RED}✗ Inference failed!${NC}\n"
    fi
}

run_attractors() {
    printf "\n${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
    printf "${GREEN} MODULE 3: ATTRACTORS ANALYSIS${NC}\n"
    printf "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n\n"
    
    printf "  ${CYAN}1)${NC} Find Attractors (and optional Perturbations)\n"
    printf "  ${CYAN}2)${NC} Simulate Trajectory Path to Attractor\n"
    printf "  ${CYAN}3)${NC} Visualize Attractors Landscape\n\n"
    
    while true; do
        read -r -p "Choose option [1-3]: " att_opt
        if [[ "$att_opt" == "1" || "$att_opt" == "2" || "$att_opt" == "3" ]]; then break; else printf "${RED}Invalid choice.${NC}\n"; fi
    done

    if [ "$att_opt" == "1" ]; then
        prompt_for_file "Path to rules table (e.g. rules_by_gene_evaluated.tsv): " rules_tsv
        read -r -p "Simulate Overexpression/Knockout? (e.g. ABF4:1,MYB44:0) [Press enter for none]: " mutation
        
        default_out=$(dirname "$rules_tsv")
        read -r -p "Enter output DIRECTORY [default: $default_out]: " out_dir
        out_dir=${out_dir:-"$default_out"}
        mkdir -p "$out_dir"
        
        printf "\n${CYAN}>> Finding Attractors...${NC}\n"
        if [ -n "$mutation" ]; then
            python3 "$ROOT_DIR/3.Attractors/1.BNI3_Attractors.py" -i "$rules_tsv" -o "$out_dir" -m "$mutation"
        else
            python3 "$ROOT_DIR/3.Attractors/1.BNI3_Attractors.py" -i "$rules_tsv" -o "$out_dir"
        fi
        [ $? -eq 0 ] && printf "\n${GREEN}✓ Attractors found successfully!${NC}\n" || printf "\n${RED}✗ Process failed!${NC}\n"

    elif [ "$att_opt" == "2" ]; then
        prompt_for_file "Path to attractors.tsv: " att_tsv
        prompt_for_file "Path to selected_rules.tsv: " sel_rules
        read -r -p "Initial Binary State (e.g. 1010010): " init_state
        
        default_out=$(dirname "$att_tsv")
        read -r -p "Enter output DIRECTORY [default: $default_out]: " out_dir
        out_dir=${out_dir:-"$default_out"}
        mkdir -p "$out_dir"
        
        printf "\n${CYAN}>> Simulating Path...${NC}\n"
        python3 "$ROOT_DIR/3.Attractors/2.BNI3_Path_to_Attractors.py" -a "$att_tsv" -r "$sel_rules" -s "$init_state" -o "$out_dir"
        [ $? -eq 0 ] && printf "\n${GREEN}✓ Simulation successful!${NC}\n" || printf "\n${RED}✗ Simulation failed!${NC}\n"

    elif [ "$att_opt" == "3" ]; then
        prompt_for_file "Path to attractors.tsv: " att_tsv
        printf "\n${CYAN}>> Rendering visualizations...${NC}\n"
        python3 "$ROOT_DIR/3.Attractors/3.BNI3_Visualize_Attractors.py" -i "$att_tsv" --heatmap --network
        [ $? -eq 0 ] && printf "\n${GREEN}✓ Visualizations rendered successfully!${NC}\n" || printf "\n${RED}✗ Visualization failed!${NC}\n"
    fi
}

# ================================================================
# MAIN MENU LOOP
# ================================================================
while true; do
    print_header
    printf "Please select the pipeline step you want to execute:\n\n"
    printf "  ${YELLOW}1)${NC} Binarization Module\n"
    printf "  ${YELLOW}2)${NC} Rules Inference Module\n"
    printf "  ${YELLOW}3)${NC} Attractors Analysis Module\n"
    printf "  ${RED}4)${NC} Exit\n\n"
    
    read -r -p "Enter option [1-4]: " option
    
    case $option in
        1) run_binarization ;;
        2) run_rules_inference ;;
        3) run_attractors ;;
        4)
            printf "\nExiting BNI3 Launcher. Happy researching!\n"
            exit 0
            ;;
        *) printf "\n${RED}Invalid option selected. Please try again.${NC}\n" ;;
    esac
    
    printf "\nPress Enter to return to the main menu..."
    read -r
done
