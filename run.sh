#!/bin/bash
set -euo pipefail

# ExpressPath — RNA-seq time course analysis pipeline launcher
# Usage: ./run.sh [snakemake options...]

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_DIR="$SCRIPT_DIR/pipeline"
DATA_DIR="$SCRIPT_DIR/data"
DESIGN_FILE="$DATA_DIR/design.yaml"
DESIGN_EXAMPLE="$DATA_DIR/design.example.yaml"
SETUP_HTML="$PIPELINE_DIR/setup_design.html"

RED=$'\033[0;31m'
GREEN=$'\033[0;32m'
YELLOW=$'\033[1;33m'
CYAN=$'\033[0;36m'
BOLD=$'\033[1m'
NC=$'\033[0m'

usage() {
    cat <<EOF
${BOLD}ExpressPath${NC} — RNA-seq time course analysis pipeline

${BOLD}Usage:${NC}
  ./run.sh [OPTIONS] [SNAKEMAKE_ARGS...]

${BOLD}Options:${NC}
  -h, --help   Show this help message

${BOLD}Examples:${NC}
  ./run.sh                       Run with 4 cores (default)
  ./run.sh -j8                   Run with 8 cores
  ./run.sh -n                    Dry-run (show what would run)
  ./run.sh --forcerun deseq2_analysis   Re-run just DESeq2 step

${BOLD}First-time setup:${NC}
  1. Open $SETUP_HTML in your browser
  2. Define cell lines, time points, treatment, replicates
  3. Browse for your TSV file → map columns
  4. Download design.yaml → save to ${DATA_DIR}/
  5. Put your TSV file in ${DATA_DIR}/
  6. Run: ./run.sh

${BOLD}Output:${NC}
  results/<timestamp>/pathway/interactive_report.html
EOF
}

# ─── Prerequisite checks ──────────────────────────────────────

check_conda() {
    if command -v conda &>/dev/null; then
        echo -e "  ${GREEN}✓${NC} conda found ($(conda --version 2>/dev/null | head -1))"
        return 0
    fi
    if command -v mamba &>/dev/null; then
        echo -e "  ${GREEN}✓${NC} mamba found"
        return 0
    fi
    echo -e "${RED}✗ conda not found${NC}"
    echo "  Install Miniconda: https://docs.conda.io/en/latest/miniconda.html"
    return 1
}

check_snakemake() {
    if command -v snakemake &>/dev/null; then
        echo -e "  ${GREEN}✓${NC} snakemake found ($(snakemake --version 2>/dev/null || true))"
        return 0
    fi
    echo -e "${RED}✗ snakemake not found${NC}"
    if command -v conda &>/dev/null; then
        echo "  Install: conda install -c conda-forge -c bioconda snakemake"
    elif command -v mamba &>/dev/null; then
        echo "  Install: mamba install -c conda-forge -c bioconda snakemake"
    fi
    return 1
}

# ─── Design file checks ───────────────────────────────────────

check_design() {
    if [ -f "$DESIGN_FILE" ]; then
        echo -e "  ${GREEN}✓${NC} design.yaml found"
        return 0
    fi

    echo -e "\n${YELLOW}design.yaml not found.${NC}"

    if [ -t 0 ]; then
        echo -e "How would you like to proceed?"
        if [ -f "$DESIGN_EXAMPLE" ]; then
            echo "  [1] Copy design.example.yaml and edit it"
        fi
        echo "  [2] Open setup_design.html in browser"
        echo "  [q] Quit"
        read -r -p "Choice: " choice
        echo

        case "$choice" in
            1)
                if [ ! -f "$DESIGN_EXAMPLE" ]; then
                    echo -e "${RED}design.example.yaml not found${NC}"
                    exit 1
                fi
                cp "$DESIGN_EXAMPLE" "$DESIGN_FILE"
                echo -e "${GREEN}Copied design.example.yaml → design.yaml${NC}"
                echo "Edit $DESIGN_FILE to match your experiment, then re-run:"
                echo "  vim $DESIGN_FILE"
                echo "  ./run.sh"
                exit 0
                ;;
            2)
                if command -v xdg-open &>/dev/null; then
                    xdg-open "$SETUP_HTML" 2>/dev/null || true
                elif command -v open &>/dev/null; then
                    open "$SETUP_HTML" 2>/dev/null || true
                else
                    echo "Open in browser: file://$SETUP_HTML"
                fi
                echo -e "\n${BOLD}Instructions:${NC}"
                echo "  1. Define cell lines, time points, treatment, replicates"
                echo "  2. Browse for your TSV file → map columns to samples"
                echo "  3. Download design.yaml → save to $DATA_DIR/"
                echo "  4. Put your TSV file in $DATA_DIR/"
                echo "  5. Re-run: ${CYAN}./run.sh${NC}"
                exit 0
                ;;
            *)
                exit 1
                ;;
        esac
    else
        echo "  To create one:"
        echo "    - Open pipeline/setup_design.html in your browser"
        echo "    - Or copy and edit data/design.example.yaml"
        exit 1
    fi
}

check_tsv() {
    local tsv_name
    tsv_name=$(python3 -c "
import yaml
with open('$DESIGN_FILE') as f:
    d = yaml.safe_load(f)
print(d.get('source_tsv', ''))
" 2>/dev/null || echo "")
    if [ -z "$tsv_name" ]; then
        echo -e "${RED}✗ design.yaml has no 'source_tsv' field${NC}"
        echo "  Re-export design.yaml from setup_design.html after browsing for your TSV file."
        exit 1
    fi
    local tsv_path="$DATA_DIR/$tsv_name"
    if [ ! -f "$tsv_path" ]; then
        echo -e "${RED}✗ TSV file not found: $tsv_path${NC}"
        echo "  Place '$tsv_name' in $DATA_DIR/"
        exit 1
    fi
    echo -e "  ${GREEN}✓${NC} $tsv_name found"
}

# ─── Main ─────────────────────────────────────────────────────

SNAKEMAKE_ARGS=(-j4)

# Parse run.sh own flags; everything else passed to snakemake
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        *)
            SNAKEMAKE_ARGS+=("$1")
            shift
            ;;
    esac
done

echo -e "\n${BOLD}${CYAN}=== ExpressPath Pipeline ===${NC}\n"

# Prerequisites
echo -e "${BOLD}Checking prerequisites:${NC}"
check_conda || exit 1
check_snakemake || exit 1

# Design
echo -e "\n${BOLD}Checking experiment design:${NC}"
check_design
check_tsv

# Run snakemake
echo -e "\n${BOLD}Running snakemake:${NC}"
echo "  cd pipeline && snakemake --use-conda ${SNAKEMAKE_ARGS[*]}"
echo

cd "$PIPELINE_DIR"
snakemake --use-conda "${SNAKEMAKE_ARGS[@]}"

# Show result path
run_dir=$(ls -dt "$SCRIPT_DIR/results"/*/ 2>/dev/null | head -1 || true)
if [ -n "$run_dir" ]; then
    run_id=$(basename "$run_dir")
    echo -e "\n${GREEN}${BOLD}✓ Pipeline complete${NC}"
    echo -e "  Report: ${CYAN}results/$run_id/pathway/interactive_report.html${NC}"
    echo -e "  Output: ${CYAN}results/$run_id/${NC}"
else
    echo -e "\n${GREEN}${BOLD}✓ Pipeline complete${NC}"
fi
