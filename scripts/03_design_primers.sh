#!/usr/bin/env bash
# =============================================================================
# 03_design_primers.sh
# Design PCR primers for InDel markers using Primer3.
#
# For each InDel in the filtered TSV:
#   1. Extract ±FLANK bp flanking sequence with samtools faidx
#   2. Write a primer3 input file
#   3. Run primer3_core
#   4. If 0 pairs returned, retry with FLANK*2, then FLANK*3
#   5. Parse output into a summary TSV
#
# Usage:
#   bash 03_design_primers.sh [options]
#
# Options:
#   --indels   FILE    Filtered InDel TSV from step 02 (required)
#   --ref      FILE    Reference FASTA, must be indexed (required)
#   --outdir   DIR     Output directory for primer3 files (required)
#   --flank    INT     Flanking sequence size bp each side (default: 200)
#   --config   FILE    Primer3 parameter file (default: config/primer3_defaults.txt)
#   -h, --help         Show this help message and exit
#
# Output:
#   <outdir>/{chrom}_{pos}.primer3.in   — Primer3 input
#   <outdir>/{chrom}_{pos}.primer3.out  — Primer3 output
#   <outdir>/primer_summary.tsv         — Summary table
#
# Requirements:
#   samtools >= 1.17, primer3_core >= 2.6 (must be in PATH)
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
INDELS=""
REF=""
OUTDIR=""
FLANK=200
CONFIG="$(dirname "$0")/../config/primer3_defaults.txt"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
usage() {
    grep '^#' "$0" | grep -v '^#!/' | sed 's/^# \{0,3\}//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --indels)  INDELS="$2";  shift 2 ;;
        --ref)     REF="$2";     shift 2 ;;
        --outdir)  OUTDIR="$2";  shift 2 ;;
        --flank)   FLANK="$2";   shift 2 ;;
        --config)  CONFIG="$2";  shift 2 ;;
        -h|--help) usage ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
for var_name in INDELS REF OUTDIR; do
    val="${!var_name}"
    if [[ -z "$val" ]]; then
        echo "[ERROR] --${var_name,,} is required." >&2
        exit 1
    fi
done

if [[ ! -f "$INDELS" ]]; then
    echo "[ERROR] InDel file not found: $INDELS" >&2; exit 1
fi
if [[ ! -f "$REF" ]]; then
    echo "[ERROR] Reference FASTA not found: $REF" >&2; exit 1
fi
if [[ ! -f "${REF}.fai" ]]; then
    echo "[ERROR] Reference index not found: ${REF}.fai — run step 01 first." >&2; exit 1
fi
if [[ ! -f "$CONFIG" ]]; then
    echo "[ERROR] Primer3 config not found: $CONFIG" >&2; exit 1
fi
if ! command -v samtools &>/dev/null; then
    echo "[ERROR] samtools not in PATH." >&2; exit 1
fi
if ! command -v primer3_core &>/dev/null; then
    echo "[ERROR] primer3_core not in PATH." >&2; exit 1
fi

mkdir -p "$OUTDIR"

# ---------------------------------------------------------------------------
# Helper: load primer3 defaults from config file (strip comments, blank lines)
# ---------------------------------------------------------------------------
load_config() {
    local cfg="$1"
    grep -v '^#' "$cfg" | grep -v '^[[:space:]]*$'
}

# ---------------------------------------------------------------------------
# Helper: run primer3 and return number of pairs found
# ---------------------------------------------------------------------------
run_primer3() {
    local input_file="$1"
    local output_file="$2"
    primer3_core < "$input_file" > "$output_file" 2>/dev/null
    local npairs
    npairs=$(grep -m1 'PRIMER_PAIR_NUM_RETURNED=' "$output_file" 2>/dev/null \
             | cut -d= -f2 | tr -d '[:space:]') || true
    echo "${npairs:-0}"
}

# ---------------------------------------------------------------------------
# Helper: write primer3 input block
# Arguments: seq_id template_seq flank_size config_block output_file
# ---------------------------------------------------------------------------
write_primer3_input() {
    local seq_id="$1"
    local template="$2"
    local flank="$3"
    local config_block="$4"
    local out="$5"

    local tlen=${#template}
    # Target: 5 bp before the indel junction, 10 bp wide
    # InDel is centered at position FLANK in the template (0-based)
    local target_start=$(( flank - 5 ))
    [[ $target_start -lt 0 ]] && target_start=0

    {
        echo "SEQUENCE_ID=${seq_id}"
        echo "SEQUENCE_TEMPLATE=${template}"
        echo "SEQUENCE_TARGET=${target_start},10"
        echo "$config_block"
        echo "="
    } > "$out"
}

# ---------------------------------------------------------------------------
# Helper: parse best pair from primer3 output
# Returns TSV line: LEFT_SEQ LEFT_TM LEFT_POS RIGHT_SEQ RIGHT_TM RIGHT_POS PRODUCT_SIZE PAIR_INDEX
# ---------------------------------------------------------------------------
parse_primer3_best() {
    local p3out="$1"
    local pair_idx="${2:-0}"

    local left right left_tm right_tm product left_pos right_pos
    left=$(grep    "PRIMER_LEFT_${pair_idx}_SEQUENCE=" "$p3out" | cut -d= -f2 | tr -d '[:space:]') || true
    right=$(grep   "PRIMER_RIGHT_${pair_idx}_SEQUENCE=" "$p3out" | cut -d= -f2 | tr -d '[:space:]') || true
    left_tm=$(grep "PRIMER_LEFT_${pair_idx}_TM=" "$p3out" | cut -d= -f2 | tr -d '[:space:]') || true
    right_tm=$(grep "PRIMER_RIGHT_${pair_idx}_TM=" "$p3out" | cut -d= -f2 | tr -d '[:space:]') || true
    product=$(grep "PRIMER_PAIR_${pair_idx}_PRODUCT_SIZE=" "$p3out" | cut -d= -f2 | tr -d '[:space:]') || true
    left_pos=$(grep "^PRIMER_LEFT_${pair_idx}=" "$p3out" | cut -d= -f2 | tr -d '[:space:]') || true
    right_pos=$(grep "^PRIMER_RIGHT_${pair_idx}=" "$p3out" | cut -d= -f2 | tr -d '[:space:]') || true

    echo "${left:-NA}	${left_tm:-NA}	${left_pos:-NA}	${right:-NA}	${right_tm:-NA}	${right_pos:-NA}	${product:-NA}	${pair_idx}"
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
SUMMARY="${OUTDIR}/primer_summary.tsv"
echo -e "CHROM\tPOS\tINDEL_SIZE\tQUAL\tALT\tLEFT_PRIMER\tLEFT_TM\tLEFT_POS\tRIGHT_PRIMER\tRIGHT_TM\tRIGHT_POS\tPRODUCT_SIZE\tFLANK\tPAIR_INDEX\tSTATUS" \
    > "$SUMMARY"

# Load config once
CONFIG_BLOCK=$(load_config "$CONFIG")

n_total=0
n_designed=0
n_failed=0

echo "[INFO] Starting primer design..." >&2
echo "[INFO] Config:  $CONFIG" >&2
echo "[INFO] InDels:  $INDELS" >&2
echo "[INFO] Ref:     $REF" >&2
echo "[INFO] Outdir:  $OUTDIR" >&2
echo "[INFO] Flank:   ${FLANK} bp (initial)" >&2

# Read filtered InDel TSV (skip header)
while IFS=$'\t' read -r chrom pos indel_size qual alt notes; do
    [[ "$chrom" == "CHROM" ]] && continue  # skip header
    [[ -z "$chrom" || "$chrom" == \#* ]] && continue

    n_total=$(( n_total + 1 ))
    echo "[INFO] Processing: ${chrom}:${pos} (size=${indel_size}, qual=${qual})" >&2

    safe_id="${chrom}_${pos}"
    p3in="${OUTDIR}/${safe_id}.primer3.in"
    p3out="${OUTDIR}/${safe_id}.primer3.out"

    npairs=0
    used_flank=$FLANK

    # Try initial flank, then 2x, then 3x
    for multiplier in 1 2 3; do
        current_flank=$(( FLANK * multiplier ))

        # Extract flanking sequence: chrom:start-end (1-based, samtools faidx is 1-based)
        extract_start=$(( pos - current_flank ))
        extract_end=$(( pos + current_flank ))
        [[ $extract_start -lt 1 ]] && extract_start=1

        region="${chrom}:${extract_start}-${extract_end}"
        template=$(samtools faidx "$REF" "$region" 2>/dev/null | grep -v '^>' | tr -d '\n' | tr '[:lower:]' '[:upper:]') || true

        if [[ -z "$template" ]]; then
            echo "[WARN]  Could not extract sequence for ${region}" >&2
            continue
        fi

        # Write primer3 input
        write_primer3_input "$safe_id" "$template" "$current_flank" "$CONFIG_BLOCK" "$p3in"

        # Run primer3
        npairs=$(run_primer3 "$p3in" "$p3out")

        if [[ "$npairs" -gt 0 ]]; then
            used_flank=$current_flank
            echo "[INFO]  Found ${npairs} pair(s) with flank=${current_flank}" >&2
            break
        else
            echo "[INFO]  No pairs with flank=${current_flank}, retrying..." >&2
        fi
    done

    if [[ "$npairs" -eq 0 ]]; then
        echo "[WARN]  No primers designed for ${safe_id}" >&2
        echo -e "${chrom}\t${pos}\t${indel_size}\t${qual}\t${alt}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t${used_flank}\tNA\tNO_PRIMERS" \
            >> "$SUMMARY"
        n_failed=$(( n_failed + 1 ))
        continue
    fi

    # Parse best pair (pair index 0)
    pair_info=$(parse_primer3_best "$p3out" 0)

    echo -e "${chrom}\t${pos}\t${indel_size}\t${qual}\t${alt}\t${pair_info}\t${used_flank}\tOK" \
        >> "$SUMMARY"
    n_designed=$(( n_designed + 1 ))

done < "$INDELS"

echo "" >&2
echo "[SUMMARY] Total InDels processed: ${n_total}" >&2
echo "[SUMMARY] Primers designed:       ${n_designed}" >&2
echo "[SUMMARY] Failed (no primers):    ${n_failed}" >&2
echo "[SUMMARY] Output summary:         ${SUMMARY}" >&2
echo "[DONE] Primer design complete." >&2
