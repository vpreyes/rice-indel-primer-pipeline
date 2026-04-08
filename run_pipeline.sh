#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh
# Master script for the Rice InDel Primer Design Pipeline.
#
# Runs all 7 steps in sequence:
#   01  Index reference genome (samtools faidx)
#   02  Filter InDel candidates (quality, size, repeat, spacing)
#   03  Design primers with Primer3
#   04  Build BLAST nucleotide database
#   05  Check primer specificity (BLAST)
#   06  Validate primer positions vs InDel
#   07  Format output (IDS naming convention)
#
# Usage:
#   bash run_pipeline.sh [options]
#
# Required options:
#   --ref      FILE   Reference genome FASTA
#   --indels   FILE   InDel TSV file (format: #CHROM POS ID REF REFLEN ALT ALTLEN QUAL)
#   --chrom    STR    Target chromosome (e.g., chr02)
#   --start    INT    Window start position (bp)
#   --end      INT    Window end position (bp)
#   --outdir   DIR    Output directory (will be created if absent)
#
# Optional:
#   --flank    INT    Flanking sequence for Primer3 (default: 200 bp)
#   --spacing  INT    Minimum spacing between markers (default: 150000 bp)
#   --threads  INT    BLAST threads (default: 4)
#   --config   FILE   Primer3 config file (default: config/primer3_defaults.txt)
#   -h, --help        Show this help message
#
# Output:
#   <outdir>/filtered_indels.tsv       Filtered InDel candidates (step 02)
#   <outdir>/primer3/primer_summary.tsv  Initial primer design (step 03)
#   <outdir>/blast_db/                 BLAST database (step 04)
#   <outdir>/primer3/final_summary.tsv  After specificity check (step 05)
#   <outdir>/primer3/validation_report.tsv  Position validation (step 06)
#   <outdir>/IDS_primers.tsv           Final IDS-named primer table (step 07)
#   <outdir>/IDS_primers.txt           Plain text primer list for ordering (step 07)
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Script location
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
REF=""
INDELS=""
CHROM=""
START=""
END=""
OUTDIR=""
FLANK=200
SPACING=150000
THREADS=4
CONFIG="${SCRIPT_DIR}/config/primer3_defaults.txt"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
usage() {
    grep '^#' "$0" | grep -v '^#!/' | sed 's/^# \{0,3\}//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)     REF="$2";     shift 2 ;;
        --indels)  INDELS="$2";  shift 2 ;;
        --chrom)   CHROM="$2";   shift 2 ;;
        --start)   START="$2";   shift 2 ;;
        --end)     END="$2";     shift 2 ;;
        --outdir)  OUTDIR="$2";  shift 2 ;;
        --flank)   FLANK="$2";   shift 2 ;;
        --spacing) SPACING="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --config)  CONFIG="$2";  shift 2 ;;
        -h|--help) usage ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate required arguments
# ---------------------------------------------------------------------------
missing=()
[[ -z "$REF"    ]] && missing+=("--ref")
[[ -z "$INDELS" ]] && missing+=("--indels")
[[ -z "$CHROM"  ]] && missing+=("--chrom")
[[ -z "$START"  ]] && missing+=("--start")
[[ -z "$END"    ]] && missing+=("--end")
[[ -z "$OUTDIR" ]] && missing+=("--outdir")

if [[ ${#missing[@]} -gt 0 ]]; then
    echo "[ERROR] Missing required arguments: ${missing[*]}" >&2
    echo "Run with --help for usage." >&2
    exit 1
fi

mkdir -p "$OUTDIR"
P3DIR="${OUTDIR}/primer3"
DBDIR="${OUTDIR}/blast_db"
mkdir -p "$P3DIR" "$DBDIR"

FILTERED="${OUTDIR}/filtered_indels.tsv"
P3SUMMARY="${P3DIR}/primer_summary.tsv"
BLASTDB="${DBDIR}/rice_ref"
FINAL_SUMMARY="${P3DIR}/final_summary.tsv"

# ---------------------------------------------------------------------------
# Log header
# ---------------------------------------------------------------------------
echo "============================================================" >&2
echo " Rice InDel Primer Design Pipeline" >&2
echo "============================================================" >&2
echo "  Reference:  $REF" >&2
echo "  InDels:     $INDELS" >&2
echo "  Chrom:      $CHROM  [$START - $END]" >&2
echo "  Outdir:     $OUTDIR" >&2
echo "  Flank:      $FLANK bp" >&2
echo "  Spacing:    $SPACING bp" >&2
echo "  Threads:    $THREADS" >&2
echo "  Config:     $CONFIG" >&2
echo "============================================================" >&2
echo "" >&2

# ---------------------------------------------------------------------------
# Step 01: Index reference
# ---------------------------------------------------------------------------
echo "[STEP 01] Indexing reference genome..." >&2
if [[ -f "${REF}.fai" ]]; then
    echo "[SKIP]   Index already exists: ${REF}.fai" >&2
else
    bash "${SCRIPT_DIR}/scripts/01_index_reference.sh" "$REF"
fi
echo "" >&2

# ---------------------------------------------------------------------------
# Step 02: Filter InDels
# ---------------------------------------------------------------------------
echo "[STEP 02] Filtering InDel candidates..." >&2
python3 "${SCRIPT_DIR}/scripts/02_filter_indels.py" \
    --indels   "$INDELS" \
    --chrom    "$CHROM" \
    --start    "$START" \
    --end      "$END" \
    --spacing  "$SPACING" \
    --output   "$FILTERED"
echo "" >&2

N_INDELS=$(tail -n +2 "$FILTERED" | wc -l | tr -d ' ')
echo "[INFO] Filtered InDels: $N_INDELS" >&2

if [[ "$N_INDELS" -eq 0 ]]; then
    echo "[ERROR] No InDels passed filtering. Adjust --start/--end or filter thresholds." >&2
    exit 1
fi
echo "" >&2

# ---------------------------------------------------------------------------
# Step 03: Design primers
# ---------------------------------------------------------------------------
echo "[STEP 03] Designing primers with Primer3..." >&2
bash "${SCRIPT_DIR}/scripts/03_design_primers.sh" \
    --indels  "$FILTERED" \
    --ref     "$REF" \
    --outdir  "$P3DIR" \
    --flank   "$FLANK" \
    --config  "$CONFIG"
echo "" >&2

# ---------------------------------------------------------------------------
# Step 04: Build BLAST database
# ---------------------------------------------------------------------------
echo "[STEP 04] Building BLAST database..." >&2
if [[ -f "${BLASTDB}.nhr" || -f "${BLASTDB}.00.nhr" ]]; then
    echo "[SKIP]   BLAST database already exists: $BLASTDB" >&2
else
    bash "${SCRIPT_DIR}/scripts/04_build_blast_db.sh" "$REF" "$DBDIR"
fi
echo "" >&2

# ---------------------------------------------------------------------------
# Step 05: Check specificity
# ---------------------------------------------------------------------------
echo "[STEP 05] Checking primer specificity via BLAST..." >&2
python3 "${SCRIPT_DIR}/scripts/05_check_specificity.py" \
    --summary    "$P3SUMMARY" \
    --outdir     "$P3DIR" \
    --db         "$BLASTDB" \
    --threads    "$THREADS" \
    --target-chr "$CHROM"
echo "" >&2

# ---------------------------------------------------------------------------
# Step 06: Validate positions
# ---------------------------------------------------------------------------
echo "[STEP 06] Validating primer positions vs InDel..." >&2
python3 "${SCRIPT_DIR}/scripts/06_validate_positions.py" \
    --summary "$FINAL_SUMMARY" \
    --outdir  "$P3DIR"
echo "" >&2

# ---------------------------------------------------------------------------
# Step 07: Format output
# ---------------------------------------------------------------------------
echo "[STEP 07] Formatting output (IDS naming)..." >&2
python3 "${SCRIPT_DIR}/scripts/07_format_output.py" \
    --summary "$FINAL_SUMMARY" \
    --output  "${OUTDIR}/IDS_primers"
echo "" >&2

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
echo "============================================================" >&2
echo " Pipeline complete!" >&2
echo "============================================================" >&2
echo "  Filtered InDels:   $FILTERED" >&2
echo "  Primer design:     $P3SUMMARY" >&2
echo "  Specificity check: ${P3DIR}/specificity_report.tsv" >&2
echo "  Position check:    ${P3DIR}/validation_report.tsv" >&2
echo "  Final primers:     ${OUTDIR}/IDS_primers.tsv" >&2
echo "  Primer order list: ${OUTDIR}/IDS_primers.txt" >&2
echo "============================================================" >&2
