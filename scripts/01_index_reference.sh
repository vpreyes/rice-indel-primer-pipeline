#!/usr/bin/env bash
# =============================================================================
# 01_index_reference.sh
# Index a reference genome FASTA file using samtools faidx.
#
# Usage:
#   bash 01_index_reference.sh <reference.fa>
#
# Output:
#   <reference.fa>.fai   (samtools FASTA index)
#
# Requirements:
#   samtools >= 1.17
# =============================================================================

set -euo pipefail

# --- Help ---
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -eq 0 ]]; then
    echo "Usage: $0 <reference.fa>"
    echo ""
    echo "Index a reference genome FASTA with samtools faidx."
    echo ""
    echo "Arguments:"
    echo "  <reference.fa>   Path to the reference genome FASTA file."
    echo ""
    echo "Output:"
    echo "  <reference.fa>.fai   Samtools FASTA index (created in same directory as input)."
    echo ""
    echo "Requirements:"
    echo "  samtools >= 1.17  (must be in PATH)"
    exit 0
fi

REF="$1"

# --- Validate input ---
if [[ ! -f "$REF" ]]; then
    echo "[ERROR] Reference file not found: $REF" >&2
    exit 1
fi

# --- Check samtools ---
if ! command -v samtools &>/dev/null; then
    echo "[ERROR] samtools is not in PATH. Install via: conda install -c bioconda samtools" >&2
    exit 1
fi

SAMTOOLS_VERSION=$(samtools --version | head -1)
echo "[INFO] Using $SAMTOOLS_VERSION" >&2

# --- Index ---
echo "[INFO] Indexing reference: $REF" >&2
samtools faidx "$REF"

INDEX="${REF}.fai"
if [[ -f "$INDEX" ]]; then
    NSEQ=$(wc -l < "$INDEX")
    echo "[INFO] Index created: $INDEX" >&2
    echo "[INFO] Number of sequences indexed: $NSEQ" >&2
else
    echo "[ERROR] Index was not created. Check that the FASTA is valid." >&2
    exit 1
fi

echo "[DONE] Reference indexed successfully." >&2
