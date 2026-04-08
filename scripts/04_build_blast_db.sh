#!/usr/bin/env bash
# =============================================================================
# 04_build_blast_db.sh
# Build a local BLAST nucleotide database from the reference FASTA.
#
# Usage:
#   bash 04_build_blast_db.sh <reference.fa> <db_dir>
#
# Arguments:
#   <reference.fa>   Reference FASTA file (must exist)
#   <db_dir>         Output directory for BLAST database files
#
# Output:
#   <db_dir>/rice_ref.*   BLAST database files (nhr, nin, nsq, etc.)
#
# Requirements:
#   BLAST+ >= 2.13 (makeblastdb must be in PATH)
# =============================================================================

set -euo pipefail

# --- Help ---
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 2 ]]; then
    echo "Usage: $0 <reference.fa> <db_dir>"
    echo ""
    echo "Build a BLAST nucleotide database from a reference FASTA."
    echo ""
    echo "Arguments:"
    echo "  <reference.fa>   Path to the reference genome FASTA file."
    echo "  <db_dir>         Directory where BLAST DB files will be written."
    echo ""
    echo "Output:"
    echo "  <db_dir>/rice_ref.*   BLAST database (use prefix: <db_dir>/rice_ref)"
    echo ""
    echo "Requirements:"
    echo "  makeblastdb >= 2.13  (must be in PATH)"
    exit 0
fi

REF="$1"
DBDIR="$2"

# --- Validate ---
if [[ ! -f "$REF" ]]; then
    echo "[ERROR] Reference file not found: $REF" >&2
    exit 1
fi

if ! command -v makeblastdb &>/dev/null; then
    echo "[ERROR] makeblastdb is not in PATH. Install via: conda install -c bioconda blast" >&2
    exit 1
fi

BLAST_VERSION=$(makeblastdb -version 2>&1 | head -1)
echo "[INFO] Using $BLAST_VERSION" >&2

mkdir -p "$DBDIR"

DBPREFIX="${DBDIR}/rice_ref"

echo "[INFO] Building BLAST nucleotide database..." >&2
echo "[INFO]   Input FASTA:  $REF" >&2
echo "[INFO]   DB prefix:    $DBPREFIX" >&2

makeblastdb \
    -in "$REF" \
    -dbtype nucl \
    -out "$DBPREFIX" \
    -parse_seqids \
    -title "rice_reference" \
    2>&1 | sed 's/^/[makeblastdb] /' >&2

# --- Verify ---
if [[ -f "${DBPREFIX}.nhr" || -f "${DBPREFIX}.00.nhr" ]]; then
    echo "[INFO] Database created successfully at: $DBPREFIX" >&2
else
    echo "[ERROR] Database files not found after makeblastdb — check input FASTA." >&2
    exit 1
fi

echo "[DONE] BLAST database ready. Use prefix: $DBPREFIX" >&2
