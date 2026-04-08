#!/usr/bin/env python3
"""
07_format_output.py
-------------------
Format the final primer summary into IDS (InDel Sequence) naming convention
and produce two output files:

  IDS_primers.tsv   Tab-separated table with IDS names, sequences, Tm, product size
  IDS_primers.txt   Plain text primer list (NAME  SEQUENCE) suitable for ordering

IDS naming convention:
  IDS{CHR}_{POS_KB}{STRAND}
  where:
    CHR     = zero-padded chromosome number (extracted from "chr02" → "02")
    POS_KB  = genomic position // 1000  (integer, e.g. 29891 for 29,891,678 bp)
    STRAND  = F (forward / LEFT) or R (reverse / RIGHT)

  Example:  chr02, pos 29891678  →  IDS02_29891F / IDS02_29891R

Input:
  --summary    final_summary.tsv from step 05 (or 06)
  --output     Output prefix (without extension); files .tsv and .txt are written

Usage:
  python 07_format_output.py \\
      --summary results/primer3/final_summary.tsv \\
      --output  results/IDS_primers
"""

import argparse
import csv
import os
import re
import sys


def ids_name(chrom: str, pos: str | int, strand: str) -> str:
    """
    Build an IDS marker name.

    chrom:  'chr02', 'chr06', etc.
    pos:    genomic position (1-based integer or string)
    strand: 'F' or 'R'
    """
    # Extract numeric chromosome number
    m = re.search(r'(\d+)', str(chrom))
    if not m:
        raise ValueError(f"Cannot extract chromosome number from: {chrom}")
    chr_num = m.group(1).zfill(2)

    pos_kb = int(pos) // 1000
    return f"IDS{chr_num}_{pos_kb}{strand}"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Format final primer summary into IDS naming convention.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--summary", required=True,
                        help="final_summary.tsv from step 05 or 06")
    parser.add_argument("--output", required=True,
                        help="Output prefix (without extension); .tsv and .txt will be created")
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.isfile(args.summary):
        print(f"[ERROR] Summary file not found: {args.summary}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    rows = []
    with open(args.summary, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(dict(row))

    print(f"[INFO] Loaded {len(rows)} entries.", file=sys.stderr)

    tsv_path = args.output + ".tsv"
    txt_path = args.output + ".txt"

    tsv_fields = [
        "NAME_F", "NAME_R",
        "SEQ_F", "SEQ_R",
        "CHROM", "POS", "INDEL_SIZE",
        "PRODUCT_SIZE", "TM_F", "TM_R",
        "FLANK", "PAIR_INDEX", "SPECIFICITY",
    ]

    n_written = 0
    n_skipped = 0

    with open(tsv_path, "w", newline="") as tsv_fh, open(txt_path, "w") as txt_fh:
        writer = csv.DictWriter(tsv_fh, fieldnames=tsv_fields, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        txt_fh.write("# Rice InDel Primer List — IDS format\n")
        txt_fh.write("# NAME\tSEQUENCE\tPRODUCT(bp)\tTm(°C)\tSPECIFICITY\n")

        for row in rows:
            spec = row.get("SPECIFICITY", "")
            if spec in ("NO_PRIMERS", "CRITICAL"):
                print(f"[SKIP] {row.get('CHROM')}:{row.get('POS')} — {spec}", file=sys.stderr)
                n_skipped += 1
                continue

            chrom = row["CHROM"]
            pos = row["POS"]
            left_seq = row.get("LEFT_PRIMER", "NA")
            right_seq = row.get("RIGHT_PRIMER", "NA")

            if left_seq == "NA" or right_seq == "NA":
                print(f"[SKIP] {chrom}:{pos} — missing primer sequence", file=sys.stderr)
                n_skipped += 1
                continue

            try:
                name_f = ids_name(chrom, pos, "F")
                name_r = ids_name(chrom, pos, "R")
            except ValueError as e:
                print(f"[ERROR] {e}", file=sys.stderr)
                n_skipped += 1
                continue

            out_row = {
                "NAME_F": name_f,
                "NAME_R": name_r,
                "SEQ_F": left_seq,
                "SEQ_R": right_seq,
                "CHROM": chrom,
                "POS": pos,
                "INDEL_SIZE": row.get("INDEL_SIZE", "NA"),
                "PRODUCT_SIZE": row.get("PRODUCT_SIZE", "NA"),
                "TM_F": row.get("LEFT_TM", "NA"),
                "TM_R": row.get("RIGHT_TM", "NA"),
                "FLANK": row.get("FLANK", "NA"),
                "PAIR_INDEX": row.get("PAIR_INDEX", "NA"),
                "SPECIFICITY": spec,
            }
            writer.writerow(out_row)

            product = row.get("PRODUCT_SIZE", "NA")
            tm_f = row.get("LEFT_TM", "NA")
            tm_r = row.get("RIGHT_TM", "NA")

            txt_fh.write(f"{name_f}\t{left_seq}\t{product}\t{tm_f}\t{spec}\n")
            txt_fh.write(f"{name_r}\t{right_seq}\t{product}\t{tm_r}\t{spec}\n")

            print(f"  {name_f}  {left_seq}", file=sys.stderr)
            print(f"  {name_r}  {right_seq}", file=sys.stderr)
            n_written += 1

    print(f"", file=sys.stderr)
    print(f"[SUMMARY] Markers written: {n_written}", file=sys.stderr)
    print(f"[SUMMARY] Markers skipped: {n_skipped}", file=sys.stderr)
    print(f"[INFO] TSV:  {tsv_path}", file=sys.stderr)
    print(f"[INFO] TXT:  {txt_path}", file=sys.stderr)
    print(f"[DONE] Output formatting complete.", file=sys.stderr)


if __name__ == "__main__":
    main()
