#!/usr/bin/env python3
"""
06_validate_positions.py
------------------------
Validate that each primer pair correctly flanks its target InDel — i.e.,
the LEFT primer ends BEFORE the InDel position in the template, and the
RIGHT primer starts AFTER it.

Position logic (using primer3 coordinate conventions, 0-based):
  LEFT  position string: "start,length"  → end_0based = start + length - 1
  RIGHT position string: "start,length"  → start is the rightmost 0-based position
                                            r_start_0based = start - length + 1
  Template: extracted as CHROM:(pos - flank)-(pos + flank), so InDel is at
            index FLANK (0-based) in the template.

  PASS:  left_end_0based  < flank   AND   r_start_0based > flank
  FAIL:  left or right primer overlaps the InDel position.

Input:
  --summary    final_summary.tsv from step 05 (must have LEFT_POS, RIGHT_POS, FLANK)

Output:
  <outdir>/validation_report.tsv    Per-marker validation result

Usage:
  python 06_validate_positions.py \\
      --summary results/primer3/final_summary.tsv \\
      --outdir  results/primer3/
"""

import argparse
import csv
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate that primers flank (do not overlap) their InDel.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--summary", required=True,
                        help="final_summary.tsv from step 05")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for validation report")
    return parser.parse_args()


def validate_row(row: dict) -> tuple[str, str]:
    """
    Return (status, message) where status is 'PASS', 'FAIL', or 'SKIP'.
    """
    if row.get("SPECIFICITY") == "NO_PRIMERS" or row.get("STATUS") == "NO_PRIMERS":
        return "SKIP", "No primers designed"

    left_pos_str = row.get("LEFT_POS", "NA")
    right_pos_str = row.get("RIGHT_POS", "NA")
    flank_str = row.get("FLANK", "NA")

    if "NA" in (left_pos_str, right_pos_str, flank_str):
        return "SKIP", f"Missing position data (LEFT_POS={left_pos_str}, RIGHT_POS={right_pos_str}, FLANK={flank_str})"

    try:
        flank = int(flank_str)

        # LEFT: "start,length" (0-based start)
        l_parts = left_pos_str.split(",")
        l_start = int(l_parts[0])
        l_len = int(l_parts[1])
        l_end = l_start + l_len - 1  # 0-based inclusive end

        # RIGHT: "start,length" where start = rightmost 0-based coordinate
        r_parts = right_pos_str.split(",")
        r_right = int(r_parts[0])   # rightmost position (0-based)
        r_len = int(r_parts[1])
        r_start = r_right - r_len + 1  # leftmost 0-based start

    except (ValueError, IndexError) as e:
        return "SKIP", f"Could not parse position: {e}"

    msgs = []
    passed = True

    if l_end >= flank:
        msgs.append(f"LEFT primer overlaps InDel (left_end={l_end} >= flank={flank})")
        passed = False

    if r_start <= flank:
        msgs.append(f"RIGHT primer overlaps InDel (right_start={r_start} <= flank={flank})")
        passed = False

    if passed:
        return "PASS", (
            f"LEFT end={l_end} < InDel@{flank}; "
            f"RIGHT start={r_start} > InDel@{flank}"
        )
    else:
        return "FAIL", "; ".join(msgs)


def main():
    args = parse_args()

    if not os.path.isfile(args.summary):
        print(f"[ERROR] Summary file not found: {args.summary}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    rows = []
    with open(args.summary, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(dict(row))

    print(f"[INFO] Validating positions for {len(rows)} entries...", file=sys.stderr)

    report_path = os.path.join(args.outdir, "validation_report.tsv")
    n_pass = 0
    n_fail = 0
    n_skip = 0

    with open(report_path, "w") as out:
        out.write("CHROM\tPOS\tFLANK\tLEFT_POS\tRIGHT_POS\tVALIDATION\tNOTES\n")
        for row in rows:
            status, msg = validate_row(row)
            chrom = row.get("CHROM", "?")
            pos = row.get("POS", "?")
            flank = row.get("FLANK", "?")
            left_pos = row.get("LEFT_POS", "NA")
            right_pos = row.get("RIGHT_POS", "NA")
            out.write(f"{chrom}\t{pos}\t{flank}\t{left_pos}\t{right_pos}\t{status}\t{msg}\n")

            if status == "PASS":
                n_pass += 1
                print(f"  [PASS] {chrom}:{pos}", file=sys.stderr)
            elif status == "FAIL":
                n_fail += 1
                print(f"  [FAIL] {chrom}:{pos} — {msg}", file=sys.stderr)
            else:
                n_skip += 1

    print(f"", file=sys.stderr)
    print(f"[SUMMARY] PASS: {n_pass}", file=sys.stderr)
    print(f"[SUMMARY] FAIL: {n_fail}", file=sys.stderr)
    print(f"[SUMMARY] SKIP: {n_skip}", file=sys.stderr)
    print(f"[INFO] Validation report: {report_path}", file=sys.stderr)

    if n_fail > 0:
        print(f"[WARN] {n_fail} marker(s) failed position validation — review before ordering primers.", file=sys.stderr)

    print(f"[DONE] Position validation complete.", file=sys.stderr)


if __name__ == "__main__":
    main()
