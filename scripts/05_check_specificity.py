#!/usr/bin/env python3
"""
05_check_specificity.py
-----------------------
Check PCR primer specificity by BLASTing all primers against the full
reference genome and classifying off-target hits.

Classification tiers (applied to hits on chromosomes OTHER than --target-chr):

  CRITICAL  pident == 100.0  AND  length / qlen >= 0.80  AND  3'-end covered
  WARN      pident >= 95.0   AND  mismatches <= 1          AND  3'-end covered
  OK        no hits meeting WARN or CRITICAL criteria

3'-end coverage: the BLAST alignment covers the last 5 bases of the primer
  (qend >= qlen - 4 in 1-based coordinates).

Rescue strategy for CRITICAL primers:
  1. Try alternate pair indices 1 and 2 from the primer3 output file.
  2. If all alternates also CRITICAL: redesign flag is set (manual step).

Input:
  --summary     Primer summary TSV from step 03 (primer_summary.tsv)
  --outdir      Directory containing *.primer3.out files from step 03
  --db          BLAST DB prefix (from step 04)
  --threads     BLAST threads (default 4)
  --target-chr  Target chromosome; off-target = all other chromosomes

Output:
  <outdir>/specificity_report.tsv    Per-primer BLAST classification
  <outdir>/final_summary.tsv         Updated summary with SPECIFICITY column

Usage:
  python 05_check_specificity.py \\
      --summary results/primer3/primer_summary.tsv \\
      --outdir  results/primer3/ \\
      --db      results/blast_db/rice_ref \\
      --target-chr chr02
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile


# ---------------------------------------------------------------------------
# BLAST helpers
# ---------------------------------------------------------------------------

BLAST_FMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
BLAST_FIELDS = ["qseqid", "sseqid", "pident", "length", "mismatch",
                "gapopen", "qstart", "qend", "sstart", "send",
                "evalue", "bitscore", "qlen"]


def run_blastn(fa_path: str, db: str, threads: int) -> list[dict]:
    """Run blastn-short and return parsed hits as list of dicts."""
    cmd = [
        "blastn",
        "-task", "blastn-short",
        "-query", fa_path,
        "-db", db,
        "-outfmt", BLAST_FMT,
        "-num_threads", str(threads),
        "-word_size", "7",
        "-evalue", "1000",
        "-perc_identity", "80",
        "-dust", "no",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] BLAST failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    hits = []
    for line in result.stdout.splitlines():
        parts = line.strip().split("\t")
        if len(parts) != len(BLAST_FIELDS):
            continue
        hit = dict(zip(BLAST_FIELDS, parts))
        hit["pident"] = float(hit["pident"])
        hit["length"] = int(hit["length"])
        hit["mismatch"] = int(hit["mismatch"])
        hit["qstart"] = int(hit["qstart"])
        hit["qend"] = int(hit["qend"])
        hit["qlen"] = int(hit["qlen"])
        hits.append(hit)
    return hits


def three_prime_covered(hit: dict) -> bool:
    """Return True if the hit covers the last 5 bases of the query (1-based)."""
    return hit["qend"] >= hit["qlen"] - 4


def classify_hit(hit: dict, target_chr: str) -> str | None:
    """
    Classify a single BLAST hit.
    Returns 'CRITICAL', 'WARN', or None (not concerning).
    Only considers off-target chromosomes.
    """
    if hit["sseqid"] == target_chr:
        return None  # on-target, expected

    cov = hit["length"] / hit["qlen"]
    tp = three_prime_covered(hit)

    if hit["pident"] == 100.0 and cov >= 0.80 and tp:
        return "CRITICAL"
    if hit["pident"] >= 95.0 and hit["mismatch"] <= 1 and tp:
        return "WARN"
    return None


def classify_primer(primer_hits: list[dict], target_chr: str) -> str:
    """Return the worst classification across all hits for a single primer."""
    worst = "OK"
    for h in primer_hits:
        c = classify_hit(h, target_chr)
        if c == "CRITICAL":
            return "CRITICAL"
        if c == "WARN":
            worst = "WARN"
    return worst


# ---------------------------------------------------------------------------
# Primer3 output parser
# ---------------------------------------------------------------------------

def parse_primer3_pair(p3out_path: str, pair_idx: int) -> dict | None:
    """
    Extract primer sequences and positions for a given pair index.
    Returns dict with keys: LEFT, RIGHT, LEFT_TM, RIGHT_TM, LEFT_POS,
    RIGHT_POS, PRODUCT_SIZE, or None if pair does not exist.
    """
    if not os.path.isfile(p3out_path):
        return None

    data = {}
    with open(p3out_path) as fh:
        for line in fh:
            line = line.strip()
            for key in [
                f"PRIMER_LEFT_{pair_idx}_SEQUENCE",
                f"PRIMER_RIGHT_{pair_idx}_SEQUENCE",
                f"PRIMER_LEFT_{pair_idx}_TM",
                f"PRIMER_RIGHT_{pair_idx}_TM",
                f"PRIMER_LEFT_{pair_idx}",
                f"PRIMER_RIGHT_{pair_idx}",
                f"PRIMER_PAIR_{pair_idx}_PRODUCT_SIZE",
            ]:
                if line.startswith(key + "="):
                    data[key] = line.split("=", 1)[1]

    left = data.get(f"PRIMER_LEFT_{pair_idx}_SEQUENCE")
    right = data.get(f"PRIMER_RIGHT_{pair_idx}_SEQUENCE")
    if not left or not right:
        return None

    return {
        "LEFT": left,
        "RIGHT": right,
        "LEFT_TM": data.get(f"PRIMER_LEFT_{pair_idx}_TM", "NA"),
        "RIGHT_TM": data.get(f"PRIMER_RIGHT_{pair_idx}_TM", "NA"),
        "LEFT_POS": data.get(f"PRIMER_LEFT_{pair_idx}", "NA"),
        "RIGHT_POS": data.get(f"PRIMER_RIGHT_{pair_idx}", "NA"),
        "PRODUCT_SIZE": data.get(f"PRIMER_PAIR_{pair_idx}_PRODUCT_SIZE", "NA"),
        "PAIR_INDEX": pair_idx,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Check PCR primer specificity via BLAST.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--summary", required=True,
                        help="Primer summary TSV from step 03 (primer_summary.tsv)")
    parser.add_argument("--outdir", required=True,
                        help="Directory containing *.primer3.out files")
    parser.add_argument("--db", required=True,
                        help="BLAST database prefix (from step 04)")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of BLAST threads")
    parser.add_argument("--target-chr", required=True,
                        help="Target chromosome (e.g., chr02); off-target = all others")
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.isfile(args.summary):
        print(f"[ERROR] Summary file not found: {args.summary}", file=sys.stderr)
        sys.exit(1)

    # --- Read primer summary ---
    rows = []
    with open(args.summary, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(dict(row))

    if not rows:
        print("[ERROR] Summary file is empty.", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Loaded {len(rows)} InDel entries from summary.", file=sys.stderr)

    # --- Build FASTA of all primers ---
    all_primers = {}  # primer_id -> sequence
    for row in rows:
        if row.get("STATUS") == "NO_PRIMERS":
            continue
        chrom = row["CHROM"]
        pos = row["POS"]
        left = row.get("LEFT_PRIMER", "NA")
        right = row.get("RIGHT_PRIMER", "NA")
        if left and left != "NA":
            pid = f"{chrom}_{pos}_LEFT"
            all_primers[pid] = left
        if right and right != "NA":
            pid = f"{chrom}_{pos}_RIGHT"
            all_primers[pid] = right

    if not all_primers:
        print("[ERROR] No primer sequences found in summary.", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Running BLAST for {len(all_primers)} primer sequences...", file=sys.stderr)

    # Write temp FASTA
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
        for pid, seq in all_primers.items():
            tmp.write(f">{pid}\n{seq}\n")
        tmp_fa = tmp.name

    # --- Run BLAST ---
    hits = run_blastn(tmp_fa, args.db, args.threads)
    os.unlink(tmp_fa)

    # Group hits by primer ID
    hits_by_primer: dict[str, list] = {pid: [] for pid in all_primers}
    for h in hits:
        qid = h["qseqid"]
        if qid in hits_by_primer:
            hits_by_primer[qid].append(h)

    # Classify each primer
    primer_class: dict[str, str] = {}
    for pid in all_primers:
        primer_class[pid] = classify_primer(hits_by_primer[pid], args.target_chr)

    # --- Write specificity report ---
    spec_report = os.path.join(args.outdir, "specificity_report.tsv")
    with open(spec_report, "w") as out:
        out.write("PRIMER_ID\tSEQUENCE\tCLASS\n")
        for pid, seq in all_primers.items():
            cls = primer_class.get(pid, "OK")
            out.write(f"{pid}\t{seq}\t{cls}\n")
    print(f"[INFO] Specificity report: {spec_report}", file=sys.stderr)

    # --- Rescue CRITICAL pairs and write final_summary.tsv ---
    final_rows = []
    n_ok = 0
    n_warn = 0
    n_critical = 0

    for row in rows:
        if row.get("STATUS") == "NO_PRIMERS":
            row["SPECIFICITY"] = "NO_PRIMERS"
            final_rows.append(row)
            continue

        chrom = row["CHROM"]
        pos = row["POS"]
        left_id = f"{chrom}_{pos}_LEFT"
        right_id = f"{chrom}_{pos}_RIGHT"

        lc = primer_class.get(left_id, "OK")
        rc = primer_class.get(right_id, "OK")

        pair_spec = "OK"
        if lc == "CRITICAL" or rc == "CRITICAL":
            pair_spec = "CRITICAL"
        elif lc == "WARN" or rc == "WARN":
            pair_spec = "WARN"

        # Rescue: try alternate pairs if CRITICAL
        if pair_spec == "CRITICAL":
            p3out = os.path.join(args.outdir, f"{chrom}_{pos}.primer3.out")
            for alt_idx in [1, 2]:
                alt_pair = parse_primer3_pair(p3out, alt_idx)
                if alt_pair is None:
                    continue

                # Write temp FASTA for alt pair
                with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
                    tmp.write(f">{left_id}_alt{alt_idx}\n{alt_pair['LEFT']}\n")
                    tmp.write(f">{right_id}_alt{alt_idx}\n{alt_pair['RIGHT']}\n")
                    tmp_fa = tmp.name

                alt_hits = run_blastn(tmp_fa, args.db, args.threads)
                os.unlink(tmp_fa)

                alt_hits_by = {f"{left_id}_alt{alt_idx}": [], f"{right_id}_alt{alt_idx}": []}
                for h in alt_hits:
                    if h["qseqid"] in alt_hits_by:
                        alt_hits_by[h["qseqid"]].append(h)

                alt_lc = classify_primer(alt_hits_by[f"{left_id}_alt{alt_idx}"], args.target_chr)
                alt_rc = classify_primer(alt_hits_by[f"{right_id}_alt{alt_idx}"], args.target_chr)

                alt_spec = "OK"
                if alt_lc == "CRITICAL" or alt_rc == "CRITICAL":
                    alt_spec = "CRITICAL"
                elif alt_lc == "WARN" or alt_rc == "WARN":
                    alt_spec = "WARN"

                if alt_spec != "CRITICAL":
                    print(f"[INFO] {chrom}:{pos} — rescued with pair index {alt_idx} ({alt_spec})", file=sys.stderr)
                    # Update row with rescued pair
                    row["LEFT_PRIMER"] = alt_pair["LEFT"]
                    row["RIGHT_PRIMER"] = alt_pair["RIGHT"]
                    row["LEFT_TM"] = alt_pair["LEFT_TM"]
                    row["RIGHT_TM"] = alt_pair["RIGHT_TM"]
                    row["LEFT_POS"] = alt_pair["LEFT_POS"]
                    row["RIGHT_POS"] = alt_pair["RIGHT_POS"]
                    row["PRODUCT_SIZE"] = alt_pair["PRODUCT_SIZE"]
                    row["PAIR_INDEX"] = str(alt_idx)
                    pair_spec = alt_spec
                    break
            else:
                print(f"[WARN] {chrom}:{pos} — all pairs CRITICAL; manual redesign needed.", file=sys.stderr)

        row["SPECIFICITY"] = pair_spec

        if pair_spec == "OK":
            n_ok += 1
        elif pair_spec == "WARN":
            n_warn += 1
        else:
            n_critical += 1

        final_rows.append(row)

    # --- Write final summary ---
    final_summary = os.path.join(args.outdir, "final_summary.tsv")
    if final_rows:
        # Build column list: original + SPECIFICITY
        fieldnames = list(final_rows[0].keys())
        if "SPECIFICITY" not in fieldnames:
            fieldnames.append("SPECIFICITY")

        with open(final_summary, "w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t",
                                    extrasaction="ignore")
            writer.writeheader()
            writer.writerows(final_rows)

    print(f"", file=sys.stderr)
    print(f"[SUMMARY] OK:       {n_ok}", file=sys.stderr)
    print(f"[SUMMARY] WARN:     {n_warn}", file=sys.stderr)
    print(f"[SUMMARY] CRITICAL: {n_critical}", file=sys.stderr)
    print(f"[INFO] Final summary: {final_summary}", file=sys.stderr)
    print(f"[DONE] Specificity check complete.", file=sys.stderr)


if __name__ == "__main__":
    main()
