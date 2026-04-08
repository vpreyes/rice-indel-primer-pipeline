#!/usr/bin/env python3
"""
02_filter_indels.py
-------------------
Filter InDel candidates from a TSV file for use in rice QTL marker design.

Filtering criteria (applied in order):
  1. Genomic window (--chrom, --start, --end)
  2. QUAL >= --qual-min  (default 500)
  3. InDel size = ALTLEN - 1, keep --size-min to --size-max bp  (default 10-30)
  4. Exclude ALT sequences with homopolymer runs >= 4 (AAAA, TTTT, CCCC, GGGG)
  5. Exclude ALT sequences with AT or GC dinucleotide repeats (ATATAT..., GCGCGC...)
  6. De-cluster: within --cluster-dist bp, keep only the higher-QUAL candidate
  7. Apply minimum inter-marker spacing (--spacing, default 150000 bp)

Input TSV columns (tab-separated, header starts with #CHROM):
  #CHROM  POS  ID  REF  REFLEN  ALT  ALTLEN  QUAL

Output TSV columns:
  CHROM  POS  INDEL_SIZE  QUAL  ALT  NOTES

Usage:
  python 02_filter_indels.py --indels indels.tsv --chrom chr02 \\
      --start 28500000 --end 31500000 --output filtered.tsv
"""

import argparse
import re
import sys


# ---------------------------------------------------------------------------
# Filter functions
# ---------------------------------------------------------------------------

def has_homopolymer(seq, min_run=4):
    """Return True if seq contains a run of >= min_run identical nucleotides."""
    pattern = r'([ACGTacgt])\1{' + str(min_run - 1) + r',}'
    return bool(re.search(pattern, seq))


def has_dinucleotide_repeat(seq, pairs=('AT', 'TA', 'GC', 'CG'), min_repeats=3):
    """Return True if seq contains a dinucleotide repeat >= min_repeats units."""
    for dinuc in pairs:
        pattern = '(' + re.escape(dinuc) + '){' + str(min_repeats) + ',}'
        if re.search(pattern, seq, re.IGNORECASE):
            return True
    return False


def decluster(candidates, cluster_dist):
    """
    Remove clustered InDels: within cluster_dist bp of another candidate,
    keep only the one with higher QUAL. Ties broken by keeping the first.

    candidates: list of dicts with keys POS (int) and QUAL (float), sorted by POS
    Returns: filtered list of dicts
    """
    if not candidates:
        return []

    kept = [candidates[0]]
    for c in candidates[1:]:
        # Check against all already-kept candidates
        conflict = False
        for k in kept:
            if abs(c['POS'] - k['POS']) < cluster_dist:
                # Within cluster distance; keep higher QUAL
                if c['QUAL'] > k['QUAL']:
                    kept.remove(k)
                    kept.append(c)
                conflict = True
                break
        if not conflict:
            kept.append(c)

    # Sort output by POS
    kept.sort(key=lambda x: x['POS'])
    return kept


def apply_spacing(candidates, spacing):
    """
    Enforce a minimum spacing between selected markers.
    Iterate in order; skip any candidate within spacing bp of the last kept.

    candidates: list of dicts sorted by POS
    Returns: filtered list of dicts
    """
    if not candidates:
        return []
    selected = [candidates[0]]
    for c in candidates[1:]:
        if c['POS'] - selected[-1]['POS'] >= spacing:
            selected.append(c)
    return selected


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description='Filter InDel candidates for rice QTL primer design.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--indels', required=True,
                        help='Input InDel TSV file (tab-separated, header #CHROM POS ID REF REFLEN ALT ALTLEN QUAL)')
    parser.add_argument('--chrom', required=True,
                        help='Chromosome to process (e.g., chr02)')
    parser.add_argument('--start', required=True, type=int,
                        help='Window start position (bp, 1-based)')
    parser.add_argument('--end', required=True, type=int,
                        help='Window end position (bp, 1-based)')
    parser.add_argument('--qual-min', type=float, default=500.0,
                        help='Minimum QUAL score to retain')
    parser.add_argument('--size-min', type=int, default=10,
                        help='Minimum InDel size (ALTLEN - 1) in bp')
    parser.add_argument('--size-max', type=int, default=30,
                        help='Maximum InDel size (ALTLEN - 1) in bp')
    parser.add_argument('--cluster-dist', type=int, default=10000,
                        help='Minimum distance (bp) between candidates; within this distance keep higher QUAL')
    parser.add_argument('--spacing', type=int, default=150000,
                        help='Minimum spacing (bp) between selected markers after clustering')
    parser.add_argument('--output', required=True,
                        help='Output filtered TSV file')
    return parser.parse_args()


def main():
    args = parse_args()

    print(f'[INFO] Reading InDel file: {args.indels}', file=sys.stderr)

    # --- Read input ---
    candidates = []
    n_total = 0
    n_chrom = 0
    n_window = 0
    n_qual = 0
    n_size = 0
    n_homopoly = 0
    n_dinuc = 0

    with open(args.indels) as fh:
        header = None
        col_idx = {}

        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue

            # Parse header
            if line.startswith('#'):
                header = line.lstrip('#').split('\t')
                # Normalize: strip whitespace and leading #
                header = [h.strip() for h in header]
                # Map column names to indices (handle #CHROM)
                for i, h in enumerate(header):
                    col_idx[h.upper()] = i
                # Validate required columns
                required = {'CHROM', 'POS', 'REF', 'REFLEN', 'ALT', 'ALTLEN', 'QUAL'}
                found = set(col_idx.keys())
                missing = required - found
                if missing:
                    print(f'[ERROR] Missing columns in header: {missing}', file=sys.stderr)
                    print(f'[ERROR] Found columns: {found}', file=sys.stderr)
                    sys.exit(1)
                continue

            if header is None:
                print('[ERROR] No header line found before data.', file=sys.stderr)
                sys.exit(1)

            parts = line.split('\t')
            if len(parts) < len(header):
                continue  # skip malformed lines

            n_total += 1

            chrom = parts[col_idx['CHROM']].strip()
            pos = int(parts[col_idx['POS']].strip())
            ref = parts[col_idx['REF']].strip()
            alt = parts[col_idx['ALT']].strip()
            try:
                altlen = int(parts[col_idx['ALTLEN']].strip())
                qual = float(parts[col_idx['QUAL']].strip())
            except ValueError:
                continue

            indel_size = altlen - 1  # insertion size relative to REF

            # Filter 1: chromosome
            if chrom != args.chrom:
                continue
            n_chrom += 1

            # Filter 2: window
            if not (args.start <= pos <= args.end):
                continue
            n_window += 1

            # Filter 3: QUAL
            if qual < args.qual_min:
                n_qual += 1
                continue

            # Filter 4: InDel size
            if not (args.size_min <= indel_size <= args.size_max):
                n_size += 1
                continue

            # Filter 5: homopolymer in ALT
            if has_homopolymer(alt, min_run=4):
                n_homopoly += 1
                continue

            # Filter 6: dinucleotide repeat in ALT
            if has_dinucleotide_repeat(alt):
                n_dinuc += 1
                continue

            candidates.append({
                'CHROM': chrom,
                'POS': pos,
                'INDEL_SIZE': indel_size,
                'QUAL': qual,
                'ALT': alt,
                'REF': ref,
            })

    print(f'[INFO] Total variants read:          {n_total}', file=sys.stderr)
    print(f'[INFO] On chromosome {args.chrom}:       {n_chrom}', file=sys.stderr)
    print(f'[INFO] In window [{args.start}-{args.end}]: {n_window}', file=sys.stderr)
    print(f'[INFO] Failed QUAL filter:           {n_qual}', file=sys.stderr)
    print(f'[INFO] Failed size filter:           {n_size}', file=sys.stderr)
    print(f'[INFO] Failed homopolymer filter:    {n_homopoly}', file=sys.stderr)
    print(f'[INFO] Failed dinucleotide filter:   {n_dinuc}', file=sys.stderr)
    print(f'[INFO] Passed all basic filters:     {len(candidates)}', file=sys.stderr)

    # Sort by position
    candidates.sort(key=lambda x: x['POS'])

    # --- De-cluster ---
    before_decluster = len(candidates)
    candidates = decluster(candidates, args.cluster_dist)
    print(f'[INFO] After de-clustering ({args.cluster_dist} bp): {len(candidates)} '
          f'(removed {before_decluster - len(candidates)})', file=sys.stderr)

    # --- Apply spacing ---
    before_spacing = len(candidates)
    candidates = apply_spacing(candidates, args.spacing)
    print(f'[INFO] After spacing filter ({args.spacing} bp): {len(candidates)} '
          f'(removed {before_spacing - len(candidates)})', file=sys.stderr)

    # --- Write output ---
    with open(args.output, 'w') as out:
        out.write('CHROM\tPOS\tINDEL_SIZE\tQUAL\tALT\tNOTES\n')
        for c in candidates:
            notes = f'size={c["INDEL_SIZE"]}bp;qual={c["QUAL"]:.1f}'
            out.write(
                f'{c["CHROM"]}\t{c["POS"]}\t{c["INDEL_SIZE"]}\t'
                f'{c["QUAL"]:.1f}\t{c["ALT"]}\t{notes}\n'
            )

    print(f'[INFO] Wrote {len(candidates)} candidates to: {args.output}', file=sys.stderr)
    print(f'[DONE] Filtering complete.', file=sys.stderr)


if __name__ == '__main__':
    main()
