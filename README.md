# Rice InDel Primer Design Pipeline

[![samtools](https://img.shields.io/badge/samtools-1.17+-blue)](http://www.htslib.org/)
[![primer3](https://img.shields.io/badge/primer3-2.6+-green)](https://primer3.org/)
[![BLAST](https://img.shields.io/badge/BLAST+-2.13+-orange)](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
[![Python](https://img.shields.io/badge/Python-3.8+-yellow)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-lightgrey)](LICENSE)

An automated pipeline for designing PCR primers targeting InDel (insertion/deletion) markers in QTL regions of rice. Primers are validated for specificity via BLAST, checked for correct flanking of the InDel, and formatted using a standardized naming convention (IDS format).

---

## Table of Contents

- [Background](#background)
- [Pipeline Overview](#pipeline-overview)
- [Directory Structure](#directory-structure)
- [Installation](#installation)
- [Input Formats](#input-formats)
- [Usage](#usage)
  - [Quick Start (run_pipeline.sh)](#quick-start)
  - [Step-by-step Usage](#step-by-step-usage)
- [Output Formats](#output-formats)
- [BLAST Specificity Thresholds](#blast-specificity-thresholds)
- [Primer Naming Convention](#primer-naming-convention)
- [Tool References](#tool-references)

---

## Background

Insertion/deletion (InDel) polymorphisms are valuable co-dominant markers for QTL mapping and marker-assisted selection in rice. This pipeline takes whole-genome resequencing-derived InDel calls (e.g., from GATK or freebayes), filters them for suitability as PCR markers, designs flanking primers with Primer3, and verifies specificity against the full rice reference genome using BLAST.

The pipeline was developed designing InDel markers.

---

## Pipeline Overview

```
rice_reference.fa   indels.tsv   qtl_peaks.tsv
        |               |
   [Step 1]         [Step 2]
  Index genome    Filter InDels
        |               |
        +-------+-------+
                |
           [Step 3]
         Design primers
          (Primer3)
                |
           [Step 4]
        Build BLAST DB
                |
           [Step 5]
      Check specificity
          (BLAST)
                |
           [Step 6]
       Validate positions
                |
           [Step 7]
        Format output
                |
        IDS primer table
```

**Key filtering criteria:**
- InDel size 10–30 bp (PCR-detectable on agarose)
- QUAL > 500
- No homopolymer runs >= 4 nt or dinucleotide repeats in ALT sequence
- No clustering within 10 kb
- Optional minimum inter-marker spacing (default 150 kb)

---

## Directory Structure

```
rice-indel-primer-pipeline/
├── README.md                    # This file
├── environment.yml              # Conda environment specification
├── run_pipeline.sh              # Master pipeline script
├── config/
│   └── primer3_defaults.txt     # Primer3 parameter defaults
├── scripts/
│   ├── 01_index_reference.sh    # Index reference genome with samtools
│   ├── 02_filter_indels.py      # Filter InDel candidates
│   ├── 03_design_primers.sh     # Design primers with Primer3
│   ├── 04_build_blast_db.sh     # Build BLAST nucleotide database
│   ├── 05_check_specificity.py  # BLAST specificity check
│   ├── 06_validate_positions.py # Validate primer-InDel overlap
│   └── 07_format_output.py      # Format final output (IDS naming)
└── example/
    ├── example_indels.tsv       # Example InDel input (10 rows)
    └── example_qtl_peaks.tsv   # Example QTL peak positions
```

---

## Installation

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Anaconda)
- [Homebrew](https://brew.sh/) (macOS, optional for primer3)

### 1. Clone the repository

```bash
git clone https://github.com/yourlab/rice-indel-primer-pipeline.git
cd rice-indel-primer-pipeline
```

### 2. Create conda environment

```bash
conda env create -f environment.yml
conda activate rice-primer
```

This installs:
- `samtools >= 1.17`
- `blast >= 2.13`
- `primer3 >= 2.6`
- `python >= 3.8`

### 3. Verify installation

```bash
samtools --version
primer3_core --version
blastn -version
python --version
```

### Alternative: Homebrew (macOS)

```bash
brew install samtools blast
brew install primer3
```

---

## Input Formats

### Reference Genome

Standard FASTA format. Must be the same assembly used to call variants.

```
>chr01
ATCGATCGATCG...
>chr02
GCTAGCTAGCTA...
```

### InDel List (TSV)

Tab-separated, with header line starting with `#CHROM`:

```
#CHROM  POS         ID   REF    REFLEN  ALT                    ALTLEN  QUAL
chr02   28512345    .    A      1       ATCGATCGATCG           12      1523.4
chr02   28698210    .    GC     2       GCATCGATCGATCGATCG     18      892.1
```

| Column  | Description                                         |
|---------|-----------------------------------------------------|
| #CHROM  | Chromosome name (must match reference FASTA header) |
| POS     | 1-based position of the variant                     |
| ID      | Variant ID (may be `.`)                             |
| REF     | Reference allele sequence                           |
| REFLEN  | Length of REF allele                                |
| ALT     | Alternate allele sequence                           |
| ALTLEN  | Length of ALT allele                                |
| QUAL    | Variant quality score                               |

### QTL Peaks (TSV)

```
TRAIT       MARKER      CHR     POS
trait1      GBS_chr02_30930592  chr02   30930592
trait2      GBS_chr02_29898387  chr02   29898387
```

---

## Usage

### Quick Start

Run the entire pipeline with default settings:

```bash
bash run_pipeline.sh \
  --ref rice_reference.fa \
  --indels indels.tsv \
  --chrom chr02 \
  --start 28500000 \
  --end 31500000 \
  --outdir results/

# With help:
bash run_pipeline.sh --help
```

### Step-by-step Usage

#### Step 1: Index reference genome

```bash
bash scripts/01_index_reference.sh rice_reference.fa
```

**Output:** `rice_reference.fa.fai`

---

#### Step 2: Filter InDels

```bash
python scripts/02_filter_indels.py \
  --indels indels.tsv \
  --chrom chr02 \
  --start 28500000 \
  --end 31500000 \
  --qual-min 500 \
  --size-min 10 \
  --size-max 30 \
  --cluster-dist 10000 \
  --spacing 150000 \
  --output results/filtered_indels.tsv
```

| Argument         | Default  | Description                                         |
|------------------|----------|-----------------------------------------------------|
| `--indels`       | required | Input InDel TSV file                                |
| `--chrom`        | required | Chromosome to filter (e.g., `chr02`)                |
| `--start`        | required | Window start position (bp)                          |
| `--end`          | required | Window end position (bp)                            |
| `--qual-min`     | 500      | Minimum QUAL score                                  |
| `--size-min`     | 10       | Minimum InDel size (ALTLEN - 1) in bp               |
| `--size-max`     | 30       | Maximum InDel size (ALTLEN - 1) in bp               |
| `--cluster-dist` | 10000    | Min distance between candidates (bp); keep higher QUAL |
| `--spacing`      | 150000   | Minimum spacing between selected markers (bp)       |
| `--output`       | required | Output filtered TSV                                 |

**Output columns:** `CHROM, POS, INDEL_SIZE, QUAL, ALT, NOTES`

---

#### Step 3: Design primers with Primer3

```bash
bash scripts/03_design_primers.sh \
  --indels results/filtered_indels.tsv \
  --ref rice_reference.fa \
  --outdir results/primer3/ \
  --flank 200 \
  --config config/primer3_defaults.txt
```

| Argument    | Default                      | Description                             |
|-------------|------------------------------|-----------------------------------------|
| `--indels`  | required                     | Filtered InDel TSV from Step 2          |
| `--ref`     | required                     | Reference FASTA (indexed)               |
| `--outdir`  | required                     | Output directory for primer3 files      |
| `--flank`   | 200                          | Flanking sequence size (bp each side)   |
| `--config`  | config/primer3_defaults.txt  | Primer3 parameter file                  |

**Output:**
- `results/primer3/{chrom}_{pos}.primer3.in` — Primer3 input files
- `results/primer3/{chrom}_{pos}.primer3.out` — Primer3 output files
- `results/primer3/primer_summary.tsv` — Summary TSV with best primer pair per InDel

---

#### Step 4: Build BLAST database

```bash
bash scripts/04_build_blast_db.sh rice_reference.fa results/blast_db/
```

**Output:** `results/blast_db/rice_ref.*` (BLAST database files)

---

#### Step 5: Check primer specificity

```bash
python scripts/05_check_specificity.py \
  --summary results/primer3/primer_summary.tsv \
  --outdir results/primer3/ \
  --db results/blast_db/rice_ref \
  --threads 4 \
  --target-chr chr02
```

| Argument       | Default | Description                                    |
|----------------|---------|------------------------------------------------|
| `--summary`    | required | Primer summary TSV from Step 3                |
| `--outdir`     | required | Directory containing primer3 output files     |
| `--db`         | required | BLAST database path (prefix)                  |
| `--threads`    | 4        | Number of BLAST threads                       |
| `--target-chr` | required | Target chromosome (off-target = other chroms) |

**Output:**
- `results/primer3/final_summary.tsv` — Final primer pairs (after specificity rescue)
- `results/primer3/specificity_report.tsv` — Per-primer BLAST classification

---

#### Step 6: Validate primer positions vs InDel

```bash
python scripts/06_validate_positions.py \
  --summary results/primer3/final_summary.tsv \
  --outdir results/primer3/
```

**Output:** `results/primer3/validation_report.tsv`

Checks that the LEFT primer end position in the template is upstream of the InDel, and the RIGHT primer start is downstream. Reports errors for any overlapping primers.

---

#### Step 7: Format output

```bash
python scripts/07_format_output.py \
  --summary results/primer3/final_summary.tsv \
  --output results/IDS_primers
```

**Output:**
- `results/IDS_primers.tsv` — Formatted table with IDS names
- `results/IDS_primers.txt` — Plain text primer list

---

## Output Formats

### final_summary.tsv

| Column          | Description                              |
|-----------------|------------------------------------------|
| CHROM           | Chromosome                               |
| POS             | Genomic position (1-based)               |
| INDEL_SIZE      | Size of InDel in bp                      |
| LEFT_PRIMER     | Forward primer sequence (5'→3')          |
| RIGHT_PRIMER    | Reverse primer sequence (5'→3')          |
| LEFT_TM         | Tm of forward primer                     |
| RIGHT_TM        | Tm of reverse primer                     |
| PRODUCT_SIZE    | Expected PCR product size                |
| FLANK           | Flanking sequence used                   |
| PAIR_INDEX      | Primer3 pair index (0, 1, or 2)          |
| SPECIFICITY     | OK / WARN / CRITICAL                     |

### IDS_primers.tsv

| Column     | Description                                   |
|------------|-----------------------------------------------|
| NAME_F     | Forward primer name (e.g., `IDS02_29891F`)    |
| NAME_R     | Reverse primer name (e.g., `IDS02_29891R`)    |
| SEQ_F      | Forward primer sequence                       |
| SEQ_R      | Reverse primer sequence                       |
| CHROM      | Chromosome                                    |
| POS        | Genomic position                              |
| PRODUCT    | Expected product size (bp)                    |
| TM_F       | Forward Tm                                    |
| TM_R       | Reverse Tm                                    |

### IDS_primers.txt

Plain text list suitable for ordering:
```
IDS02_29891F  ATCGATCGATCGATCGATCG
IDS02_29891R  TAGCTAGCTAGCTAGCTAGC
IDS02_30451F  GCATGCATGCATGCATGCAT
IDS02_30451R  CATGCATGCATGCATGCATG
```

---

## BLAST Specificity Thresholds

Primers are classified using the following rules applied to off-target BLAST hits (hits on chromosomes other than the target chromosome):

| Class    | Criteria                                                                                  |
|----------|-------------------------------------------------------------------------------------------|
| CRITICAL | Hit with 100% identity, length coverage >= 80% of primer length, and 3' end is covered   |
| WARN     | Hit with >= 95% identity, <= 1 mismatch, and 3' end is covered                           |
| OK       | No hits meeting WARN or CRITICAL criteria                                                 |

**3' end coverage** is defined as the hit covering the last 5 bases of the primer (positions qlen-4 to qlen on the query).

**Rescue strategy for CRITICAL primers:**
1. Try alternate primer pairs (pair index 1, then 2) from Primer3 output
2. If all alternates fail: redesign with stricter parameters (`PRIMER_MIN_GC=45`, `PRIMER_MAX_POLY_X=3`, wider flank)
3. If only one primer in a pair fails: pin the good primer with `PRIMER_LEFT_INPUT` or `PRIMER_RIGHT_INPUT` and redesign the other

---

## Primer Naming Convention

Primers follow the IDS (InDel Sequence) naming convention:

```
IDS{CHR}_{POS_KB}{STRAND}
```

- `CHR`: Zero-padded chromosome number (e.g., `02`, `06`)
- `POS_KB`: Genomic position divided by 1000 (integer, e.g., `29891` for position 29,891,000 bp)
- `STRAND`: `F` (forward/LEFT) or `R` (reverse/RIGHT)

**Examples:**
- `IDS02_29891F` — chr02, position ~29.891 Mb, forward primer
- `IDS02_30451R` — chr02, position ~30.451 Mb, reverse primer

---

## Tool References

- **samtools / htslib**: Danecek P, et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008

- **Primer3**: Untergasser A, et al. (2012). Primer3 — new capabilities and interfaces. *Nucleic Acids Research*, 40(15), e115. https://doi.org/10.1093/nar/gks596

- **BLAST+**: Camacho C, et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. https://doi.org/10.1186/1471-2105-10-421

- **Rice reference genome (IRGSP-1.0)**: International Rice Genome Sequencing Project (2005). The map-based sequence of the rice genome. *Nature*, 436, 793–800. https://doi.org/10.1038/nature03895

---

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contact

For questions, please open an issue on GitHub or contact Vincent Pamugas Reyes, Graduate School of Bioagricultural Sciences, Nagoya University.
