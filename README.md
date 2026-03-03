# PCR Deduper

## Overview
This tool is a specific implementation of a **Reference Based PCR Duplicate Remover**. It is designed to process **sorted SAM files** containing single-end sequencing data. The script identifies and removes PCR duplicates by accounting for unique molecular identifiers (UMIs) and correcting for alignment artifacts like soft-clipping.

The goal is to retain only a single copy of each unique read (based on alignment position, strand, and UMI) to reduce bias in downstream analysis.

## Repository Structure
* **`israel_deduper.py`**: The main Python script that performs the deduplication.
* **`STL96.txt`**: A text file containing the list of 96 known valid UMIs used for filtering.
* **`pseudocode_revised.txt`**: A detailed outline of the algorithm and logic used in the tool.
* **`summary.txt`**: Contains summary statistics or output details from test runs.
* **`test_files/`**: Directory containing sample data for testing the pipeline.
* **`test.sam`**: A small SAM file used for development and testing.

## Requirements
* **Python 3.12+**
* **Standard Libraries**: `argparse`, `re` (Regex)

## Usage
The script is run from the command line. It requires a sorted SAM file and a file containing known UMIs.

### General Syntax
```bash
python israel_deduper.py -f <input_sorted.sam> -o <output_deduped.sam> -u <umi_file.txt>
```

### Assumptions about the input

- **Single-end data**: The tool assumes reads are single-end (no paired-end logic).
- **UMI location**: The UMI is expected to be the **last colon-separated field** in the read name (QNAME), e.g.  
  `NS500451:154:HWKTMBGXX:1:11101:15364:1139:GAACAGGT` (UMI = `GAACAGGT`).
- **Known UMIs only**: Only reads whose UMI appears in the UMI list file (`-u/--umi`, e.g. `STL96.txt`) are considered.  
  Reads with UMIs **not** in this list are skipped and counted as "Reads with wrong UMIs" in the summary.

### Output and report

The script writes a deduplicated SAM file to the path given by `-o/--outfile`.  
In addition, it writes a simple text report alongside the SAM file:

- **Report path**: `<output_deduped.sam>.report.txt`

The report and stdout summary both include:

- Number of header lines
- Total alignments processed
- Reads with wrong UMIs
- Duplicates removed
- Unique reads written
- Unique reads per chromosome

## Test Files

The `test_files/` directory contains small SAM files for validating behavior:

- `test_files/duped.sam`: Contains intentionally duplicated reads (including cases that differ by UMI, chromosome, strand, and CIGAR) to exercise the deduplication logic.
- `test_files/dedupped.sam`: Expected output produced by running `israel_deduper.py` on `test_files/duped.sam` with `STL96.txt` as the UMI list.

You can rerun the test manually with:

```bash
python israel_deduper.py -f test_files/duped.sam -o test_files/dedupped.sam -u STL96.txt
```
