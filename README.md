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
