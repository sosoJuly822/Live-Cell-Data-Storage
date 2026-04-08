# Live-Cell-Data-Storage

## 1 Sequence Matching Module

This directory contains the pipeline for reference-assisted analysis of DNA storage sequencing data, including read matching, error profiling, and downstream result processing.

### Overview

The sequence matching workflow is designed to:

Assign sequencing reads (FASTQ) to reference sequences in a predefined library
Perform multi-stage alignment (primer → index → payload/unique region)
Generate base-level error statistics (insertion, deletion, substitution, match)
Support large-scale processing via multiprocessing and cluster execution

This pipeline is intended for analysis purposes (e.g., error profiling and sequence recovery evaluation) and assumes access to the reference library.

### File Structure

```text
sequence_matching/
├── analysis_result.ipynb
├── analysis_result.py
├── cpu-fanlab-cluster_library_read1.slurm
├── FINAL.xlsx
├── result_output.sh
└── seq_match_multiprocessing.py
```



#### Core Scripts

`FINAL.xlsx` Reference library containing all designed DNA sequences used for read assignment.

`seq_match_multiprocessing.py` The main pipeline for sequence matching and error counting:

- Reads FASTQ files
- Filters reads based on quality
- Performs primer localization and sequence matching
- Asigns reads to reference sequences
- Outputs per-position error statistics as `.pkl` files

`analysis_result.py/ipynb` Scripts for analyzing the generated error dictionaries (`.pkl`), including: 

- Error rate calculation
- Distribution analysis
- Visualization

#### Execute Script

Use the SLURM script to process large datasets:

```text
sbatch cpu-fanlab-cluster_library_read1.slurm
```

Convert `.pkl` results into CSV and compute error rates:

```text
bash result_output.sh
```