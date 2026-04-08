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
1_Sequence_matching/
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

## 2 Single Select Module

This directory contains scripts for evaluating sequence recovery efficiency under different sequencing depths via random downsampling and read traversal analysis.


### File Structure
```text
1_Single_select/
├── run_downsampling_cellpool.sh
├── run_downsampling_ID1.sh
├── run_downsampling_ID5000.sh
├── run_downsampling_ID10000.sh
├── run_downsampling_ID20000.sh
├── run_downsampling_ID30000.sh
├── single_select_count_reads_1Indexs.py
├── single_select_count_reads_10Indexs.py
└── single_select_findTrueIndex_multiprocessing.py
```

#### Workflow

1. Random Downsampling
   
    Run one of the provided shell scripts to generate downsampled FASTQ files:

    ```text
    bash run_downsampling_ID10000.sh
    ```

    Each script produces a downsampled FASTQ dataset corresponding to a specific experimental condition or sequencing depth.

2. Parameter Configuration
   
    Before running the recovery analysis, the script
`single_select_count_reads_10Indexs.py` must be manually configured according to the downsampled dataset.

Specifically, users need to modify:

fastq_name → name of the downsampled FASTQ file
Single_ids_list → list of target reference indices
output_file_name → output file name

An example configuration is shown below:

| Dataset       | FASTQ file               | Single_ids_list                                   | Output file                                                |
|---------------|--------------------------|---------------------------------------------------|------------------------------------------------------------|
| Cell pool     | PB-222_1_depth_30x       | [1-1, 5000-1, 10000-1, 20000-1, 30000-1]         | single_select_cellpool_count_reads_for_perfect_recovery.pkl |
| ID1           | BG192_1_depth_30x        | [1-1]                                             | single_select_ID1_count_reads_for_perfect_recovery.pkl     |
| ID5000        | BG-233_1_depth_30x       | [5000-1]                                          | single_select_ID5000_count_reads_for_perfect_recovery.pkl  |
| ID10000       | BG252_1_depth_30x        | [10000-1]                                         | single_select_ID10000_count_reads_for_perfect_recovery.pkl |
| ID20000       | BG284_1_depth_30x        | [20000-1]                                         | single_select_ID20000_count_reads_for_perfect_recovery.pkl |
| ID30000       | BG-JX29_1_depth_30x      | [30000-1]                                         | single_select_ID30000_count_reads_for_perfect_recovery.pkl |

3. Recovery Analysis

    After configuring the parameters, run the appropriate script depending on the number of target sequences:

    Multiple sequences (e.g., cell pool)

    ```bash
    python single_select_count_reads_10Indexs.py
    ```

    Single sequence:

    ```text
    python single_select_count_reads_1Indexs.py \
    --fastq_name PB-222_1_depth_30x \
    --single_ids_list 0 \
    --output_file_name single_select_ID1_count_reads_for_perfect_recovery.pkl
    ```