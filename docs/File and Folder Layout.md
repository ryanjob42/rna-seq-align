# File and Folder Layout
This document describes (at a high level) what files and folders exist here.
This is copy/pasted from an email with Reema, and needs additional work to be considered complete.
For now though, it is still useful information.

## config.yaml
The "config.yaml" file is what you can edit to point Snakemake at the files you want to process and control where the results get output to.
These must be on Slipstick.

## start-workflow.sh
The "start-workflow.sh" file is a script you can run to start the pipeline.

## Snakefile
The "Snakefile" file is what controls the tasks that need to be done. The high-level steps are:

1. "download_fasta" rule:
   1. Download the "faSplit" tool if needed.
2. "split_fasta" rule:
   1. If the user wants to split (i.e., the "split count" is set to a number greater than 1), the FASTA file for the genome will be split into that many pieces.
   2. If the "split count" is set to 1, then this will be skipped.
3. "generate_single_index" rule and script:
   1. Uses STAR to generate a genome index folder for each of the FASTA splits.
   2. If there is no splitting (the "split count" is 1), then it will use the full FASTA file provided.
4. "generate_all_indexes" checkpoint:
   1. This is just to make sure that all the indexes are generated.
   2. It also helps me handle the case where a genome index doesn't generate because of the "no valid exon lines" message we discussed previously.
5. "perform_single_alignment" rule and script:
   1. This takes a single input FASTQ file and a single split genome index, then runs STAR on them.
   2. Currently, it only handles single-ended reads, but I'll fix that soon to handle pair-ended as well.
6. "perform_all_alignments" rule:
   1. This is just to make sure that all the alignments are completed.
7. "compute_split_read_counts" rule and script:
   1. This will compute the read counts for all aligned .bam files for a single split of the genome index.
8. "compute_all_read_counts" rule:
   1. This makes sure all the read counts are computed.
9.  "all" rule:
    1.  This makes sure everything gets done.

## scripts/
The "scripts" folder holds all of the Python and R scripts for the more complicated processes:

- Generating a genome index from a FASTA file.
- Aligning a single-ended FASTQ file using a single genome index (I'll update this soon to allow pair-ended reads as well).
- Computing the read counts for a set of BAM files.

## profiles/default/config.yaml
The "profiles/default/config.yaml" file provides settings for how to run this on Slipstick. It should work on Riviera as well, but I haven't had time to test that for sure yet.
