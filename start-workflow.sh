#!/bin/bash

source ~/.bashrc
conda activate rna-seq-align
module load singularity
module load r/4.1.0
snakemake --scheduler greedy
