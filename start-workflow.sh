#!/bin/bash

source ~/.bashrc
conda activate rna-seq-align
module load singularity
snakemake --scheduler greedy --groups generate_single_index=mygroup --group-components mygroup=3
