#!/bin/bash

source ~/.bashrc
conda activate rna-seq-align
module load singularity
snakemake --scheduler greedy
