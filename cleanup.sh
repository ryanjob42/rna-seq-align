#!/bin/bash

to_clean=(
    .snakemake/
    alignment/
    flags/
    read_counts/
    split_fasta/
    split_genome_index/
    tmp
)

for d in ${to_clean[@]}; do
    if [ -d "$d" ]; then
        rm -rf "$d"
    fi
done
