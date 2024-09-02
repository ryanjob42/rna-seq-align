with open(snakemake.output[0], 'w') as out_file:
    out_file.writelines(idx + '\n' for idx in snakemake.input.indexes)
