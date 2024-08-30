import subprocess

command = [
    'STAR',
    '--runMode', 'genomeGenerate',
    '--genomeDir', f'{snakemake.params.genome_index_dir}/split_{snakemake.wildcards.number}.fa',
    '--genomeFastaFiles', snakemake.input.fasta,
    '--sjdbGTFfile', snakemake.input.gtf,
    '--sjdbOverhang', '99',
    '--genomeSAindexNbases', '12'
]

subprocess.run(command)
exit(0)
