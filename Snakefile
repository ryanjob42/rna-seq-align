genome_fasta = 'data/Mus_musculus.GRCm39.dna_sm.toplevel.fa'
genome_gtf = 'data/Mus_musculus.GRCm39.110.gtf'

rule all:
    input:
        split_genome = rules.perform_split.output

rule download_faSplit:
    output: 'tools/faSplit'
    shell: 'rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faSplit ./faSplit'

rule perform_split:
    input:
        fasplit = 'tools/faSplit',
        fasta = genome_fasta
    output:
        directory('split_genome')
    shell: '{input.fasplit} byname {input.fasta:q} {output:q}'
