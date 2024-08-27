genome_fasta = 'data/Mus_musculus.GRCm39.dna_sm.toplevel.fa'
genome_gtf = 'data/Mus_musculus.GRCm39.110.gtf'

rule all:
    input:
        split_genome = 'split_genome'

rule download_faSplit:
    output: 'tools/faSplit'
    shell: 'rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faSplit {output}'

# This is a checkpoint rule so we can find all the files generated
# in the output directory and create a genome index from each
# as a separate rule execution.
checkpoint perform_split:
    input:
        fasplit = 'tools/faSplit',
        fasta = genome_fasta
    output:
        directory('split_genome')
    shell: '''
        mkdir -p "{output}"
        {input.fasplit} byname {input.fasta:q} "{output}/"
    '''

def get_all_split_indices(wildcards):
    # Get the directory that the "perform_split" checkpoint created.
    split_dir = checkpoints.perform_split.get(**wildcards).output[0]

    # Create a string that Snakemake understands for finding all the
    # fasta (.fa) files in that directory.
    fa_search_string = os.path.join(split_dir, '{name}.fa')

    # Use Snakemake's "glob_wildcards" to find the names of all the
    # fasta files. This just returns the "name" part in the search string.
    fa_names = glob_wildcards(fa_search_string).name

    # Use Snakemake's "expand" to return a list of all the fasta files.
    return expand(fa_search_string, name=fa_names)

rule generate_all_indexes:
    input:
        get_all_split_indices
    output:
        directory('genome_index/')

rule generate_single_index:
    input:
        split_fasta = 'split_genome/{name}.fa',
        gtf_file = genome_gtf
    output:
        directory('genome_index/{name}')
    shell: '''
        STAR \
        --runThreadN 40 \
        --runMode genomeGenerate \
        --genomeDir "GenomeIndex/{name}" \
        --genomeFastaFiles {input.split_fasta:q} \
        --sjdbGTFfile {input.gtf_file:q} \
        --sjdbOverhang 99 \
        --genomeSAindexNbases 12
    '''
