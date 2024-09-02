GENOME_FASTA = 'data/Mus_musculus.GRCm39.dna_sm.toplevel.fa'
GENOME_GTF = 'data/Mus_musculus.GRCm39.110.gtf'
SPLIT_COUNT = 20
SPLIT_FASTA_DIR = 'split_fasta'
SPLIT_GENOME_INDEX_DIR = 'split_genome_index'
FASTQ_DIR = 'fastq'
ALIGNMENT_DIR = 'alignment'

# STAR nubmers each split 0, 1, 2, etc. However, they wil be 0-padded to the longest length.
# For example, if there are 11 splits (i.e., 0 to 10), they will be 00, 01, 02, ..., 09, 10.
SPLIT_NUMBERS = [str(i).zfill(len(str(SPLIT_COUNT-1))) for i in range(SPLIT_COUNT)]

THREAD_COUNT = 128
RAM_AMOUNT = 60000000000

if SPLIT_COUNT < 1:
    print(f'The split count must be at least 1, but was: {SPLIT_COUNT}')
    exit(1)

# Runs everything.
rule all:
    input:
        'split_index_list.txt'

# Downloads the faSplit tool.
rule download_faSplit:
    output: 'tools/faSplit'
    shell: 'rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faSplit {output}'

# Splits single FASTA file into pieces.
# This will be skipped if the user asked for only a single split.
rule split_fasta:
    input:
        fasplit = 'tools/faSplit',
        fasta = GENOME_FASTA
    output:
        expand(f'{SPLIT_FASTA_DIR}/split_{{number}}.fa', number=SPLIT_NUMBERS)
    shell: f'''
        mkdir -p "{SPLIT_FASTA_DIR}"
        "{{input.fasplit}}" sequence "{{input.fasta:}}" {SPLIT_COUNT} "{{SPLIT_FASTA_DIR}}/split_"
    '''

# Depending on how the splits are performed, some of the generated FASTA files
# may fail to generate an index. Per Reema, this is OK and we can just ignore those.
# This doesn't directly work with Snakemake however, as it requires the outputs actually
# get created. To work around this, trying to generate an index will output a "flag" file
# that indicates that it tried to generate the index. This rule will gather all those,
# and is set up as a "checkpoint" so Snakemake will detect which indices are here.
# Also, we use "ensure" to indicate that the directory can't be empty.
checkpoint generate_all_indexes:
    input: expand(f'{SPLIT_GENOME_INDEX_DIR}/split_index_{{number}}.done', number=SPLIT_NUMBERS)
    output: temp(touch('all_generated.txt'))

# If there's only supposed to be a single split, there's no need to actually do the split.
# In that case, we just use the full FASTA file the user gave.
# Since generating an index may fail (which is OK per Reema), this rule outputs a "flag" file
# indicating it tried to generate the index. Additionally, the command ends with "exit 0" to
# make sure STAR's exit code won't cause Snakemake to fail the execution.
rule generate_single_index:
    input:
        fasta = branch(SPLIT_COUNT > 1, f'split_fasta/split_{{number}}.fa', GENOME_FASTA),
        gtf = GENOME_GTF
    params:
        genome_index_dir = SPLIT_GENOME_INDEX_DIR
    output: touch(f'{SPLIT_GENOME_INDEX_DIR}/split_index_{{number}}.done')
    script: 'scripts/generate_single_index.py'

# This function will find all of the genome indexes which were actually created
# by all the "generate_single_index" rule executions. To make sure this behaves
# with Snakemake's planning, this uses their suggested approach for finding what
# got created. See the link below.
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution
def all_generated_indexes(wildcards):
    # Get the specific execution of the "generate_all_indexes" rule.
    generate_indexes_rule_execution = checkpoints.generate_all_indexes.get(**wildcards)

    # Get the output from that rule. This will be the genome index directory.
    index_dir = generate_indexes_rule_execution.output[0]

    # Use Snakemake's "glob_wildcards" to get the set of split numbers which
    # successfully created a genome index.
    search_string = os.path.join(index_dir, "split_{number}.fa")
    split_numbers = glob_wildcards(search_string).number

    # Use Snakemake's "expand" to create the list of discovered indexes.
    return expand(search_string, number=split_numbers)

rule split_index_list:
    input:
        flag = 'all_generated.txt',
        indexes = all_generated_indexes
    output: 'split_index_list.txt'
    script: 'scripts/split_index_list.py'
