##############################################
# Configuration
##############################################

# Read the configuration file and create some shorthand variables
# to make it easier to use, and so there is only one place to make
# changes if the names of the configurations change.
configfile: 'config.yaml'

GENOME_FASTA = config['genome-fasta-file']
GENOME_GTF = config['genome-gtf-file']
FASTQ_INPUT_FILE_PATH = config['fastq-input-file']
SPLIT_COUNT = config['split-count']
SPLIT_FASTA_DIR = config['split-fasta-output-directory']
SPLIT_GENOME_INDEX_DIR = config['genome-index-output-directory']
STAR_GENOME_INDEX_THREADS = config['star-genome-index-thread-count']
ALIGNMENT_DIR = config['alignment-output-directory']
STAR_ALIGNMENT_THREADS = config['star-alignment-thread-count']
STAR_ALIGNMENT_RAM = config['star-alignment-bam-sort-ram']
READ_COUNTS_DIR = config['read-counts-output-directory']

# STAR nubmers each split 0, 1, 2, etc. However, they wil be 0-padded to the longest length.
# For example, if there are 11 splits (i.e., 0 to 10), they will be 00, 01, 02, ..., 09, 10.
SPLIT_NUMBERS = [str(i).zfill(len(str(SPLIT_COUNT-1))) for i in range(SPLIT_COUNT)]

if SPLIT_COUNT < 1:
    print(f'The split count must be at least 1, but was: {SPLIT_COUNT}')
    exit(1)

# Parse the FASTQ inputs file, which is a tab-separated file with two columns.
# The first column is called "Name", and contains a unique ID for a sample.
# The second column is "Fastq", and contains a comma-separated list of FASTQ files
# for the same sample. If there are multiple, we assume it's a pair-ended read.
# Otherwise, we assume it's a single-ended read.
# From this data, we want a dictionary mapping each sample ID to its list of FASTQ files
# and separate lists of the pair-ended and single-ended sample IDs.
FASTQ_FILES = {}
PAIR_ENDED_SAMPLES = []
SINGLE_ENDED_SAMPLES = []

import csv
with open(FASTQ_INPUT_FILE_PATH) as fastq_input_file:
    reader = csv.DictReader(fastq_input_file, delimiter='\t')
    for line in reader:
        sample_id = line['Name']
        sample_fastq_files = line['Fastq'].split(',')
        FASTQ_FILES[sample_id] = sample_fastq_files
        if len(sample_fastq_files) > 1:
            PAIR_ENDED_SAMPLES.append(sample_id)
        else:
            SINGLE_ENDED_SAMPLES.append(sample_id)

# Read counts must be computed separately for pair-ended and single-ended alignments.
# However, we're not guaranteed to have both.
# To make it easier for us to define the rule that computes all the counts,
# we'll create a list which contains "PE" if we have any pair-ended reads
# and "SE" if we have any single-ended reads (or both if we have both).
# This way, we can use Snakemake's "expand" function to make the rules easier to read.
READ_COUNT_ENDEDNESS = []
if len(PAIR_ENDED_SAMPLES) > 0:
    READ_COUNT_ENDEDNESS.append('PE')
if len(SINGLE_ENDED_SAMPLES) > 0:
    READ_COUNT_ENDEDNESS.append('SE')


##############################################
# Top-Level Rule
##############################################

# Runs everything.
rule all:
    input:
        'all_counts.txt'


##############################################
# FASTA Splitting
##############################################

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


##############################################
# Genome Index Creation
##############################################

# Depending on how the splits are performed, some of the generated FASTA files
# may fail to generate an index. Per Reema, this is OK and we can just ignore those.
# This doesn't directly work with Snakemake however, as it requires the outputs actually
# get created. To work around this, trying to generate an index will output a "flag" file
# that indicates that it tried to generate the index. This rule will gather all those,
# and is set up as a "checkpoint" so Snakemake will detect which indices are here.
# Also, we use "ensure" to indicate that the directory can't be empty.
checkpoint generate_all_indexes:
    input: expand(f'flags/split_index_{{number}}.done', number=SPLIT_NUMBERS)
    output:
        flag = temp(touch('all_generated.txt')),
        out_dir = directory(SPLIT_GENOME_INDEX_DIR)

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
        star_threads = STAR_GENOME_INDEX_THREADS,
        genome_index_path = os.path.join(SPLIT_GENOME_INDEX_DIR, f'split_{{number}}.fa')
    output: touch(f'flags/split_index_{{number}}.done')
    script: 'scripts/generate_single_index.py'


##############################################
# Sequence Alignment
##############################################

# This function will find all of the genome indexes which were actually created
# by all the "generate_single_index" rule executions. To make sure this behaves
# with Snakemake's planning, this uses their suggested approach for finding what
# got created. See the link below.
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution
def successful_index_generation_numbers(wildcards):
    search_string = os.path.join(SPLIT_GENOME_INDEX_DIR, f'split_{{number}}.fa')
    return glob_wildcards(search_string).number

# Ensures all of the alignments are completed.
# We need one BAM file for all combinations of split numbers and sample IDs.
rule perform_all_alignments:
    input:
        generation_complete_flag = ancient('all_generated.txt'),
        pair_ended_alignments = expand(
            f'{ALIGNMENT_DIR}/split_{{number}}/{{sample_id}}_PE_Aligned.sortedByCoord.out.bam',
            number=successful_index_generation_numbers,
            sample_id=PAIR_ENDED_SAMPLES),
        single_ended_alignments = expand(
            f'{ALIGNMENT_DIR}/split_{{number}}/{{sample_id}}_SE_Aligned.sortedByCoord.out.bam',
            number=successful_index_generation_numbers,
            sample_id=SINGLE_ENDED_SAMPLES),
    output: temp(touch('all_aligned.txt'))

# Performs a single alignment.
# The list of FASTQ files needed is determined by looking up all the file names
# from the "FASTQ_FILES" dictionary for the associated sample ID.
# The way we call the script doesn't change if it's a pair-ended or single-ended read,
# just the name of the output file, which we'll use a wildcard to represent.
rule perform_single_alignment:
    input:
        genome_index = f'{SPLIT_GENOME_INDEX_DIR}/split_{{number}}.fa',
        fastq = lookup(dpath=f'{{sample_id}}', within=FASTQ_FILES)
    params:
        star_threads = STAR_ALIGNMENT_THREADS,
        star_memory = STAR_ALIGNMENT_RAM
    output:
        f'{ALIGNMENT_DIR}/split_{{number}}/{{sample_id}}_{{endedness}}_Aligned.sortedByCoord.out.bam'
    script: 'scripts/perform_single_alignment.py'


##############################################
# Read Counts and Stats
##############################################

# Computes all read counts for all splits.
# This needs to be done separately for each split of our genome index,
# and separately for pair-ended and single-ended reads sper split.
rule compute_all_read_counts:
    input:
        # Since we need to know what all the split numbers are to determine
        # which read counts to compute, we need to have the "all_generated.txt"
        # flag as an input. This way, Snakemake will update this rule's inputs
        # once the "generate_all_indexes" checkpoint is complete.
        ancient('all_generated.txt'),

        # Read counts are computed for each successful split of our genome index.
        # Furthermore, it's done separately for pair-ended and single-ended reads.
        # The read count script produces two files: one for the read counts
        # and one for the read count statistics.
        expand(
            f'{READ_COUNTS_DIR}/split_{{number}}_{{endedness}}_{{count_type}}.txt',
            number=successful_index_generation_numbers,
            endedness=READ_COUNT_ENDEDNESS,
            count_type=['read_counts', 'read_count_stats'])

    output: temp(touch('all_counts.txt'))

# Computes the counts for pair-ended reads.
# We have separate rules for pair-ended and single-ended reads
# just to make the rules easier to read.
rule compute_pe_read_counts:
    input:
        annotations = GENOME_GTF,

        # Only compute read counts over the pair-ended alignments.
        # We want the split number to be populated from the wildcard in our output,
        # so we need to tell "expand" that it's OK that it won't be populating all
        # of the wildcards (i.e., the "allow_missing=True" argument).
        bam_files = expand(
            f'{ALIGNMENT_DIR}/split_{{number}}/{{sample_id}}_PE_Aligned.sortedByCoord.out.bam',
            sample_id=PAIR_ENDED_SAMPLES,
            allow_missing=True)
    params:
        pair_ended = True
    output:
        read_counts = f'{READ_COUNTS_DIR}/split_{{number}}_PE_read_counts.txt',
        read_stats = f'{READ_COUNTS_DIR}/split_{{number}}_PE_read_count_stats.txt'
    script: 'scripts/read_counts.R'

# Computes the counts for pair-ended reads.
# We have separate rules for pair-ended and single-ended reads
# just to make the rules easier to read.
rule compute_se_read_counts:
    input:
        annotations = GENOME_GTF,

        # Only compute read counts over the single-ended alignments.
        # We want the split number to be populated from the wildcard in our output,
        # so we need to tell "expand" that it's OK that it won't be populating all
        # of the wildcards (i.e., the "allow_missing=True" argument).
        bam_files = expand(
            f'{ALIGNMENT_DIR}/split_{{number}}/{{sample_id}}_SE_Aligned.sortedByCoord.out.bam',
            sample_id=SINGLE_ENDED_SAMPLES,
            allow_missing=True)
    params:
        pair_ended = False
    output:
        read_counts = f'{READ_COUNTS_DIR}/split_{{number}}_SE_read_counts.txt',
        read_stats = f'{READ_COUNTS_DIR}/split_{{number}}_SE_read_count_stats.txt'
    script: 'scripts/read_counts.R'
