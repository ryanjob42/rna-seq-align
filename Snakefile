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

# Read the FASTQ inputs file into a dictionary of lists.
# The file is a tab-separated file with two columns.
# The first is the "sample_id": the unique name given to the sample.
# The second is a comma-separated list of FASTQ read files associated with it.
# All FASTQ files in the same "sample_id" are put into STAR together.
# If there is only a single FASTQ file, it's treated as a single-ended read.
# Otherwise, if there are multiple, it's treated as a pair-ended read.
# While parsing, make sure there are no duplicate names.
FASTQ_FILES = {}
import csv
with open(FASTQ_INPUT_FILE_PATH) as fastq_input_file:
    reader = csv.DictReader(fastq_input_file, delimiter='\t')
    for line in reader:
        if line['Name'] in FASTQ_FILES:
            print(f'The FASTQ input file contained a duplicate name: {line['Name']}')
            exit(1)
        FASTQ_FILES[line['Name']] = line['Fastq'].split(',')

# This function checks if a specific sample ID (per the "Name" column
# of the FASTQ input file) is a pair-ended read (i.e., True) or a single-ended
# read (i.e., False). If there is only one file, it's assumed to be single-ended.
# Otherwise, it's assumed to be pair-ended.
# This function takes in Snakemake wildcards and checks the "sample_id" wildcard
# against the "Name" column of the FASTQ input file.
def is_pair_ended(wildcards):
    input_files = FASTQ_FILES[wildcards.sample_id]
    return isinstance(input_files, list) and (len(input_files) > 1)

IS_PAIR_ENDED = {
    sample_id: (isinstance(fastq_files, list) and (len(fastq_files) > 1))
    for sample_id, fastq_files in FASTQ_FILES.items()
}

# The samples must all be single-ended or must all be pair-ended.
# We can't allow a mix of the two currently.
if (True in IS_PAIR_ENDED.values()) and (False in IS_PAIR_ENDED.values()):
    print('Currently, we require all samples to be single-ended or all samples to be pair-ended, but you provided some of both.')
    exit(1)
ALL_PAIR_ENDED = (True in IS_PAIR_ENDED.values())


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
        alignments = expand(f'{ALIGNMENT_DIR}/split_{{number}}/{{sample_id}}Aligned.sortedByCoord.out.bam',
            number=successful_index_generation_numbers,
            sample_id=FASTQ_FILES.keys())
    output: temp(touch('all_aligned.txt'))

# Performs a single alignment.
# The list of FASTQ files needed is determined by looking up all the file names
# from the "FASTQ_FILES" dictionary for the associated sample ID.
rule perform_single_alignment:
    input:
        genome_index = f'{SPLIT_GENOME_INDEX_DIR}/split_{{number}}.fa',
        fastq = lookup(dpath=f'{{sample_id}}', within=FASTQ_FILES)
    params:
        star_threads = STAR_ALIGNMENT_THREADS,
        star_memory = STAR_ALIGNMENT_RAM
    output:
        f'{ALIGNMENT_DIR}/split_{{number}}/{{sample_id}}Aligned.sortedByCoord.out.bam'
    script: 'scripts/perform_single_alignment.py'


##############################################
# Read Counts and Stats
##############################################

# Computes all read counts for all splits.
rule compute_all_read_counts:
    input:
        generation_complete_flag = ancient('all_generated.txt'),
        alignment_complete_flag = ancient('all_aligned.txt'),
        read_counts = expand(f'{READ_COUNTS_DIR}/split_{{number}}_read_counts.txt', number=successful_index_generation_numbers),
        read_stats = expand(f'{READ_COUNTS_DIR}/split_{{number}}_read_count_stats.txt', number=successful_index_generation_numbers)
    output: temp(touch('all_counts.txt'))

# Computes the read counts for a single split.
rule compute_split_read_counts:
    input:
        annotations = GENOME_GTF,
        bam_files = expand(
            f'{ALIGNMENT_DIR}/split_{{number}}/{{sample_id}}Aligned.sortedByCoord.out.bam',
            sample_id=FASTQ_FILES.keys(),
            allow_missing=True)
    params:
        pair_ended = ALL_PAIR_ENDED
    output:
        read_counts = f'{READ_COUNTS_DIR}/split_{{number}}_read_counts.txt',
        read_stats = f'{READ_COUNTS_DIR}/split_{{number}}_read_count_stats.txt'
    script:
        'scripts/read_counts.R'
