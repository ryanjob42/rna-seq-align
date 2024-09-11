import logging
import os
import subprocess

from argparse import ArgumentParser
from tempfile import TemporaryDirectory
from typing import Union

# This script can either be run directly, or be run by Snakemake.
# The "if" statement at the very end detects which case we're in
# and either runs the "direct_main" or the "snakemake_main" function.
# These functions grab the input and output paths,
# then use the "perform_alignment" function to actually run STAR.

# If you're experiencing issues with this script, uncomment the line below
# to turn on debug logging.
# logging.basicConfig(level=logging.DEBUG)

def direct_main() -> None:
    '''Runs STAR to perform alignment using data provided via command-line arguments.
    The first argument is the number of threads that STAR can use. This is STAR's "--runThreadN" argument.
    The second argument is the amount of memory that STAR can use. This is STAR's "--limitBAMsortRAM" argument.
    The third argument is the path to the genome index. This is STAR's "--genomeDir" argument.
    The fourth argument is the path to the FASTQ file. This is STAR's "--readFilesIn" argument.
    The fifth argument is the prefix STAR uses for all output files. This is STAR's "--outFileNamePrefix" argument.
    '''
    parser = ArgumentParser(
        prog='STAR Alignment',
        description='Runs STAR to perform a genome alignment.')
    parser.add_argument('star_threads', help='''The number of threads that STAR can use. This is STAR's "--runThreadN" argument.''')
    parser.add_argument('star_memory', help='''The amount of memory that STAR can use. This is STAR's "--limitBAMsortRAM" argument.''')
    parser.add_argument('genome_index_path', help='''The path to the genome index. This is STAR's "--genomeDir" argument.''')
    parser.add_argument('fastq_path', help='''The path to the FASTQ file. This is STAR's "--readFilesIn" argument.''')
    parser.add_argument('output_prefix', help='''The prefix STAR uses for all output files. This is STAR's "--outFileNamePrefix" argument.''')
    args = parser.parse_args()
    perform_alignment(args.star_threads, args.star_memory, args.genome_index_path, args.fastq_path, args.output_prefix)

def snakemake_main() -> None:
    '''Runs STAR to perform alignment using data provided by Snakemake.'''
    # Snakemake will give us the entire name of the output file we want.
    # However, since this is given to STAR as the "--outputPrefix" argument,
    # it will append "Aligned.sortedByCoord.out.bam" to whatever we give it.
    # Thus, we need to strip that prefix off first so it doesn't appear twice.
    output_prefix = snakemake.output[0].removesuffix('Aligned.sortedByCoord.out.bam')
    perform_alignment(
        snakemake.params.star_threads,
        snakemake.params.star_memory,
        snakemake.input.genome_index,
        snakemake.input.fastq,
        output_prefix
    )

def perform_alignment(threads: Union[str,int], memory: Union[str,int], genome_index_path: str, fastq_path: str, output_prefix: str) -> None:
    '''Runs STAR to perform a genome alignment.
    @param genome_index_path The path to the genome index.
    @param fastq_path The path to the FASTQ file.
    @param output_path The output path.
    '''
    logging.debug('STAR thread count: %s', str(threads))
    logging.debug('STAR RAM limit: %s', str(memory))
    logging.debug('Genome index: %s', genome_index_path)
    logging.debug('FASTQ file: %s', fastq_path)
    logging.debug('Output prefix: %s', output_prefix)

    # We want to make sure the "./tmp" directory exists,
    # as all instances of STAR will use it for temporary storage.
    # This is preferable to the "/tmp" directory, as on Slipstick
    # or Riviera, it may not be large enough.
    TEMP_DIR_PARENT = 'tmp'
    os.makedirs(TEMP_DIR_PARENT, exist_ok=True)

    # Inside the "./tmp" directory, we want a subdirectory
    # dedicated for this specific instance of STAR.
    with TemporaryDirectory(dir=TEMP_DIR_PARENT) as temp_dir:
        # STAR doesn't want the directory to exist yet, however.
        # Because of that, we'll tell star to use the "star_tmp"
        # directory inside this instance's temporary directory.
        # It's slightly roundabout, but works well enough.
        star_temp_dir = os.path.join(temp_dir, 'star_tmp')

        command = [
            'STAR',
            '--runThreadN', str(threads),
            '--genomeDir', genome_index_path,
            '--readFilesIn', fastq_path,
            '--readFilesCommand', 'gunzip -c',
            '--outSAMtype', 'BAM', 'SortedByCoordinate',
            '--limitBAMsortRAM', str(memory),
            '--outBAMsortingThreadN', '128',
            '--outFileNamePrefix', output_prefix,
            '--genomeLoad', 'NoSharedMemory',
            '--outSAMunmapped', 'Within',
            '--outSAMstrandField', 'intronMotif',
            '--outSJtype', 'Standard',
            '--outTmpDir', star_temp_dir,
        ]
        logging.debug('STAR command: %s', ' '.join(command))
        subprocess.run(command, check=True)

if __name__ == '__main__':
    # If the "snakemake" variable is defined, assume this is being run by Snakemake.
    # Otherwise, assume this is being run directly.
    if 'snakemake' in locals():
        logging.debug('This script was run by Snakemake.')
        snakemake_main()
    else:
        logging.debug('This script was not run by Snakemake.')
        direct_main()
