import logging
import os
import subprocess

from argparse import ArgumentParser
from tempfile import TemporaryDirectory

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
    The first argument is the path to the genome index.
    The second argument is the path to the FASTQ file.
    The third argument is the output path.
    '''
    parser = ArgumentParser(
        prog='STAR Alignment',
        description='Runs STAR to perform a genome alignment.')
    parser.add_argument('genome_index_path', help='The path to the genome index.')
    parser.add_argument('fastq_path', help='The path to the FASTQ file.')
    parser.add_argument('output_path', help='The output path.')
    args = parser.parse_args()
    perform_alignment(args.genome_index_path, args.fastq_path, args.output_path)

def snakemake_main() -> None:
    '''Runs STAR to perform alignment using data provided by Snakemake.'''
    output_path = snakemake.output[0].removesuffix('Aligned.sortedByCoord.out.bam')
    perform_alignment(
        snakemake.input.genome_index,
        snakemake.input.fastq,
        output_path
    )

def perform_alignment(genome_index_path: str, fastq_path: str, output_path: str) -> None:
    '''Runs STAR to perform a genome alignment.
    @param genome_index_path The path to the genome index.
    @param fastq_path The path to the FASTQ file.
    @param output_path The output path.
    '''
    TEMP_DIR_PARENT = 'tmp'
    os.makedirs(TEMP_DIR_PARENT, exist_ok=True)

    with TemporaryDirectory(dir=TEMP_DIR_PARENT) as temp_dir:
        star_temp_dir = os.path.join(temp_dir, 'star_tmp')
        command = [
            'STAR',
            '--runThreadN', '128',
            '--genomeDir', genome_index_path,
            '--readFilesIn', fastq_path,
            '--readFilesCommand', 'gunzip -c',
            '--outSAMtype', 'BAM', 'SortedByCoordinate',
            '--limitBAMsortRAM', '60000000000',
            '--outBAMsortingThreadN', '128',
            '--outFileNamePrefix', output_path,
            '--genomeLoad', 'NoSharedMemory',
            '--outSAMunmapped', 'Within',
            '--outSAMstrandField', 'intronMotif',
            '--outSJtype', 'Standard',
            '--outTmpDir', star_temp_dir,
        ]
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
