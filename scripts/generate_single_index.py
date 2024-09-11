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
# then use the "generate_index" function to actually run STAR.

# If you're experiencing issues with this script, uncomment the line below
# to turn on debug logging.
# logging.basicConfig(level=logging.DEBUG)

def direct_main() -> None:
    '''Runs STAR to try generating a genome index using data provided by command-line arguments.
    The first argument is the number of threads that STAR can use. This is STAR's "--runThreadN" argument.
    The second argument is the directory for the genome index to be created. This is STAR's "--genomeDir" argument.
    The third argument is the FASTA file to generate from. This is STAR's "--genomeFastaFiles" argument.
    The fourth argument is the GTF file to generate from. This is STAR's "--sjdbGTFfile" argument.
    '''
    parser = ArgumentParser(
        prog='STAR Genome Index Generation',
        description='Runs STAR to try generating a genome index.')
    parser.add_argument('star_threads', help='''The number of threads that STAR can use. This is STAR's "--runThreadN" argument.''')
    parser.add_argument('genome_index_path', help='''The directory for the genome index to be created. This is STAR's "--genomeDir" argument.''')
    parser.add_argument('fasta_path', help='''The FASTA file to generate from. This is STAR's "--genomeFastaFiles" argument.''')
    parser.add_argument('gtf_path', help='''The GTF file to generate from. This is STAR's "--sjdbGTFfile" argument.''')
    args = parser.parse_args()
    generate_index(args.star_threads, args.genome_index_path, args.fasta_path, args.gtf_path)

def snakemake_main() -> None:
    '''Runs STAR to try generating a genome index using data provided by Snakemake.'''
    generate_index(
        snakemake.params.star_threads,
        snakemake.params.genome_index_path,
        snakemake.input.fasta,
        snakemake.input.gtf
    )

def generate_index(threads: Union[str,int], genome_index_path: str, fasta_path: str, gtf_path) -> None:
    '''Runs STAR to generate a genome index.
    @param star_threads The number of threads that STAR can use. This is STAR's "--runThreadN" argument
    @param genome_index_path The directory for the genome index to be created. This is STAR's "--genomeDir" argument
    @param fasta_path The FASTA file to generate from. This is STAR's "--genomeFastaFiles" argument
    @param gtf_path The GTF file to generate from. This is STAR's "--sjdbGTFfile" argument
    '''
    logging.debug('STAR thread count: %s', str(threads))
    logging.debug('Genome index: %s', genome_index_path)
    logging.debug('FASTA file: %s', fasta_path)
    logging.debug('GTF file: %s', gtf_path)

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

        # We want to try generating the genome index in this location.
        # The location can't be provided as an output however, as an index
        # won't be created if the given FASTA file doesn't contain
        # useful information for the given GTF file.
        # This is OK, and we want to skip it in that case.
        # However, if we put the directory as an output of this rule
        # and it doesn't get created, Snakemake will fail the pipeline.
        # As of Sep 2024, there is no way to make an output optional.
        command = [
            'STAR',
            '--runThreadN', str(threads),
            '--runMode', 'genomeGenerate',
            '--genomeDir', genome_index_path,
            '--genomeFastaFiles', fasta_path,
            '--sjdbGTFfile', gtf_path,
            '--sjdbOverhang', '99',
            '--genomeSAindexNbases', '12',
            '--outTmpDir', star_temp_dir,
            '--outFileNamePrefix', star_temp_dir,
        ]
        logging.debug('STAR command: %s', ' '.join(command))
        subprocess.run(command, check=False)

if __name__ == '__main__':
    # If the "snakemake" variable is defined, assume this is being run by Snakemake.
    # Otherwise, assume this is being run directly.
    if 'snakemake' in locals():
        logging.debug('This script was run by Snakemake.')
        snakemake_main()
    else:
        logging.debug('This script was not run by Snakemake.')
        direct_main()
