import os
import subprocess

from tempfile import TemporaryDirectory

TEMP_DIR_PARENT = 'tmp'
os.makedirs(TEMP_DIR_PARENT, exist_ok=True)

with TemporaryDirectory(dir=TEMP_DIR_PARENT) as temp_dir:
    star_temp_dir = os.path.join(temp_dir, 'star_tmp')
    command = [
        'STAR',
        '--runThreadN', '128',
        '--genomeDir', snakemake.input.genome_index,
        '--readFilesIn', snakemake.input.fastq,
        '--readFilesCommand', 'gunzip -c',
        '--outSAMtype', 'BAM SortedByCoordinate',
        '--limitBAMsortRAM', '60000000000',
        '--outBAMsortingThreadN', '128',
        '--outFileNamePrefix', snakemake.output[0],
        '--genomeLoad', 'NoSharedMemory',
        '--outSAMunmapped', 'Within',
        '--outSAMstrandField', 'intronMotif',
        '--outSJtype', 'Standard'
    ]
    subprocess.run(command, check=True)
