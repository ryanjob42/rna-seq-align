import os
import subprocess

from tempfile import TemporaryDirectory

TEMP_DIR_PARENT = 'tmp'
os.makedirs(TEMP_DIR_PARENT, exist_ok=True)

with TemporaryDirectory(dir=TEMP_DIR_PARENT) as temp_dir:
    star_temp_dir = os.path.join(temp_dir, 'star_tmp')
    command = [
        'STAR',
        '--runMode', 'genomeGenerate',
        '--genomeDir', f'{snakemake.params.genome_index_dir}/split_{snakemake.wildcards.number}.fa',
        '--genomeFastaFiles', snakemake.input.fasta,
        '--sjdbGTFfile', snakemake.input.gtf,
        '--sjdbOverhang', '99',
        '--genomeSAindexNbases', '12',
        '--outTmpDir', star_temp_dir,
        '--outFileNamePrefix', temp_dir,
    ]
    subprocess.run(command)

exit(0)
