# rna-seq-align
A workflow for RNA sequence alignment.

## Notes - Testing Splitting

Downloaded `faSplit`.

```shell
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faSplit ./
```

Running with no arguments will display the help info.
Using `--help` will not.

Downloaded the example files Reema used in her email.

```shell
scp racecar:/data4/singh/DATA/Mus_musculus.GRCm39.dna_sm.toplevel.fa ./
scp racecar:/data4/singh/DATA/Mus_musculus.GRCm39.110.gtf ./
```

Split the fasta file how Reema suggested.
This created 61 different fasta files.

```shell
mkdir -p Genome
./faSplit byname ./Mus_musculus.GRCm39.dna_sm.toplevel.fa Genome/
```

This is the loop Reema gave (with slight modifications to work for me).

```shell
do
    STAR --runThreadN 40 \
        --runMode genomeGenerate \
        --genomeDir GenomeIndex/$file \
        --genomeFastaFiles $file \
        --sjdbGTFfile ./Mus_musculus.GRCm39.110.gtf \
        --sjdbOverhang 99 \
        --genomeSAindexNbases 12
done
```

Ran a single iteration of it by manually assigning `export file=GenomeIndex/10.fa`.
It took about 2.5 minutes to run that one iteration.
STAR created the folder `GenomeIndex/Genome/10.fa/` with a bunch of files in it.
Reema says she "ran the above script two times", but I don't see how that would affect anything.

At this point, I'm going to start creating the Snakemake workflow.
Then I can run this more properly and see why she might want to run it twice.
