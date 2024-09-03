# This script computes the read counts for a .bam file
# using Rsubread's "featureCounts" function.

# This script can either be run directly, or be run by Snakemake.
# The "if" statement at the top detects which case we're in
# and captures the script arguments appropriately.
# Then, "featureCounts" is run and results are saved to output files.

# When running via the command line, the arguments are as follows:
# 1. annotation_file: The path to the annotations file (e.g., a GTF file).
#        This is given to the "annot.ext" argument of featureCounts.
# 2. bam_directory: The path to the directory containing the BAM files to use.
#        All .bam files found are given to the "files" argument of featureCounts.
# 3. pair_ended: Whether or not the reads are pair ended. Provide either "True" or "False".
#        This is given to the "isPairEnded" argument of featureCounts.
# 4. read_count_file: The file path to write the read counts to.
# 5. read_stats_file: The file path to write the read count statistics to.
# 6. thread_count: The number of CPU threads to use.

# Load the "Rsubread" library so we can use the "featureCounts" function.
library(Rsubread)

# Either pull the arguments from Snakemake or the command line arguments.
# Note: Snakemake provides a list of BAM files directly, while the command line arguments
# accept a directory to search through for where the BAM files can be found.
if (exists("snakemake")) {
    annotation_file <- snakemake@input[['annotations']]
    bam_files <- snakemake@input[['bam_files']]
    pair_ended <- snakemake@params[['pair_ended']]
    read_count_file <- snakemake@output[['read_counts']]
    read_stats_file <- snakemake@output[['read_stats']]
    thread_count <- 128
} else {
    args = commandArgs(trailingOnly=TRUE)
    if (length(args) != 6) {
        stop("Exactly 6 arguments are required. Read the script for details.")
    }

    annotation_file <- args[1]
    bam_directory <- args[2]
    pair_ended <- args[3]
    read_count_file <- args[4]
    read_stats_file <- args[5]
    thread_count <- args[6]

    # Get a list of the BAM files in the BAM directory.
    bam_files = list.files(path=bam_directory, patterm="*.bam$")
}

# Use featureCounts to get the read counts and their statistics.
Sample1 <- featureCounts(
    files=bam_files,
    annot.ext=annotation_file,
    isGTFAnnotationFile=TRUE,
    GTF.featureType="exon",
    GTF.attrType="gene_id",
    isPairedEnd=pair_ended,
    nthreads=thread_count)

# Write the read counts to the output file.
write.table(
    Sample1$counts,
    file=read_count_file,
    sep="\t",
    row.names=TRUE,
    col.names=TRUE,
    quote=FALSE)

# Write the read count statistics to the output file.
write.table(
    Sample1$stat,
    file=read_stats_file,
    sep="\t",
    row.names=TRUE,
    col.names=TRUE,
    quote=FALSE)
