# RNA Sequence Alignment
A Snakemake workflow for RNA sequence alignment.

## Setup
There is a one-time setup process for this workflow.
See the [First Time Setup](./docs/First%20Time%20Setup.md) instructions for more details.

## Usage
There are a few important pieces to using this workflow.

First, you will update the `config.yaml` configuration file.
This is a text file that you can use to tell the workflow what you want to do.
The file is commented as to what each setting does.

Second, you'll upload your genome's FASTA file and annotations file (the .gtf file).
The default place for them is in the `data` folder here, but you can change the folder using the `config.yaml` file.
Note: you may need to change the names of these files in the `config.yaml` file.

Third, you'll upload your FASTQ files you want to align.
The default place for them is in the `fastq` folder here, but you can change that using the `config.yaml` file.
Currently, the workflow will assume these are all single-ended files, and that you want to use all of them.
Ryan is working to change this, and will have an update soon.

Finally, once everything is ready, you can simply run the `start-workflow.sh` script.
This contains the couple of commands needed to run Snakemake.
While it's running, you MUST leave the terminal window open, otherwise it will stop the workflow.

If you'd like to be able to close the terminal window, here is what you can do:
1. Run the `tmux` command. After a moment, the terminal will refresh, and you'll see a (potentially green) status bar at the bottom.
   1. The `tmux` program lets you start a command, leave, then come back.
2. Run the `start-workflow.sh` command. Give it a minute or two to make sure it's running.
3. When you want to leave, press `Ctrl+b`, then tap `d` to "detatch" from your tmux session.
   1. It's the same on all platforms: Windows, Linux, or Mac.
4. To re-attach later, simply run the command `tmux attach`.
