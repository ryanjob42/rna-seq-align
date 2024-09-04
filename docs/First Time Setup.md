# First Time Setup
There are a few first-time setup steps to run before you can use this pipeline.

## Download the Pipeline
First, you will need to download the pipeline.

1. Log in to Slipstick
2. Use the `cd` command to go to a location where you want to download this pipeline to.
3. Run the command below to download the pipeline. It will be saved as the `rna-seq-align` folder.
   1. `git clone https://github.com/ryanjob42/rna-seq-align.git`

## Install STAR
Note: even if you have STAR already, you may need to do the last two steps still.
For the pipeline to work, STAR must be on your `PATH`, meaning you can run it simply by typing `STAR` from any location.

Next, you will need to download STAR and add it to your `PATH` environment variable.
If you are downloading it for the first time, you can do the following:

1. In your browser, navigate to the link below.
   1. https://github.com/alexdobin/STAR/releases/latest
2. Copy the link to the `STAR_*.zip` file (where `*` is the version number).
3. On Slipstick, use the `cd` command to go to a location where you want to download STAR to.
4. Use the `wget` command followed by the link to the Zip file to download it.
   1. For example: `wget https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip`
5. Use the `unzip` command to unzip STAR.
   1. For example: `unzip STAR_*.zip`
6. Use the `cd` command to go into the unzipped folder, then into the `Linux_x86_64_static` folder.
   1. For example: `cd STAR_*/Linux_x86_64_static`
7. Run the command below to add the STAR executable to your `PATH` environment variable. Note: the single vs. double quotes are necessary, so copy/paste it please!
   1. `echo 'export PATH="'$(pwd)':$PATH"' >> ~/.bashrc && source ~/.bashrc`
   2. This will add a line to your `.bashrc` file, which is a file that loads all your settings when you log in.
   3. This line specifically will make it so you can run the `STAR` command from anywhere, and it will always point to that version of STAR.
   4. Note: if you need to update this in the future, you can simply edit your `.bashrc` file using `vim` (or your preferred text editor).

If you've downloaded STAR before, but it's not on your `PATH` (i.e., you can't enter `STAR` on the command line and have it work), do the following:

## Install Python and the Conda Environment
Note: while Slipstick has modules for Python and Miniconda, I recommend you follow these instructions anyways.
This will have you install your own version of Miniconda (which contains Python), that way you can keep it up to date and use the latest package versions.

Run the commands below to install the latest version of Miniconda:

```shell
mkdir -p "$HOME/miniconda3"
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O "$HOME/miniconda3/miniconda.sh"
bash "$HOME/miniconda3/miniconda.sh" -b -u -p "$HOME/miniconda3"
```

Run the commands below to do some cleanup and make sure Conda is initialized.

```shell
rm "$HOME/miniconda3/miniconda.sh"
source "$HOME/.bashrc"
"$HOME/miniconda3/bin/conda init bash"
```

From inside the `rna-seq-align` folder, run the command below to install the Conda environment for this workflow:

```shell
conda env create -f environment.yml
```

## Install the Rsubread R Package
Note: even if you have done this before, you should still follow these instructions to make sure the pipeline can access it.

Computing the read counts is done using the Rsubread package in R.
We need to make sure it is installed and accessible from within the pipeline.
To do so, follow these steps:

1. On Slipstick, run the `module load r/4.1.0` command to load the current version of R:
2. Run the `R` command to launch R in an interactive mode.
3. Enter the command below to install the BiocManager package, which will help us install Rsubread.
   1. `install.packages("BiocManager")`
   2. If it asks you if you'd like to make a personal library, enter `yes`. You may need to agree to multiple similar prompts.
   3. If it asks you which CRAN mirror to use, enter the number of the one that is closest to you. `66` worked for me (not all of them work, though).
   4. If there is an error, try continuing on anyways. If it still doesn't work, use `q()` to exit, then `n` to skip saving the workspace image. Restart from step 2, selecting a different CRAN mirror.
4. Enter the command below to install the Rsubread package.
   1. `BiocManager::install("Rsubread")`
   2. This may take a minute or two.
5. Enter the command below to see if it installed correctly. If it is installed correctly, it will just continue on with no message.
   1. `library(Rsubread)`
6. Enter the `q()` command to exit, then `n` to skip saving the workspace image.
