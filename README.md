# Bioinf_JIGCD
Bioinformatics group project

## Notes from Lecture for README

- Include steps for running our pipeline

- Have an introduction section and a dependencies section

- Include an install section (?)

- If using conda environments, you can upload the `config.yml`, which details all the packages in the environment:
  
  ```bash
  conda env export > environment_drop.yml

## Introduction
This pipeline compares ATAC-seq data for two species and two tissues to assess:
1. Is transcriptional regulatory element activity more conserved across tissues or species? How does this transcriptional regulatory element conservation and code differ between tissues and species?
2. To what extent does the transcriptional regulatory code differ between enhancers and promoters? Does this depend on tissue or species?
3. To what extent are the biological processes upregulated in tissue conserved across species?

## Dependencies
There are a couple tools that are required to run this pipeline: HALPER, halLiftover, bedtools, GREAT, and the MEME suite. Other than GREAT, which is utilized outside of the command line, dependency set-up can be referenced with <.yml file>. 

*(We gotta make sure that setting up a conda environment with the tools works for running the eventual pipeline script file)*

## Installation Instructions
Once the system has been configured as specified in **Dependencies**, you should be good to go to run the pipeline as specified in **Running the Pipeline**. 

## Inputs
* *(do we include the ENCODE data in here?)*
* .bed or .narrowPeak ATAC-seq files for two species and two tissues (resulting in four files total)
* Cactus whole genome alignment .hal file, including the two species of interest
* Names of the two species of interest in the Cactus alignment, specified as the source species and the target species
* *TSS peak data*
* *(am I forgetting anything else?)*

## Running the Pipeline
After you have configured your system as specifed in **Installation Instructions**, use the command below to run the pipeline on a slurm cluster:

``` bash
sbatch <name of script> <peakfile1> <peakfile2> <peakfile3> <peakfile4> <wholegenomealignment> <speciesName1> <speciesName2> <TSSdata>
```

This pipeline can take a while to run, depending on the supercomputer - up to a day or so.

*(do we intend on including instructions for not running it on slurm cluster?)*

## Outputs
### Step 1:
* *(do we include this in here?)*
### Step 2:
* *filenames*: the first species' open chromatin regions mapped to the second species' genome, and the second species' open chromatin regions mapped to the first species' genome. 
* *filenames*: open chromatin regions in each species whose orthologs in the other species are open and those whose orthologs in the other species are closed.
### Step 3:
* *filenames*: open chromatin regions in each species that are open in both tissues and those that are open in only one tissue.
* *filenames*: quantification of which species/tissue combination have more open chromatin regions.
### Step 4:
* *filenames*: candidate biological processes regulated by open chromatin regions in each species/tissue combination, as well as open chromatin regions shared across tissues, specific to each tissue, shared across species, and specific to each species.
### Step 5:
* *filenames*: quantification of how the percentage of promotors and enhancers in the open chromatin data in each species/tissue combination compares across tissues and across species.
### Step 6:
* *filenames*: sequence patterns that occur more than expected by chance, in the full set of peaks for each species/tissue combination, enhancers and promoters for each species/tissue combination, and enhancers that are shared across tissues, specific to each tissue, shared across species, and specific to each species.

## Citations
* ENCODE
* HALPER/halLiftover
* GREAT
* bedtools
* MEME suite

## Contributors
* Darrick Lo
* Jinqi Hou
* Grace Du
* Isabelle D'Amico
