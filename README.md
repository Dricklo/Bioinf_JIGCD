## 03-713: Bioinformatics Data Integration Practicum (Spring 2025) Group Project
<!--
## Notes from Lecture for README

- Include steps for running our pipeline

- Have an introduction section and a dependencies section

- Include an install section (?)

- If using conda environments, you can upload the `config.yml`, which details all the packages in the environment:
  
  ```bash
  conda env export > environment_drop.yml
-->
## Introduction
This pipeline compares ATAC-seq data for two species and two tissues to assess:
1. Is transcriptional regulatory element activity more conserved across tissues or species? How does this transcriptional regulatory element conservation and code differ between tissues and species?
2. To what extent does the transcriptional regulatory code differ between enhancers and promoters? Does this depend on tissue or species?
3. To what extent are the biological processes upregulated in tissue conserved across species?

## Dependencies
There are a couple tools that are required to run this pipeline: HALPER, halLiftover, bedtools, GREAT, and the MEME suite.

**Bedtools**
Note that the PSC cluster comes with bedtools already installed, but must be loaded by entering this command into the bash command line:
```
module load bedtools/2.30.0
```
If using the PSC cluster to run all the scripts below, make sure this module is always loaded beforehand!
 
If the cluster does not have bedtools installed or you are running everything locally, you can also install bedtools with two methods: 

**Method 1: Environment installation**  
To install bedtools to a previous environment. You first need to activate your environment
```
conda activate environment_name 
```
You can also create a new environment, with a specific name instead of `<my-env>`
```
conda create --name <my-env>
conda activate <my-env>
``` 
Then, install with:
```
conda install -c bioconda bedtools 
```

**Method 2: Local installation**  
If installing to the local machine for operation outside of the environment or for a SLURM cluster, follow the instructions at the given website:
https://bedtools.readthedocs.io/en/latest/content/installation.html
For clarification sake, here are some simple steps to install it running in bash: 
Step 1: make a directory for bedtools from the home directory and navigate to it.
```
mkdir -p ~/tools
cd ~/tools
```

Step 2: clone the github repo for it, and then navigate into it
```
git clone https://github.com/arq5x/bedtools2.git
cd bedtools2
``` 

Step 3: Compile bedtools 
```
make
```
There should be files now located in ~/tools/bedtools2/bin/bedtools 

Step 4: Verify bedtools is now installed:
```
./bin/bedtools --version
```
There should be an output of bedtools and its version.

Step 5: Add it permanently to bash profile, and then reload. 
```
export PATH=$HOME/tools/bedtools2/bin:$PATH
source ~/.bashrc
```


**HALPER/HalLiftOver installation**

To be able to produce the outputs for step 2. You need to make sure that HalLiftover/HALPER is installed. The main instructions are detailed in the README.md from:

https://github.com/pfenninglab/halLiftover-postprocessing](https://github.com/pfenninglab/halLiftover-postprocessing/blob/master/hal_install_instructions.md

Note that PSC does not come with HALPER/HalLiftover installed, and thus must be followed manually. There are several sub-dependencies in this installation process itself (from the installation README):
* Python version 3.6 or 3.7
* Python libraries matplotlib and numpy
  * numpy (http://www.numpy.org/)
    * HALPER has been tested using numpy versions 1.14.3, 1.16.0, 1.16.4, 1.16.6, and 1.18.2
  * matplotlib (https://matplotlib.org/downloads.html)
    * HALPER has been tested using matplotlib versions 1.5.1, 2.2.3, 3.2.1
* HALPER has been tested on Linux (CentOS 6, CentOS 7, and Ubuntu 18.04.4), Windows (Windows 10), and Mac

**Conda Environment**

Although the Conda Environment is properly setup in the HALPER/halLiftover installation methods, I have also included it in the environment folder as a hal.yml file. It can be found in the github repo under the folder _environment_ as "hal.yml". 
To create a new environment with the hal.yml file, you would enter the following into bash:
```
conda env create -f hal.yml
```
To update a previous environment that is already activated with conda activate, you would use:
```
conda env update -f hal.yml
``` 
**MEME suite**

MEME suite is installed following the steps from this website: 
https://meme-suite.org/meme/doc/install.html?man_type=web  

We will also be using the Bed2FASTA utility to convert .bed files to fasta files for MEME suite.  
<!--
//include eventual section about using bed2Fasta to convert .bed to fasta files.  
//include installation instructions for MEME.
-->

## Installation Instructions
Once the system has been configured as specified in **Dependencies**, you need to clone this GitHub repository onto your machine. The inputs specified below are included in the repository, and the script is hard-coded to pull these inputs. 

## Inputs
* .bed or .narrowPeak ATAC-seq files for two species and two tissues (resulting in four files total)
* Quality control data from running ENCODE on the ATAC-seq data for each species and tissue
* Cactus whole genome alignment .hal file, including the two species of interest
* Names of the two species of interest in the Cactus alignment, specified as the source species and the target species
* .bed file of TSS peak enrichment data for the two species of interest

## Running the Pipeline
After you have configured your system as specifed in **Installation Instructions**, please cd into the scripts folder and use the command below to run the pipeline on a slurm cluster:

``` bash
sbatch wrapperScript.sh
```

This pipeline can take a while to run, depending on the supercomputer - up to a day or so.

## Common Issues
If you run into issues running the above script due to your system not locating the halLiftover script, make sure you run the below code from steps 12 and 13 of the [hal/HALPER instruction manual]([url](https://github.com/pfenninglab/halLiftover-postprocessing/blob/master/hal_install_instructions.md)). 

``` bash
export PATH=[repos dir]/hal/bin:${PATH}
export PYTHONPATH=[repos dir]/halLiftover-postprocessing:${PYTHONPATH}
source ~/.bash_profile
```

## Outputs
### Step 1:
* *data/Human_qc/qc_TISSUE_human.html, data/Mouse_qc/qc_TISSUE_mouse.html*: the ENCODE data for each species' tissue ATAC-seq data.
### Step 2:
* *results/step2_HALPER_results/idr.conservative_peak.TISSUE_SPECIES1ToSPECIES2.HALPER.narrowPeak*: the first species' open chromatin regions mapped to the second species' genome, and the second species' open chromatin regions mapped to the first species' genome. 
* *results/step2a_bedtools_results/SPECIES1_TISSUE_peaks_with_open_SPECIES2_orthologs.bed, results/step2a_bedtools_results/SPECIES1_TISSUE_peaks_with_closed_SPECIES2_orthologs.bed*: open chromatin regions in each species whose orthologs in the other species are open and those whose orthologs in the other species are closed.
### Step 3:
* *result/step3_comparison/cross_tissue/SPECIES_shared_across_tissues.bed, result/step3_comparison/cross_tissue/SPECIES_TISSUE_<specific>_OCR.bed, result/step3_comparison/cross_species/SPECIES1_TISSUE_shared_with_SPECIES2_TISSUE.bed*: open chromatin regions in each species that are open in both tissues and those that are open in only one tissue.
* *result/step3_comparison/cross_tissue/cross_tissue_percentages.txt, result/step3_comparison/cross_species/cross_species_percentages.txt*: quantification of which species/tissue combination have more open chromatin regions.
### Step 4:
* *filenames*: candidate biological processes regulated by open chromatin regions in each species/tissue combination, as well as open chromatin regions shared across tissues, specific to each tissue, shared across species, and specific to each species.
### Step 5:
* *filenames*: quantification of how the percentage of promotors and enhancers in the open chromatin data in each species/tissue combination compares across tissues and across species.
### Step 6:
* *filenames*: sequence patterns that occur more than expected by chance, in the full set of peaks for each species/tissue combination, enhancers and promoters for each species/tissue combination, and enhancers that are shared across tissues, specific to each tissue, shared across species, and specific to each species.

## Citations

- **ENCODE**  
  ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome. *Nature*. 2012 Sep 6;489(7414):57-74. [PMID: 22955616]

- **HALPER / halLiftover**  
  Zhang Y, Li T, Preissl S, Amaral ML, Grinstein JD, Farah EN, et al. HALPER facilitates the identification of regulatory element orthologs across species. *Nucleic Acids Res*. 2020 Jul 9;48(12):e72. [https://pubmed.ncbi.nlm.nih.gov/32407523/]  
  GitHub: https://github.com/pfenninglab/halLiftover-postprocessing

- **GREAT**  
  McLean CY, Bristor D, Hiller M, Clarke SL, Schaar BT, Lowe CB, et al. GREAT improves functional interpretation of cis-regulatory regions. *Nat Biotechnol*. 2010 May;28(5):495-501. [https://pubmed.ncbi.nlm.nih.gov/20436461/]  
  Website: https://great.stanford.edu/great/public/html/

- **bedtools**  
  Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*. 2010 Mar 15;26(6):841-2. [https://pubmed.ncbi.nlm.nih.gov/20110278/]  
  Docs: https://bedtools.readthedocs.io/en/latest/

- **MEME Suite**  
  Bailey TL, Boden M, Buske FA, Frith M, Grant CE, Clementi L, et al. MEME Suite: tools for motif discovery and searching. *Nucleic Acids Res*. 2009 Jul;37(Web Server issue):W202-8. [https://pubmed.ncbi.nlm.nih.gov/21486936/]  
  Website: https://meme-suite.org/meme/  
  Docs: https://meme-suite.org/meme/doc/meme-chip.html


## Contributors
* Darrick Lo
* Jinqi Hou
* Grace Du
* Isabelle D'Amico
