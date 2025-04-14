#!/bin/bash

#SBATCH --job-name=bedtools
#SBATCH --output=bedtools.out
#SBATCH --error=bedtools.err
#SBATCH --partition=RM-shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

#SBATCH --mail-type=ALL
#SBATCH --mail-user=idamico@andrew.cmu.edu

module load bedtools/2.30.0

bedtools intersect -a mouse_liver_conservative_peak.narrowPeak.gz -b mouse_adrenal_conservative_peak.narrowPeak.gz > mouse_shared_across_tissues.bed

bedtools intersect -v -a mouse_liver_conservative_peak.narrowPeak.gz -b mouse_adrenal_conservative_peak.narrowPeak.gz > mouse_liver_specific_OCR.bed

bedtools intersect -v -a mouse_adrenal_conservative_peak.narrowPeak.gz -b mouse_liver_conservative_peak.narrowPeak.gz > mouse_adrenal_specific_OCR.bed

