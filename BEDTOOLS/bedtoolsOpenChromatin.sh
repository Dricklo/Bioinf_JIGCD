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

bedtools intersect -a /ocean/projects/bio230007p/ikaplow/HumanAtac/AdrenalGland/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz -b /ocean/projects/bio230007p/idamico/results/idr.conservative_peak.HumanToMouse.HALPER.narrowPeak -wa > human_adrenal_peaks_with_open_mouse_orthologs.bed

bedtools intersect -a /ocean/projects/bio230007p/ikaplow/HumanAtac/AdrenalGland/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz -b /ocean/projects/bio230007p/idamico/results/idr.conservative_peak.HumanToMouse.HALPER.narrowPeak -v > human_adrenal_peaks_with_closed_mouse_orthologs.bed

