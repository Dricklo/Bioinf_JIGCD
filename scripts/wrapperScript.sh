#!/bin/bash

#SBATCH --job-name=halper
#SBATCH --output=halper.out
#SBATCH --error=halper.err
#SBATCH --partition=RM-shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

module load anaconda3/2024.10-1

conda init
conda activate hal

# Step 2

# Human Adrenal to Mouse
bash halper_map_peak_orthologs.sh \
 -b ../data/OCR_idr_conservative/Human_adrenal_peak/idr.conservative_peak.narrowPeak \
 -o ../results/step2_HALPER_results/ \
 -s Human \
 -t Mouse \
 -c ../data/Alignment/10plusway-master.hal

# Human Liver to Mouse
bash halper_map_peak_orthologs.sh \
 -b ../data/OCR_idr_conservative/Human_liver_peak/idr.conservative_peak.narrowPeak \
 -o ../results/step2_HALPER_results/ \
 -s Human \
 -t Mouse \
 -c ../data/Alignment/10plusway-master.hal

# Mouse Adrenal to Human
bash halper_map_peak_orthologs.sh \
 -b ../data/OCR_idr_conservative/Mouse_adrenal_peak/idr.conservative_peak.narrowPeak \
 -o ../results/step2_HALPER_results/ \
 -s Mouse \
 -t Human \
 -c ../data/Alignment/10plusway-master.hal

# Mouse Liver to Human
bash halper_map_peak_orthologs.sh \
 -b ../data/OCR_idr_conservative/Mouse_liver_peak/idr.conservative_peak.narrowPeak \
 -o ../results/step2_HALPER_results/ \
 -s Mouse \
 -t Human \
 -c ../data/Alignment/10plusway-master.hal

# Step 2a

module load bedtools/2.30.0

# Intersect human adrenal and mouse
# open orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Adrenal_HumanToMouse.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Human_adrenal_peak/idr.conservative_peak.narrowPeak \
 -wa > ../results/step2a_bedtools_results/human_adrenal_peaks_with_open_mouse_orthologs.bed

# closed orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Adrenal_HumanToMouse.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Human_adrenal_peak/idr.conservative_peak.narrowPeak \
 -v > ../results/step2a_bedtools_results/human_adrenal_peaks_with_closed_mouse_orthologs.bed

# Intersect human liver and mouse
# open orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Liver_HumanToMouse.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Human_liver_peak/idr.conservative_peak.narrowPeak \
 -wa > ../results/step2a_bedtools_results/human_liver_peaks_with_open_mouse_orthologs.bed

# closed orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Liver_HumanToMouse.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Human_liver_peak/idr.conservative_peak.narrowPeak \
 -v > ../results/step2a_bedtools_results/human_liver_peaks_with_closed_mouse_orthologs.bed

# Intersect mouse adrenal and human
# open orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Adrenal_MouseToHuman.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Mouse_adrenal_peak/idr.conservative_peak.narrowPeak \
 -wa > ../results/step2a_bedtools_results/Mouse_adrenal_peaks_with_open_human_orthologs.bed

# closed orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Adrenal_MouseToHuman.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Mouse_adrenal_peak/idr.conservative_peak.narrowPeak \
 -v > ../results/step2a_bedtools_results/Mouse_adrenal_peaks_with_closed_human_orthologs.bed

# Intersect mouse liver and human
# open orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Liver_MouseToHuman.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Mouse_liver_peak/idr.conservative_peak.narrowPeak \
 -wa > ../results/step2a_bedtools_results/Mouse_liver_peaks_with_open_human_orthologs.bed

# closed orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Liver_MouseToHuman.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Mouse_liver_peak/idr.conservative_peak.narrowPeak \
 -v > ../results/step2a_bedtools_results/Mouse_liver_peaks_with_closed_human_orthologs.bed

# Step 3
# cross species
bedtools intersect \ 
 -a ../results/step2_HALPER_results/idr.conservative_peak.Adrenal_MouseToHuman.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Human_adrenal_peak/idr.conservative_peak.narrowPeak \
 -wa -u > ../results/step3_comparison/cross_species/mouse_adrenal_shared_with_human_adrenal.bed

bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Liver_MouseToHuman.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Human_liver_peak/idr.conservative_peak.narrowPeak \
 -wa -u > ../results/step3_comparison/cross_species/mouse_liver_shared_with_human_liver.bed

# calculate percentages


# cross tissues
bedtools intersect \
 -a ../data/OCR_idr_conservative/Mouse_liver_peak/idr.conservative_peak.narrowPeak \
 -b ../data/OCR_idr_conservative/Mouse_adrenal_peal/idr.conservative_peak.narrowPeak \
 > ../results/step3_comparison/cross_tissues/mouse_shared_across_tissues.bed

bedtools intersect -v \
 -a ../data/OCR_idr_conservative/Mouse_liver_peak/idr.conservative_peak.narrowPeak \
 -b ../data/OCR_idr_conservative/Mouse_adrenal_peak/idr.conservative_peak.narrowPeak \
 > ../results/step3_comparison/cross_tissues/mouse_liver_specific_OCR.bed

bedtools intersect -v \
 -a ../data/OCR_idr_conservative/Mouse_adrenal_peak/idr.conservative_peak.narrowPeak \
 -b ../data/OCR_idr_conservative/Mouse_liver_peak/idr.conservative_peak.narrowPeak \
 > ../results/step3_comparison/cross_tissues/mouse_adrenal_specific_OCR.bed

bedtools intersect \
 -a ../data/OCR_idr_conservative/Human_liver_peak/idr.conservative_peak.narrowPeak \
 -b ../data/OCR_idr_conservative/Human_adrenal_peal/idr.conservative_peak.narrowPeak \
 > ../results/step3_comparison/cross_tissues/human_shared_across_tissues.bed

bedtools intersect -v \
 -a ../data/OCR_idr_conservative/Human_liver_peak/idr.conservative_peak.narrowPeak \
 -b ../data/OCR_idr_conservative/Human_adrenal_peak/idr.conservative_peak.narrowPeak \
 > ../results/step3_comparison/cross_tissues/human_liver_specific_OCR.bed

bedtools intersect -v \
 -a ../data/OCR_idr_conservative/Human_adrenal_peak/idr.conservative_peak.narrowPeak \
 -b ../data/OCR_idr_conservative/Human_liver_peak/idr.conservative_peak.narrowPeak \
 > ../results/step3_comparison/cross_tissues/human_adrenal_specific_OCR.bed

# calculate percentages

# Step 4
# run outside of the pipeline with the GREAT webserver

# Step 5
# run outside of the pipeline with results from the GREAT webserver

# Step 6
# skipping for the rough draft submission
