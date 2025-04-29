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
# Find sequence patterns that occur more than expected by chance

# Define paths for step 6
STEP6_DIR="${RESULTS_DIR}/step6_motif_analysis"

# Create output directories
mkdir -p ${STEP6_DIR}/6b
mkdir -p ${STEP6_DIR}/6c
mkdir -p ${STEP6_DIR}/6d
mkdir -p ${STEP6_DIR}/6e
mkdir -p ${STEP6_DIR}/6f
mkdir -p ${STEP6_DIR}/6g

# Define genome references
HUMAN_GENOME="${DATA_DIR}/genomes/hg38.fa"
MOUSE_GENOME="${DATA_DIR}/genomes/mm10.fa"

echo "Starting motif analysis pipeline for step 6..."

# Step 6b: Enhancers for each species, tissue combination
echo "Processing Step 6b: Enhancers for each species-tissue combination"

# Human Adrenal enhancers
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/promoters/human_adrenal_enhancers.bed" \
    "${HUMAN_GENOME}" \
    "${STEP6_DIR}/6b/human_adrenal_enhancers"

# Human Liver enhancers
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/enhancers/human_liver_enhancers.bed" \
    "${HUMAN_GENOME}" \
    "${STEP6_DIR}/6b/human_liver_enhancers"

# Mouse Adrenal enhancers
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/enhancers/mouse_adrenal_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "${STEP6_DIR}/6b/mouse_adrenal_enhancers"

# Mouse Liver enhancers
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/enhancers/mouse_liver_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "${STEP6_DIR}/6b/mouse_liver_enhancers"

# Step 6c: Promoters for each species, tissue combination
echo "Processing Step 6c: Promoters for each species-tissue combination"

# Human Adrenal promoters
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/promoters/human_adrenal_promoters.bed" \
    "${HUMAN_GENOME}" \
    "${STEP6_DIR}/6c/human_adrenal_promoters"

# Human Liver promoters
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/promoters/human_liver_promoters.bed" \
    "${HUMAN_GENOME}" \
    "${STEP6_DIR}/6c/human_liver_promoters"

# Mouse Adrenal promoters
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/promoters/mouse_adrenal_promoters.bed" \
    "${MOUSE_GENOME}" \
    "${STEP6_DIR}/6c/mouse_adrenal_promoters"

# Mouse Liver promoters
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/promoters/mouse_liver_promoters.bed" \
    "${MOUSE_GENOME}" \
    "${STEP6_DIR}/6c/mouse_liver_promoters"

# Step 6d: Enhancers that are shared across tissues in each species
echo "Processing Step 6d: Enhancers shared across tissues in each species"

# Human shared enhancers
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/enhancers/human_shared_enhancers.bed" \
    "${HUMAN_GENOME}" \
    "${STEP6_DIR}/6d/human_shared_enhancers"

# Mouse shared enhancers
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/enhancers/mouse_shared_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "${STEP6_DIR}/6d/mouse_shared_enhancers"

# Step 6e: Enhancers that are specific to each tissue in each species
echo "Processing Step 6e: Tissue-specific enhancers in each species"

# Human adrenal-specific enhancers
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/enhancers/human_adrenal_specific_enhancers.bed" \
    "${HUMAN_GENOME}" \
    "${STEP6_DIR}/6e/human_adrenal_specific_enhancers"

# Human liver-specific enhancers
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/enhancers/human_liver_specific_enhancers.bed" \
    "${HUMAN_GENOME}" \
    "${STEP6_DIR}/6e/human_liver_specific_enhancers"

# Mouse adrenal-specific enhancers
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/enhancers/mouse_adrenal_specific_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "${STEP6_DIR}/6e/mouse_adrenal_specific_enhancers"

# Mouse liver-specific enhancers
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/enhancers/mouse_liver_specific_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "${STEP6_DIR}/6e/mouse_liver_specific_enhancers"

# Step 6f: Enhancers that are shared across species for each tissue
echo "Processing Step 6f: Enhancers shared across species for each tissue"

# Adrenal shared enhancers (mouse-to-human mapping)
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/enhancers/adrenal_shared_enhancers.bed" \
    "${HUMAN_GENOME}" \
    "${STEP6_DIR}/6f/adrenal_shared_enhancers"

# Liver shared enhancers (mouse-to-human mapping)
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/enhancers/liver_shared_enhancers.bed" \
    "${HUMAN_GENOME}" \
    "${STEP6_DIR}/6f/liver_shared_enhancers"

# Step 6g: Enhancers that are specific to each species for each tissue
echo "Processing Step 6g: Species-specific enhancers for each tissue"

# Human-specific adrenal enhancers
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/enhancers/human_specific_adrenal_enhancers.bed" \
    "${HUMAN_GENOME}" \
    "${STEP6_DIR}/6g/human_specific_adrenal_enhancers"

# Human-specific liver enhancers
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/enhancers/human_specific_liver_enhancers.bed" \
    "${HUMAN_GENOME}" \
    "${STEP6_DIR}/6g/human_specific_liver_enhancers"

# Mouse-specific adrenal enhancers
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/enhancers/mouse_specific_adrenal_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "${STEP6_DIR}/6g/mouse_specific_adrenal_enhancers"

# Mouse-specific liver enhancers
bash prep_fasta_and_memechip.sh \
    "${RESULTS_DIR}/enhancers/mouse_specific_liver_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "${STEP6_DIR}/6g/mouse_specific_liver_enhancers"

echo "Step 6 motif analysis complete!"
