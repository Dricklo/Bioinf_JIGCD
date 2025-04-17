#!/bin/bash
#SBATCH --job-name=mouse2human_liver_intersect
#SBATCH --output=mouse2human_liver_intersect.out
#SBATCH --error=mouse2human_liver_intersect.err
#SBATCH --time=01:00:00
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=2000M

# Load bedtools module
module load bedtools/2.30.0

# Define input files
# Human peak is used here because MAPPED_MOUSE is in human coordinates
HUMAN_PEAK="/ocean/projects/bio230007p/jhou3/repos/raw_data/Human_liver_peak/idr.conservative_peak.narrowPeak.gz"
MAPPED_MOUSE="/ocean/projects/bio230007p/jhou3/repos/results/halper_liver_m2h_out/idr.conservative_peak.MouseToHuman.HALPER.narrowPeak.gz"

# Define output directory
OUTPUT_DIR="/ocean/projects/bio230007p/jhou3/repos/results/step2a"
mkdir -p $OUTPUT_DIR

# Open orthologs: mouse OCRs whose human orthologs are open
bedtools intersect -a $MAPPED_MOUSE -b $HUMAN_PEAK -wa \
  > $OUTPUT_DIR/mouse_liver_peaks_with_open_human_orthologs.bed

# Closed orthologs: mouse OCRs whose human orthologs are not open
bedtools intersect -a $MAPPED_MOUSE -b $HUMAN_PEAK -v \
  > $OUTPUT_DIR/mouse_liver_peaks_with_closed_human_orthologs.bed


