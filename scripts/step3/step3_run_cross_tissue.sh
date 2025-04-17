#!/bin/bash
#SBATCH --job-name=step3_cross_tissue
#SBATCH --output=step3_cross_tissue.out
#SBATCH --error=step3_cross_tissue.err
#SBATCH --time=01:00:00
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=2000M

# Load bedtools module
module load bedtools/2.30.0

# Set output directory
OUTPUT_DIR="/ocean/projects/bio230007p/jhou3/repos/results/step3_comparison/cross_tissue"
mkdir -p $OUTPUT_DIR

echo "Starting Step 3: Cross-Tissue Comparison"
echo "Output will be saved to: $OUTPUT_DIR"
echo


# Mouse: liver vs adrenal

bedtools intersect \
  -a /ocean/projects/bio230007p/jhou3/repos/raw_data/Mouse_liver_peak/idr.conservative_peak.narrowPeak.gz \
  -b /ocean/projects/bio230007p/jhou3/repos/raw_data/Mouse_adrenal_peak/idr.conservative_peak.narrowPeak.gz \
  > $OUTPUT_DIR/mouse_shared_across_tissues.bed

bedtools intersect -v \
  -a /ocean/projects/bio230007p/jhou3/repos/raw_data/Mouse_liver_peak/idr.conservative_peak.narrowPeak.gz \
  -b /ocean/projects/bio230007p/jhou3/repos/raw_data/Mouse_adrenal_peak/idr.conservative_peak.narrowPeak.gz \
  > $OUTPUT_DIR/mouse_liver_specific_OCR.bed

bedtools intersect -v \
  -a /ocean/projects/bio230007p/jhou3/repos/raw_data/Mouse_adrenal_peak/idr.conservative_peak.narrowPeak.gz \
  -b /ocean/projects/bio230007p/jhou3/repos/raw_data/Mouse_liver_peak/idr.conservative_peak.narrowPeak.gz \
  > $OUTPUT_DIR/mouse_adrenal_specific_OCR.bed


# Human: liver vs adrenal

bedtools intersect \
  -a /ocean/projects/bio230007p/jhou3/repos/raw_data/Human_liver_peak/idr.conservative_peak.narrowPeak.gz \
  -b /ocean/projects/bio230007p/jhou3/repos/raw_data/Human_adrenal_peak/idr.conservative_peak.narrowPeak.gz \
  > $OUTPUT_DIR/human_shared_across_tissues.bed

bedtools intersect -v \
  -a /ocean/projects/bio230007p/jhou3/repos/raw_data/Human_liver_peak/idr.conservative_peak.narrowPeak.gz \
  -b /ocean/projects/bio230007p/jhou3/repos/raw_data/Human_adrenal_peak/idr.conservative_peak.narrowPeak.gz \
  > $OUTPUT_DIR/human_liver_specific_OCR.bed

bedtools intersect -v \
  -a /ocean/projects/bio230007p/jhou3/repos/raw_data/Human_adrenal_peak/idr.conservative_peak.narrowPeak.gz \
  -b /ocean/projects/bio230007p/jhou3/repos/raw_data/Human_liver_peak/idr.conservative_peak.narrowPeak.gz \
  > $OUTPUT_DIR/human_adrenal_specific_OCR.bed


# Count summary

echo "Cross-tissue OCR counts:" > $OUTPUT_DIR/cross_tissue_counts.txt
for file in $OUTPUT_DIR/*.bed; do
    echo -n "$(basename $file): " >> $OUTPUT_DIR/cross_tissue_counts.txt
    wc -l < "$file" >> $OUTPUT_DIR/cross_tissue_counts.txt
done

echo "Step 3 cross-tissue comparison completed."


# Compute percentage summaries and save to file
echo "" >> $OUTPUT_DIR/cross_tissue_counts.txt
PERCENT_FILE="$OUTPUT_DIR/cross_tissue_percentages.txt"
echo "Cross-tissue OCR percentages:" > "$PERCENT_FILE"

# ========== MOUSE ==========
mouse_liver_total=$(zcat /ocean/projects/bio230007p/jhou3/repos/raw_data/Mouse_liver_peak/idr.conservative_peak.narrowPeak.gz | wc -l)
mouse_adrenal_total=$(zcat /ocean/projects/bio230007p/jhou3/repos/raw_data/Mouse_adrenal_peak/idr.conservative_peak.narrowPeak.gz | wc -l)
mouse_total=$((mouse_liver_total + mouse_adrenal_total))

mouse_shared=$(wc -l < $OUTPUT_DIR/mouse_shared_across_tissues.bed)
mouse_liver_spec=$(wc -l < $OUTPUT_DIR/mouse_liver_specific_OCR.bed)
mouse_adrenal_spec=$(wc -l < $OUTPUT_DIR/mouse_adrenal_specific_OCR.bed)

echo "Mouse Shared: $(echo "scale=2; 100*$mouse_shared/$mouse_total" | bc)% ($mouse_shared/$mouse_total)" >> "$PERCENT_FILE"
echo "Mouse Liver-Specific: $(echo "scale=2; 100*$mouse_liver_spec/$mouse_total" | bc)% ($mouse_liver_spec/$mouse_total)" >> "$PERCENT_FILE"
echo "Mouse Adrenal-Specific: $(echo "scale=2; 100*$mouse_adrenal_spec/$mouse_total" | bc)% ($mouse_adrenal_spec/$mouse_total)" >> "$PERCENT_FILE"

# ========== HUMAN ==========
human_liver_total=$(zcat /ocean/projects/bio230007p/jhou3/repos/raw_data/Human_liver_peak/idr.conservative_peak.narrowPeak.gz | wc -l)
human_adrenal_total=$(zcat /ocean/projects/bio230007p/jhou3/repos/raw_data/Human_adrenal_peak/idr.conservative_peak.narrowPeak.gz | wc -l)
human_total=$((human_liver_total + human_adrenal_total))

human_shared=$(wc -l < $OUTPUT_DIR/human_shared_across_tissues.bed)
human_liver_spec=$(wc -l < $OUTPUT_DIR/human_liver_specific_OCR.bed)
human_adrenal_spec=$(wc -l < $OUTPUT_DIR/human_adrenal_specific_OCR.bed)

echo "Human Shared: $(echo "scale=2; 100*$human_shared/$human_total" | bc)% ($human_shared/$human_total)" >> "$PERCENT_FILE"
echo "Human Liver-Specific: $(echo "scale=2; 100*$human_liver_spec/$human_total" | bc)% ($human_liver_spec/$human_total)" >> "$PERCENT_FILE"
echo "Human Adrenal-Specific: $(echo "scale=2; 100*$human_adrenal_spec/$human_total" | bc)% ($human_adrenal_spec/$human_total)" >> "$PERCENT_FILE"

