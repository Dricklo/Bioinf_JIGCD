#!/bin/bash -l
#SBATCH --job-name=step3_cross_species
#SBATCH --output=step3_cross_species.out
#SBATCH --error=step3_cross_species.err
#SBATCH --time=01:00:00
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=2000M

# Load bedtools
module load bedtools/2.30.0

# ==== Input Files ====

# Mouse liver mapped to human coordinates
MOUSE_LIVER_M2H="/ocean/projects/bio230007p/jhou3/repos/results/halper_liver_m2h_out/idr.conservative_peak.MouseToHuman.HALPER.narrowPeak"

# Mouse adrenal mapped to human coordinates
MOUSE_ADRENAL_M2H="/ocean/projects/bio230007p/jhou3/repos/results/halper_adrenal_m2h_out/idr.conservative_peak.MouseToHuman.HALPER.narrowPeak"

# Human liver original peaks
HUMAN_LIVER="/ocean/projects/bio230007p/jhou3/repos/raw_data/Human_liver_peak/idr.conservative_peak.narrowPeak"

# Human adrenal original peaks
HUMAN_ADRENAL="/ocean/projects/bio230007p/jhou3/repos/raw_data/Human_adrenal_peak/idr.conservative_peak.narrowPeak"

# ==== Output directory ====
OUTDIR="/ocean/projects/bio230007p/jhou3/repos/results/step3_comparison/cross_species"
mkdir -p $OUTDIR

# ==== Intersections ====

# Mouse liver vs Human liver
bedtools intersect -a $MOUSE_LIVER_M2H -b $HUMAN_LIVER -wa -u > $OUTDIR/mouse_liver_shared_with_human_liver.bed

# Mouse adrenal vs Human adrenal
bedtools intersect -a $MOUSE_ADRENAL_M2H -b $HUMAN_ADRENAL -wa -u > $OUTDIR/mouse_adrenal_shared_with_human_adrenal.bed

# ==== Count shared peaks and save ====

COUNT_FILE="$OUTDIR/cross_species_counts.txt"
echo "[RESULT] Shared peak counts (cross-species)" > $COUNT_FILE

echo -n "Mouse liver vs Human liver: " >> $COUNT_FILE
wc -l < $OUTDIR/mouse_liver_shared_with_human_liver.bed >> $COUNT_FILE

echo -n "Mouse adrenal vs Human adrenal: " >> $COUNT_FILE
wc -l < $OUTDIR/mouse_adrenal_shared_with_human_adrenal.bed >> $COUNT_FILE


# ==== Compute percentage of shared peaks and save ====
PERCENT_FILE="$OUTDIR/cross_species_percentages.txt"
echo "[RESULT] Shared peak percentages (cross-species)" > "$PERCENT_FILE"

# Mouse liver mapped to human vs Human liver total
human_liver_total=$(zcat /ocean/projects/bio230007p/jhou3/repos/raw_data/Human_liver_peak/idr.conservative_peak.narrowPeak.gz | wc -l)
mouse_liver_shared=$(wc -l < $OUTDIR/mouse_liver_shared_with_human_liver.bed)

percent_liver=$(echo "scale=2; 100 * $mouse_liver_shared / $human_liver_total" | bc)
echo "Mouse liver vs Human liver: $percent_liver% ($mouse_liver_shared / $human_liver_total)" >> "$PERCENT_FILE"

# Mouse adrenal mapped to human vs Human adrenal total
human_adrenal_total=$(zcat /ocean/projects/bio230007p/jhou3/repos/raw_data/Human_adrenal_peak/idr.conservative_peak.narrowPeak.gz | wc -l)
mouse_adrenal_shared=$(wc -l < $OUTDIR/mouse_adrenal_shared_with_human_adrenal.bed)

percent_adrenal=$(echo "scale=2; 100 * $mouse_adrenal_shared / $human_adrenal_total" | bc)
echo "Mouse adrenal vs Human adrenal: $percent_adrenal% ($mouse_adrenal_shared / $human_adrenal_total)" >> "$PERCENT_FILE"

