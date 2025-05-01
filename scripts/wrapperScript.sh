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

# Intersect human adrenal and mouse - human adrenal orthologs in the mouse are open & closed
# open orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Adrenal_HumanToMouse.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Mouse_adrenal_peak/idr.conservative_peak.narrowPeak \
 -wa > ../results/step2a_bedtools_results/human_adrenal_peaks_with_open_mouse_orthologs.bed

# closed orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Adrenal_HumanToMouse.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Human_adrenal_peak/idr.conservative_peak.narrowPeak \
 -v > ../results/step2a_bedtools_results/Mouse_adrenal_peaks_with_closed_mouse_orthologs.bed

# Intersect human liver and mouse - human liver orthologs in the mouse are open & closed
# open orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Liver_HumanToMouse.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Mouse_liver_peak/idr.conservative_peak.narrowPeak \
 -wa > ../results/step2a_bedtools_results/mouse_liver_peaks_with_open_mouse_orthologs.bed

# closed orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Liver_HumanToMouse.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Mouse_liver_peak/idr.conservative_peak.narrowPeak \
 -v > ../results/step2a_bedtools_results/mouse_liver_peaks_with_closed_mouse_orthologs.bed

# Intersect mouse adrenal and human - mouse adrenal orthologs in the human are open & closed
# open orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Adrenal_MouseToHuman.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Human_adrenal_peak/idr.conservative_peak.narrowPeak \
 -wa > ../results/step2a_bedtools_results/human_adrenal_peaks_with_open_human_orthologs.bed

# closed orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Adrenal_MouseToHuman.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Human_adrenal_peak/idr.conservative_peak.narrowPeak \
 -v > ../results/step2a_bedtools_results/human_adrenal_peaks_with_closed_human_orthologs.bed

# Intersect mouse liver and human - mouse liver orthologs in the human are open & closed
# open orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Liver_MouseToHuman.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Human_liver_peak/idr.conservative_peak.narrowPeak \
 -wa > ../results/step2a_bedtools_results/human_liver_peaks_with_open_human_orthologs.bed

# closed orthologs
bedtools intersect \
 -a ../results/step2_HALPER_results/idr.conservative_peak.Liver_MouseToHuman.HALPER.narrowPeak \
 -b ../data/OCR_idr_conservative/Human_liver_peak/idr.conservative_peak.narrowPeak \
 -v > ../results/step2a_bedtools_results/human_liver_peaks_with_closed_human_orthologs.bed

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

# Step 3 - calculate cross-species shared peak percentages

# set the path 
adrenal_shared=../results/step3_comparison/cross_species/mouse_adrenal_shared_with_human_adrenal.bed
adrenal_total=../results/step2_HALPER_results/idr.conservative_peak.Adrenal_MouseToHuman.HALPER.narrowPeak

liver_shared=../results/step3_comparison/cross_species/mouse_liver_shared_with_human_liver.bed
liver_total=../results/step2_HALPER_results/idr.conservative_peak.Liver_MouseToHuman.HALPER.narrowPeak

# Count number of peaks
adrenal_shared_count=$(wc -l < "$adrenal_shared")
adrenal_total_count=$(wc -l < "$adrenal_total")
liver_shared_count=$(wc -l < "$liver_shared")
liver_total_count=$(wc -l < "$liver_total")

# calculate percentages
adrenal_pct=$(awk -v a="$adrenal_shared_count" -v b="$adrenal_total_count" 'BEGIN {printf "%.2f", a / b * 100}')
liver_pct=$(awk -v a="$liver_shared_count" -v b="$liver_total_count" 'BEGIN {printf "%.2f", a / b * 100}')

# out put the file
output_file=../results/step3_comparison/cross_species/cross_species_percentages.txt
mkdir -p "$(dirname "$output_file")"
cat <<EOF > "$output_file"
[RESULT] Shared peak percentages (cross-species)
Mouse liver vs Human liver: ${liver_pct}% (${liver_shared_count} / ${liver_total_count})
Mouse adrenal vs Human adrenal: ${adrenal_pct}% (${adrenal_shared_count} / ${adrenal_total_count})
EOF

echo "[Done] Percentages written to $output_file"



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

# Define file paths
mouse_shared=../results/step3_comparison/cross_tissues/mouse_shared_across_tissues.bed
mouse_liver_specific=../results/step3_comparison/cross_tissues/mouse_liver_specific_OCR.bed
mouse_adrenal_specific=../results/step3_comparison/cross_tissues/mouse_adrenal_specific_OCR.bed

human_shared=../results/step3_comparison/cross_tissues/human_shared_across_tissues.bed
human_liver_specific=../results/step3_comparison/cross_tissues/human_liver_specific_OCR.bed
human_adrenal_specific=../results/step3_comparison/cross_tissues/human_adrenal_specific_OCR.bed

# Define original total peak files
mouse_total=../data/OCR_idr_conservative/Mouse_liver_peak/idr.conservative_peak.narrowPeak
human_total=../data/OCR_idr_conservative/Human_liver_peak/idr.conservative_peak.narrowPeak

# Count number of peaks
ms=$(wc -l < "$mouse_shared")
ml=$(wc -l < "$mouse_liver_specific")
ma=$(wc -l < "$mouse_adrenal_specific")
mt=$(wc -l < "$mouse_total")  # Assume total = Mouse liver peak count

hs=$(wc -l < "$human_shared")
hl=$(wc -l < "$human_liver_specific")
ha=$(wc -l < "$human_adrenal_specific")
ht=$(wc -l < "$human_total")  # Assume total = Human liver peak count

# Function to calculate percentage
pct() { awk -v a=$1 -v b=$2 'BEGIN { printf "%.2f", a / b * 100 }'; }

# Calculate percentages
msp=$(pct $ms $mt)
mlp=$(pct $ml $mt)
map=$(pct $ma $mt)

hsp=$(pct $hs $ht)
hlp=$(pct $hl $ht)
hap=$(pct $ha $ht)

# Write output
output_file=../results/step3_comparison/cross_tissues/cross_tissue_percentages.txt
mkdir -p "$(dirname "$output_file")"

cat <<EOF > "$output_file"
Cross-tissue OCR percentages:
Mouse Shared: ${msp}% (${ms}/${mt})
Mouse Liver-Specific: ${mlp}% (${ml}/${mt})
Mouse Adrenal-Specific: ${map}% (${ma}/${mt})
Human Shared: ${hsp}% (${hs}/${ht})
Human Liver-Specific: ${hlp}% (${hl}/${ht})
Human Adrenal-Specific: ${hap}% (${ha}/${ht})
EOF

echo "[Done] Cross-tissue OCR percentages written to $output_file"







# Step 4
# run outside of the pipeline with the GREAT webserver

# Step 5
# run outside of the pipeline with results from the GREAT webserver

# Step 6
# obtain the bed file for enhancer & promotor
bedtools closest \
-a human_adrenal_specific_OCR.bed \
-b gencode.Human.v27.annotation._TSSsWithStrand_sorted.bed \
-d | awk '$10>5000' > enhancers_human_adrenal.bed

bedtools closest \
-a human_liver_specific_OCR.bed \
-b gencode.Human.v27.annotation._TSSsWithStrand_sorted.bed \
-d | awk '$10>5000' > enhancers_human_liver.bed

bedtools closest \
-a mouse_adrenal_specific_OCR.bed \
-b gencode.Mouse.v27.annotation._TSSsWithStrand_sorted.bed \
-d | awk '$10>5000' > enhancers_mouse_adrenal.bed


bedtools closest \
-a mouse_liver_specific_OCR.bed \
-b gencode.Mouse.v27.annotation._TSSsWithStrand_sorted.bed \
-d | awk '$10>5000' > enhancers_mouse_liver.bed
