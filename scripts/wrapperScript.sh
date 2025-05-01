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

# Step 3, cross tissues
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

# Step 5: Divide data into enhancers and promoters
# For this step, we will divide the previously produced .bed files into enhancers and promoters using a 5kb threshold
# These steps will tie into Step 6

echo "Starting Step 5: Identifying enhancers for motif analysis..."

# Create output directories
mkdir -p ../results/step5/enhancers_promoters

# Step 5 (continued)
# Let's start with preparation for step 6b, we do need to perform an intiial data cleansing step to sort/curtail some files.
# Process Human Adrenal peaks
echo "Processing Human Adrenal peaks..."
cat ../data/OCR_idr_conservative/Human_adrenal_peak/idr.conservative_peak.narrowPeak | \
  cut -f1-3 | sort -k1,1 -k2,2n -u > ../results/step5/human_adrenal_sorted.bed

# Find enhancers (>5kb from TSS)
bedtools closest -a ../results/step5/human_adrenal_sorted.bed \
  -b ../data/TSS_sorted/gencode.Human.v27.annotation._TSSsWithStrand_sorted.bed -d | \
  awk '$10>5000' > ../results/step5/enhancers_promoters/human_adrenal_enhancers.bed

# Process Human Liver peaks
echo "Processing Human Liver peaks..."
cat ../data/OCR_idr_conservative/Human_liver_peak/idr.conservative_peak.narrowPeak | \
  cut -f1-3 | sort -k1,1 -k2,2n -u > ../results/step5/human_liver_sorted.bed

# Find enhancers (>5kb from TSS)
bedtools closest -a ../results/step5/human_liver_sorted.bed \
  -b ../data/TSS_sorted/gencode.Human.v27.annotation._TSSsWithStrand_sorted.bed -d | \
  awk '$10>5000' > ../results/step5/enhancers_promoters/human_liver_enhancers.bed

# Process Mouse Adrenal peaks
echo "Processing Mouse Adrenal peaks..."
cat ../data/OCR_idr_conservative/Mouse_adrenal_peak/idr.conservative_peak.narrowPeak | \
  cut -f1-3 | sort -k1,1 -k2,2n -u > ../results/step5/mouse_adrenal_sorted.bed

# Find enhancers (>5kb from TSS)
bedtools closest -a ../results/step5/mouse_adrenal_sorted.bed \
  -b ../data/TSS_sorted/gencode.Mouse.vM15.annotation_TSSWithStrand_sorted.bed -d | \
  awk '$10>5000' > ../results/step5/enhancers_promoters/mouse_adrenal_enhancers.bed

# Process Mouse Liver peaks
echo "Processing Mouse Liver peaks..."
cat ../data/OCR_idr_conservative/Mouse_liver_peak/idr.conservative_peak.narrowPeak | \
  cut -f1-3 | sort -k1,1 -k2,2n -u > ../results/step5/mouse_liver_sorted.bed

# Find enhancers (>5kb from TSS)
bedtools closest -a ../results/step5/mouse_liver_sorted.bed \
  -b ../data/TSS_sorted/gencode.Mouse.vM15.annotation_TSSWithStrand_sorted.bed -d | \
  awk '$10>5000' > ../results/step5/enhancers_promoters/mouse_liver_enhancers.bed


# STep 5 (continued)
# Preparation for step 6c, similar to before except threshold for awk is reversed:
# Now we are finding promoters. 

echo "Identifying promoters for Step 6c..."

# Process Human Adrenal peaks for promoters
echo "Processing Human Adrenal promoters..."
bedtools closest -a ../results/step5/human_adrenal_sorted.bed \
  -b ../data/TSS_sorted/gencode.Human.v27.annotation._TSSsWithStrand_sorted.bed -d | \
  awk '$10<=5000' > ../results/step5/enhancers_promoters/human_adrenal_promoters.bed

# Process Human Liver peaks for promoters
echo "Processing Human Liver promoters..."
bedtools closest -a ../results/step5/human_liver_sorted.bed \
  -b ../data/TSS_sorted/gencode.Human.v27.annotation._TSSsWithStrand_sorted.bed -d | \
  awk '$10<=5000' > ../results/step5/enhancers_promoters/human_liver_promoters.bed

# Process Mouse Adrenal peaks for promoters
echo "Processing Mouse Adrenal promoters..."
bedtools closest -a ../results/step5/mouse_adrenal_sorted.bed \
  -b ../data/TSS_sorted/gencode.Mouse.vM15.annotation_TSSWithStrand_sorted.bed -d | \
  awk '$10<=5000' > ../results/step5/enhancers_promoters/mouse_adrenal_promoters.bed

# Process Mouse Liver peaks for promoters
echo "Processing Mouse Liver promoters..."
bedtools closest -a ../results/step5/mouse_liver_sorted.bed \
  -b ../data/TSS_sorted/gencode.Mouse.vM15.annotation_TSSWithStrand_sorted.bed -d | \
  awk '$10<=5000' > ../results/step5/enhancers_promoters/mouse_liver_promoters.bed

echo "Step 5 promoter identification complete."

# Step 5 (continued)
# Now, preparation for Step 6d: Identify enhancers shared across tissues in each species
echo "Processing Step 6d: Enhancers shared across tissues in each species..."

# Process human shared peaks across tissues
echo "Processing human shared enhancers across tissues..."
cat ../results/step3_comparison/cross_tissues/human_shared_across_tissues.bed | \
  cut -f1-3 | sort -k1,1 -k2,2n -u > ../results/step5/human_shared_across_tissues_sorted.bed

# Find enhancers shared across tissues
bedtools closest -a ../results/step5/human_shared_across_tissues_sorted.bed \
  -b ../data/TSS_sorted/gencode.Human.v27.annotation._TSSsWithStrand_sorted.bed -d | \
  awk '$10>5000' > ../results/step5/enhancers_promoters/human_shared_enhancers.bed

# Process mouse shared peaks across tissues
echo "Processing mouse shared enhancers across tissues..."
cat ../results/step3_comparison/cross_tissues/mouse_shared_across_tissues.bed | \
  cut -f1-3 | sort -k1,1 -k2,2n -u > ../results/step5/mouse_shared_across_tissues_sorted.bed

# Find enhancers shared across tissues
bedtools closest -a ../results/step5/mouse_shared_across_tissues_sorted.bed \
  -b ../data/TSS_sorted/gencode.Mouse.vM15.annotation_TSSWithStrand_sorted.bed -d | \
  awk '$10>5000' > ../results/step5/enhancers_promoters/mouse_shared_enhancers.bed

echo "Step 6d preparation complete: Files for enhancers shared across tissues created"

# Step 5 (continued)
# Now, preparation for step 6e: Enhancers that are specific to each tissue in each species

echo "Processing Step 6e: Tissue-specific enhancers in each species..."

# Process human adrenal-specific peaks
echo "Processing human adrenal-specific enhancers..."
cat ../results/step3_comparison/cross_tissues/human_adrenal_specific_OCR.bed | \
  cut -f1-3 | sort -k1,1 -k2,2n -u > ../results/step5/human_adrenal_specific_sorted.bed

# Find human adrenal-specific enhancers
bedtools closest -a ../results/step5/human_adrenal_specific_sorted.bed \
  -b ../data/TSS_sorted/gencode.Human.v27.annotation._TSSsWithStrand_sorted.bed -d | \
  awk '$10>5000' > ../results/step5/enhancers_promoters/human_adrenal_specific_enhancers.bed

# Process human liver-specific peaks
echo "Processing human liver-specific enhancers..."
cat ../results/step3_comparison/cross_tissues/human_liver_specific_OCR.bed | \
  cut -f1-3 | sort -k1,1 -k2,2n -u > ../results/step5/human_liver_specific_sorted.bed

# Find human liver-specific enhancers
bedtools closest -a ../results/step5/human_liver_specific_sorted.bed \
  -b ../data/TSS_sorted/gencode.Human.v27.annotation._TSSsWithStrand_sorted.bed -d | \
  awk '$10>5000' > ../results/step5/enhancers_promoters/human_liver_specific_enhancers.bed

# Process mouse adrenal-specific peaks
echo "Processing mouse adrenal-specific enhancers..."
cat ../results/step3_comparison/cross_tissues/mouse_adrenal_specific_OCR.bed | \
  cut -f1-3 | sort -k1,1 -k2,2n -u > ../results/step5/mouse_adrenal_specific_sorted.bed

# Find mouse adrenal-specific enhancers
bedtools closest -a ../results/step5/mouse_adrenal_specific_sorted.bed \
  -b ../data/TSS_sorted/gencode.Mouse.vM15.annotation_TSSWithStrand_sorted.bed -d | \
  awk '$10>5000' > ../results/step5/enhancers_promoters/mouse_adrenal_specific_enhancers.bed

# Process mouse liver-specific peaks
echo "Processing mouse liver-specific enhancers..."
cat ../results/step3_comparison/cross_tissues/mouse_liver_specific_OCR.bed | \
  cut -f1-3 | sort -k1,1 -k2,2n -u > ../results/step5/mouse_liver_specific_sorted.bed

# Find mouse liver-specific enhancers
bedtools closest -a ../results/step5/mouse_liver_specific_sorted.bed \
  -b ../data/TSS_sorted/gencode.Mouse.vM15.annotation_TSSWithStrand_sorted.bed -d | \
  awk '$10>5000' > ../results/step5/enhancers_promoters/mouse_liver_specific_enhancers.bed

echo "Step 6e preparation complete: Tissue-specific enhancer files created"

# Step 5 (continued)
# Now, preparation for step 6f: Enhancers that are shared across species for each tissue

echo "Processing Step 6f: Enhancers shared across species for each tissue..."

# Process shared adrenal peaks across species
echo "Processing adrenal enhancers shared across species..."
cat ../results/step3_comparison/cross_species/mouse_adrenal_shared_with_human_adrenal.bed | \
  cut -f1-3 | sort -k1,1 -k2,2n -u > ../results/step5/adrenal_shared_across_species_sorted.bed

# Find enhancers shared across species for adrenal tissue
bedtools closest -a ../results/step5/adrenal_shared_across_species_sorted.bed \
  -b ../data/TSS_sorted/gencode.Mouse.vM15.annotation_TSSWithStrand_sorted.bed -d | \
  awk '$10>5000' > ../results/step5/enhancers_promoters/adrenal_shared_enhancers.bed

# Process shared liver peaks across species
echo "Processing liver enhancers shared across species..."
cat ../results/step3_comparison/cross_species/mouse_liver_shared_with_human_liver.bed | \
  cut -f1-3 | sort -k1,1 -k2,2n -u > ../results/step5/liver_shared_across_species_sorted.bed

# Find enhancers shared across species for liver tissue
bedtools closest -a ../results/step5/liver_shared_across_species_sorted.bed \
  -b ../data/TSS_sorted/gencode.Mouse.vM15.annotation_TSSWithStrand_sorted.bed -d | \
  awk '$10>5000' > ../results/step5/enhancers_promoters/liver_shared_enhancers.bed

echo "Step 6f preparation complete: Files for enhancers shared across species created"

# Step 5 (continued)
# Now, preparation for step 6g: Enhancers that are shared across species for each tissue

echo "Processing Step 6g: Species-specific enhancers for each tissue..."

# Process mouse adrenal-specific peaks
echo "Processing mouse-specific adrenal enhancers..."
cat ../results/step3_comparison/cross_species_different/mouse_adrenal_different_than_human_small.bed | \
  cut -f1-3 | sort -k1,1 -k2,2n -u > ../results/step5/mouse_specific_adrenal_sorted.bed

# Find mouse-specific adrenal enhancers
bedtools closest -a ../results/step5/mouse_specific_adrenal_sorted.bed \
  -b ../data/TSS_sorted/gencode.Mouse.vM15.annotation_TSSWithStrand_sorted.bed -d | \
  awk '$10>5000' > ../results/step5/enhancers_promoters/mouse_specific_adrenal_enhancers.bed

# Process mouse liver-specific peaks
echo "Processing mouse-specific liver enhancers..."
cat ../results/step3_comparison/cross_species_different/mouse_liver_different_than_human_small.bed | \
  cut -f1-3 | sort -k1,1 -k2,2n -u > ../results/step5/mouse_specific_liver_sorted.bed

# Find mouse-specific liver enhancers
bedtools closest -a ../results/step5/mouse_specific_liver_sorted.bed \
  -b ../data/TSS_sorted/gencode.Mouse.vM15.annotation_TSSWithStrand_sorted.bed -d | \
  awk '$10>5000' > ../results/step5/enhancers_promoters/mouse_specific_liver_enhancers.bed

echo "Step 6g preparation complete: Species-specific enhancer files created"


# Finally, with preparation for step 6 complete, we can used the resulting bed files and move to step 6:
# Step 6: Find sequence patterns that occur more than expected by chance
echo "Starting motif analysis pipeline for step 6..."

# Define paths for genome reference files, previously downloaded and detailed in readme. 
HUMAN_GENOME="../fa_reference/hg38.fa"
MOUSE_GENOME="../fa_reference/mm10.fa"

# Create output directories
mkdir -p ../results/step6_motifs/6b
mkdir -p ../results/step6_motifs/6c
mkdir -p ../results/step6_motifs/6d
mkdir -p ../results/step6_motifs/6e
mkdir -p ../results/step6_motifs/6f

# Step 6b: Enhancers for each species, tissue combination
echo "Processing Step 6b: Enhancers for each species-tissue combination"

# Human Adrenal enhancers
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/human_adrenal_enhancers.bed" \
    "${HUMAN_GENOME}" \
    "../results/step6_motifs/6b/human_adrenal_enhancers"

# Human Liver enhancers
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/human_liver_enhancers.bed" \
    "${HUMAN_GENOME}" \
    "../results/step6_motifs/6b/human_liver_enhancers"

# Mouse Adrenal enhancers
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/mouse_adrenal_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "../results/step6_motifs/6b/mouse_adrenal_enhancers"

# Mouse Liver enhancers
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/mouse_liver_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "../results/step6_motifs/6b/mouse_liver_enhancers"

# Step 6c: Promoters for each species, tissue combination
echo "Processing Step 6c: Promoters for each species-tissue combination"

# Human Adrenal promoters
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/human_adrenal_promoters.bed" \
    "${HUMAN_GENOME}" \
    "../results/step6_motifs/6c/human_adrenal_promoters"

# Human Liver promoters
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/human_liver_promoters.bed" \
    "${HUMAN_GENOME}" \
    "../results/step6_motifs/6c/human_liver_promoters"

# Mouse Adrenal promoters
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/mouse_adrenal_promoters.bed" \
    "${MOUSE_GENOME}" \
    "../results/step6_motifs/6c/mouse_adrenal_promoters"

# Mouse Liver promoters
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/mouse_liver_promoters.bed" \
    "${MOUSE_GENOME}" \
    "../results/step6_motifs/6c/mouse_liver_promoters"

# Step 6d: Enhancers that are shared across tissues in each species
echo "Processing Step 6d: Enhancers shared across tissues in each species"

# Human shared enhancers
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/human_shared_enhancers.bed" \
    "${HUMAN_GENOME}" \
    "../results/step6_motifs/6d/human_shared_enhancers"

# Mouse shared enhancers
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/mouse_shared_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "../results/step6_motifs/6d/mouse_shared_enhancers"

# Step 6e: Enhancers that are specific to each tissue in each species
echo "Processing Step 6e: Tissue-specific enhancers in each species"

# Human adrenal-specific enhancers
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/human_adrenal_specific_enhancers.bed" \
    "${HUMAN_GENOME}" \
    "../results/step6_motifs/6e/human_adrenal_specific_enhancers"

# Human liver-specific enhancers
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/human_liver_specific_enhancers.bed" \
    "${HUMAN_GENOME}" \
    "../results/step6_motifs/6e/human_liver_specific_enhancers"

# Mouse adrenal-specific enhancers
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/mouse_adrenal_specific_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "../results/step6_motifs/6e/mouse_adrenal_specific_enhancers"

# Mouse liver-specific enhancers
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/mouse_liver_specific_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "../results/step6_motifs/6e/mouse_liver_specific_enhancers"

# Step 6f: Enhancers that are shared across species for each tissue
echo "Processing Step 6f: Enhancers shared across species for each tissue"

# Adrenal shared enhancers
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/adrenal_shared_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "../results/step6_motifs/6f/adrenal_shared_enhancers"

# Liver shared enhancers
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/liver_shared_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "../results/step6_motifs/6f/liver_shared_enhancers"
	
# Step 6g: Enhancers that are specific to each species for each tissue
echo "Processing Step 6g: Species-specific enhancers for each tissue"

# Create output directory if it doesn't exist
mkdir -p ../results/step6_motifs/6g

# Mouse-specific adrenal enhancers
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/mouse_specific_adrenal_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "../results/step6_motifs/6g/mouse_specific_adrenal_enhancers"

# Mouse-specific liver enhancers
bash ../scripts/step6/prep_fasta_and_memechip.sh \
    "../results/step5/enhancers_promoters/mouse_specific_liver_enhancers.bed" \
    "${MOUSE_GENOME}" \
    "../results/step6_motifs/6g/mouse_specific_liver_enhancers"

echo "Step 6 motif analysis complete!"