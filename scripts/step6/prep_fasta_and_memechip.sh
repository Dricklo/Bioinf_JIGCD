#!/bin/bash

# Usage: bash prep_fasta_and_memechip.sh <BEDFILE> <GENOMEFA> <OUT_PREFIX>

# Capture input arguments
BEDFILE="$1"
GENOMEFA="$2"
OUT_PREFIX="$3"

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: bash prepare_fasta_for_memechip.sh <BEDFILE> <GENOMEFA> <OUT_PREFIX>"
    exit 1
fi

# Make sure the output directory exists
OUT_DIR=$(dirname "$OUT_PREFIX")
mkdir -p "$OUT_DIR"

# Step 1: Filter out regions smaller than 6bp (if needed)
awk '($3-$2) >= 6' "$BEDFILE" > "${OUT_PREFIX}_filtered.bed"

# Step 2: Resize regions to 256bp centered at midpoint
awk '{
    mid = int(($2 + $3)/2);
    start = mid - 128;
    end = mid + 128;
    if (start < 0) start = 0;  # prevent negative coordinates
    print $1 "\t" start "\t" end
}' "${OUT_PREFIX}_filtered.bed" > "${OUT_PREFIX}_256bp.bed"

# Step 3: Extract fasta sequences
bedtools getfasta -fi "$GENOMEFA" -bed "${OUT_PREFIX}_256bp.bed" -fo "${OUT_PREFIX}_256bp.fasta"

# Step 4: Run MEME-CHIP, load module first
module load MEME-suite/5.4.1  # (if not already loaded)

# Create a MEME-CHIP output directory
MEMECHIP_OUTDIR="${OUT_PREFIX}_memechip"

# Run meme-chip
meme-chip "${OUT_PREFIX}_256bp.fasta" -oc "$MEMECHIP_OUTDIR"

# Print completion message
echo "MEME-ChIP analysis for $BEDFILE complete. Results in $MEMECHIP_OUTDIR"