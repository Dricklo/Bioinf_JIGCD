#!/bin/bash

#SBATCH -p RM-shared
#SBATCH --cpus-per-task=1
#SBATCH -o halper_%j.out
#SBATCH -e halper_%j.err
#SBATCH -t 12:00:00

# Set up environment, have to SPECIFICALLY specify the export path and python path otherwise i get the halLiftover/ortholog command not found error
export PATH=/jet/home/lod/repos/hal/bin:${PATH}

export PYTHONPATH=/jet/home/lod/repos/halLiftover-postprocessing:${PYTHONPATH}

source $(conda info --base)/etc/profile.d/conda.sh

conda activate hal



# Run the script with this command, adapated from Isabelle and Jinqi

bash repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh \
-b Human_peak/idr.conservative_peak.narrowPeak.gz \
-o result/ \
-s Human \
-t Mouse \
-c /ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal
