#!/bin/bash

#SBATCH -p RM-shared
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --ntasks=16
#SBATCH --account bio230007p
#SBATCH -J halper_liver_m2h
#SBATCH -e /ocean/projects/bio230007p/jhou3/repos/results/halper_liver_m2h.err
#SBATCH -o /ocean/projects/bio230007p/jhou3/repos/results/halper_liver_m2h.out
#SBATCH --export=ALL

set -x

bash /ocean/projects/bio230007p/jhou3/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh \
  -b /ocean/projects/bio230007p/jhou3/repos/raw_data/idr.conservative_peak.narrowPeak \
  -o /ocean/projects/bio230007p/jhou3/repos/results/halper_liver_m2h_out \
  -s Mouse \
  -t Human \
  -c /ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal

