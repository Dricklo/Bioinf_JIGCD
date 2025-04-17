This file contains the results of running bedtools of a source species against
its converted orthologs from the other species. It contains both the open (overlapping) and 
closed (non-overlapping) chromatin regions between the two species. 

Input: Halper file of converted species (e.g, human_to_mouse), and original ATAC-seq data 
(mouse in this case). 
Output: bed file containing the overlapping (-wa) or non-overlapping (-v) regions.