This file should contain the results of running bedtools intersect between two datasets,
whether it is across species or across tissue. 

Depending on the flag used to run it, it will find either overlapping regions or regions
that are specific to one species. (-wa or -v respectively).

Input: two genomic .bed datasets (ATAC-Seq data for example)
Output: One dataset representing overlapping or exclusive regions.
