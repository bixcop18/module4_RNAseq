#!/bin/bash -e


samples="Sample_118A Sample_119B Sample_120A DF5 DF19 Sample_103A Sample_104B Sample_105B spike_Z32_rep1 spike_Z32_rep2"

sample_folder="/Users/ramirezr/Dropbox/JIC/Workshops/2018-05 Beca/reads"

mkdir -p kallisto_out

for s in $samples; do
	mkdir -p "kallisto_out/$s"
	kallisto quant -i ../Trinity-mini-test-CDHit90.fasta.k31  -o kallisto_out/$s "$sample_folder/$s/$s.r1.fastq.gz" "$sample_folder/$s/$s.r2.fastq.gz"
done