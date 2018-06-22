#!/bin/bash

for f in  1529553255.fas.0 1529553255.fas.1 1529553255.fas.2 1529553255.fas.3; do
	contigs=`grep -c ">" $f`
	echo "$f: $contigs"
	gmap --min-intronlength=20 --format=gff3_gene --npaths=1 --ordered --min-identity=0.95 -d transformed_coordinates -D ../../../References/gmap $f > $f.gff
	gt gff3 -tidy yes -retainids yes  -sort yes $f.gff > $f.sorted.gff
	gt seqtranslate $f > $f.pep.fas
	gffread -x $f.gff.sorted.nuc.fasta -y $f.gff.sorted.pep.fasta -F -g ../../../References/transformed_coordinates.fasta $f.sorted.gff
	gt sketch -format pdf -seqid chr1B -start 320232 -end 345985 $f.pdf $f.sorted.gff
done
