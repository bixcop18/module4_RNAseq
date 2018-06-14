### Practical 5: Variant Calling with RNA-Seq from Arabidopsis Accessions


**Background to the original study**

* [The original paper](https://www.nature.com/articles/nature10414)

* [Supplementary data from the original paper](https://media.nature.com/original/nature-assets/nature/journal/v477/n7365/extref/nature10414-s1.pdf)

* [Website with the data from the study](http://mtweb.cs.ucl.ac.uk/mus/www/19genomes/index.html)

* This study formed part of the [Arabidopsis 1000 genomes Project](http://1001genomes.org/index.html)


**Objectives**

* Align RNA-Seq reads from seedling material of two Arabidopsis ecotypes, Columbia and Landsberg erecta with the STAR splice-aware aligner

* Call SNPs with Freebayes and filter the SNPs

* Annotate SNP calls with snpEff to indicate whether the SNP changes are functionally significant e.g. produce a stop codon

<!-- Prepare a working directory for this project -->
<br>

**Inputs**

* Location of fastq reads: /var/scratch/baileyp/Practical\_5\_Arabidopsis\_accns/fastq_files

* fastq read file for Columbia (Col_0): R38\_L1\_Col\_0\_sdlg\_R1\_nss\_barcode\_GTA.fastq
(78 base pair reads)

* fastq read file for Landberg erecta (Ler_0): R38\_L3\_Ler\_0\_sdlg\_R1\_nss\_barcode\_TGC.fastq
 (78 base pair reads)

* Location of the Arabidopsis genome reference and GTF file:
	/var/scratch/baileyp/Practical\_5\_Arabidopsis_accns/genome\_reference

* Reference: Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

* GTF file: Arabidopsis_thaliana.TAIR10.39.gtf

<br>

**Packages to load**

```sh
module load star/2.5.2b
module load trimmomatic/0.38
module load samtools/1.8
```
<br>

**Step 1 - Index the Arabidopsis reference with STAR**

* [The STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

* First, copy over the Arabidopsis reference to your working directory:

```sh
cp /var/scratch/baileyp/Practical\_2\_Arabidopsis_accns/genome\_reference/ .
```

* Now index the reference. Before running, adjust the code so that STAR uses four cpus instead of one. Also adjust the sjdbOverhang parameter according to the manual instructions. 100 will not  produce an index with this data set!

```sh
STAR --runThreadN 1 \
--runMode  genomeGenerate \
--genomeDir  ./ \
--genomeFastaFiles  Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
--sjdbGTFfile  Arabidopsis_thaliana.TAIR10.39.gtf \
--sjdbOverhang 100
```

<br>

**Step 2 - Prepare the reads for aligning to the  reference**

* **Note:** The FreeBayes SNP caller requires that the Phred sequence quality scores are Phred+33 but in this data set these reads were processed with the Illumnina v1.6 pipeline which outputs reads with Phred+64 scores.

* For more information on phred scores, see [the fastq format Wikipedia page](https://en.wikipedia.org/wiki/FASTQ_format)

* Therefore we need to convert these fastq records from Phred+64 to Phred+33.

* Amazinghly, the Wiki page tells you how to convert the scores, like this:
 
```sh
sed -e '4~4y/@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghi/!"#$%&'\''()*+,-.\/0123456789:;<=>?@ABCDEFGHIJ/' \
/var/scratch/baileyp/Practical_2_Arabidopsis_accns/fastq_files/R38_L1_Col_0_sdlg_R1_nss_barcode_GTA.fastq \
> Col_0_R1_Phred33.fastq
```

* Use the above command to convert the phred scores in both the Col\_0 and Ler\_0 fastq files and give the output files a short sensible name as shown.

* Now we are ready to trim the read ends to remove bases with low quality. We will use the Trimmomatic program:

```sh
At_accn=Col_0 	# Col_0 Ler_0
cpu=4
trimmomatic SE \
-threads $cpu \
-phred33 \
-trimlog ${At_accn}_R1_trimmomatic.log \
Col_0_R1_Phred33.fastq\
${At_accn}_R1_trimmomatic.fq.gz \
LEADING:1 \
TRAILING:1 \
SLIDINGWINDOW:5:18 \
MINLEN:36 \
MINLEN:36 \
> ${At_accn}_R1_trimmomatic.log 2>&1 &
```

* The output file is a compressed file with name ending ".fq.gz" which can be used as input to the STAR aligner.

<br>

**Step 3 - Align prepared reads to the reference with STAR**

* Align separately the trimmed Col\_0, then Ler\_0 reads to the Arabidopsis genome reference: 

```sh
At_accn=Col_0      # Col_0 Ler_0
STAR --runThreadN 4 \
--runMode alignReads \
--genomeDir 	[your_genome_reference_directory] \
--readFilesCommand zcat \
--readFilesIn {$At_accn}_R1_trimmomatic.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outSAMattrRGline ID:$At_accn SM:$At_accn \
--outFileNamePrefix ${At_accn}
```

* **Note:** We need to insert the read group information into the bam file for FreeBayes. It can be added as part of the STAR command using the outSAMattrRGline flag. It avoids having to add it manually afterwards, as we did for Practical 1.

<br>

**Step 4 - Remove duplicate reads from the bamfiles**

* It is important to remove duplicate reads before calling SNPs. We will use Samtools to do this step. It is a four step process and uses exactly the same commands used in the same steps in Practical 1.

* Run each of the commands after exchanging the asterix for the input or output file name for each sample:

```sh
samtools sort -n -o *_st_namesort.bam *.bam
```
* Add ms and MC tags to the bam records (used by the markdup command):

```sh
samtools fixmate -m *_st_namesort.bam *_st_fixmate.bam
```

* Now sort the file again by position (field 4 of the sam record) 

```sh
samtools sort -o *_st_positionsort.bam *_st_fixmate.bam
```

* Finally, we are ready to mark and remove duplicate reads:

```sh
samtools markdup -s -r *_st_positionsort.bam  *_st_rmdup.bam
```

* **Question:** From the output printed to the screen, what is the percent of duplicate reads?

<br>

**Step 5 - Merge the bam files**

* Merge the bam files for both samples. For this practical, we do not need to use the -r flag to attach the read group information because it has already been added:

```sh
samtools merge -f merged.bam Col_0_st_rmdup.bam Ler_0_st_rmdup.bam
```


* Qickly check that the read groups are present and the SM flag contains the abbreviations for the Arabidopsis accessions:

```sh
samtools view -H Col_0_Ler_0_merged.bam
```

* Index the bam file:

```sh
samtools  index  Col_0_Ler_0_merged.bam
```

<br>

**Step 6 - Call SNPs with FreeBayes**

* Now we are ready to call SNPs. Calling SNPs on the whole genome takes quite a long time but we can restrict the analysis a particular chromosome arm with the -r flag. We are interested in chromosome 4:

```sh
module load freebayes/1.0.2
freebayes \
--fasta-reference Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
Col_0_Ler_0_merged.bam \
-r 4 \
> merged_chr4_ONLY.vcf &
```
<!--At1G64970 - the SNP is in an intron!!!!-->

<br>

**Step 7 - Filter the chromosome 4 SNP calls**

* We will filter exactly as for the previous practical:

```sh
module load vcftools/0.1.15
vcftools --vcf merged_chr4_ONLY.vcf \
--recode --recode-INFO-all \
--minQ 20 \
--minDP 6 \
--remove-indels \
--max-missing 1 \
--out chr4_ONLY_minQ20_minDP6_noindels_maxm1
```

<br>

**Step 8 - Annotate snpEFF Annotation**

* [The snpEff Manual](http://snpeff.sourceforge.net/SnpEff_manual.html)

* There are programs to annotate SNPs according to whether the SNP falls within a gene to produce a deleterious mutation e.g. stop codon or potentially a non-synomymous substitution. We will use the snpEff tool:

```sh
module load snpeff/4.1g
snpEff \
-c /var/scratch/baileyp/Practical_5_Arabidopsis_accns/snpEff/snpEff.config \
athalianaTair10 \
chr4_ONLY_minQ20_minDP6_noindels_maxm1.recode.vcf \
> chr4_ONLY_filtered_snpeff.vcf
```

* The output comes with a summary html file containing information on the number of genes with SNP effects or consequences.

* A complete list of consequences can be found at Ensembl [here] (https://www.ensembl.org/info/genome/variation/predicted_data.html#consequences).

* How would you quickly count the number of stop codons and non-synomymous substitutions effects in the new snpeff.vcf file. How many of each are there?
It is better to use a filtering tool however e.g. snpSift (see via the snpEff manual above)

* Another SNP effect annotation tool to look at is the  
the Ensembl [Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html) which also has its own [filtering tool](https://www.ensembl.org/info/docs/tools/vep/script/vep_filter.html).

<!-- **Answers to the questions** -->








