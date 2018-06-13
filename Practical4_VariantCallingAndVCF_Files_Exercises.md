### Practical 4: An Exercise in Variant Calling and Interpreting VCF Files


**Objectives**



* Prepare bam files for calling SNPs.

	In a previous practical, you used Kallisto to align the RNA-Seq reads from spike material of two wheat cultivars, Chinese Spring and Azhurnaya and obtained a bam file. You can use these files.

* Learn how to call SNPs with Freebayes

* Become familiar with VCF format files

* Learn some criteria for filtering  VCF files

<br>

**Step 1 - Remove duplicate reads from the bamfiles with Samtools**

* [Samtools manual](http://www.htslib.org/doc/samtools.html)

* [SAM format rules](https://samtools.github.io/hts-specs/SAMv1.pdf)


**Inputs**

* Location: /var/scratch/baileyp/Practical\_4\_small\_wheat\_dataset/sorted\_bam\_from\_HISAT2\_mapping 

* bam files for Chinese Spring:&nbsp; DF5 &nbsp;and &nbsp;DF19&nbsp; (two replicate samples)

* bam files for Azhurnaya: &nbsp;119B &nbsp;and &nbsp;120A &nbsp;(two replicate samples

* Create a working directory for the practical

<br>

**Package to load**

```sh
module load samtools/1.8
```

* It is important to remove duplicate reads before calling SNPs if PCR was used in the preparation of your sequencing library. We will use Samtools to do this step. It is a four step process:

* Take the DF5 bam file and sort it  in read name order (this step can be omitted if the file is already name ordered, but it is not in this case!):

```sh
samtools sort -n -o CS_DF5_st_namesort.bam /var/scratch/baileyp/Practical_4_small_wheat_dataset/sorted_bam_from_HISAT2_mapping/sorted_DF5_P.bam
```
* Add MC and ms tags to the bam records (used by the markdup command):

```sh
samtools fixmate -m CS_DF5_st_namesort.bam CS_DF5_st_fixmate.bam
```

* Now sort the file again by position (field 4 of the sam record)

```sh
samtools sort -o CS_DF5_st_positionsort.bam CS_DF5_st_fixmate.bam
```

* Finally, we are ready to mark and remove duplicate reads:

```sh
samtools markdup -s -r CS_DF5_st_positionsort.bam  CS_DF5_st_rmdup.bam
```

* **Question:** From the output printed to the screen or your log files, what is the percent of duplicate reads?

* View the output with samtools view:

```sh
samtools view CS_DF5_st_rmdup.bam | head
```
- can you see ther information that has been added to the end of bam records (ms and MC fields)

* These steps need to be repeated for the other three samples. To make this slightly easier to do, we can set two variables that contain the sample name, then run the first command again using those variables that are set up to process the next sample, DF19:

```sh
sample=DF19				# Other sample names: 119B, 120A
samplePrefix=CS_DF19	# Other sample prefixes: Az_119B, Az_120A
samtools sort -n -o ${samplePrefix}_st_namesort.bam  /var/scratch/baileyp/Practical_4_small_wheat_dataset/sorted_bam_from_HISAT2_mapping/sorted_${samplePrefix}_P.bam
```

* After running the above command, reset the variables to process the other two samples with the same command. Use the up arrow key to retrieve the command.

* **Note:** If you get confused using the variables, you can perform each step manually for each sample, as described that the top of this step for the first sample, CS\_DF5.
*  **Note:** You can place an ampersand character ('&') at the end of any command to make it run in the background. Then you can carry on using the command line for something else. When the program has finished you will get a message printed to the screen - e.g.:

		[1]+  Done   name_of_command  


* Now perform the reamining steps of de-duplication process. An example is given below for the DF19 sample and running the fixmate step. Remember, you need to set the variable for each sample and retrieve the command you want to run using the up arrow key. Suggested sample prefix names: CS\_DF19, Az\_119B, Az\_120A 

```sh
samplePrefix=CS_119B
samtools fixmate -m ${samplePrefix}_st_namesort.bam ${samplePrefix}_st_fixmate.bam
```
<br>
<br>

**Step 2 - Merge the bam files and add read group information**

* Merge the bamfiles for all samples. The -r flag attaches an identifier (read group) to the bam record, based on the file name:

```sh
samtools merge -r merged.bam CS_DF5_st_rmdup.bam CS_DF19_st_rmdup.bam Az_119B_st_rmdup.bam Az_120A_st_rmdup.bam
```

* Next, you need to add a short recognisable name for each sample into the read group information. To do this you first need to extract the header from the bam file:

```sh
samtools view -H merged.bam > merged_header.txt 
```

* Now add the read group info at the bottom of the header file with your favourite editor or use vi and type:

```sh
vi merged_header.txt
```

* Go to the last line in the file.

* Type a 'G'.

* Press 'i' to insert text.

* Copy and paste these lines at the bottom of the file (the spaces need to be tabs!). Check that the "ID" tag is the same as the bam file name (minus the .bam ending). The SM (sample) tag is used as the sample name by Freebayes:

```sh 
@RG		ID:CS_DF5_st_rmdup		SM:CS_DF5	LB:LIB1	PL:Illumina
@RG		ID:CS_DF19_st_rmdup		SM:CS_DF19	LB:LIB2	PL:Illumina
@RG		ID:Az_119B_st_rmdup		SM:Az_119B	LB:LIB3	PL:Illumina
@RG		ID:Az_120A_st_rmdup		SM:Az_120A	LB:LIB4	PL:Illumina
```

* Then press 'shift :', then type wq to save and quit the vi editor.



* **Tip:** Other useful commands in the vi editor escape mode:
	dd = delete current line; gg goes to the top of the file again

* Next reheader the bam file:

```sh
samtools reheader merged_header.txt merged.bam > merged_readgroups.bam
```

* **Question:** How else could read group info be added to the samples?


* Check that the header was applied properly:

```sh
samtools view -H merged_readgroups.bam
```

* Now index the bam file. The index can be used for quick retrieval of records by downstream applications:

```sh
samtools  index  merged_readgroups.bam
```

<br>

**Step 3: Call SNPs with FreeBayes**

* [FreeBayes info and manual](https://github.com/ekg/freebayes)

* [Paper reference](https://arxiv.org/abs/1207.3907) - it's very technical 🙈!

**Inputs**

* merged_readgroups.bam prepared above in step 2.

* Genome reference:
	/var/scratch/baileyp/Practical_1_RicardosExample/References/transformed_coordinates.fasta
	
<br>

**Package to load**

```sh
module load freebayes/1.0.2
```

* We can now run FreeBayes on the alignments of all 4 samples simultaneously and generate a VCF file:


```sh
freebayes \
--fasta-reference /var/scratch/baileyp/Practical_1_RicardosExample/References/transformed_coordinates.fasta \
merged_readgroups.bam \
> all_snps.vcf &
```

* Note the ' \ ' characters. This allows the command to straddle across multiple lines. This makes commands clearer to visualise when many flags are being used.

* Freebayes will take about 30 - 35 minutes to run. To speed up the analysis, it is possible to call SNPs on a particular part of the reference. The program will then run more quickly. To demonstrate this just run the program for wheat chromosome 1A using the -r flag. In the reference fasta file, the name of the chromosome is "chr1A":     


```sh
freebayes \
-r chr1A \
--fasta-reference /var/scratch/baileyp/Practical_1_RicardosExample/References/transformed_coordinates.fasta \
merged_readgroups.bam \
> chr1A_snps.vcf &
```

* If you had a large number of samples you could analyse each chromosome arm separately but run them in parallel at the same time. This could save a lot of time.

* **Question:** How would you specify to only call SNPs on particular parts of a chromosome and also multiple parts of multiple chromosomes? Take a Look in the command line help and see if you can spot the flags to use:

```sh
freebayes -h | more
```

* **Question:** How would you find out the name/identifier of the chromsomes to specify in these flags?

* There are some filtering parameters avaialable in FreeBayes but it might be better to filter the SNPs after collecting all the raw calls - examples:
	* -C --min-alternate-count (default = 2) &nbsp;&nbsp;&nbsp;&nbsp; 6 might be a better minimum to tolerate 
	* -i --no-indel
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Ignore insertion and deletion alleles
<!-- * -p --ploidy (default = 2, diploid)
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FreeBayes can consider polyploids when calling SNPs. -->

<br>

**Step 4: Understand the VCF format and filter the VCF file**

* The VCF format looks very complicated! Here are three links that explain the format. (You may find one easier to understand than the others): 

	[https://en.wikipedia.org/wiki/Variant\_Call\_Format](https://en.wikipedia.org/wiki/Variant_Call_Format)
	
	[What is a VCF and how should I interpret it?](https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it)

	[Another explanation of VCF format](http://www.internationalgenome.org/wiki/Analysis/vcf4.0)

	[The formal specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
	
* The aim here is to introduce you to the essential syntax and fields. Let's do some filtering with VCFtools:

* [VCFtools usage examples](https://vcftools.github.io/documentation.html)

* [VCFtools manual](https://vcftools.github.io/man_latest.html)

* [VCFtools utilities](https://vcftools.github.io/perl_module.html)

* **Package to load**

```sh
module load vcftools/0.1.15
```

* We will filter the SNPs on the following basic criteria (vcftools flags in brackets):
	* minimum quality (--minQ)	20
	* depth for an individual sample (--minDP)	6
	* remove indels (--remove-indels)
	* allow no missing genotypes (max-missing 1)

* Add these flags to the command below, one after the other and keep a record of the number of SNPs you find after adding each additional filter. Also alter the output filenames accordingly.


```sh
vcftools --vcf all_snps.vcf \
[ add your filters here ]

```

* Try also adding these parameters to obtain the calculations for the transition (Ts)/tranversion ratio (Tv):
--TsTv-summary
--TsTv-by-count

* You will probably get an error but it's very simple to correct.

* Is the ratio Ts/Tv as you would expect or too high, too low?

* Finally, print out the filtered SNPs to a file. You need to add the recode flag(s) before a filtered vcf file will be produced.

```sh
vcftools --vcf merged_readgroups.vcf \
--recode --recode-INFO-all \
--minQ 20 \
--minDP 6 \
--remove-indels \
--max-missing 1 \
--out all_minQ20_minDP6_noindels_maxm1
```

* Try adding up the values you find in the INFO and FORMAT fields and see whether they add up to a field that represents the total for the number of alleles or depth. Examples to try:

<br>

**List of tools you could also use to carry some of the above tasks**

* For marking duplicates: [Picard](https://broadinstitute.github.io/picard/) MarkDuplicates program

* Other SNP callers:
	* GATK - another haplotype caller like FreeBayes (but FreeBayes is easier to use)
	* [Index of GATK programs](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/index) (see under Variant Discovery tools)
	* Using GATK is a 2 step process. For Calling germline SNPs you would use HaplotypeCaller followed by GenotypeGVCFs
	* [Samtools](http://samtools.sourceforge.net/mpileup.shtml)

* Other VCF filtering tools:
	* bcftools
	* [SnpSift](http://snpeff.sourceforge.net/SnpSift.html)
	* [VEP](https://www.ensembl.org/info/docs/tools/vep/script/vep_filter.html)


<br>
<br>
<br>
<br>
**Answers to the questions**




Question: How else could read group info be added to the samples?
Read groups can be added when setting up many read aligners e.g. BWA mem - go to the manual and look at the mem command and the -R option. ([BWA manual](http://bio-bwa.sourceforge.net/bwa.shtml))












