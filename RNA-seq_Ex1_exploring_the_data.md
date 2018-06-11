# RNA-seq: exercise 1
## Exploring the data:

## Objectives:

* Explore FASTA and FASTQ files
* Run FASTQC on all 4 samples
* Trim reads from all 4 samples
* Prepare reference indices to run the alignments

## Organising a project

NGS analysis can lead to have many temporary files. For that reason, it is recommended to have subfolders with each step on the analysis. In that way, it is easy to remove and repeat failed attempts and it is easier to find the results and the steps you followed in the future. 

### Inputs

* **References/** You can symlink to the actual reference, to avoid coping big files in the cluster
* **Reads/** Depending on how you get the reads, it may be wise to have sub-folders per sample, or give the reads an standard name. For this tutorial, the name of the library is used

To explore the content of the files you can use ```less``` and ```gunzip -c``` if the fastq file is compressed 

We will be carrying out two sets of mapping. One using the transcriptome reference and one using the genome reference.

Have a peek at the genome fasta reference file

```sh
less References/transformed_coordinates.fasta
```

* The arrows to go up and down the file
* ```Ng``` can be used to go to a particular line (```700g``` goes to line 700)
* ```g``` and ```G``` go to the first or last line 
* ```-N``` displays the line number
* ```/``` can be used to search a pattern and ```n``` goes to the next result
* ```&``` displays only the lines with the pattern pattern (```&>``` displays the fasta headers)
* ```h``` shows all the commands you can use in less


If you want to know how many scaffolds the file has, you can use ```grep```.

```sh
grep -c ">" References/transformed_coordinates.fasta
```

* It is very important to use the speech marks as you will overwrite your input file if you omit them. 

You can also look at the transcriptome reference we will be using.
```sh
less References/transformed_coordinates.fasta
```
Similarly, if you want to know how many chromosomes the file has, you can again use ```grep```.

```sh
grep -c ">" References/selected_refseq1.0.fasta
```

To inspect the ```fastq.gz``` files you can use less. However, if you want to do any other analysis with a text processing tool, you need to uncompress it and pass it to your command. For example, to count the number of lines with ```wc -l```, you do the following:

```sh
gunzip -c reads/104B/Sample_104B.r1.fastq.gz | wc -l
``` 

#### Questions:
* Which are the contigs in the reference file?
* How many reads does each fastq have? (Hint: Each read is four lines in Illumina ```fastq``` files)

### Outputs
For this tutorial, create folders for the outputs of each program. 

```sh
mkdir -p fastqc
mkdir -p trim
mkdir -p gtf_splice
mkdir -p hisat_index
mkdir -p kall_index
```
These are the tools we will use in the tutorial:

```sh
fastqc/0.11.5
trimmomatic/0.38
cufflinks/2.2.1
python/3.6.2
hisat2/2.0.5
kallisto/0.43.0
``` 

* Keep the record of which version of the software you use. Sometimes updates may change the output. To be consistent in the analysis, **use always the same version** within a experiment. 
* You can write scripts with all your commands to run the same analysis again.


## Verifying the quality of the reads

To review the quality of the fastq, fastqc can be used to plot and summarise the 

```sh
fastqc/0.11.5
fastqc -o fastqc/ -f fastq reads/Sample_104B.r1.fastq.gz
```
After fastqc is done, open the report in ```fastqc/Sample_104B.r1_fastqc.html```

Repeat this for each of the remaining samples (r1 and r2 for each):
```sh
/reads/Sample_104B.r2.fastq.gz
/reads/Sample_105B.r1.fastq.gz
/reads/Sample_105B.r2.fastq.gz
/reads/Sample_119B.r1.fastq.gz
/reads/Sample_119B.r2.fastq.gz
/reads/Sample_120A.r1.fastq.gz
/reads/Sample_120A.r2.fastq.gz
```
#### Question
* Is the quality of the reads good?

## Trimming the reads
We will use Trimmomatic to trim the PE reads for quality. For example:

```sh
module load trimmomatic/0.38

trimmomatic PE -threads 2 \
reads/Sample_104B/Sample_104B.r1.fastq.gz \
reads/Sample_104B/Sample_104B.r2.fastq.gz \
-baseout /trim/104B.fastq \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:2 MINLEN:60 AVGQUAL:20
```
Run Trimmomatic on all the samples as PE reads.

These are fairly relaxed trimming parameters. Look up the options in the documenation to find out what they mean: 
http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

For these samples, the majority of reads will survive trimming as PE. It is these that we will use for mapping.

```sh
samplename_1P.fastq and samplename_2P.fastq
```
Any reads that only survived as SE will be named: 

```sh
samplename_1U.fastq and samplename_2U.fastq
```
You can examine the results of the trimming, and find out how many reads passed the quality filters by reading the .err file.

## Indexing the references
Reference files tend to be big, so sequence aligners require an index to be able to run. We will make indices for both our HISAT2 mapping and our kallisto mapping.

### HISAT2 index
```sh

module load hisat2/2.0.5

hisat2-build -p 2 \
References/transformed_coordinates.fasta \
hisat_index/transformed_coordinates_fasta_indx
```
You can provide HISAT2 with a list of known splice sites, which it then uses when aligning reads with small anchors. We will use this file tomorrow.

You can create such a list using the python script:

```sh
hisat2_extract_splice_sites.py
```
Where ```hisat2_extract_splice_sites.py``` is included in the HISAT2 package.

First we need to convert our gff file to a gtf, using a utility found in the Cufflinks package.

```sh
module load cufflinks/2.2.1

gffread References/transformed_coordinates.gff -T -o \
gtf_splice/transformed_coordinates.gtf

```
Then we can run the python script to extract the splice sites:

```sh
module load python/3.6.2

extract_splice_sites.py \
gtf_splice/transformed_coordinates.gtf \
> gtf_splice/transformed_coordinates.gff.txt
```

### kallisto index

```sh
module load kallisto/0.43.0

kallisto index -i kall_index/kall_selected_refseq1 \
References/selected_refseq1.0.fasta
```

#### Questions:

* How many new files do you have? Which ones correspond to which aligners?

## Useful links

* **FastQ format**: https://en.wikipedia.org/wiki/FASTQ_format
* **Trimmomatic manual**:
    http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
* **HISAT2 manual**:
    https://ccb.jhu.edu/software/hisat2/manual.shtml
* **kallisto manual**:
    https://pachterlab.github.io/kallisto/manual
