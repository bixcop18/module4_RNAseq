# RNA-seq: exercise 3
## Calling differential expression using Ballgown and sleuth

## Objectives:
* Import the raw count data into the R packages: Ballgown and sleuth
* Filter and carry out between sample normalization
* Call DE between spike and root samples
* Visualise the data

* As before, Keep the record of which version of the software you use. Sometimes updates may change the output. To be consistent in the analysis, **use always the same version** within an experiment. 
* You can write scripts with all your commands to run the same analysis again.

For this exercise we will be working in R studio, on your local machine. First copy the expression data we will be using from the module4_RNA_seq folder on github, into your R working directory.
You will need the following two folders:

* For Ballgown:

```sh
alignments_HISAT2/
```

This folder contains another folder ```ballgown_str_mrg``` that contains 4 folders containing the expression data quantified by StringTie. 


* For sleuth:

```sh
/alignments_kallisto
```

Which contains the expression data for the same four samples, quantified by kallisto.

## Ballgown
## Load the relevant packages in R

These include the Ballgown package that you will use for performing most of the analyses, as well as a few other packages that help with these analyses, 
specifically RSkittleBrewer (for setting up colours), genefilter (for fast calculation of means and variances), dplyr (for sorting and arranging results)
and devtools (for reproducibility and installing packages):

```sh
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
library(ggplot2)
library(cowplot)
```
For any that are missing, use:

```sh
install.packages('package')
```
Next we need to load the phenotype data for the samples. But first you will have to create this file yourself as a ```.csv```. Do this using the text editor of your choice, and save it in the ```alignments_HISAT2/``` directory as ```spike_root_pheno_data```. It contains information about the RNA-seq samples. Each sample should be described on one row of the file, and each column should contain one variable. The contents of the file should look like below. 

```sh
sample,tissue
104B,root
105B,root
119B,spike
120A,spike
```

Next, to read this file into R, we use the command ```read.csv```. In this file, the values are separated by commas, but if the file were tab-delimited you could use the function ```read.table```.

```sh
pheno_data = read.csv("spike_root_pheno_data.csv")
```

Read in the expression data that was calculated by StringTie. 
To do this, we use the ballgown command with the following three parameters: 
	the directory in which the data are stored ```dataDir```
	a pattern that appears in the sample names ```samplePattern```
	phenotypic information that we loaded in the previous step ```pData```.

```sh
bg = ballgown(dataDir = "ballgown_strg_mrg/", samplePattern = "1", pData=pheno_data)
```

Next we filter to remove low-abundance genes. One common issue with RNA-seq data is that genes often have very few or zero counts. A common step is to filter out some of these. Another approach that has been used for gene expression analysis is to apply a variance filter. Here we remove all transcripts with a variance across samples less than one:

```sh
bg_filt = subset(bg,"rowVars(texpr(bg)) >1",genomesubset=TRUE)
```

We can find out how many transcripts are in the Ballgown object by typing:

```sh
bg
bg_filt
```

We can make a table containing the information from the filtered dataset, such as the transcript id ```t_id```, transcript name ```t_name```, count information, genomic coordinates and gene id. This table will be useful when identifying specific genes and transcripts later in the analysis.

```sh
bg_filt_table = texpr(bg_filt, 'all')
```

Use the ```head()``` function to look at the contents of this table.

Next, using the function ```stattest``` on the filtered results, we can identify transcripts that show statistically significant differences between root and spike.

```sh
results_transcripts <- stattest(bg_filt, feature="transcript", covariate="tissue", getFC=T, meas = "FPKM")
```

Then we can identify genes that show statistically significant differences between groups. 

```sh
results_genes <- stattest(bg_filt, feature="gene", covariate="tissue", getFC=T, meas = "FPKM")
```
The parameters we have specified here, in both cases are:

* ```feature``` the type of genomic feature to be tested for DE
* ```covariate``` this is the covariate of interest and must correpond to the nsame of a column in ```pheno_data```
* ```getFC``` returns the estimated fold change
* ```meas``` specifies the expression measurement to be used for statistical tests


Next, add in the gene names and ids to ```results_transcripts```

```sh
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filt),
                                 geneIDs=ballgown::geneIDs(bg_filt), results_transcripts)
```

You will notice, if you look at ```results_transcripts```, that only some of the transcripts have been associated with gene names. One of the 'features' of StringTie ```--merge``` is that it renames most of the transcripts and genes and some of this information is left behind in the merged ```gtf``` file. It is possible to retrieve this information from the ```gtf``` with a custom script and merge it with the ```results_transcripts``` and ```results_genes``` tables in R.


Sort from the results from smallest pval to largest:

```sh
results_transcripts <- arrange(results_transcripts, pval)
results_genes <- arrange(results_genes, pval)
```

We can count the number of transcripts with a p-value <= 0.05.

```sh
table(results_transcripts$pval <= 0.05)
```
Do this for the ```results_genes``` table.

Write the results to a ```.csv``` file

```sh
write.csv(results_transcripts, "spike_root_ballgown_results_trans.csv", row.names=F, quote=F)
write.csv(results_genes, "spike_root_ballgown_results_genes.csv", row.names=F, quote=F)
```

Identify transcripts and genes with a p-value <=0.05:

```sh
significant_results_transcripts <- subset(results_transcripts,results_transcripts$pval<=0.05)
significant_results_genes <- subset(results_genes,results_genes$pval<=0.05)
```
#### Questions:
* How many 'low abundance' genes with variance in expression less than 1 FPKM, were removed?
* How many transcripts were identified as showing statistically significant differences between spike and root?
* And how many of these transcripts had a p-value <=0.05?
* How many genes were identified as showing statistically significant differences between spike and root?
* And how many of these genes had a p-value <=0.05?

You can use Ballgown to visualize RNA-seq results in a variety of ways. The plots produced by Ballgown make the results easier to view and compare expression data. First, we specify the colour palette.

```sh
tropical=c('darkorange', 'dodgerblue','hotpink', 'limegreen', 'yellow')
palette(tropical)
```

Next we can plot the distribution of gene abundances (measured as FPKM values) across samples, colored by tissue type. 

Ballgown stores a variety of measurements that can be compared and visualized. We will compare the FPKM measurements for all the transcripts, but you can also obtain measurements for each gene in the data set. 
This command below extracts the expression measurements from the ballgown object ```bg_filt```.

```sh
fpkm = texpr(bg_filt, meas="FPKM")
```

The plots will be easier to visualize if you first transform the FPKM data. We use a log2 transformation that adds one to all FPKM values because log2(0) is undefined.
 

```sh
fpkm = log2(fpkm+1)
```

Then we create the plot: 

```sh
boxplot(fpkm, col=as.numeric(pheno_data$tissue), las=2,ylab='log2(FPKM+1)')
```
If you see the following error message in the Console window:

```sh
Error in plot.new() : figure margins too large
```
Try increasing the size of the plot window in R studio and run the command again.


We can also make plots of individual transcripts across samples. For example, here we show how to create a plot for the 1948th transcript in the dataset. The first command below shows the name of the transcript:

```sh
ballgown::transcriptNames(bg_filt)[1948]

plot(fpkm[1948,] ~ pheno_data$tissue, border=c(1,2), main=ballgown::transcriptNames(bg)[1948], 
	pch=19, xlab="tissue", ylab='log2(FPKM+1)')

points(fpkm[1948,] ~ jitter(as.numeric(pheno_data$tissue)), col=as.numeric(pheno_data$tissue))

```

We can also plot the average expression levels for all transcripts of a gene within different groups using the plotMeans function. We need to specify which gene to plot, which Ballgown object to use and which variable to group by.

As an example, we will plot the transcripts associated with the example we looked at above ```mrg.1188``` by using the following command. The plot shows the average expression, between samples, coloured according to expression level.

```sh
plotMeans('mrg.1188', bg_filt, groupvar='tissue', legend=T)
```

We can also create an MA plot, with DE expressed genes (p<=0.05) highlighted in red.

```sh
results_transcripts$mean <- rowMeans(texpr(bg_filt))

ggplot(results_transcripts, aes(log2(mean), log2(fc), colour = pval<=0.05)) +
  scale_colour_manual(values=c("#999999", "#FF0000")) +
  geom_point() +
  geom_hline(yintercept=0)
```

Ballgown has many other capabilities that are worth exploring, though it may take a little time and forum reading to follow some of them!

## sleuth
Key features of sleuth include:

The ability to perform both transcript-level and gene-level analysis. Compatibility with kallisto enabling a fast and accurate workflow from reads to results.
The use of boostraps to ascertain and correct for technical variation in experiments.

To use sleuth, RNA-Seq data must first be quantified with kallisto, which is a program for very fast RNA-seq quantification based on pseudo-alignment. An important feature of kallisto is that it outputs bootstraps along with the estimates of transcript abundances.

These can serve as proxies for technical replicates, allowing for an ascertainment of the variability in estimates due to the random processes underlying RNA-seq as well as the statistical procedure of read assignment.

sleuth has been designed to facilitate the exploration of RNA-Seq data by utilizing the Shiny web application framework by RStudio.


## Load the relevant packages in R

```sh
library("sleuth")
library(dplyr)
library(MASS)
library(shiny)
```

The first step in a sleuth analysis is to specify where the kallisto results are stored. 
A variable is created for this purpose.

```sh
sample_id <- dir(file.path(".", "alignments_kallisto/"))
```
Check you have entered the correct file path to your kallisto files by typing:

```sh
sample_id
```

A list of paths to the kallisto results indexed by the sample IDs is collated with:

```sh
kal_dirs <- file.path(".", "alignments_kallisto", sample_id)
```
Check you have specified the right file paths:

```sh
kal_dirs
```

The next step is to create and load an auxillary table that describes the experimental design and the relationship between the kallisto directories and the samples.
Make this table and save it as ```experiment_data``` to the ```alignments_kallisto/``` directory.

```sh
sample condition
104B	root
105B	root
119B	spike
120A	spike

```
Now the directories must be appended in a new column to the table describing the experiment. This column must be labelled path, otherwise sleuth will report an error. This is to ensure that samples can be associated with kallisto quantifications.

```sh
s2c <- read.table(file.path(".", "experiment_data"), header=TRUE, stringsAsFactors = FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
```
It is important to check that the pairings are correct:

```sh
s2c
```

Next, the “sleuth object (so)” can be constructed. This object will store not only the information about the experiment, but also details of the model to be used for differential testing, and the results.
It is prepared and used with a set of commands that load the kallisto processed data into the object, estimate parameters for the sleuth response error measurement (full) model, estimate parameters for the sleuth reduced model, and perform differential analysis (testing) using the likelihood ratio test or the wald test. On a laptop the four steps should take about a few minutes altogether.

The sleuth object must first be initialized with

```sh
so <- sleuth_prep(s2c, ~condition)
```

Then the models are fitted with:

```sh
so <-sleuth_fit(so)
so <-sleuth_fit(so, ~1, 'reduced')
so <-sleuth_lrt(so, 'reduced', 'full')
```

We can collate the results onto a table:

```sh
results_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
```

If we want to create a table of significantly differentially expressed transcripts:

```sh
results_table_significant <- dplyr::filter(results_table, pval <= 0.05)
```

#### Questions:
* How many transcripts are significantly differentially expressed?

We can use the ```shiny``` app to visualize the results.

```sh
sleuth_live(so)
```
This opens a webpage in your browser. In your console window in R studio you will see:
```sh
Listening on http://xxxxxxxxxx
```
When you want to return to working in R studio, you will need to quit this command to return to the ```>``` prompt, before you can continue.

Navigate around the webpage and explore the different tabs.
From the ```maps``` tab. You can visualize a ```sample heatmap``` and ```PCA``` of the samples.
Look at the ```test table``` under the ```analyses``` tab.
Sort the table by qval in descending order to identify the most significantly differentially expressed gene.

Use ```transcript view``` under the ```analyses``` tab to visualize how expression of these transcripts changes across samples.

#### Questions:
* How would you go about finding out more about the most significantly DE transcripts?

If we want to visualise an MA plot, we need to run the 'Wald test'.

```sh
so <- sleuth_wt(so, 'conditionspike')
```

This will allow us to then look at an MA plot of our results.

```sh
sleuth_live(so)
```
Navigate to ```MA plot``` under the ```analyses``` tab on the webpage.

Take time to explore the results, and the options available using the ```shiny``` app.

## Useful links
* **Ballgown manual**:
    https://bioconductor.org/packages/release/bioc/manuals/ballgown/man/ballgown.pdf
* Ballgown tutorial adapted from:
	Kim, D. et al. Transcript-level expression analysis of RNA- seq experiments with HISAT , StringTie and Transcript-level expression analysis
	of RNA-seq experiments with HISAT , StringTie and Ballgown. Nat. Protoc. 11, 1650–1667 (2016).
* **sleuth manual**:
	https://pachterlab.github.io/sleuth/manual
* sleuth tutorial adapted from: 
	https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html
