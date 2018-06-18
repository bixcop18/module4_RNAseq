# De Novo Transcriptome assembly with Trinity. 


We are going to divide the group in 4 groups. 

### Group 1 

Samples to assemble: 

```
Sample_104B
Sample_105B
Sample_119B
Sample_120A
```

Trimming: Enabled

### Group 2

```
Sample_104B
Sample_105B
Sample_119B
Sample_120A
```

Trimming: Dissabled 

### Group 3

```
DF19
DF5
Sample_103A
Sample_104B
Sample_105B
Sample_118A
Sample_119B
Sample_120A
spike_Z32_rep1
spike_Z32_rep2
```

Trimming: Enabled

### Group 4

```
DF19
DF5
Sample_103A
Sample_104B
Sample_105B
Sample_118A
Sample_119B
Sample_120A
spike_Z32_rep1
spike_Z32_rep2
```

Trimming: Dissabled 


## Inputs

The inputs are in:

```
/var/scratch/Practical4_Example/reads
```


## Bash template 

Trinity requires    

```sh
#!/bin/bash -e
#SBATCH -p batch
#SBATCH -w taurus
#SBATCH -o log/list_files.%N.%j.out
#SBATCH -n 8
#SBATCH --mem=35G

module load trinity/v2.6.6

left="sample1.r1.fastq.gz,sample2.r1.fastq.gz,sample3.r1.fastq.gz"
right="sample1.r2.fastq.gz,sample2.r2.fastq.gz,sample3.r2.fastq.gz"

Trinity --trimmomatic --seqType fq --max_memory 35G --left $left --right $right --CPU 8 --output my_trinity_output --full_cleanup
echo "DONE"
```


using the ```full_cleanup``` to remove at the end all the temporary files for trinity.  The parameter ```CPU``` is consistent with ```-n``` in ```#SBATCH```. Likewise, the memory should be consistent.  

* What are the steps that trinity is doing? What do you think each step does? 
* Can you try the same argumetns in trimmomatic as used on the previous excercises? 




## Comparing the annotations. 

The samples are coming from the same regions of the reference in 