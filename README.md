# Nanopore workflow on UW's WildIris cluster

# UNDER ACTIVE DEVELOPMENT to adapt to WildIris!!!!
- Notes for conversion:
	- Still to run:
		- Pilon runs out of memory
	- Canu - check options Joe uses
	- don't need all scripts in the git repo - they're not setup for WildIris

	- Minimap running now with more memory




Various commands for handling Nanopore data.
![alt text](sequencing-animated.gif)


# Table of Contents

* [Overview](https://github.com/seanharrington256/WildIris-Nanopore-Workflow#Overview)  
* [Basecalling](https://github.com/seanharrington256/WildIris-Nanopore-Workflow#Basecalling)
* [Read Processing](https://github.com/seanharrington256/WildIris-Nanopore-Workflow#Assessing-and-filtering-Nanopore-Data)  
    * [Adapter Trimming](https://github.com/seanharrington256/WildIris-Nanopore-Workflow#Trim-adapters-with-porechop) 
    * [Filter Reads](https://github.com/seanharrington256/WildIris-Nanopore-Workflow#Accessing-and-filtering-Nanopore-Data)  
    * [Assessment of reads](https://github.com/seanharrington256/WildIris-Nanopore-Workflow#Read-Assessment-with-Nanoplot)  
* [Nanopore-only Assembly](https://github.com/seanharrington256/WildIris-Nanopore-Workflow#Assessing-and-filtering-Nanopore-Data)  
    * [Canu](https://github.com/seanharrington256/WildIris-Nanopore-Workflow#Accessing-and-filtering-Nanopore-Data)  
    * [Miniasm](https://github.com/seanharrington256/WildIris-Nanopore-Workflow#Trim-adapters-with-porechop) 
* [Illumina and Nanopore Hybrid Assembly](https://github.com/seanharrington256/WildIris-Nanopore-Workflow#Illumina-and-Nanopore-Hybrid-Assembly)
* [Assembly Polishing](https://github.com/seanharrington256/WildIris-Nanopore-Workflow#Assessing-and-filtering-Nanopore-Data)

# Overview

How data was produced etc.  
We will use Illumina data in some instances throughout this tutorial (hybrid-assembly, polishing, and assessment). You'll want to run adapter trimming on the illumina reads prior to using them in any instance. In addition, you will need an illumina-only assembly for the quality assessment of the nanopore assemblies. FOllow my main genome assembly tutorial to produce a high quality assembly and to determine information such as insert size and estimated genome size.  

## Sample Preperation and Sequencing

Recommended/standardized workflow - 
  
alternative methods
  
multiplexing 24 max.



## Get the data and set up our session

The data for this exercise are already on WildIris in the wy\_t3_2022 project directory. If you are not already logged in, do so.

Make a new directory called "nanopore" in your home directory and copy the data there:

```
mkdir ~/nanopore/
cp -r /project/wy_t3_2022/example-nanopore-illumina/* ~/nanopore/
```

Note: if you are not in this workshop and want to run through this protocol, please contact me.


Now we'll need to set up an interactive session so that we're not doing anything intensive on the login node. We'll also load up a module that contains a lot of the tools we'll use:

```
salloc -A wy_t3_2022 -t 0-05:00 --mem=10G --cpus-per-task=2
module load gentools/1.0.0
```

Note that `gentools` is not a single piece of software, it is an environment on WildIris that contains several different programs that we are using in this tutorial. If you are running through this on a different cluster or your own computer, you'll have to either load up multiple modules or install all of these programs into your own conda environment.



# Read processing
![read processing flowchart for nanopore and illumina data](read_processing_workflow.png)

## Basecalling

Basecalling is typically run automatically on the sequencing instrument. For example the GirdIOn will run the Guppy basecaller as soon as the fast5s are produced. When using a MinION sequencer plugged into a standard Mac desktop, we have found basecalling to be very slow, and therefore would recommend running guppy on WildIris. This is still somewhat slow, but you can at least speed it up by giving it a large number of threads. The fastest way to run Guppy is on a GPU rather than a CPU, but WildIris does not currently have any GPU nodes (Teton, UW's other cluster does). 


Our example data here are already basecalled, and so we will skip this step, but it is documented here for future reference.

 
Manual: https://denbi-nanopore-training-course.readthedocs.io/en/latest/basecalling/basecalling.html
Alternative Tools: Albacore, DeepNano-blitz, minKNOW, Chiron, Bonito
Comparison of base callers: https://github.com/rrwick/Basecalling-comparison

**Again, we're skipping this step** because our data are already base called.

```
# basecalling with guppy
guppy_basecaller -i <inputdir> -s <output_dir> --flowcell FLO-MIN106 --kit SQK-LSK109 –fast5_out -r -t 15
```


### Repair corrupted read files produced with guppy


**The fastq files need to be decompressed for both methods**. For Joe's custom workflow, the column for the sort program must match up with the order the fastq was produced (the numbers in the filename). Be sure to test that this works properly.

#### Joes' custom method.
```
cd ~/nanopore
# ls *.fastq | sort -t'_' -k2 -n | xargs cat - > raw_reads_prepped.fastq  # This line is only needed when your reads are in multiple files - make sure this only includes nanopore reads if you run it
python /project/wy_t3_2022/fix_guppy_fastq.py raw_reads_prepped.fastq > fixed_raw_reads.fastq.gz
```

#### Nanopore community way
https://community.nanoporetech.com/posts/fastq-errors-on-gridion-an
```
pip3 install ont-fastq-deconcatenate
fix_concatenated_fastqs -i <path_to_folder_of_fastqs>
```

## Trim adapters with porechop

The porechop "check\_reads" option removes the need to specify adapters. It will automatically check and detemrine which ones to remove.

```
porechop --check_reads 1000 -i Kphil49844-ONT-1.fastq.gz -o adapter_trimmed.fastq
```


## Filter reads with filtlong
This step is optional. I did not run it through my first attempts.

```
filtlong --min_mean_q 80 --min_length 2000 adapter_trimmed.fastq > filtered.fq
```

## Read Assessment with Nanoplot
https://github.com/wdecoster/NanoPlot
I usually run this on the raw reads and after any adapter/quality trimming. Run time ~ 2 hrs per 10 GB


```
NanoPlot --fastq filtered.fq --threads 2 -o nanoplot_out
```

# Assemblies overview

![Assembly methodology for microbial genomes with Nanopore and/or Illumina data as a base](assemblies_workflow.png)


# Nanopore only assembly

## Canu

Reads > 1kb an genome size of 130MB

De Bruijn graph contigs were generated with Platanus


Canu will take a little while, so let's create a new slurm script to run this. Either use `nano canu.slurm` to create and edit a file for this script, or do this by using Cyberduck to create a new file and then edit it in your tet editor of choice. Put the below into this slurm script -- replace `YOUR_EMAIL@EMAIL.com` with your actual email if you want to get notifications.

```
#!/bin/bash

#SBATCH --job-name canu
#SBATCH -A wy_t3_2022
#SBATCH -t 0-08:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL@EMAIL.com
#SBATCH -e err_canu_%A.err
#SBATCH -o std_canu_%A.out

source ~/.bash_profile


# load modules necessary
module load gentools/1.0.0

# Set working directory
cd ~/nanopore/

# Run canu
canu -d canu-assembly -p filt genomeSize=3.5m maxMemory=32 maxThreads=24 gnuplotTested=true useGrid=false -nanopore-raw filtered.fq
```

Once you've created and saved that file, submit the job:

`sbatch canu.slurm`

Check the status of the job once it's running:

`squeue -u YOUR_USERNAME`




## Miniasm


Make a new directory and enter it.

```
mkdir miniasm_assembly
cd miniasm_assembly
```

Then run minimap2 and miniasm:

# THIS CURRENTLY DOESN'T WORK

```
#!/bin/bash

#SBATCH --job-name miniasm
#SBATCH -A wy_t3_2022
#SBATCH -t 0-08:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL@EMAIL.com
#SBATCH -e err_miniasm_%A.err
#SBATCH -o std_miniasm_%A.out

source ~/.bash_profile


# load modules necessary
module load gentools/1.0.0

# Set working directory
cd ~/nanopore/miniasm_assembly

# Run miniasm


minimap2 -x ava-ont ../adapter_trimmed.fastq ../adapter_trimmed.fastq | gzip -1 > mapped.paf.gz
miniasm -f ../adapter_trimmed.fastq mapped.paf.gz > miniasm.gfa
awk '/^S/{print ">"$2"\n"$3}' miniasm.gfa > miniasm.fasta
```


# Illumina and Nanopore Hybrid Assembly

### Hybrid Assembly w/ spades
https://www.ncbi.nlm.nih.gov/pubmed/26589280
Run this as a normal spades job but specify the --nanopore reads. Note that with low coverage data I found that the nanopore data does not significantly improve the assembly.


Make yet another SLURM script, `spades_hyb.slurm`


```
#!/bin/bash

#SBATCH --job-name spades_hyb
#SBATCH -A wy_t3_2022
#SBATCH -t 0-08:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL@EMAIL.com
#SBATCH -e err_spades_hyb_%A.err
#SBATCH -o std_spades_hyb_%A.out

cd ~/nanopore

module load gcc/11.2.0 spades/3.15.3

forward='KphilA_CGCTCATT-AGGCGAAG_L001_R1_001.fastq.gz'
reverse='KphilA_CGCTCATT-AGGCGAAG_L001_R2_001.fastq.gz'
nanopore='Kphil49844-ONT-1.fastq.gz'

spades.py -t 24 -1 $forward -2 $reverse --nanopore $nanopore  -o spades-hybrid-only -m 30 --only-assembler
```

Submit the job and check on it:

```
sbatch spades_hyb.slurm.slurm
squeue -u YOUR_USERNAME
```


### Hybrid Assembly w/ Masurca

Masurca is another assembler that we can use for hyrbid assembly. Masurca requires a config file that specifies the parameters to run it. Full documentation is here: [https://github.com/alekseyzimin/masurca](https://github.com/alekseyzimin/masurca).

Let's start by making a config file called `masurca_config.txt`. **Note that here you have to use the full, absolute path to the data files and cannot use ~, you will need to edit this to your specific paths**

```
DATA
PE= pe 500 50  /FULL/PATH_TO/nanopore/KphilA_CGCTCATT-AGGCGAAG_L001_R1_001.fastq.gz /FULL/PATH_TO/nanopore/KphilA_CGCTCATT-AGGCGAAG_L001_R2_001.fastq.gz
NANOPORE=/FULL/PATH_TO/nanopore/Kphil49844-ONT-1.fastq.gz
END

PARAMETERS
NUM_THREADS=16
JF_SIZE=2000000000
USE_LINKING_MATES=0
GRAPH_KMER_SIZE=auto
SOAP_ASSEMBLY = 0
END
```


Then, run the Masurca script, which itself generates an executable script to run the actual assembly:

```
module load gentools
masurca masurca_config.txt
```

Make a new directory to run that file from and move it there:

```
mkdir masurca_out
mv assemble.sh masurca_out
```


Then we can make a slurm script called `masurca.slurm` to submit that executable script as a job:


```
#!/bin/bash

#SBATCH --job-name masurca
#SBATCH -A wy_t3_2022
#SBATCH -t 0-08:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL@EMAIL.com
#SBATCH -e err_masurca_%A.err
#SBATCH -o std_masurca_%A.out

cd ~/nanopore/masurca_out

module load gentools

./assemble.sh
```

Submit it and check that it's running:

```
sbatch masurca.slurm
squeue -u YOUR_USERNAME
```



# Genome Assembly Polishing

There are many routes you can take, and many programs to choose. We will be using three. An initial polishing of a nanopore assembly with nanopore reads with racon. Further polishing using nanopore reads with makon, and a final polishing of the assembly using illumina reads with pilon. You can do any combination of these tools, for example skip right to pilon. Try them all and compare the results!


Racon
Minpolish
nanopolish
nextpolish
pilon

![Pipeline for genome assembly assessment, polishing, and merging](polishing_merging_assessment_workflow.png)

## Racon and Medaka

Manuscript:  
Tutoral: https://denbi-nanopore-training-course.readthedocs.io/en/latest/polishing/medaka/racon.html  

Racon can be used as a polishing tool after the assembly with either Illumina data or data produced by third generation of sequencing. The type of data inputed is automatically detected.  
  
Racon takes as input only three files: contigs in FASTA/FASTQ format, reads in FASTA/FASTQ format and overlaps/alignments between the reads and the contigs in MHAP/PAF/SAM format. Output is a set of polished contigs in FASTA format printed to stdout. All input files can be compressed with gzip (which will have impact on parsing time).
  
The medaka documentation advises to do four rounds with racon before polishing with medaka, since medaka has been trained with racon polished assemblies. We need to iterate all the steps four times.

Create a slurm script called `racon.slurm`:

```
#!/bin/bash

#SBATCH --job-name racon
#SBATCH -A wy_t3_2022
#SBATCH -t 0-08:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL@EMAIL.com
#SBATCH -e err_racon_%A.err
#SBATCH -o std_racon_%A.out

# set the working directory
cd ~/nanopore

# Load up modules
module load gentools/1.0.0 gcc/11.2.0 bwa/0.7.17

# starting data
genome=canu-assembly/filt.contigs.fasta
nanopore_reads=filtered.fq

# index genome
bwa index $genome

# map ont reads to assembly
bwa mem -t 24 -x ont2d $genome $nanopore_reads > mapping-filteredONT.sam

# polish with racon and produce new consensus sequence
racon -m 8 -x -6 -g -8 -w 500 -t 24 $nanopore_reads mapping-filteredONT.sam $genome > racon.fasta

# repeat (2)
bwa index racon.fasta
bwa mem -t 24 -x ont2d racon.fasta $nanopore_reads > mapping-filteredONT2.sam
racon -m 8 -x -6 -g -8 -w 500 -t 24 $nanopore_reads mapping-filteredONT2.sam racon.fasta > racon_round2.fasta

# repeat (3)
bwa index racon_round2.fasta
bwa mem -t 24 -x ont2d racon_round2.fasta $nanopore_reads > mapping-filteredONT3.sam
racon -m 8 -x -6 -g -8 -w 500 -t 24 $nanopore_reads mapping-filteredONT3.sam racon_round2.fasta > racon_round3.fasta

# repeat (4)
bwa index racon_round3.fasta
bwa mem -t 24 -x ont2d racon_round3.fasta $nanopore_reads > mapping-filteredONT4.sam
racon -m 8 -x -6 -g -8 -w 500 -t 24 $nanopore_reads mapping-filteredONT4.sam racon_round3.fasta > racon-round4.fasta

```

Submit it and check that it's running:

```
sbatch racon.slurm
squeue -u YOUR_USERNAME
```


## Medaka

When this finishes, we can then feed the final assembly into Medaka, in a slurm script called `medaka.slurm`:

```
#!/bin/bash

#SBATCH --job-name medaka
#SBATCH -A wy_t3_2022
#SBATCH -t 0-08:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL@EMAIL.com
#SBATCH -e err_medaka_%A.err
#SBATCH -o std_medaka_%A.out
#SBATCH --partition=wildiris-phys

# load modules
module load gentools samtools bcftools

# move to the right directory
cd ~/nanopore

# run medaka
medaka_consensus -i Kphil49844-ONT-1.fastq.gz -d racon-round4.fasta -o medaka_out -t 24 -m r941_min_high_g303
```
* Note that we had to add `--partition=wildiris-phys` to specify that medaka should be run on a physical and not virtual node. Don't worry too much about what this means, it has to do with specifics of medaka and WildIris that won't typically apply.


Submit it and check that it's running:

```
sbatch medaka.slurm
squeue -u YOUR_USERNAME
```



## Pilon

Here we will polish an assembly using Illumina reads. You can do this right away, or after the racon/medaka polishing. We'll demonstrate this on the unpolished assembly so that we can start this up now:

  
#### Step 1: Map Illumina reads to Assembly/Read FASTA. This should look familiar from the Illumina-only tutorial.

```
# load modules
module load gcc bwa
# index the canu assembly with bwa
bwa index canu-assembly/filt.contigs.fasta 
# map the illumina reads to the canu assembly
bwa mem -t 2 canu-assembly/filt.contigs.fasta ~/nanopore/KphilA_CGCTCATT-AGGCGAAG_L001_R1_001.fastq.gz ~/nanopore/KphilA_CGCTCATT-AGGCGAAG_L001_R2_001.fastq.gz > mapped_to_canu.sam
# Remove sequencing reads that did not match to the assembly and convert the SAM to a BAM.
samtools view -@ 2 -Sb  mapped_to_canu.sam  | samtools sort -@ 2 - sorted_mapped_canu
# index the new bam file
samtools index sorted_mapped_canu.bam
```


#### Step 2: Set up the slurm script to submit the job:

Make a slurm script `pilon.slurm`:

```
#!/bin/bash

#SBATCH --job-name pilon
#SBATCH -A wy_t3_2022
#SBATCH -t 0-08:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL@EMAIL.com
#SBATCH -e err_pilon_%A.err
#SBATCH -o std_pilon_%A.out

# load modules
module load gentools

# move to the right directory
cd ~/nanopore

# run pilon
pilon --genome canu-assembly/filt.contigs.fasta --bam sorted_mapped_canu.bam --outdir pilon_out
```


#### Step 3: Rinse and repeat
Repeat the entire process on the newly polished genome. Then again and agin, until you're happy.


This may require splitting the job up into smaller chunks or running on a system with more memory.



## Scaffolding w/ LINKS
Note that this process also requires a ton of memory. Maybe use a reduced set of reads or run on a high memory node on a different system if this doesn't run.


LINKS takes a "file of filenames" (.fof) file as the input specifying the Nanopore reads file. Let's make this real quick:

```
# this again has to be a full path that does not include "~"
readlink -f ~/nanopore/Kphil49844-ONT-1.fastq.gz > nanopore_reads.fof
```

Then use that to run LINKS, again with a SLURM script, this time called `links.slurm`:

```
#!/bin/bash

#SBATCH --job-name links
#SBATCH -A wy_t3_2022
#SBATCH -t 0-08:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL@EMAIL.com
#SBATCH -e err_links_%A.err
#SBATCH -o std_links_%A.out

# load modules
module load gentools

# move to the right directory
cd ~/nanopore

# run links
LINKS -f canu-assembly/filt.contigs.fasta -s nanopore_reads.fof -b out_links
```


## Assembly Assessment Scripts


Quast for contiguity, BUSCO for completness, BWA for quality and correctness.

```
quast.py miniasm.fasta
~/scripts/quality_check__genome_busco.sh <assembly.fasta>
```




step 1.) map the reads to the assembly





