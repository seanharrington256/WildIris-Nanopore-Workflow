# Nanopore workflow on UW's WildIris cluster

# UNDER ACTIVE DEVELOPMENT to adapt to WildIris!!!!


Various commands for handling Nanopore data.
![alt text](sequencing-animated.gif)


# Table of Contents I NEED TO EDIT THIS

* [Overview](https://github.com/Joseph7e/Nanopore-Workflow#Overview)  
    * [Basecalling](https://github.com/Joseph7e/Nanopore-Workflow#Basecalling)
    * [Read Processing](https://github.com/Joseph7e/Nanopore-Workflow#Assessing-and-filtering-Nanopore-Data)  
       * [Adapter Trimming](https://github.com/Joseph7e/Nanopore-Workflow#Trim-adapters-with-porechop) 
       * [Filter Reads](https://github.com/Joseph7e/Nanopore-Workflow#Accessing-and-filtering-Nanopore-Data)  
       * [Assessment of reads](https://github.com/Joseph7e/Nanopore-Workflow#Read-Assessment-with-Nanoplot)  
    * [Nanopore-only Assembly](https://github.com/Joseph7e/Nanopore-Workflow#Assessing-and-filtering-Nanopore-Data)  
       * [Canu](https://github.com/Joseph7e/Nanopore-Workflow#Accessing-and-filtering-Nanopore-Data)  
       * [Miniasm](https://github.com/Joseph7e/Nanopore-Workflow#Trim-adapters-with-porechop) 
    * [Assembly Polishing](https://github.com/Joseph7e/Nanopore-Workflow#Assessing-and-filtering-Nanopore-Data)

# Overview

How data was produced etc.  
We will use Illumina data in some instances throughout this tutorial (hybrid-assembly, polishing, and assessment). You'll want to run adapter trimming on the illumina reads prior to using them in any instance. IN addition, you will need an illumina-only assembly for the quality assessment of the nanopore assemblies. FOllow my main genome assembly tutorial to produce a high quality assembly and to determine information such as insert size and estimated genome size.  

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

```
# basecalling with guppy
guppy_basecaller -i <inputdir> -s <output_dir> --flowcell FLO-MIN106 --kit SQK-LSK109 –fast5_out -r -t 15
```


### Repair corrupted read files produced with guppy


The fastq files need to be decompressed for both methods. For my custom workflow, the column for the sort program must match up with the order the fastq was produced (the numbers in the filename). Be sure to test that this works properly.

#### Joes' custom method.
```
cd <fastq_directory>
ls *.fastq | sort -t'_' -k2 -n | xargs cat - > ../raw_reads.fastq
/mnt/lustre/hcgs/joseph7e/scripts/nanopore_fix_fastq.py <fastq_from_above> > <fixed.fastq>
```

#### Nanopore community way
https://community.nanoporetech.com/posts/fastq-errors-on-gridion-an
```
pip install ont-fastq-deconcatenate
apt-get update && apt-get install python3-pip
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
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
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
canu -d canu-assembly -p filt genomeSize=3.5m gnuplotTested=true useGrid=false -nanopore-raw filtered.fq
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
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
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


Make yet another SLURM script


```
#!/bin/bash

#SBATCH --job-name spades_hyb
#SBATCH -A wy_t3_2022
#SBATCH -t 0-08:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL@EMAIL.com
#SBATCH -e err_spades_hyb_%A.err
#SBATCH -o std_spades_hyb_%A.out

cd ~/nanopore

module load gcc/11.2.0 spades/3.15.3

forward='KphilA_CGCTCATT-AGGCGAAG_L001_R1_001.fastq.gz'
reverse='KphilA_CGCTCATT-AGGCGAAG_L001_R2_001.fastq.gz'
nanopore='Kphil49844-ONT-1.fastq.gz'

spades.py -t 24 -1 $forward -2 $reverse --nanopore $nanopore  -o spades-hybrid-only -m 10 --only-assembler
```



### Hybrid Assembly w/ Masurca


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

```
# starting data
genome=canu-assembly.fasta 
nanopore_reads=../../filtered_1000_80.fastq

# index genome
bwa index $genome

# map ont reads to assembly
bwa mem -t 24 -x ont2d $genome $nanopore_reads > mapping-filteredONT.sam

# polish with racon and produce new consensus sequence
racon -m 8 -x -6 -g -8 -w 500 -t 24 $nanopore_reads mapping-filteredONT.sam <genome.fasta> > racon.fasta

# repeat (2)
bwa index racon.fasta
bwa mem -t 24 -x ont2d racon.fasta $nanopore_reads > mapping-filteredONT2.sam
racon -m 8 -x -6 -g -8 -w 500 -t 24 $nanopore_reads mapping-filteredONT2.sam racon.fasta > racon_round2.fasta

# repeat (3)
bwa index racon.fasta
bwa mem -t 24 -x ont2d racon2.fasta $nanopore_reads > mapping-filteredONT3.sam
racon -m 8 -x -6 -g -8 -w 500 -t 24 $nanopore_reads mapping-filteredONT3.sam racon2.fasta > racon_round3.fasta

# repeat (4)
bwa index racon.fasta
bwa mem -t 24 -x ont2d racon3.fasta $nanopore_reads > mapping-filteredONT4.sam
racon -m 8 -x -6 -g -8 -w 500 -t 24 $nanopore_reads mapping-filteredONT4.sam racon3.fasta > racon-round4.fasta

# clean up all the iterations
```

## Medaka


## Pilon
Here we will polish an assembly using Illumina reads. You can do this right away, or after the racon/medaka polishing.
  
#### Step 1: Map Illumina reads to Assembly/Read FASTA.
The script below outputs a lot of extra (useful) data. We only need the mapping file.
```
sbatch ~/scripts/bwa_index_and_mapV2.sh <reference.fasta> <forward_reads> <reverse_reads> <sample_name>
mkdir assembly_ploshing && cd assembly_polishing
mv ../bwa_mapping*/sorted_mapped.bam* ./
```

#### Step 2: Prepare genome chunks for array job.
```
grep ">" <assembly.fasta> \
| tr -d ">" \
| shuf \
| split -d -l 200 - genomechunk.
rename genomechunk.0 genomechunk. genomechunk.0*
mkdir pilon
mv genomechunk* pilon/
```

#### Step 3: Edit pilon script and run pilon
At this point you should be in a directory with four files. One mapping file with its index, a directory named pilon (which contains all genome chunk data), and a copy or symlink of your reference assembly.
```
# copy over generic slurm script
cp ~/scripts/nanopore_pilon.slurm pilon.slurm
```
You need to edit a few things in this script.
#SBATCH --array=0-X%8 (change X to equal the number of genome chunks found in pilon dir.
genome="YOUR GENOME FILE HERE"
--frags "NAME OF MAPPING FILE HERE"

```
sbatch ./pilon.slurm
```

#### Step 4: View logs and concatenate polished genome
```
mkdir logs_pilon
mv *.log logs_pilon
cat pilon/*.fasta > polished_genome.fasta
```

#### Step 5: Rinse and repeat
Repeat the entire process on the newly polished genome. Then again and agin, until you're happy.






## Scaffolding w/ LINKS
Note that this process requires a ton of memory. Maybe use a high memory node or used a reduced set of reads.
```
LINKS -f hybrid_assembly_fixed.fasta  -s <txt_file_with_nanopore_read_paths> -b <output_base>
```
```
sbatch ~/nanopore_LINKS.sh <assembly> <nanopore_reads>
```
## Assembly Assessment Scripts
Quast for contiguity, BUSCO for completness, BWA for quality and correctness.
```
quast.py miniasm.fasta
~/scripts/quality_check__genome_busco.sh <assembly.fasta>
```




step 1.) map the reads to the assembly


- Got up to Racon & Medaka (mostly)

- come back to repairing reads


- Canu takes long and is set up in a dumb way to submit jobs
	- check options Joe uses - genome size
- for WildIris intro tutorial:
	- add in editing of text files using bbedit/notepad++
	- worker nodes don't inherit the loaded modules from the login node
	- I use '.slurm' for my slurm scripts to make them easily identifiable
	- useful shortcuts like ctrl+a, cd -, etc.	
- prob don't need the scripts in the git repo - they're not setup for WildIris
- Minimap returns empty file - what's UP???
- Masurca has nothing for it - add or skip????



