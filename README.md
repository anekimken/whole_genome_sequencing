# whole_genome_sequencing
 Whole genome sequencing workflow. Based on code from Chloe Girard and Melissa Pickett
 
## File structure in the cloud
Data folder: /labs/PI_SUNetID/your_name/your_experiment

Script folder: /labs/PI_SUNetID/your_name/scripts


## Load files to the cluster
1. Load files to SCG Genomics Cluster using sftp
 - Make sure your .fastq files are unzipped.
 - Put the files in a directory with your name and a subdirectory with your experiment (eg your_name/your_experiment)

2. Upload reference sequence in the same directory
 - This should be a .tar.gz file.
 
## Prepare reference sequence
This can go in the folder for your experiment
1. Unzip it with the command: tar xvzf file.tar.gz
2. Convert it to a dict file
```
module add picard
picard CreateSequenceDictionary R=genome.fa O=genome.dict
module load samtools
samtools faidx genome.fa
```

## Upload analysis script to your scripts folder
/labs/PI_SUNetID/your_name/scripts

## Run job
Submit batch to slurm:
sbatch analysis_script_location read_1.fastq read_2.fastq output_name  reference_genome_location

For example: 
sbatch ../scripts/generate_hq-vcf.sh Adam-TAAGGCGA-ATAGAGAG_S8_R1_001.fastq Adam-TAAGGCGA-ATAGAGAG_S8_R2_001.fastq CB1267  ./Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/Bowtie2Index/genome
	- note there is no file extension at the end of the reference genome filename

Wait for the script to run on the cluster!

If you're doing a mutagenesis screen, you'll have to add some steps here to remove mutations present in your parent strain.

## Export data



Add modules to your .bashrc on SCG
	- enter command cat>~/.bashrc	
	- add modules to your account by pasting these lines after starting cat:
		module add bowtie/2.0.5
		module add samtools
		module add gatk/3.6
		module add picard
		module add python
		module add fastqc
		module add bedtools
		module add vcftools/0.1.13
		export PICARD_DIR='/srv/gsfs0/software/picard'
		module add java/8u40
		module add snpeff/4.3i
	- ctrl-c to quit out of cat

Add scripts for analysis
	- upload bash shell script that actually runs the analysis
 
