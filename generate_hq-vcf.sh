#!/bin/bash

#SBATCH --job-name=genomeindex
#SBATCH --output=genomeindex_%j.out
#SBATCH --error=genomeindex_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-user=anekimke@stanford.edu
#SBATCH --partition=interactive
#SBATCH --account=default

#arguments:
#1= read pair 1
#2= read pair 2
#3= output name prefix for all files (e.g. sample name)
#4= path to genome, e.g. /labs/feldmanj/Melissa/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/Bowtie2Index/genome
#path to fastq files: /labs/feldmanj/Melissa/12032019spar1/MP1/MP1-TAGGCATG-TCTACTCT_S11_R1_001.fastq
#/labs/feldmanj/Melissa/12032019spar1/MP1/MP1-TAGGCATG-TCTACTCT_S11_R2_001.fastq
source ~/.bashrc

#your .bashrc should contain at least
#module add bowtie/2.0.5
# module add samtools
# module add gatk/3.6
# module add picard
# module add python
# module add fastqc
# module add bedtools
# module add vcftools/0.1.13
# export PICARD_DIR='/srv/gsfs0/software/picard'
# module add java/8u40
# module add snpeff/4.3i

#Mapping of reads using bowtie2
bowtie2 -q --very-sensitive -p 8 -x $4 -1 $1 -2 $2 -S $3_mapped.sam 

#where
#-q  fastq files as entry files
#--very-sensitive  longer but maps better
#-p 8 = nbr of threads
#-x: genome index file
#-1: read1_val.fastq
#-2: read2_val.fastq
#-S: output file name

#SAM to BAM
samtools view -b $3_mapped.sam  -o $3_mapped.bam  

#Where
#-b convert to BAM
#-o output file name

#Sort BAM by coordinates and create index
module purge
module load java/8u40
module load picard
module load samtools
picard SortSam \
   INPUT=$3_mapped.bam \
   OUTPUT=$3_mapped-sorted.bam \
   SORT_ORDER=coordinate \
   CREATE_INDEX=true


#Mark duplicates
picard MarkDuplicates \
   INPUT=$3_mapped-sorted.bam \
   OUTPUT=$3_mapped-sorted-markdup.bam \
   METRICS_FILE=$3_metrics_sorting.txt

#where
#INPUT is your aligned reads .bam file (output.bam from previous step)
#METRICS_FILE is a text (.txt) file with the metrics of the process

#Add ReadGroups
#source ~/.bashrc
module purge
module load java/8u40
module load picard
module load samtools
picard AddOrReplaceReadGroups \
     I=$3_mapped-sorted-markdup.bam \
     O=$3_mapped-sorted-markdup-RG.bam \
     RGID=ID_$3 \
     RGLB=LB_$3 \
     RGPL=illumina \
     RGPU=PU_$3 \
     RGSM=Sample_$3 \
 	 SORT_ORDER=coordinate \
   	 CREATE_INDEX=true

#Read groups have to be added for the next steps, but in our cases we don't really care what they are. For more info, read GATK :)


#Call variants

module load java/latest
java -Xmx3G -jar /srv/gsfs0/software/gatk/gatk-3.6/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R $4.fa \
-I $3_mapped-sorted-markdup-RG.bam \
--genotyping_mode DISCOVERY \
-stand_emit_conf 10 \
-stand_call_conf 30 \
-o $3_raw-variants.vcf

module purge

#where
#R= reference genome.fa (path)
#I= bam file with treated reads
#stand_emit_conf= minimum confidence threshold (phred-scaled) at which the program should emit sites that appear to be possibly variant. Typically 10
#stand_call_conf = minimum confidence threshold (phred-scaled) at which the program should emit variant sites as called. Typically 30.

#Select for high quality (HQ) SNPs and INDELs 
source ~/.bashrc
module load java/latest
#select SNPs
java -Xmx3G -jar /srv/gsfs0/software/gatk/gatk-3.6/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R $4.fa \
   -V $3_raw-variants.vcf \
   -selectType SNP \
   -o $3_raw-variants_SNPs.vcf

#where
#R genome.fa
#V input.vcf
#o output.vcf

#Mark HQ SNPs
java -Xmx3G -jar /srv/gsfs0/software/gatk/gatk-3.6/GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R $4.fa \
  -V $3_raw-variants_SNPs.vcf \
  --filterExpression "QD < 2.00 || FS > 60.00 || MQ < 40.00" \
  --filterName "SNP-FILTER" \
  -o $3_SNPs-HQ.vcf

#where
#R genome.fa
#V raw_snp.vcf
#o output.vcf


#Select INDELs
java -Xmx3G -jar /srv/gsfs0/software/gatk/gatk-3.6/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R $4.fa \
   -V $3_raw-variants.vcf\
   -selectType INDEL \
   -o $3_raw-variants_INDELs.vcf


#where
#1 genome.fa
#2 combined_variants.vcf
#3 output.vcf

#Mark HQ INDELs
java -Xmx3G -jar /srv/gsfs0/software/gatk/gatk-3.6/GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R $4.fa \
  -V $3_raw-variants_INDELs.vcf \
  --filterExpression "QD < 2.00 || FS > 200.00" \
  --filterName "INDEL-FILTER" \
  -o $3_INDELs-HQ.vcf

#where
#1 genome.fa
#2 raw_INDELs.vcf
#3 output.vcf

#Combine  HQ SNPs and INDELs vcf
java -Xmx3G -jar /srv/gsfs0/software/gatk/gatk-3.6/GenomeAnalysisTK.jar \
  -T CombineVariants \
  -R $4.fa \
  --variant $3_SNPs-HQ.vcf \
  --variant $3_INDELs-HQ.vcf \
  -genotypeMergeOptions UNIQUIFY \
  -o $3_HQ-variants.vcf 

module purge

#Annotate VCF
source ~/.bashrc
module load java/latest

java -Xmx4G -jar /srv/gsfs0/software/snpeff/4.3i/snpEff/snpEff.jar \
eff \
WBcel235.86 \
$3_HQ-variants.vcf  >  $3_HQ-variants.ann.vcf


