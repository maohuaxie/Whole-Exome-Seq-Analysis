#!/usr/bin/bash
#Prepare the tools for emory cbioportal automate process pipeline.
# Exit this script on any error.
set -euo pipefail
cd ~/project/bcl2

time bwa mem -aM -t 12 -R "@RG\tID:sample17\tSM:Seq01\tPL:ILLUMINA\tPI:330" ~/ref/m10/m10noMT.fa index17_GTAGAG_L001-L002_R1_001.fastq index17_GTAGAG_L001-L002_R2_001.fastq > sample17.sam
time samtools view -Sb sample17.sam > sample17.bam 
time samtools sort sample17.bam > sample17.sorted.bam 

java -jar ~/src/picard/build/libs/picard.jar MarkDuplicates INPUT=sample17.sorted.bam \
MAX_RECORDS_IN_RAM=2000000 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true \
METRICS_FILE=sample17.sorted.dups OUTPUT=sample17.sortedDeDup.bam

samtools index sample17.sortedDeDup.bam

java -jar ~/src/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ~/ref/m10/m10noMT.fa \
-I sample17.sortedDeDup.bam \
-known ~/ref/m10/mm10.INDELS.dbSNP142.sorted.vcf \
-o sample17.ForIndelRealigner.intervals

java -jar ~/src/gatk/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R ~/ref/m10/m10noMT.fa \
-I sample17.sortedDeDup.bam \
-known ~/ref/m10/mm10.INDELS.dbSNP142.sorted.vcf \
-targetIntervals sample17.ForIndelRealigner.intervals \
--filter_bases_not_stored \
-log sample17.realigned.log \
-o sample17.realigned.bam 

samtools mpileup -f ~/ref/m10/m10noMT.fa \
sample17.realigned.bam > sample17.realigned.bam.mpileup.txt

samtools mpileup -C 0 -A -B -d 10000 -v -u -f ~/ref/m10/m10noMT.fa \
sample17.realigned.bam | bcftools call -O v -v -c -n 0.05 -p 1 -A -o sample17.bcftools.vcf

java -jar ~/src/snpEff/snpEff.jar  mm10 sample17.bcftools.vcf >sample17.bcftools.snps.vcf

samtools mpileup -f /home/emorycpl/ref/m10/m10noMT.fa sample17.sortedDeDup.bam >sample17.mpileup.txt
samtools mpileup -C 0 -A -B -d 10000 -v -u -f ~/ref/m10/m10noMT.fa \
sample17.sortedDeDup.bam | bcftools call -O v -v -c -n 0.05 -p 1 -A -o sample17m.bcftools.vcf
java -jar ~/src/snpEff/snpEff.jar  mm10 sample17m.bcftools.vcf >sample17m.bcftools.snps.vcf
