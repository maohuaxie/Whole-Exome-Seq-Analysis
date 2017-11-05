please check with the following web site: https://github.com/clfougner/MouseExomeSequencing/tree/master/Install
Basically, My pipeline based on the pipeline provied by the above website.

# install R 
sudo add-apt-repository 'deb https://mirror.ibcp.fr/pub/CRAN/bin/linux/ubuntu xenial/' 
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo apt-get update
sudo apt-get install libgdal1-dev libproj-dev libgeos-dev
sudo apt-get install r-base-core

#Format known INDELS file be compatible with mm10.fa
#Unzip known del file
gunzip C57BL_6NJ.mgp.v5.indels.dbSNP142.normed.vcf.gz  

gawk '{ if($0 !~ /^#/) print "chr"$0; else if(match($0,/(##contig=<ID=)(.*)/,m)) print m[1]"chr"m[2]; else print $0 }' C57BL_6NJ.mgp.v5.indels.dbSNP142.normed.vcf > mm10.INDELS.vcf

#Format known SNP file be compatible with mm10.fa
#Unzip known SNP file
gunzip C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz

gawk '{ if($0 !~ /^#/) print "chr"$0; else if(match($0,/(##contig=<ID=)(.*)/,m)) print m[1]"chr"m[2]; else print $0 }' C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf > mm10.SNPs.unordered.vcf

java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar \
SortVcf \
I=mm10.SNPs.unordered.vcf \
O=mm10.C57.SNPs.vcf

# index mm10
bwa index m10.fa

# CreateSequenceDictionary

time java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar \
CreateSequenceDictionary \
REFERENCE=m10.fa \
OUTPUT=mm10.dict

# BWA Alignment
bcl2_N_Fe_R1.fastq.gz  bcl2_N_X_R1.fastq.gz  bcl2_P_Fe_R1.fastq.gz  bcl2_P_X_R1.fastq.gz  ontrol_R2.fastq.gz  
bcl2_N_Fe_R2.fastq.gz  bcl2_N_X_R2.fastq.gz  bcl2_P_Fe_R2.fastq.gz  bcl2_P_X_R2.fastq.gz  contorl_R1.fastq.gz  

time bwa mem -t 12 -M  -R "@RG\tID:bcl2_N_Fe\tSM:bcl2_N_Fe\tLB:WES\tPL:Illumina" /home/emorycpl/ref/mm10/m10.fa bcl2_N_Fe_R1.fastq.gz bcl2_N_Fe_R2.fastq.gz 1>bcl2_N_Fe.sam 2>/dev/null
time bwa mem -t 12 -M  -R "@RG\tID:bcl2_P_Fe\tSM:bcl2_P_Fe\tLB:WES\tPL:Illumina" /home/emorycpl/ref/mm10/m10.fa bcl2_P_Fe_R1.fastq.gz bcl2_P_Fe_R2.fastq.gz 1>bcl2_P_Fe.sam 2>/dev/null
time bwa mem -t 12 -M  -R "@RG\tID:bcl2_N_X\tSM:bcl2_N_X\tLB:WES\tPL:Illumina" /home/emorycpl/ref/mm10/m10.fa bcl2_N_X_R1.fastq.gz bcl2_N_X_R2.fastq.gz 1>bcl2_N_X.sam 2>/dev/null
time bwa mem -t 12 -M  -R "@RG\tID:bcl2_P_X\tSM:bcl2_P_X\tLB:WES\tPL:Illumina" /home/emorycpl/ref/mm10/m10.fa bcl2_P_X_R1.fastq.gz bcl2_P_X_R2.fastq.gz 1>bcl2_P_X.sam 2>/dev/null
time bwa mem -t 12 -M  -R "@RG\tID:Control\tSM:control\tLB:WES\tPL:Illumina" /home/emorycpl/ref/mm10/m10.fa control_R1.fastq.gz control_R2.fastq.gz 1>control.sam 2>/dev/null
#NOTE: USE REAL SAMPLE NAMES FOR RG! MUTECT DOES NOT WORK IF SAMPLE NAMES ARE THE SAME FOR
#TUMOR AND NORMAL SAMPLES

# Sort mapped reads by coordinate

samtools faidx m10.fa 

time java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar \
SortSam \
SORT_ORDER=coordinate \
INPUT=bcl2_N_Fe.sam \
OUTPUT=bcl2_N_Fe.reorder.bam

time java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar \
SortSam \
SORT_ORDER=coordinate \
INPUT=bcl2_P_Fe.sam \
OUTPUT=bcl2_P_Fe.reorder.bam 

time java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar \
SortSam \
SORT_ORDER=coordinate \
INPUT=bcl2_N_X.sam \
OUTPUT=bcl2_N_X.reorder.bam

time java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar \
SortSam \
SORT_ORDER=coordinate \
INPUT=bcl2_P_X.sam \
OUTPUT=bcl2_P_X.reorder.bam

time java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar \
SortSam \
SORT_ORDER=coordinate \
INPUT=control.sam \
OUTPUT=control.reorder.bam

#Mark duplicates
nohup time java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar \
MarkDuplicates \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT \
REMOVE_DUPLICATES=false \
I=bcl2_N_Fe.reorder.bam \
O=bcl2_N_Fe.reorder.dedup.bam \
M=bcl2_N_Fe.metric

nohup time java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar \
MarkDuplicates \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT \
REMOVE_DUPLICATES=false \
I=bcl2_P_Fe.reorder.bam \
O=bcl2_P_Fe.reorder.dedup.bam \
M=bcl2_P_Fe.metric


nohup time java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar \
MarkDuplicates \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT \
REMOVE_DUPLICATES=false \
I=bcl2_N_X.reorder.bam \
O=bcl2_N_X.reorder.dedup.bam \
M=bcl2_N_X.metric


nohup time java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar \
MarkDuplicates \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT \
REMOVE_DUPLICATES=false \
I=bcl2_P_X.reorder.bam \
O=bcl2_P_X.reorder.dedup.bam \
M=bcl2_P_X.metric

nohup time java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar \
MarkDuplicates \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT \
REMOVE_DUPLICATES=false \
I=conrol.reorder.bam \
O=control.reorder.dedup.bam \
M=control.metric

# index marked bam files:

time samtools index bcl2_N_Fe.reorder.dedup.bam
time samtools index bcl2_P_Fe.reorder.dedup.bam
time samtools index bcl2_N_X.reorder.dedup.bam
time samtools index bcl2_P_X.reorder.dedup.bam
time samtools index control.reorder.dedup.bam

# RealignerTargetCreator

java -jar ~/src/gatk/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /home/emorycpl/ref/mm10/m10.fa \
-I bcl2_N_Fe.reorder.dedup.bam \
-known /home/emorycpl/ref/mm10/mm10.INDELS.vcf \
-o bcl2_N_Fe.reorder.dedup.realign.intervals


java -jar ~/src/gatk/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /home/emorycpl/ref/mm10/m10.fa \
-I bcl2_P_Fe.reorder.dedup.bam \
-known /home/emorycpl/ref/mm10/mm10.INDELS.vcf \
-o bcl2_P_Fe.reorder.dedup.realign.intervals

java -jar ~/src/gatk/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /home/emorycpl/ref/mm10/m10.fa \
-I bcl2_N_X.reorder.dedup.bam \
-known /home/emorycpl/ref/mm10/mm10.INDELS.vcf \
-o bcl2_N_X.reorder.dedup.realign.intervals


java -jar ~/src/gatk/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /home/emorycpl/ref/mm10/m10.fa \
-I bcl2_P_X.reorder.dedup.bam \
-known /home/emorycpl/ref/mm10/mm10.INDELS.vcf \
-o bcl2_P_X.reorder.dedup.realign.intervals


java -jar ~/src/gatk/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /home/emorycpl/ref/mm10/m10.fa \
-I control.reorder.dedup.bam \
-known /home/emorycpl/ref/mm10/mm10.INDELS.vcf \
-o control.reorder.dedup.realign.intervals

# 




