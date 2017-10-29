# paper；https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4812767/a whole-exome sequencing (WES) study was performed in 161 NPC cases and 895 controls of Southern Chinese descent
set -ueo pipefail
for id in  SRR1139956 SRR1139958 SRR1139966 SRR1139973 SRR1139999 SRR1140007 SRR1140015 SRR1140023 SRR1140044 SRR1140045
do echo $id
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP035/SRP035573/$id/$id.sra
done

# do it by for loop, download the sra files and convert to fastq files
for id in SRR1139956 SRR1139966 SRR1139973 SRR1139999 SRR1140007 SRR1140015 SRR1140023 SRR1140044 SRR1140045
do
nohup fastq-dump --gzip --split-files $id
done
# do it in a hard way
set -ueo pipefail
nohup fastq-dump --gzip --split-files SRR1139956
nohup fastq-dump --gzip --split-files SRR1139966
nohup fastq-dump --gzip --split-files SRR1139973
nohup fastq-dump --gzip --split-files SRR1139999
nohup fastq-dump --gzip --split-files SRR1140007
nohup fastq-dump --gzip --split-files SRR1140015
nohup fastq-dump --gzip --split-files SRR1140023
nohup fastq-dump --gzip --split-files SRR1140044
nohup fastq-dump --gzip --split-files SRR1140045

# for demo, I will try 3 cases 
for id in SRR1139956 SRR1139958 SRR1139966 SRR1139973 
do
nohup fastq-dump --gzip --split-files $id
done

# please go to  ftp://ftp.broadinstitute.org/bundle/hg19/ download the following files:
ucsc.hg19.fasta
ucsc.hg19.fasta.fai
1000G_omni2.5.hg19.vcf
1000G_omni2.5.hg19.vcf.idx
1000G_phase1.indels.hg19.vcf
1000G_phase1.indels.hg19.vcf.idx
dbsnp_137.hg19.vcf
dbsnp_137.hg19.vcf.idx
hapmap_3.3.hg19.vcf
hapmap_3.3.hg19.vcf.idx
Mills_and_1000G_gold_standard.indels.hg19.vcf
Mills_and_1000G_gold_standard.indels.hg19.vcf.idx

# let's have a briefly look at the quality of the sequencing data
ls *.gz |xargs fastqc -o ./ -t 5
###############################
### pipe for SNP-calling ######
###############################
# let's try GATK golden standard pipeline

#####################################################
################ Step 1 : Alignment #################
#####################################################
# preparing a reference for use with BWA and GATK
# Generate the BWA index
bwa index ucsc.hg19.fasta
# you can do for loop or do it in a hard way
for SAMPLE in SRR1139966 SRR1139973
do
R1=${SAMPLE}_1.fastq.gz
R2=${SAMPLE}_2.fastq.gz
nohup time bwa mem -t 12 -M  -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tLB:WES\tPL:Illumina" /home/maohuaxie/ref/h19/ucsc.hg19.fasta $R1 $R2 1>$SAMPLE.sam 2>/dev/null
done

time bwa mem -t 12 -M  -R "@RG\tID:SRR1139966\tSM:SRR1139966\tLB:WES\tPL:Illumina" /home/emorycpl/ref/h19/ucsc.hg19.fasta \
/home/emorycpl/project/exome/SRR1139966_1.fastq.gz SRR1139966_2.fastq.gz >SRR1139966.sam 

time bwa mem -t 12 -M  -R "@RG\tID:SRR1139973\tSM:SRR1139973\tLB:WES\tPL:Illumina" /home/emorycpl/ref/h19/ucsc.hg19.fasta \
/home/emorycpl/project/exome/SRR1139973_1.fastq.gz SRR1139973_2.fastq.gz >SRR1139973.sam 

#####################################################
################ Step 2: Sort and Index #############
#####################################################

java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar SortSam SORT_ORDER=coordinate INPUT=SRR1139966.sam OUTPUT=SRR1139966.bam
java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar SortSam SORT_ORDER=coordinate INPUT=SRR1139973.sam OUTPUT=SRR1139973.bam

nohup time samtools index SRR1139966.bam
nohup time samtools index SRR1139973.bam

rm -rf SRR1139966.sam SRR1139973.sam

#####################################################
################ Step 3: Basic Statistics ###########
#####################################################

samtools flagstat SRR1139966.bam 
samtools stats  SRR1139966.bam 

#####################################################
####### Step 4: multiple filtering for bam files ####
#####################################################

###MarkDuplicates###

java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar MarkDuplicates \
INPUT=SRR1139966.bam OUTPUT=SRR1139966_marked.bam METRICS_FILE=SRR1139966.metrics
rm SRR1139966.bam  

java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar ~/src/picard/build/libs/picard.jar MarkDuplicates \
INPUT=SRR1139973.bam OUTPUT=SRR1139973_marked.bam METRICS_FILE=SRR1139973.metrics
rm SRR1139973.bam  
samtools index SRR1139966_marked.bam
samtools index SRR1139973_marked.bam

##
java  -Xmx10g -jar -jar ~/src/picard/build/libs/picard.jar MarkDuplicates.jar \
I=SRR1139966.bam \
O=SRR1139966_marked.bam \
VALIDATION_STRINGENCY=LENIENT \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
REMOVE_DUPLICATES= false \
M=SRR1139966.metrics

### we also can try the following code from http://www.jianshu.com/p/938d362fc48d
java -jar ~/src/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES= false MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 INPUT=SRR1139973.bam OUTPUT=SRR1139973_marked.bam METRICS_FILE=SRR1139973.metrics
java -jar ~/src/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES= false MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 INPUT=SRR1139966.bam OUTPUT=SRR1139966_marked.bam METRICS_FILE=SRR1139966.metrics
###FixMateInfo### some old version GATK need do this step, new version does not require do this step, so you may skip this step

java -Djava.io.tmpdir=/home/emorycpl/tmp  -Xmx40g -jar ~/src/picard/build/libs/picard.jar FixMateInformation \
INPUT=SRR1139966_marked.bam OUTPUT=SRR1139966_marked_fixed.bam SO=coordinate
samtools index SRR1139966_marked_fixed.bam

java -Djava.io.tmpdir=/home/emorycpl/tmp  -Xmx40g -jar ~/src/picard/build/libs/picard.jar FixMateInformation \
INPUT=SRR1139973_marked.bam OUTPUT=SRR1139973_marked_fixed.bam SO=coordinate
samtools index SRR1139966_marked_fixed.bam

#####################################################
####### Step 5: gatk process bam files ##############
#####################################################

### SplitNCigar ###

java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T SplitNCigarReads \
-R /home/emorycpl/ref/h19/ucsc.hg19.fasta  -I SRR1139966_marked.bam  -o SRR1139966_marked_fixed_split.bam \
-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS --fix_misencoded_quality_scores 
## --fix_misencoded_quality_scores only if phred 64 

java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T SplitNCigarReads \
-R /home/emorycpl/ref/h19/ucsc.hg19.fasta  -I SRR1139973_marked.bam  -o SRR1139973_marked_fixed_split.bam \
-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS --fix_misencoded_quality_scores 


## for GATK snp calling, you have to create fa.dict file
java -jar ~/src/picard/build/libs/picard.jar CreateSequenceDictionary REFERENCE=/home/emorycpl/ref/h19/ucsc.hg19.fasta OUTPUT=/home/emorycpl/ref/h19/ucsc.hg19.fasta.dict

###RealignerTargetCreator###

java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator \
-I SRR1139966_marked.bam -R /home/emorycpl/ref/h19/ucsc.hg19.fasta -o SRR1139966_marked_target.intervals \
-known /home/emorycpl/ref/h19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
-known /home/emorycpl/ref/h19/1000G_phase1.indels.hg19.sites.vcf -nt 5 --fix_misencoded_quality_scores

java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator \
-I SRR1139973_marked.bam -R /home/emorycpl/ref/h19/ucsc.hg19.fasta -o SRR1139973_marked_target.intervals \
-known /home/emorycpl/ref/h19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
-known /home/emorycpl/ref/h19/1000G_phase1.indels.hg19.sites.vcf -nt 5 --fix_misencoded_quality_scores


###IndelRealigner###

java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T IndelRealigner \
-I SRR1139966_marked.bam  -R /home/emorycpl/ref/h19/ucsc.hg19.fasta -targetIntervals SRR1139966_marked_target.intervals \
-o SRR1139966_marked_realigned.bam \
-known /home/emorycpl/ref/h19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
-known /home/emorycpl/ref/h19/1000G_phase1.indels.hg19.sites.vcf --fix_misencoded_quality_scores

###BaseRecalibrator###

java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T BaseRecalibrator \
-I SRR1139966_marked_realigned.bam -R /home/emorycpl/ref/h19/ucsc.hg19.fasta -o SRR1139966_marked_realigned_temp.table -knownSites /home/emorycpl/ref/h19/dbsnp_138.hg19.vcf
samtools index SRR1139966_marked_realigned.bam

###PrintReads###

java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T PrintReads \
-R /home/emorycpl/ref/h19/ucsc.hg19.fasta -I SRR1139966_marked_realigned.bam -o SRR1139966_marked_realigned.recal.bam -BQSR SRR1139966_marked_realigned_temp.table

samtools index SRR1139966_marked_realigned.recal.bam

#####################################################
############## Step 6: gatk call snp/indel##########
#####################################################


###  call variants with GATK's UnifiedGenotyper (faster, lots of raw calls, must run filter after)
time java -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -R /home/emorycpl/ref/h19/ucsc.hg19.fasta \
-T UnifiedGenotyper \
-glm BOTH \
-D /home/emorycpl/ref/h19/dbsnp_138.hg19.vcf \
-metrics SRR1139966_marked_realigned.recal.bam.gatk.snps.metrics \
-stand_call_conf 50.0 -stand_emit_conf 10.0 \
-o SRR1139966_marked_realigned.recal.bam.gatk.vcf \
-I SRR1139966_marked_realigned.recal.bam

​
### 
java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T HaplotypeCaller  \
-R /home/emorycpl/ref/h19/ucsc.hg19.fasta -I SRR1139966_marked_realigned.recal.bam --dbsnp /home/emorycpl/ref/h19/dbsnp_138.hg19.vcf  \
-o  SRR1139966_marked_realigned.recal_raw.snps.indels.vcf

java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T SelectVariants  -R /home/emorycpl/ref/h19/ucsc.hg19.fasta  \
-selectType SNP \
-V SRR1139966_marked_realigned.recal_raw.snps.indels.vcf -o SRR1139966_marked_realigned.recal_raw_snps.vcf
​
java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T SelectVariants  -R /home/emorycpl/ref/h19/ucsc.hg19.fasta  \
-selectType INDEL  \
-V SRR1139966_marked_realigned.recal_raw.snps.indels.vcf -o SRR1139966_marked_realigned.recal_raw_indels.vcf

## for SNP
​
java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/emorycpl/ref/h19/ucsc.hg19.fasta  \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"  \
--filterName "my_snp_filter" \
-V SRR1139966_marked_realigned.recal_raw_snps.vcf  -o SRR1139966_marked_realigned.recal_raw_filtered_snps.vcf   
​
java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/emorycpl/ref/h19/ucsc.hg19.fasta  \
--excludeFiltered \
-V SRR1139966_marked_realigned.recal_raw_filtered_snps.vcf  -o  SRR1139966_marked_realigned.recal_raw_filtered.passed.snps.vcf
​
​java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/emorycpl/ref/h19/ucsc.hg19.fasta  \
-eval SRR1139966_marked_realigned.recal_raw_filtered.passed.snps.vcf -o  SRR1139966_marked_realigned.recal_raw_filtered.passed.snps.vcf.eval
​
​
## for  INDEL
​
​java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/emorycpl/ref/h19/ucsc.hg19.fasta  \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"  \
--filterName "my_indel_filter" \
-V SRR1139966_marked_realigned.recal_raw_indels.vcf  -o SRR1139966_marked_realigned.recal_raw.filtered_indels.vcf   
​
java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/emorycpl/ref/h19/ucsc.hg19.fasta  \
--excludeFiltered \
-V SRR1139966_marked_realigned.recal_raw_indels.vcf  -o  SRR1139966_marked_realigned.recal_raw.PASS.indels.vcf
​
java -Djava.io.tmpdir=/home/emorycpl/tmp -Xmx40g -jar /home/emorycpl/src/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/emorycpl/ref/h19/ucsc.hg19.fasta  \
-eval SRR1139966_marked_realigned.recal_raw_indels.vcf -o  SRR1139966_marked_realigned.recal_raw.indels.vcf.eval
​

 

 
