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
​
###FixMateInfo###

java -Djava.io.tmpdir=/home/emorycpl/tmp  -Xmx40g -jar ~/src/picard/build/libs/picard.jar FixMateInformation \
INPUT=SRR1139966_marked.bam OUTPUT=SRR1139966_marked_fixed.bam SO=coordinate
samtools index SRR1139966_marked_fixed.bam

java -Djava.io.tmpdir=/home/emorycpl/tmp  -Xmx40g -jar ~/src/picard/build/libs/picard.jar FixMateInformation \
INPUT=SRR1139973_marked.bam OUTPUT=SRR1139973_marked_fixed.bam SO=coordinate
samtools index SRR1139966_marked_fixed.bam



