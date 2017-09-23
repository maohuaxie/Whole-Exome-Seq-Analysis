# Reproduction the result for the muttional landscapes of genetic and chemical models of Kras-driven lung cancer
# Pull down the data from ENA database 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR171/ERR171411/ERR171411_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR171/ERR171411/ERR171411_1.fastq.gz

# preparing a reference for use with BWA and GATK
# Generate the BWA index
bwa index ~/ref/m10/m10noMT.fa
# Generate the fastafile index
samtools faidx ~/ref/m10/m10noMT.fa
# Generate the sequence dictionary
java -jar ~/src/picard/build/libs/picard.jar CreateSequenceDictionary REFERENCE=/home/emorycpl/ref/m10/m10noMT.fa OUTPUT=/home/emorycpl/ref/m10/m10noMT.dict
# Mapping sequencing reads against a reference genome
bwa mem -aM -t 12 -R "@RG\tID:Seq01p\tSM:Seq01\tPL:ILLUMINA\tPI:330" ~/ref/m10/m10noMT.fa ERR171411_1.fastq ERR171411_2.fastq > Seq01pairs.sam
# Converting a SAM file to a BAM file and sorting a BAM file by coordinates
samtools view -Sb Seq01pairs.sam > Seq01pairs.bam 
samtools sort Seq01pairs.bam > Seq01pairs.sorted.bam 
# Examining aligned records in the BAM file to locate duplicate reads
java -jar ~/src/picard/build/libs/picard.jar MarkDuplicates INPUT=Seq01pairs.sorted.bam \
MAX_RECORDS_IN_RAM=2000000 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true \
METRICS_FILE=Seq01pairs.sorted.dups OUTPUT=Seq01pairs.sortedDeDup.bam
# You may get the error: Failed to write core dump. Core dumps have been disabled. To enable core dumping, try ulimit -c unlimited before starting Java again
ulimit -c unlimited
ulimit -c -l
# Indexing sorted alignment for fast random access
samtools index Seq01pairs.sortedDeDup.bam
# Performing local realignment around indelsto correct mapping-related artifacts
# Create a target list of intervals to be realigned
java -jar ~/src/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ~/ref/m10/m10.fa \
-I Seq01pairs.sortedDeDup.bam \
-known ~/ref/m10/mm10c67indel.vcf \
-o Seq01pairs.ForIndelRealigner.intervals
# Get the indel.dbsnp file
wget -O ./mm10.INDELS.dbSNP137.vcf.gz ftp://ftp-mouse.sanger.ac.uk/REL-1303-SNPs_Indels-GRCm38/mgp.v3.indels.rsIDdbSNPv137.vcf.gz
gunzip mm10.INDELS.dbSNP137.vcf.gz 
sudo nano indelvcf.py 
# make indelvcf.py file
import re,sys,fileinput
Argument = []
Argument = sys.argv[1:] 
Filepath = Argument[0]
Strain_column = int(Argument[1])
Outpath = Argument[2]
newfile = open(str(Outpath),"w")
for line in fileinput.input([Filepath]):
        if line.startswith("#"):
                if line.startswith("##"):
                        newfile.write(str(line))
                        continue
                else:
                        header = ""
                        header = ("\t".join(line.split("\t")[0:9]))+"\t"+line.split("\t")[Strain_column]+"\n"
                        newfile.write(str(header))

        rowlist = []
        rowlist = line.split("\t")

        genotype = []
        genotype = rowlist[Strain_column].split(":")

        if genotype[-1] == "1":
                if genotype[0] == "1/1" or genotype[0] == "2/2" or genotype[0] == "3/3" or genotype[0] == "4/4" or genotype[0] == "5/5" or genotype[0] == "6/6" or genotype[0] == "7/7":
                        newline = ""
                        newline = "chr"+("\t".join(rowlist[0:9]))+"\t"+rowlist[Strain_column]+"\n"
                        newfile.write(str(newline))

newfile.close()
# running the indelvcf.py code
python indelvcf.py mm10.INDELS.dbSNP137.vcf 16 
# you may will get the error message: MESSAGE: Input files /home/emorycpl/ref/m10/mm10c67indel.vcf and reference have incompatible contigs. Please see https://software.broadinstitute.org/gatk/documentation/article?id=63for more information. Error details: No overlapping contigs found.
# ERROR   /home/emorycpl/ref/m10/mm10c67indel.vcf contigs = [chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chrX]
# ERROR   reference contigs = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 1, 2, 3, 4, 5, 6, 7, 8, 9, MT, X, Y]
# To fix this, download the indel.dbsnp again
wget -O ./mm10.INDELS.dbSNP142.vcf.gz ftp://ftp-mouse.sanger.ac.uk/current_indels/strain_specific_vcfs/C57BL_6NJ.mgp.v5.indels.dbSNP142.normed.vcf.gz
gunzip ./mm10.INDELS.dbSNP142.vcf.gz

java -jar ~/src/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ~/ref/m10/m10noMT.fa \
-I Seq01pairs.sortedDeDup.bam \
-known ~/ref/m10/mm10.INDELS.dbSNP142.sorted.vcf \
-o Seq01pairs.ForIndelRealigner.intervals

java -jar ~/src/picard/build/libs/picard.jar SortVcf \
I=mm10.INDELS.dbSNP142.vcf \
O=mm10.INDELS.dbSNP142.sorted.vcf \
# You may get error as below:
# ERROR   /home/emorycpl/ref/m10/mm10.INDELS.dbSNP142.vcf contigs = [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 3, 4, 5, 6, 7, 8, 9, X, Y]
# ERROR   reference contigs = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 1, 2, 3, 4, 5, 6, 7, 8, 9, MT, X, Y]
# To fix this, I will do the following steps:
wget ftp://ftp.ensembl.org/pub/release-88/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.{1..19}.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-88/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.{X,Y}.fa.gz
gunzip -c Mus_musculus.GRCm38.dna.chromosome.* > m10noMT.fa

# Performing local realignment around indelsto correct mapping-related artifacts
# Performing realignment of the target intervals

java -jar ~/src/gatk/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R ~/ref/m10/m10noMT.fa \
-I Seq01pairs.sortedDeDup.bam \
-known ~/ref/m10/mm10.INDELS.dbSNP142.sorted.vcf \
-targetIntervals Seq01pairs.ForIndelRealigner.intervals \
--filter_bases_not_stored \
-log Seq01pairs.realigned.log \
-o Seq01pairs.realigned.bam
# Recalibrating base quality scores in order to correct sequencing errors and other experimental artifacts
# Download the mm10.dbSNP137.vcf.gz file
wget -O ./mm10.dbSNP137.vcf.gz ftp://ftp-mouse.sanger.ac.uk/REL-1303-SNPs_Indels-GRCm38/mgp.v3.snps.rsIDdbSNPv137.vcf.gz
gunzip ./mm10.dbSNP137.vcf.gz
python dbsnpvcf.py mm10.dbSNP137.vcf 16 mm10.dbSNP137.C57.vcf 
# Sorting mm10.dbSNP137.C57.vcf, Notice, The chromoe order will match the config order. You may change vcf file order by changing the config order as you want. 
java -jar ~/src/picard/build/libs/picard.jar SortVcf \
I=mm10.dbSNP137.C57.vcf \
O=mm10.dbSNP137.C57.sorted.vcf \
# Recalibrating base quality scores in order to correct sequencing errors and other experimental artifacts
java -jar ~/src/gatk/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R ~/ref/m10/m10noMT.fa \
-I Seq01pairs.realigned.bam \
--default_platform ILLUMINA \
--force_platform ILLUMINA \
-knownSites ~/ref/m10/mm10.dbSNP137.C57.sorted.vcf \
-knownSites ~/ref/m10/mm10.INDELS.dbSNP142.sorted.vcf \
-o output.GATKrealigned.recal_data.table

java -jar ~/src/gatk/GenomeAnalysisTK.jar \
-T PrintReads \
-R ~/ref/m10/m10noMT.fa \
-I Seq01pairs.realigned.bam \
-BQSR output.GATKrealigned.recal_data.table \
-log output.BQnewQual.log
–o Seq01pairs.GATKrealigned.Recal.bam

# Generating a plot report to assess the quality of a recalibration
java -jar ~/src/gatk/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R ~/ref/m10/m10noMT.fa \
-I Seq01pairs.GATKrealigned.Recal.bam \
--default_platform ILLUMINA \
--force_platform ILLUMINA \
-knownSites ~/ref/m10/mm10.dbSNP137.C57.sorted.vcf \
-knownSites ~/ref/m10/mm10.INDELS.dbSNP142.sorted.vcf \
-BQSR output.GATKrealigned.recal_data.table \
-log output.BQRecal.After.log
-o output.GATKrealigned.recal_data_after.table

java -jar ~/src/gatk/GenomeAnalysisTK.jar \
-T AnalyzeCovariates \
-R ~/ref/m10/m10noMT.fa \
-before output.GATKrealigned.recal_data.table \
-after output.GATKrealigned.recal_data_after.table \
-plots seq01pairs.plots.pdf
-csv seq01pairs.plots.csv
# Calling SNVs and indelssimultaneously via local de-novo assembly of haplotypes
java -jar ~/src/gatk/GenomeAnalysisTK.jar \
-R ~/ref/m10/m10noMT.fa \
-T HaplotypeCaller \
-I Seq01pairs.GATKrealigned.Recal.bam \
--dbsnp ~/ref/m10/mm10.dbSNP137.C57.sorted.vcf \
-stand_call_conf 30 \
-L Seq01pairs.ForIndelRealigner.intervals \
-o Seq01pairs.raw.snps.indels.vcf

# Tips:-stand_call_conf: Qualscore at which to call the variant
# -stand_emit_conf: Qualscore at which to emit the variant as filtered
# -minPruning: Amount of pruning to do in the deBruijngraph
# Raw variant files are often very large and full of false positive variant calls
# Calling SNVs and indelssimultaneously via a Bayesian genotype likelihood model
java -jar ~/src/gatk/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R ~/ref/m10/m10noMT.fa \
-I Seq01pairs.realigned.bam \
--dbsnp ~/ref/m10/mm10.dbSNP137.C57.sorted.vcf \
-stand_call_conf 30 \
-o Seq01pairs.gatk.UG.vcf

# Bam2mpipeup
samtools mpileup -f ~/ref/m10/m10noMT.fa \
Seq01pairs.realigned.bam > Seq01pairs.realigned.bam.mpileup.txt
# Samtools mpileup to bcftools to call variants in VCF file
samtools mpileup -C 0 -A -B -d 10000 -v -u -f ~/ref/m10/m10noMT.fa \
Seq01pairs.realigned.bam | bcftools call -O v -v -c -n 0.05 -p 1 -A -o Seq01pairs.bcftools.vcf
# SNP-calling by freebayes
freebayes -f ~/ref/m10/m10noMT.fa Seq01pairs.realigned.bam \
1> Seq01pairs.bayes.vcf 2> Seq01pairs.bayes.vcf.log
# SNP-calling by varscan
java -jar ~/src/Varscan/VarScan.jar mpileup2snp \
~/project/bcl2/Seq01pairs.realigned.bam.mpileup.txt --output-vcf 1 >~/project/bcl2/Seq01pairs.varscan.vcf

java -jar ~/src/Varscan/VarScan.jar mpileup2indel \
~/project/bcl2/Seq01pairs.realigned.bam.mpileup.txt --output-vcf 1 >~/project/bcl2/Seq01pairs.varscan.indel.vcf
# SNP-calling by unified Genotyper
java -jar ~/src/gatk/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R ~/ref/m10/m10noMT.fa \
-I Seq01pairs.realigned.bam \
--dbsnp ~/ref/m10/mm10.dbSNP137.C57.sorted.vcf \
-stand_call_conf 30 \
-o Seq01pairs.gatk.UG.vcf

# filter the snp
perl -alne '{next if $F[5]<20;/DP=(\d+)/;next if $1<10;next if /INDEL/;/(DP4=.*?);/;print "$F[0]\t$F[1]\t$F[3]\t$F[4]:$1"}' Seq01pairs.bcftools.vcf >Seq01pairs.bcftools.vcf.filter
perl -alne '{next if $F[5]<20;/DP=(\d+)/;next if $1<10;next unless /TYPE=snp/;@tmp=split/:/,$F[9];print "$F[0]\t$F[1]\t$F[3]\t$F[4]:$tmp[0]:$tmp[1]"}'  Seq01pairs.bayes.vcf >  Seq01pairs.bayes.vcf.filter
perl -alne '{@tmp=split/:/,$F[9];next if $tmp[3]<10;print "$F[0]\t$F[1]\t$F[3]\t$F[4]:$tmp[0]:$tmp[3]"}' Seq01pairs.varscan.vcf >Seq01pairs.varscan.vcf.filter
perl -alne '{next if $F[5]<20;/DP=(\d+)/;next if $1<10;next if length($F[3]) >1;next if length($F[4]) >1;@tmp=split/:/,$F[9];print "$F[0]\t$F[1]\t$F[3]\t$F[4]:$tmp[0]:$tmp[1]:$tmp[2]"}'  Seq01pairs.gatk.UG.vcf >Seq01pairs.gatk.UG.vcf.filter
# annotation snp with Snpef 
java -jar ~/src/snpEff/snpEff.jar  mm10 Seq01pairs.varscan.vcf >Seq01pairs.varscan.snps.vcf

# contents in ~/project/bcl2 folder
emorycpl@emorycpl-Precision-Tower-5810:~/project/bcl2$ ls -al
total 89391312
drwxrwxr-x 2 emorycpl emorycpl       12288 Sep  3 14:44 .
drwxrwxr-x 4 emorycpl emorycpl        4096 Aug 27 18:26 ..
-rw-r--r-- 1 root     root             709 Sep  1 09:34 clean_swap.sh
-rw------- 1 emorycpl emorycpl 22381461504 Sep  1 09:29 core
-rw-rw-r-- 1 emorycpl emorycpl  5204085380 Aug 31 22:01 ERR171411_1.fastq
-rw-rw-r-- 1 emorycpl emorycpl  5204085380 Aug 31 22:13 ERR171411_2.fastq
-rw-rw-r-- 1 emorycpl emorycpl 10862717096 Aug 27 17:48 index17_GTAGAG_L001-L002_R1_001.fastq
-rw-rw-r-- 1 emorycpl emorycpl 10862717096 Aug 27 17:49 index17_GTAGAG_L001-L002_R2_001.fastq
-rw-rw-r-- 1 emorycpl emorycpl       44142 Sep  2 20:06 output.BQnewQual.log
-rw-rw-r-- 1 emorycpl emorycpl      669214 Sep  2 13:40 output.GATKrealigned.recal_data.table
-rw-rw-r-- 1 emorycpl emorycpl  4366905552 Sep  1 21:31 Seq01pairs.bam
-rw-rw-r-- 1 emorycpl emorycpl   110855642 Sep  3 11:13 Seq01pairs.bayes.vcf
-rw-rw-r-- 1 emorycpl emorycpl     2185895 Sep  3 14:41 Seq01pairs.bayes.vcf.filter
-rw-rw-r-- 1 emorycpl emorycpl           0 Sep  3 09:33 Seq01pairs.bayes.vcf.log
-rw-rw-r-- 1 emorycpl emorycpl   295519480 Sep  2 23:35 Seq01pairs.bcftools.vcf
-rw-rw-r-- 1 emorycpl emorycpl     3367271 Sep  3 14:41 Seq01pairs.bcftools.vcf.filter
-rw-rw-r-- 1 emorycpl emorycpl     2813529 Sep  2 08:59 Seq01pairs.ForIndelRealigner.intervals
-rw-rw-r-- 1 emorycpl emorycpl    45960433 Sep  3 12:46 Seq01pairs.gatk.UG.vcf
-rw-rw-r-- 1 emorycpl emorycpl     3448956 Sep  3 14:44 Seq01pairs.gatk.UG.vcf.filter
-rw-rw-r-- 1 emorycpl emorycpl     5131944 Sep  3 12:46 Seq01pairs.gatk.UG.vcf.idx
-rw-rw-r-- 1 emorycpl emorycpl     5394184 Sep  2 10:02 Seq01pairs.realigned.bai
-rw-rw-r-- 1 emorycpl emorycpl  3159684628 Sep  2 10:02 Seq01pairs.realigned.bam
-rw-rw-r-- 1 emorycpl emorycpl  9663595470 Sep  2 22:48 Seq01pairs.realigned.bam.mpileup.txt
-rw-rw-r-- 1 emorycpl emorycpl        6442 Sep  2 10:02 Seq01pairs.realigned.log
-rw-rw-r-- 1 emorycpl emorycpl 13006707560 Sep  1 21:13 Seq01pairs.sam
-rw-rw-r-- 1 emorycpl emorycpl  3117393730 Sep  1 21:53 Seq01pairs.sorted.bam
-rw-rw-r-- 1 emorycpl emorycpl  3206933608 Sep  1 22:08 Seq01pairs.sortedDeDup.bam
-rw-rw-r-- 1 emorycpl emorycpl     3911768 Sep  1 22:18 Seq01pairs.sortedDeDup.bam.bai
-rw-rw-r-- 1 emorycpl emorycpl        2811 Sep  1 22:08 Seq01pairs.sorted.dups
-rw-rw-r-- 1 emorycpl emorycpl     1548658 Sep  3 13:58 Seq01pairs.varscan.indel.vcf
-rw-rw-r-- 1 emorycpl emorycpl    17178094 Sep  3 13:29 Seq01pairs.varscan.vcf
-rw-rw-r-- 1 emorycpl emorycpl     2248817 Sep  3 14:37 Seq01pairs.varscan.vcf.filter

# contents in ~/ref/m10 folder
emorycpl@emorycpl-Precision-Tower-5810:~/ref/m10$ ls -al
total 61319344
drwxrwxr-x 3 emorycpl emorycpl        4096 Sep  2 11:31 .
drwxrwxr-x 4 emorycpl emorycpl        4096 Aug 27 21:55 ..
drwxrwxr-x 3 emorycpl emorycpl        4096 Aug 30 22:28 data
-rw-r--r-- 1 root     root            1145 Sep  2 10:33 dbsnpvcf.py
-rw-r--r-- 1 root     root            1147 Sep  1 17:15 indelvcf.py
-rw-rw-r-- 1 emorycpl emorycpl        2125 Aug 29 15:10 m10.dict
-rw-rw-r-- 1 emorycpl emorycpl  2770964555 Aug 27 22:07 m10.fa
-rw-rw-r-- 1 emorycpl emorycpl        9189 Aug 28 00:18 m10.fa.amb
-rw-rw-r-- 1 emorycpl emorycpl        1779 Aug 28 00:18 m10.fa.ann
-rw-rw-r-- 1 emorycpl emorycpl  2725537772 Aug 28 00:18 m10.fa.bwt
-rw-rw-r-- 1 emorycpl emorycpl         624 Aug 29 14:51 m10.fa.fai
-rw-rw-r-- 1 emorycpl emorycpl   681384419 Aug 28 00:18 m10.fa.pac
-rw-rw-r-- 1 emorycpl emorycpl  1362768888 Aug 28 00:29 m10.fa.sa
-rw-rw-r-- 1 emorycpl emorycpl   930312808 Aug 27 22:25 m10.gtf
-rw-rw-r-- 1 emorycpl emorycpl        2116 Sep  1 20:56 m10noMT.dict
-rw-rw-r-- 1 emorycpl emorycpl  2770947930 Sep  1 19:55 m10noMT.fa
-rw-rw-r-- 1 emorycpl emorycpl        9189 Sep  1 20:38 m10noMT.fa.amb
-rw-rw-r-- 1 emorycpl emorycpl        1705 Sep  1 20:38 m10noMT.fa.ann
-rw-rw-r-- 1 emorycpl emorycpl  2725521464 Sep  1 20:37 m10noMT.fa.bwt
-rw-rw-r-- 1 emorycpl emorycpl         598 Sep  1 20:54 m10noMT.fa.fai
-rw-rw-r-- 1 emorycpl emorycpl   681380344 Sep  1 20:38 m10noMT.fa.pac
-rw-rw-r-- 1 emorycpl emorycpl  1362760736 Sep  1 20:48 m10noMT.fa.sa
-rw-rw-r-- 1 emorycpl emorycpl     9488149 Sep  1 17:18 mm10c67indel.vcf
-rw-rw-r-- 1 emorycpl emorycpl     5041001 Sep  2 11:39 mm10.dbSNP137.C57.sorted.vcf
-rw-rw-r-- 1 emorycpl emorycpl      118177 Sep  2 11:39 mm10.dbSNP137.C57.sorted.vcf.idx
-rw-rw-r-- 1 emorycpl emorycpl     5045727 Sep  2 11:31 mm10.dbSNP137.C57.vcf
-rw-rw-r-- 1 emorycpl emorycpl 42196785097 Sep  2 10:46 mm10.dbSNP137.vcf
-rw-r--r-- 1 root     root            1024 Sep  2 11:26 .mm10.dbSNP137.vcf.swp
-rw-rw-r-- 1 emorycpl emorycpl  4508269676 Sep  1 17:06 mm10.INDELS.dbSNP137.vcf
-rw-rw-r-- 1 emorycpl emorycpl     8839422 Sep  1 23:31 mm10.INDELS.dbSNP142.sorted.vcf
-rw-rw-r-- 1 emorycpl emorycpl       46466 Sep  1 23:31 mm10.INDELS.dbSNP142.sorted.vcf.idx
-rw-rw-r-- 1 emorycpl emorycpl     8853970 Sep  1 23:23 mm10.INDELS.dbSNP142.vcf
-rw-rw-r-- 1 emorycpl emorycpl    36782152 Aug 30 22:19 snpEff_v4_3_mm10.zip
