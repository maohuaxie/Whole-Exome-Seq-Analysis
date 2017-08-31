# Whole-Exome-Seq-Analysis
# Before we are learning how to perform bionformatic data analysis, we need to set up the ubuntu 16.04 LTS linux system.
# one step I want to point out is that it is much better for the linux beginner to partition the disk as below.
#For simplicity, we partition 3 disks:
/ # one is for boot, if you have enough disk space, please allocate more than 100 GB for boot.
swap # one is for swap, basically, It is allocated twice size with your RAM, if you RAM is 32 GB, please allocate 64 GB for swap
/data # one is for data, all the left disk space are allocated to data disk. 

# After installtion of your ubuntu operation system, you may need to set up static IP address:
# First thing you need to do is to enable SSH in ubuntu 16.04.
sudo apt-get install openssh-server
sudo vim /etc/ssh/sshd_config   // you may need install vim, to do this, please run sudo apt-get install vim or you can run sudo nano /etc/ssh/sshd_config
change Permit RootLogin to yes
# then go to /etc/network/interfaces folder to set up the static IP address
sudo vim /etc/network/interfaces or sudo nano /etc/network/interfaces

# interfaces(5) file used by ifup(8) and ifdown(8)
auto lo
iface lo inet loopback
# change to 

auto eth0
iface eth0 inet static
address 192.168.71.160
gateway 192.168.71.2 // you can get the infor from the properties of your network
netmask 255.255.255.0
dns-nameservers 8.8.8.8

#and then, go to /etc/NetworkManager/NetworkManager.conf folder
sudo vim /etc/NetworkManager/NetworkManager.conf  or sudo nano /etc/NetworkManager/NetworkManager.conf
[if updown] managed = false  // change false to true
 
#Now, the Ubuntu is ready to install bioinformtics tools
#Considering some tools running in python2.7, we are going to install tools via anaconda2. 
#Download anaconda2.sh file
bash anacoda2.sh //please note, it may not exactly the same.
#and then export the path
export PATH=/home/emorycpl/anaconda2/bin:$PATH >> ~/.bashrc
source ~/.bashrc
# There is an easy way to install the bioinformtics tools, you may reference with the biostarhandbook, and give the credit to the authors. 
# Now download and run the doctor.py script from a terminal:
mkdir -p ~/bin
curl http://data.biostarhandbook.com/install/doctor.py > ~/bin/doctor.py
chmod +x ~/bin/doctor.py
# Run it from the terminal and 'doctor' will tell you what (if anything) is out of order:
~/bin/doctor.py
# Later, at any time, when you have a problem run the doctor.py command again to check your settings.
# The doctor.py can also help your system "get better":
doctor.py --fixme
The following steps need to do before we are going to install the bioinformatic tools. 
1. Update the Linux system.
2. Install the required Linux libraries and Java.
3. Install the conda package manager.
sudo apt-get update && sudo apt-get upgrade -y
sudo apt-get install -y curl unzip build-essential ncurses-dev
sudo apt-get install -y byacc zlib1g-dev python-dev git cmake
sudo apt-get install -y python-pip libhtml-parser-perl libwww-perl
sudo apt-get install -y default-jdk ant
# you also can run the script as below. It will install most of the libraries which are required for bioinformatic data analysis. 
curl http://data.biostarhandbook.com/install/aptget.txt | xargs sudo apt-get -y install
#Then, add r, conda-forge and bioconda channels 
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels bioconda
# After you run the following script. It will insall most of the tools for bioinformatic data analysis. 
# and you will be noticed. 
# Doctor! Doctor! Give me the news.
# Checking 13 symptoms...
# Optional program not found: wonderdump
# Optional program not found: global-align.sh
# Optional program not found: local-align.sh
# You are doing well!
curl http://data.biostarhandbook.com/install/conda.txt | xargs conda install -y
export PATH=$PATH:$HOME/bin
# one more program you need to install is R, How to install R in linux?
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 
sudo add-apt-repository ppa:marutter/rdev
sudo apt-get update
sudo apt-get upgrade 
sudo apt-get install r-base r-base-dev 
# if you get error, you can fix the problem by doing:
sudo rm -f /etc/apt/sources.list.d/*
sudo dpkg -i /var/cache/apt/archives/*.deb
sudo dpkg --configure -a

# Install Blast
cd ~/src
curl -O ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.5.0+-x64-linu
x.tar.gz
tar zxvf ncbi-blast-2.5.0+-x64-linux.tar.gz
export PATH=~/src/ncbi-blast-2.5.0+/bin:$PATH
# See the programs that are available
ls ~/src/ncbi-blast-2.5.0+/bin/

# Install GATK
mkdir -p ~/src/gatk
mv ~/Downloads/GenomeAnalysisTK-3.4-46.tar.bz2 ~/src/gatk
cd ~/src/gatk
tar jxvf GenomeAnalysisTK-3.4-46.tar.bz2

# Test the installation
java -jar ~/src/gatk/GenomeAnalysisTK.jar -h

# Create a script that launches gatk :
echo '#!/bin/bash' > ~/bin/gatk
echo 'java -jar ~/src/gatk/GenomeAnalysisTK.jar $@' >> ~/bin/gatk
# Make the script executable
chmod +x ~/bin/gatk

# Install Varscan
mkdir -p ~/src/VarScan
cd ~/src/VarScan
# Test the installation
java -jar ~/src/Varscan/VarScan.v2.3.9.jar 
# Create a script that launches Varscan :
echo '#!/bin/bash' > ~/bin/Varscan
echo 'java -jar ~/src/Varscan/VarScan.v2.3.9.jar $@' >> ~/bin/Varscan
# Make the script executable
chmod +x ~/bin/Varscan

# Install cufflinks
cd ~/src
curl -kL http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar xzv
# Enable the path.
echo 'export PATH=~/src/cufflinks-2.2.1.Linux_x86_64:$PATH' >> ~/.bash_profile
source ~/.bash_profile

# Install hisat2
# Find and download the Linux version (also applies as an alternative on Mac OSX)
cd ~/src
curl -OL ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.3-beta-Linux_x86_64.zip
unzip hisat2-2.0.3-beta-Linux_x86_64.zip
# The Linux link the executables:
ln -s ~/src/hisat2-2.0.3-beta/hisat2-build ~/bin
ln -s ~/src/hisat2-2.0.3-beta/hisat2 ~/bin
# Test that the programs work with
hisat2

# Install picard 
cd ~/src
# Obtain the picard distribution.
git clone https://github.com/broadinstitute/picard.git
cd picard/
./gradlew shadowJar
# Test that the installation succeeded:
java -jar ~/src/picard/build/libs/picard.jar
Create a script that launches picard:
echo '#!/bin/bash' > ~/bin/picard
echo 'java -jar ~/src/picard/build/libs/picard.jar $@' >> ~/bin/picard
# Make the script executable
chmod +x ~/bin/picard
# Test that the script works with:
picard

# This is helpful when one wants to understand what type of files come with fastqc
cd ~/src
curl -O http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
# Link the fastqc executable to the ~/bin folder that
# you have already added to the path.
ln -sf ~/src/FastQC/fastqc ~/bin/fastqc
# Due to what seems a packaging error
# the executable flag on the fastqc program is not set.
# We need to set it ourselves.
chmod +x ~/bin/fastqc
# Test installation by running:
fastqc -h

# Install trimmomatic 
cd ~/src
curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip
# The program can be invoked via
java -jar ~/src/Trimmomatic-0.36/trimmomatic-0.36.jar
# The ~/src/Trimmomatic-0.36/adapters/ directory contains
# Illumina specific adapter sequences.
ls ~/src/Trimmomatic-0.36/adapters/
# How to run trimmomatic
# Unfortunately running trimmomatic is as user unfriendly as it gets. To run it we "simply" type:
java -jar ~/src/Trimmomatic-0.36/trimmomatic-0.36.jar
# That gets old very quickly. To simplify the invocation create a script in the ~/bin folder:
echo '#!/bin/bash' > ~/bin/trimmomatic
echo 'java -jar ~/src/Trimmo
# test installation by running
trimmomatic

# Install tophat
cd ~/src
curl -kO https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
tar xzvf tophat-2.1.1.Linux_x86_64.tar.gz
# This is a suite of tools thus is best to
# add the entire directory to the search path.
export PATH=~/src/tophat-2.1.1.Linux_x86_64:$PATH
# See the progams included with TopHat.
ls ~/src/tophat-2.1.1.Linux_x86_64
# test installation by running
tophat -h

# Install SubRead
# This package contains an aligner and feature counter. The programs that will use are:
*subread-align - short read aligner
*featureCounts - feature counter
Website: http://bioinf.wehi.edu.au/subread-package/
# For all other platforms Obtain the source code:
mkdir -p ~/src ~/bin
cd ~/src
curl -OkL http://sourceforge.net/projects/subread/files/subread-1.5.1/subread-1.5.1-source.tar.gz
tar zxvf subread-1.5.1-source.tar.gz
cd subread-1.5.1-source/src
# Compile on Linux
make -f Makefile.Linux
# Link the executables into the ~bin folder:
ln -fs ~/src/subread-1.5.1-source/bin/featureCounts ~/bin
ln -fs ~/src/subread-1.5.1-source/bin/subread-align ~/bin
ln -fs ~/src/subread-1.5.1-source/bin/subread-buildindex ~/bin
# Check that the programs work
featureCounts
subread-align

# Install kallisto
# Make the source directory if they are not exist.
mkdir -p ~/src ~/bin
# The URL to the download.
cd ~/src
curl -kLO https://github.com/pachterlab/kallisto/releases/download/v0.43.0/kallisto_linux-v0.43.0.tar.gz | tar zvz -C ~/src
# Link kallisto to the path
ln -s kallisto_linux-v0.43.0/kallisto ~/bin
# Install snpEff
# A separate directory for downloads, source and executables
mkdir -p ~/down ~/src ~/bin
# Download the source.
curl -kL http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip > ~/down/snpEff_latest_core.zip
# Unpack into the source directory
unzip ~/down/snpEff_latest_core.zip -d ~/src
Create a starter script to snpEff
echo '#!/bin/bash' > ~/bin/snpEff
echo 'java -jar ~/src/snpEff/snpEff.jar $@' >> ~/bin/snpEff
chmod +x ~/bin/snpEff

***Bioinformatic project for job searching***
# Refer to the publication titled "Mapping RNA-seq Reads with STAR"
# As to install star software, if you have anaconda installed in your os, you can install it as below:
conda install star //please be note that you have to put "export PATH=/home/emorycpl/anaconda2/bin:$PATH" on your .bashrc file in advance, otherwise you will get error message.
# create new folder as project, Download project data and unzip it:
mkdir ~/project
# now we need download the most recent GRcm38 file.
# please go to your home directory, here my home directory is "/home/emorycpl", and run the script as below.
mkdir -p ~/ref/m10
cd ~/ref/m10
wget ftp://ftp.ensembl.org/pub/release-88/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.{1..19}.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-88/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.{MT,X,Y}.fa.gz
gunzip -c Mus_musculus.GRCm38.dna.chromosome.* > m10.fa

# then we need download the annotation file:
cd ~/ref/m10
wget ftp://ftp.ensembl.org/pub/release-88/gtf/mus_musculus/Mus_musculus.GRCm38.88.gtf.gz
gunzip -c Mus_musculus.GRCm38.88.gtf.gz > m10.gtf

# Until now, we have all the staff for generating the star indice. All right, let's do it.
index15_ATGTCA_L001-L002_R1_001.fastq-009.gz  index17_GTAGAG_L001-L002_R2_001.fastq.gz      index19_GTGAAA_L001-L002_R1_001.fastq-005.gz  index20_GTGGCC_L001-L002_R2_001.fastq.gz  
index15_ATGTCA_L001-L002_R2_001.fastq-008.gz  index18_GTCCGC_L001-L002_R1_001.fastq.gz      index19_GTGAAA_L001-L002_R2_001.fastq.gz                                    
index17_GTAGAG_L001-L002_R1_001.fastq.gz      index18_GTCCGC_L001-L002_R2_001.fastq-007.gz  index20_GTGGCC_L001-L002_R1_001.fastq-007.gz  
# Align the paired reads to reference genome using bwa mem.
time bwa index /home/emorycpl/ref/m10/m10.fa 
time bwa mem -t 12 -M /home/emorycpl/ref/m10/m10.fa \
/home/emorycpl/project/bcl2/index15_ATGTCA_L001-L002_R1_001.fastq-009.gz \
/home/emorycpl/project/bcl2/index15_ATGTCA_L001-L002_R2_001.fastq-008.gz \
1>/home/emorycpl/project/bcl2/Sample15.sam 2>Sample15.bwa.log 

nohup time bwa mem -t 12 -M /home/emorycpl/ref/m10/m10.fa \
/home/emorycpl/project/bcl2/index17_GTAGAG_L001-L002_R1_001.fastq.gz \
/home/emorycpl/project/bcl2/index17_GTAGAG_L001-L002_R2_001.fastq.gz \
1>/home/emorycpl/project/bcl2/Sample17.sam 2>Sample17.bwa.log &

nohup time bwa mem -t 12 -M /home/emorycpl/ref/m10/m10.fa \
/home/emorycpl/project/bcl2/index18_GTCCGC_L001-L002_R1_001.fastq.gz \
/home/emorycpl/project/bcl2/index18_GTCCGC_L001-L002_R2_001.fastq-007.gz \
1> /home/emorycpl/project/bcl2/Sample18.sam 2>Sample18.bwa.log &

nohup time bwa mem -t 12 -M /home/emorycpl/ref/m10/m10.fa \
/home/emorycpl/project/bcl2/index19_GTGAAA_L001-L002_R1_001.fastq-005.gz \
/home/emorycpl/project/bcl2/index19_GTGAAA_L001-L002_R2_001.fastq.gz \
1> /home/emorycpl/project/bcl2/Sample19.sam 2>Sample19.bwa.log &

nohup time bwa mem -t 12 -M /home/emorycpl/ref/m10/m10.fa \
/home/emorycpl/project/bcl2/index20_GTGGCC_L001-L002_R1_001.fastq-007.gz \
/home/emorycpl/project/bcl2/index20_GTGGCC_L001-L002_R2_001.fastq.gz \
1> /home/emorycpl/project/bcl2/Sample20.sam 2>Sample20.bwa.log 

# Now, we can use the samtools view command to convert the BAM to SAM so we mere mortals can read it.
reference=/home/emorycpl/ref/m10/m10.fa
for i in *sam 
do 
echo $i
echo ${i%.*}.sorted.bam
samtools view -Sb $i > ${i%.*}.sorted.bam >${i%.*}.AddOrReplaceReadGroups.log 
done
cat Sample15.sam | head
cat Sample15.sam | grep -v @ | cut -f 1 | head -5
samtools view Sample15.sam | cut -f 1 | head -5
samtools view Sample15.sam | cut -f 1,3,4 | head -5
samtools view Sample15.sam | cut -f 1,3,4,5,6 | head -5
samtools view Sample15.sam | cut -f 1,2,4,6,7,8,9 | head -1
samtools view -b Sample15.sam | bedtools bamtobed | head -1
samtools view Sample15.sam | cut -f 10,11 | head -1
samtools view Sample15.sam | cut -f 12,13,14,15 | head -1
samtools view Sample15.bam | head -5
# sort the bam file
samtools sort Sample15.bam -o Sample15.sorted.bam
samtools view Sample15.sorted.bam | head -5
# Index the bam file
samtools index Sample15.sorted.bam
samtools view Sample15.sorted.bam 1:33000000-34000000 | wc -l
samtools view -c -f 4 Sample15.sorted.bam //unmatch
samtools view -c -F 4 Sample15.sorted.bam //match
# How to get an overview of the alignments in a BAM file?
samtools flagstat Sample15.sorted.bam
samtools idxstats Sample15.sorted.bam
bamtools stats -in Sample15.sorted.bam
# How can I find out the depth of coverage?
samtools depth Sample15.sorted.bam | head
bamtools coverage -in Sample15.sorted.bam | head
# To find the most covered reagion
samtools depth Sample15.sorted.bam | sort -k 3 -rn | head
# What is a "SECONDARY" alignment?
# A secondary alignment refers to a read that produces multiple alignments in the genome
samtools flags SECONDARY
samtools view -c -F 4 -f 256 Sample15.sorted.bam
# Picard re-sort as preparation step for marking duplicates
java -jar ~/src/picard/build/libs/picard.jar MarkDuplicates \
I=Sample15.sorted.bam \
O=Sample15_marked_duplicates.bam \
M=marked_dup_metrics.txt
# Create a target list of intervals to be realigned
wget -O ./mm10.INDELS.dbSNP142.vcf.gz ftp://ftp-mouse.sanger.ac.uk/current_indels/strain_specific_vcfs/C57BL_6NJ.mgp.v5.indels.dbSNP142.normed.vcf.gz
gunzip ./mm10.INDELS.dbSNP142.vcf.gz
awk '{ if($0 !~ /^#/) print "chr"$0; else if(match($0,/(##contig=<ID=)(.*)/,m)) print m[1]"chr"m[2]; else print $0 }' mm10.INDELS.dbSNP142.vcf > mm10.INDELS.vcf

java -jar ~/src/gatk/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /home/emorycpl/ref/m10/m10.fa \
-I Sample15_marked_duplicates.bam \
-known mm10.INDELS.dbSNP142.vcf \
-o forIndelRealigner.intervals

java -jar ~/src/picard/build/libs/picard.jar CreateSequenceDictionary REFERENCE=/home/emorycpl/ref/m10/m10.fa OUTPUT=/home/emorycpl/ref/m10/m10.dict

# bam2mpileup
Samtools mpileup -f /home/emorycpl/ref/m10/m10.fa Sample15_marked_duplicates.bam 1>Sample15.mpileup.txt 2>Sample15.mpileup.log
# Samtools mpileup to bcftools to call variants in VCF file
samtools mpileup -C 0 -A -B -d 10000 -v -u -f /home/emorycpl/ref/m10/m10.fa \
Sample15_marked_duplicates.bam | bcftools call -O v -v -c -n 0.05 -p 1 -A -o Sample15.vcf
# SNP-calling by freebayes
freebayes -f /home/emorycpl/ref/m10/m10.fa Sample15_marked_duplicates.bam \
1> Sample15.bayes.vcf 2> Sample15.bayes.vcf.log
# SNP-calling by varscan
java -jar ~/src/Varscan/VarScan.jar mpileup2snp \
~/project/bcl2/Sample15.mpileup.txt --output-vcf 1 >~/project/bcl2/Sample15.varscan.vcf

java -jar ~/src/Varscan/VarScan.jar mpileup2indel \
~/project/bcl2/Sample15.mpileup.txt --output-vcf 1 >~/project/bcl2/Sample15.varscan.indel.vcf
# filter the snp
perl -alne '{next if $F[5]<20;/DP=(\d+)/;next if $1<10;next if /INDEL/;/(DP4=.*?);/;print "$F[0]\t$F[1]\t$F[3]\t$F[4]:$1"}' Sample15.vcf >Sample15.bcftools.vcf.filter
perl -alne '{next if $F[5]<20;/DP=(\d+)/;next if $1<10;next unless /TYPE=snp/;@tmp=split/:/,$F[9];print "$F[0]\t$F[1]\t$F[3]\t$F[4]:$tmp[0]:$tmp[1]"}'  Sample15.bayes.vcf > Sample15.bayes.vcf.filter
perl -alne '{@tmp=split/:/,$F[9];next if $tmp[3]<10;print "$F[0]\t$F[1]\t$F[3]\t$F[4]:$tmp[0]:$tmp[3]"}' Sample15.varscan.vcf >Sample15.varscan.vcf.filter
# annotation snp with Annovar
./convert2annovar.pl -format vcf4 --includeinfo -coverage 20 -fraction 0.05 ~/project/bcl2/Sample15.varscan.vcf > ~/project/bcl2/Tumor_cleaned.avinput
./table_annovar.pl ~/project/bcl2/Tumor_cleaned.avinput mousedb/ -buildver mm10 -out Tumor_cleaned -remove -protocol nonsyn_splicing,genomicSuperDups,SC_MOUSE_GENOMES.genotype.vcf,dominant -operation g,rr,f,m -nastring NA
# annotation snp with Snpef 
java -jar ~/src/snpEff/snpEff.jar  mm10 Sample15.varscan.vcf >Sample15.varscan.snps.vcf
