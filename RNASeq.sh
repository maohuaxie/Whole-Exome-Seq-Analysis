mkdir -p Rseq/star/GRCh38/star_indices_overhang100
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ~/Rseq/star/GRCh38/star_indices_overhang100 --genomeFastaFiles ~/Rseq/star/GRCh38/sequence/GRCh38_r88.all.fa --sjdbGTFfile ~/Rseq/star/GRCh38/annotation/Homo_sapiens.GRCh38.88.gtf --sjdbOverhang 100
# Please be note that the splice-junction-data-base-overhang parameter should have a value of read length – 1. 
# In this case, we have set it to 100, because our read length is 101. 
# Furthermore, the runThreadN parameter determines how many threads you want to ask for. E.g.  Dell T5610 with a 2.6 GHz Intel xeon e5 2640 processor. 
# has 16 “cores” (16 “real” and 16“virtual”) via Hyper Threading.

STAR --runThreadN 16 \
--genomeDir ~/Rseq/star/GRCh38/star_indices_overhang100 \
--readFilesIn ~/Rseq/fastq/SRR1039508_1.fastq ~/Rseq/fastq/SRR1039508_2.fastq

Generating the genome:
Note the sjdbOverhang is used for constructing the splice junction database. It should be set to (read length - 1), and according to the manual a general value of 100 will work as well.
For this limited demonstration, I am only going to align to the genes on chromosome 1, so I subset the GTF file:
grep -P '^1\t' Homo_sapiens.GRCh38.79.gtf > Homo_sapiens.GRCh38.79.chrom1.gtf
We then moved files in subdirectories, and created one for the STAR genome index:
mkdir gtf
mkdir genome
mv *.gtf gtf
mv *.fa genome
mkdir GRCh38.79.chrom1
The STAR command to generate the genome index:
STAR --runMode genomeGenerate \
--genomeDir GRCh38.79.chrom1 \
--genomeFastaFiles genome/Homo_sapiens.GRCh38.dna.chromosome.1.fa \
--sjdbGTFfile gtf/Homo_sapiens.GRCh38.79.chrom1.gtf \
--sjdbOverhang 62 
Mapping the reads:
STAR --runThreadN 12 \
--genomeDir GRCh38.79.chrom1 \
--readFilesIn fastq/SRR1039508_1.fastq fastq/SRR1039508_2.fastq
Downloading from the ENA:
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_2.fastq.gz
Alias for ls:
alias ll='ls -lGh'
Unzipping:
gunzip *.fastq.gz
Looking at the FASTQ files:
less SRR1039508_1.fastq
wc -l SRR1039508_1.fastq
Quality control with fastqc
head -n 8 SRR1039508_1.fastq
head -n 8 SRR1039508_2.fastq
# I'm going to call fastqc the first file, the second file, and I'll say no extract. 
# and I'll say time just so we see how long this takes.
time fastqc --noextract SRR1039508_1.fastq SRR1039508_2.fastq 
RSEM expects a GTF file with only exons, which are each assigned to a transcript_id.
awk '$3 == "exon"' gtf/Homo_sapiens.GRCh38.88.gtf > gtf/Homo_sapiens.GRCh38.88.exons.gtf
RSEM will then prepare a reference transcriptome against which to align reads.
mkdir rsemGenome
time rsem-prepare-reference --gtf gtf/Homo_sapiens.GRCh38.88.exons.gtf genome/GRCh38_r88.all.fa rsemGenome/GRCh38.88.all
time bowtie-build -f rsemGenome/GRCh38.88.all.idx.fa rsemGenome/GRCh38.88.all
# please note that you may can try bowtie2-build -f ...
time rsem-calculate-expression -p 16 --paired-end fastq/SRR1039508_1.fastq fastq/SRR1039508_2.fastq rsemGenome/GRCh38.88.all SRR1039508
Once you have installed R, running the following lines in your console will install Bioconductor:
source("http://bioconductor.org/biocLite.R")
biocLite()
Make sure to hit [a] to update all packages. This is important so that your answers will match the answers accepted by the grading bot.
To install specific packages from Bioconductor use, for example:
biocLite(c("pasilla","DEXSeq"))
We will provide a list of all packages we will use here.
If you want to see what version of Biocondutor you are using and whether your packages are up to date:
library(BiocInstaller)
biocVersion()
biocValid()
library(knitr)
opts_chunk$set(fig.path=paste0("figure/", sub("(.*).Rmd","\\1",basename(knitr:::knit_concord$get('infile'))), "-"))
The DEXSeq package offers differential testing of exon usage within each gene. Here we will explore the R code used in a DEXSeq analysis. We omit the python calls for preparing the annotation and count tables, but these can be found in the vignette at the above link. The python calls are generally along the lines of:
python dexseq_prepare_annotation.py gtffile.gtf dexseq.gff
python dexseq_count.py dexseq.gff sample1.sam sample1.txt
Once we have repeated the dexseq_count script for each sample, we can read the data into R using the code chunks below. As we are working with pre-prepared data, we first point to these files which live within the pasilla package.
The pasilla package contains counts from an experiment by Brooks et al
We will run DEXSeq on a subset of the genes, for demonstration purposes.
library("pasilla")
inDir = system.file("extdata", package="pasilla", mustWork=TRUE)
countFiles = list.files(inDir, pattern="fb.txt$", full.names=TRUE)
flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
genesForSubset = read.table(file.path(inDir, "geneIDsinsubset.txt"),
  stringsAsFactors=FALSE)[[1]]
As in DESeq2 we use a sampleTable to define the samples:
sampleTable = data.frame(
  row.names = c( "treated1", "treated2", "treated3",
    "untreated1", "untreated2", "untreated3", "untreated4" ),
  condition = c("knockdown", "knockdown", "knockdown",
    "control", "control", "control", "control" ),
  libType = c( "single-end", "paired-end", "paired-end",
    "single-end", "single-end", "paired-end", "paired-end" ) )
sampleTable
We now read the data into a DEXSeqDataSet object:
library("DEXSeq")
dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )
Subset the genes, for demonstration purposes:
dxd = dxd[geneIDs( dxd ) %in% genesForSubset,]
Now we run the estimation and testing functions:

dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )
dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
The following code extracts a results table, makes an MA-plot, and draws the expression levels over the exons to highlight differential exon usage:

dxr = DEXSeqResults( dxd )
plotMA( dxr, cex=0.8 )
plotDEXSeq( dxr, "FBgn0010909", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
Again, drawing the expression levels, now showing the annotated transcripts below:

plotDEXSeq( dxr, "FBgn0010909", displayTranscripts=TRUE, legend=TRUE,
              cex.axis=1.2, cex=1.3, lwd=2 )
For more details on the DEXSeq software, see the vignette and the paper, which is linked from the vignette page:
browseVignettes("DEXSeq")
We conclude by adding the session information:
sessionInfo()


# Rsem software installation: 
# go to the directory with RSEM software
make 
make install DESTDIR=/home/emorycpl prefix=/anaconda2

##
time STAR --runThreadN 16 \
--genomeDir ~/Rseq/star/GRCh38/star_indices_overhang100 \
--readFilesIn ~/emory/fastq/ENCFF001RFH.fastq.gz ~/emory/fastq/ENCFF001RFG.fastq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate Unsorted
ls -shl
##
time STAR --runMode inputAlignmentsFromBAM --inputBAMfile Aligned.sortedByCoord.out.bam --outWigType bedGraph --outWigStrand Stranded

##
mkdir rsemGenome
time rsem-prepare-reference --gtf ~/emory/gtf/Homo_sapiens.GRCh38.88.exons.gtf ~/emory/genome/GRCh38_r88.all.fa ~/emory/rsemGenome/GRCh38.88.all
time bowtie2-build -f rsemGenome/GRCh38.88.all.idx.fa rsemGenome/GRCh38.88.all
## please note that you may can try bowtie2-build -f ...
time rsem-calculate-expression -p 16 --paired-end ~/emory/fastq/ENCFF001RFG.fastq ~/emory/fastq/ENCFF001RFH.fastq rsemGenome/GRCh38.88.all ENCFF001RF

##
sudo apt-get install your_pkgname
sudo apt-cache search probably_pkgname
wget http://somewhere/your_debname.deb
sudo dpkg -i your_debname.deb
wget http://somewhere/your_software.tar.gz
tar xzf your_software.tar.gz -C ~/app
echo "export PATH=$PATH:$HOME/app/your_software_path" >> ~/.bashrc
source ~/.bashrc
git clone https://github.com/some_rep
cd some_rep && make
sudo make install
## or you can do:

wget http://somewhere/software_sourcecode.tar.gz
tar xzf software_sourcecode.tar.gz -C ~/app
cd ~/app/software_sourcecode_path && make
sudo make install


time STAR --runThreadN 16 \
--genomeDir ~/Rseq/idex \
--readFilesIn ~/emory/fastq/WT_1_1.fq ~/emory/fastq/WT_1_2.fq \
--quantMode TranscriptomeSAM

###--outSAMtype BAM SortedByCoordinate 


awk '$3 == "exon"' ~/Rseq/gtfm/Mus_musculus.GRCm38.89.gtf > ~/Rseq/gtfm/Mus_musculus.GRCm38.89.exons.gtf

mkdir rsem_ref

rsem-prepare-reference --gtf ~/Rseq/gtfm/Mus_musculus.GRCm38.89.exons.gtf \
~/Rseq/ mreference/GRCm38.all.fa ~/Rseq/rsem_ref/ref

time bowtie2-build -f ~/Rseq/rsem_ref/ref.idx.fa ~/Rseq/rsem_ref/ref
time rsem-calculate-expression -p 16 --paired-end --bowtie2 ~/emory/fastq/MT_2_1.fq ~/emory/fastq/MT_2_2.fq ~/Rseq/rsem_ref/ref MT
