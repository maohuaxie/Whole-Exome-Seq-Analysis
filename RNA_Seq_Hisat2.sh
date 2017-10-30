# some RNA-seq pipeline, please check the following link : https://wiki2.org/en/List_of_RNA-Seq_bioinformatics_tools
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50177
# download the sra files
for ((i=677;i<=680;i++)) ;
do wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP029/SRP029245/SRR957$i/SRR957$i.sra;
done
# convert to fastq files
ls *sra |while read id; 
do fastq-dump --split-3 $id;
done
# check the data quality with fastqc
ls *fastq |xargs fastqc -t 10
# mapping read to h19 reference genome
reference=/home/emorycpl/ref/h19/ucsc.hg19.fasta
IDX=/home/emorycpl/ref/h19/ucsc.hg19.fasta
# Build a hisat2 index for the genome.
hisat2-build $reference $IDX

hisat2 -p 12 -x $IDX -U SRR957677.fastq -S control_1.sam 2>control_1.log
hisat2 -p 12 -x $IDX -U SRR957678.fastq -S control_2.sam 2>control_2.log
hisat2 -p 5 -x $IDX -U SRR957679.fastq -S siSUZ12_1.sam 2>siSUZ12_1.log
hisat2 -p 5 -x $IDX -U SRR957680.fastq -S siSUZ12_2.sam 2>siSUZ12_2.log
# convert sam files to bam files
ls *sam |while read id;
do (nohup samtools sort -n -@ 12 -o ${id%%.*}.Nsort.bam $id &);
done
# using htseq-counts to quantify the gene experssion level
ls *.Nsort.bam |while read id;
do (nohup samtools view $id | nohup htseq-count -f sam -s no -i gene_name - ~/reference/gtf/gencode/gencode.v25lift37.annotation.gtf 1>${id%%.*}.geneCounts 2>${id%%.*}.HTseq.log&);
done

In the case of HISAT2 , in addition to the reference genome index, other annotation files are also present.
A pre-built index for GRCh38 assembly can be obtained from the command line:
# make a directory to store the references.
mkdir -p ./refs
# This is the URL for the index.
URL=ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz
# Download and unpack the index.
curl $URL | tar zxv
# Move the index directory to the reference directory
mv grch38 ./refs
# The index has the prefix of "genome".
IDX=refs/grch38/genome

How do I quantify the HISAT2 alignments?
First, you need a file that contains the genomic coordinates of annotated genes:
# Set up some shortcuts.
mkdir -p refs
GTF=refs/GRCh38.gtf
URL=ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz
#Download feature file from Ensembl and unzip before saving.
curl $URL | gunzip -c > $GTF

Run featureCounts on a single dataset to count reads on genes:
featureCounts -a $GTF -o counts.txt bam/SRR3191542.bam
