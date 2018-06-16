#!/bin/sh
set -euo pipefail
if [ ! -d /home/mxie/Projects/miRNA/test ]; 
then mkdir /home/mxie/Projects/miRNA/test;
fi
if [ ! -d /home/mxie/Projects/miRNA/test/counted ]; 
then mkdir /home/mxie/Projects/miRNA/test/counted;
fi

if [ ! -d /home/mxie/Projects/miRNA/test/QC_mappedlength ]; 
then mkdir /home/mxie/Projects/miRNA/test/QC_mappedlength;
fi
if [ ! -d /home/mxie/Projects/miRNA/test/trimmed ]; 
then mkdir /home/mxie/Projects/miRNA/test/trimmed;
fi
if [ ! -d /home/mxie/Projects/miRNA/test/QC_readlength ]; 
then mkdir /home/mxie/Projects/miRNA/test/QC_readlength;
fi

if [ ! -d /home/mxie/Projects/miRNA/test/unmapped ]; 
then mkdir /home/mxie/Projects/miRNA/test/unmapped;
fi
if [ ! -d /home/mxie/Projects/miRNA/test/mapped ]; 
then mkdir /home/mxie/Projects/miRNA/test/mapped;
fi

ls  *.fastq.gz| while read id
do
echo $id
file="${id:0:8}"
zcat $file* | pigz > /home/mxie/Projects/miRNA/test/$file.fastq.gz; 
done

cd /home/mxie/Projects/miRNA/test
ls *fastq.gz |while read id
do
echo $id
#/sw/hgcc/Pkgs/Anaconda2/4.2.0/bin/cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o /home/mxie/Projects/miRNA/test/trimmed/${id%%.*}_miRNAtrimmed.fastq.gz $id > /home/mxie/Projects/miRNA/test/trimmed/${id%%.*}_adaprm.log
/sw/hgcc/Pkgs/Anaconda2/4.2.0/bin/cutadapt -m 15 -a AGATCGGAAGAGCACACGTCT $id > /home/mxie/Projects/miRNA/test/trimmed/${id%%.*}_miRNAtrimmed.fastq.gz
#java -Xmx2g -jar /sw/hgcc/Pkgs/Trimmomatic/0.36/trimmomatic-0.36.jar SE -phred33 $id  /home/mxie/Projects/miRNA/test/trimmed/${id%%.*}_miRNAtrimmed.fastq.gz ILLUMINACLIP:/home/mxie/Projects/miRNA/src/illumina_RGA_3p_adapter.fa:2:30:10;
done     

cd /home/mxie/Projects/miRNA/test/trimmed
ls *_miRNAtrimmed.fastq.gz |while read id
do
echo $id
gunzip -c $id | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c >/home/mxie/Projects/miRNA/test/QC_readlength/${id%%.*}_length.txt;
done

cd /sw/hgcc/Pkgs/SHRiMP/2.2.3
export SHRIMP_FOLDER=$PWD
cd /home/mxie/Projects/miRNA/test/trimmed
ls *miRNAtrimmed.fastq.gz |while read id
do
echo $id
$SHRIMP_FOLDER/bin/gmapper-ls --strata -o 1 --qv-offset 33 -Q -L /home/mxie/Projects/miRNA/src/miRbase/mature.human-ls $id \
-N 14 --mode mirna -w 170% -E --un /home/mxie/Projects/miRNA/test/unmapped/${id}_unmapped-miRNA.fastq \
 >/home/mxie/Projects/miRNA/test/mapped/${id%%.*}_mapping.sam 2>/home/mxie/Projects/miRNA/test/mapped/${id%%.*}_maping.log;
done

cd /home/mxie/Projects/miRNA/test/mapped

ls  *.sam| while read id
do
echo $id
file="${id:0:8}"
/sw/hgcc/Pkgs/samtools/1.5/bin/samtools view -Su $id > /home/mxie/Projects/miRNA/test/mapped/${file}_miRNAtrimmed_mapping.bam;
/sw/hgcc/Pkgs/samtools/1.5/bin/samtools sort /home/mxie/Projects/miRNA/test/mapped/${file}_miRNAtrimmed_mapping.bam > /home/mxie/Projects/miRNA/test/mapped/${file}_miRNAtrimmed_to_MB_mt_sorted.bam
/sw/hgcc/Pkgs/samtools/1.5/bin/samtools index /home/mxie/Projects/miRNA/test/mapped/${file}_miRNAtrimmed_to_MB_mt_sorted.bam
done

cd /home/mxie/Projects/miRNA/test/mapped
for tbf in *_miRNAtrimmed_to_MB_mt_sorted.bam; 
do sid=${tbf%%_*}; 
if [ ! -e /home/mxie/Projects/miRNA/test/counted/${sid}_mt_subcount.txt ]; 
then echo /home/mxie/Projects/miRNA/test/counted/${sid}_mt_subcount.txt; 
for mirna in `/sw/hgcc/Pkgs/samtools/1.5/bin/samtools view -h ${tbf} | gawk '$1 == "@SQ"{printf("%s\n",substr($2,4))}' | sort`; do cnt=`/sw/hgcc/Pkgs/samtools/1.5/bin/samtools view -c ${tbf} ${mirna}`; echo "${mirna}	${cnt}"; 
done > /home/mxie/Projects/miRNA/test/counted/${sid}_mt_subcount.txt; fi; done

cd /home/mxie/Projects/miRNA/test/mapped
ls  *.sam| while read id
do
echo $id
file="${id:0:8}"
/sw/hgcc/Pkgs/samtools/1.5/bin/samtools view $id | gawk '{mlarr[length($10)]++}END{for (len = 1; len <=101; len++) printf("%d\t%d\n",len,mlarr[len])}' | sort -k 1n,1n  > /home/mxie/Projects/miRNA/test/QC_mappedlength/${file}_mapped_length_distribution_data.txt;
done
