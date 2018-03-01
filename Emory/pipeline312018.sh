
/home/twingo/Projects/AG056533/Amp-AD/RosMap/WGS
/home/twingo/Projects/AG056533/Amp-AD/Mayo/WGS
/home/twingo/Projects/AG056533/Amp-AD/MSBB/WGS

ls directory | wc -l

module avail
http://wingolab.org/hgcc/
qrsh  # change node00 to node07

zcat file | more to see *.gz file 

module load vcftools/0.1.13

 1  ls -s
    2  ls -al
    3  cd /home/twingo/Projects/AG056533/Amp-AD/RosMap/WGS
    4  cd ~
    5  ls -al
    6  module avail 
    7  cd /home/twingo/Projects/AG056533/Amp-AD/RosMap/WGS
    8  ls -wc
    9  ls -c
   10  wc -l
   11  ls | wc -l
   12  module avail
   13  module load samtools/1.5
   14  more DEJ_11898_B01_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz
   15  cd /home/twingo/Projects/AG056533/Amp-AD/MSBB/WGS
   16  ls | wc -l
   17  clear
   18  cat SCH_11923_B02_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz
   19  cd /home/twingo/Projects/AG056533/Amp-AD/RosMap/WGS
   20  tar zxvf DEJ_11898_B01_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz
   21  gunzip .gz
   22  tar: This does not look like a tar archive
   23  gunzip DEJ_11898_B01_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz
   24  zcat DEJ_11898_B01_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz
   25  zcat DEJ_11898_B01_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz | more
   26  moudule avail
   27  module avail
   28  cd /home/twingo/Projects/AG056533/Amp-AD/RosMap/WGS
   29  module avail 
   30  module load vcftools/0.1.13
   31  vcf-merge DEJ_11898_B01_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz DEJ_11898_B01_GRM_WGS_2017-05-15_11.recalibrated_variants.vcf.gz > text.vcf.gz
   32  sudo vcf-merge DEJ_11898_B01_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz DEJ_11898_B01_GRM_WGS_2017-05-15_11.recalibrated_variants.vcf.gz > text.vcf.gz
   33  cd ..
   34  cd twingo
   35  cd projects
   36  cd Projects/
   37  cd AG056533/
   38  cd Amp-AD/
   39  cd Mayo/
   40  cd WGS
   41  ls -ish
   42  module avail
   43  module load vcftools/0.1.13
   44  vcf-merge SCH_11923_B01_GRM_WGS_2017-04-27_Y.recalibrated_variants.annotated.vcf.gz SCH_11923_B01_GRM_WGS_2017-04-27_X.recalibrated_variants.annotated.vcf.gz > test.vcf.gz
   45  pwd
   46  cd 
   47  mkdir Projects
   48  cd Projects/
   49  mkdir WGS
   50  cd WGS/
   51  mkdir data src doc analysis
   52  ll
   53  cd data/
   54  ln -sv /home/twingo/Projects/AG056533/Amp-AD/Mayo/WGS/*.vcf.gz* .
   55  l
   56  ls -lah
   57  cd ..
   58  cd analysis/
   59  l
   60  history |grep bcf
   61  history |grep vcf
   62  cd ..
   63  vim README.md
   64  git init
   65  git status
   66  git add -A.
   67  git add -A .
   68  git status
   69  l
   70  git commit -m 'initial commit'
   71  git log
   72  l
   73  cd data/
   74  l
   75  cd ..
   76  pwd
   77  cd -
   78  mkdir Mayo
   79  mv * Mayo/
   80  l
   81  mkdir MSBB RosMap
   82  cd Mayo/
   83  ls -ish
   84  cd ..
   85  history 


bcftools concat SCH_11923_B02_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz SCH_11923_B02_GRM_WGS_2017-05-15_11.recalibrated_variants.vcf.gz > ../../analysis/test.vcf
bcftools concat *.vcf.gz > ../../analysis/Mayo.vcf.gz


   1  ls -s
    2  ls -al
    3  cd /home/twingo/Projects/AG056533/Amp-AD/RosMap/WGS
    4  cd ~
    5  ls -al
    6  module avail 
    7  cd /home/twingo/Projects/AG056533/Amp-AD/RosMap/WGS
    8  ls -wc
    9  ls -c
   10  wc -l
   11  ls | wc -l
   12  module avail
   13  module load samtools/1.5
   14  more DEJ_11898_B01_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz
   15  cd /home/twingo/Projects/AG056533/Amp-AD/MSBB/WGS
   16  ls | wc -l
   17  clear
   18  cat SCH_11923_B02_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz
   19  cd /home/twingo/Projects/AG056533/Amp-AD/RosMap/WGS
   20  tar zxvf DEJ_11898_B01_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz
   21  gunzip .gz
   22  tar: This does not look like a tar archive
   23  gunzip DEJ_11898_B01_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz
   24  zcat DEJ_11898_B01_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz
   25  zcat DEJ_11898_B01_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz | more
   26  moudule avail
   27  module avail
   28  cd /home/twingo/Projects/AG056533/Amp-AD/RosMap/WGS
   29  qstat
   30  qstat –u ‘\*’
   31  qstat -r
   32  qrsh
   33  /home/mxie
   34  cd mxie
   35  pwd
   36  cd projects
   37  cd Projects/
   38  cd WGS
   39  cd analysis/
   40  zcat test.vcf.gz | more
   41  zcat test.vcf.gz 
   42  rm -rf test.vcf.gz 
   43  pwd
   44  cat test.vcf 
   45  cat test.vcf
   46  rm -rf test.vcf
   47  clear
   48  rm -rf test.vcf 
   49  more test.vcf
   50  clear
   51  ls -ish
   52  rm -rf text.vcf
   53  rm -rf test.vcf
   54  ls -ish
   55  rm -rf MSBB.vcf
   56  clear
   57  ls -ish
   58  cd /home/twingo/Projects/AG056533/Amp-AD/RosMap/WGS
   59  module avail 
   60  module load vcftools/0.1.13
   61  vcf-merge DEJ_11898_B01_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz DEJ_11898_B01_GRM_WGS_2017-05-15_11.recalibrated_variants.vcf.gz > text.vcf.gz
   62  sudo vcf-merge DEJ_11898_B01_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz DEJ_11898_B01_GRM_WGS_2017-05-15_11.recalibrated_variants.vcf.gz > text.vcf.gz
   63  cd ..
   64  cd twingo
   65  cd projects
   66  cd Projects/
   67  cd AG056533/
   68  cd Amp-AD/
   69  cd Mayo/
   70  cd WGS
   71  ls -ish
   72  module avail
   73  module load vcftools/0.1.13
   74  vcf-merge SCH_11923_B01_GRM_WGS_2017-04-27_Y.recalibrated_variants.annotated.vcf.gz SCH_11923_B01_GRM_WGS_2017-04-27_X.recalibrated_variants.annotated.vcf.gz > test.vcf.gz
   75  pwd
   76  cd 
   77  mkdir Projects
   78  cd Projects/
   79  mkdir WGS
   80  cd WGS/
   81  mkdir data src doc analysis
   82  ll
   83  cd data/
   84  ln -sv /home/twingo/Projects/AG056533/Amp-AD/Mayo/WGS/*.vcf.gz* .
   85  l
   86  ls -lah
   87  cd ..
   88  cd analysis/
   89  l
   90  history |grep bcf
   91  history |grep vcf
   92  cd ..
   93  vim README.md
   94  git init
   95  git status
   96  git add -A.
   97  git add -A .
   98  git status
   99  l
  100  git commit -m 'initial commit'
  101  git log
  102  l
  103  cd data/
  104  l
  105  cd ..
  106  pwd
  107  cd -
  108  mkdir Mayo
  109  mv * Mayo/
  110  l
  111  mkdir MSBB RosMap
  112  cd Mayo/
  113  ls -ish
  114  cd ..
  115  history 
  116  cd RosMap/
  117  ln -sv /home/twingo/Projects/AG056533/Amp-AD/RosMap/WGS/*.vcf.gz* .
  118  cd ..
  119  ln -sv /home/twingo/Projects/AG056533/Amp-AD/MSBB/WGS/*.vcf.gz* .
  120  rm -rf *.vcf.gz*
  121  cd MSBB
  122  ln -sv /home/twingo/Projects/AG056533/Amp-AD/MSBB/WGS/*.vcf.gz* .
  123  ls -ihs
  124  ls -ish
  125  clear
  126  history
  127  vcf-merge *variants.vcf.gz > SCH_11923_B02_GRM_WGS_2017-05-15_recalibrated_variants.vcf.gz
  128  module avail
  129  module load bcftools/1.5
  130  bcftools merge *_variants.vcf.gz > ../../analysis/test.vcf.gz
  131  bcftools merge *_variants.vcf.gz > ../../analysis/test.vcf
  132  bcftools merge SCH_11923_B02_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz SCH_11923_B02_GRM_WGS_2017-05-15_11.recalibrated_variants.vcf.gz > ../../analysis/test.vcf
  133  bcftools concatenate VCF/BCF SCH_11923_B02_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz SCH_11923_B02_GRM_WGS_2017-05-15_11.recalibrated_variants.vcf.gz > ../../analysis/test.vcf
  134  bcftools concat SCH_11923_B02_GRM_WGS_2017-05-15_10.recalibrated_variants.vcf.gz SCH_11923_B02_GRM_WGS_2017-05-15_11.recalibrated_variants.vcf.gz > ../../analysis/test.vcf
  135  rm -rf *X*.gz
  136  rm -rf *Y*.gz
  137  rm -rf *.tbi
  138  rm -rf *annotated.vcf.gz
  139  rm -rf SCH_11923_B02_GRM_WGS_2017-05-15_others.recalibrated_variants.vcf.gz SCH_11923_B02_GRM_WGS_2017-05-15_recalibrated_variants.vcf.gz
  140  bcftools concat *.vcf.gz > ../../analysis/MSBB.vcf
  141  clear
  142  nohup bcftools concat *.vcf.gz > ../../analysis/MSBB.vcf.gz
  143  cd Projects/WGS/
  144  cd data
  145  cd MSBB/
  146  rm -rf *.gz
  147  ln -sv /home/twingo/Projects/AG056533/Amp-AD/MSBB/WGS/*.vcf.gz* .
  148  cd ..
  149  cd analysis/
  150  cd ..
  151  nano README.md 
  152  cd src
  153  cd ..
  154  clear
  155  doc
  156  cd doc
  157  cd ..
  158  nano README.md 
  159  clear
  160  cd data
  161  cd ..
  162  cp README.md ./data/Mayo/
  163  cp README.md ./data/MSBB/
  164  cp README.md ./data/RosMap/
  165  cd data
  166  cd Mayo/
  167  clear
  168  cd ..
  169  cd MSBB/
  170  nano README.md 
  171  clear
  172  cd ..
  173  cd RosMap/
  174  cd ..
  175  cp README.md ./data/RosMap/
  176  cd data
  177  cd RosMap/
  178  nano README.md 
  179  clear
  180  cd ..
  181  cd Mayo
  182  rm -rf *.annotated.*
  183  rm -rf *.tbi
  184  rm -rf *X*
  185  rm -rf *Y*
  186  rm -rf SCH_11923_B01_GRM_WGS_2017-04-27_others.recalibrated_variants.vcf.gz
  187  nohup bcftools concat *.vcf.gz > ../../analysis/Mayo.vcf.gz
  188  history
  
nohup bcftools concat *.vcf.gz > ../../analysis/RosMap.vcf.gz

git add analysis
git commit -m "adding analysis directory"


chmod –R 755 /home/mxie/Projects/WGS/analysis


Yep, and it would be faster to do in pigz if you have it installed (or install it if not; multithreaded, will matter for 2.2TB):

Make a command.sh
#!/bin/sh
for file in *.vcf.gz; do pigz $file; done

qsub -cwd -pe smp 24 command.sh

It would be even better to do that in parallel, submitting 1 qsub job per file. So each qsub would be in form of

nthCommand.sh
#!/bin/sh
pigz fileN.vcf.gz

qsub -cwd -pe smp 24 nthCommand.sh

Thanks!
Alex


grep PASS <(pigz -d -c *.vcf.gz )

cut -f7 file.vcf | grep “.” | wc -l

grep -P “\#|PASS” file.vcf | pigz -c - > file_pass.vcf.gz

pigz -d -c Mayo.vcf.gz | cut -f7 | grep -P "^\.$" | wc -l


vcf2snp.sh
#!/bin/sh
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE[\t%IUPACGT]\n' MSBB.vcf.gz | vcfToSnp > MSBB.snp

qsub -cwd -pe smp 24 vcf2snp.sh



Hi Maohua,

1.	The goal is to use the MBSS.vcf.gz WGS to generate personal protein databases for https://github.com/wingolab-org/GenPro. 
2.	We need to know the TMT batch (proteomics batches) of the IDs so that each sample that was run in the batch can be included together into a single database.
3.	We only really care to make proteomic databases for the individuals that had proteomics done (see column “"Proteomics Barcode BM-10#" in the Emory_MSS_MASTER_03252015.xls file.
4.	My helper scripts and the binary database - /mnt/22qFSS/data2/GenPro

I will get you the lookup table for the Mt. Sinai (e.g., “MSBB”). 

Best wishes,

Thomas

find -name 


mergeVcf.sh

#!/bin/sh

if [ $# != 1 ]; then
  echo "$0 <out name>"
  exit
fi

PRJDIR="/mnt/22qFSS/data2/GenPro/1000g/phase3_all"

if [ -e /bin/mktemp ]; then
  TMPDIR=`/bin/mktemp -d -p /scratch/` || exit
elif [ -e /usr/bin/mktemp ]; then
  TMPDIR=`/usr/bin/mktemp -d –p /scratch/` || exit
else
  echo “Error. Cannot find mktemp to create tmp directory”
  exit
fi

/home/twingo/local/bin/bcftools concat -a -d none -O z -o ${TMPDIR}/${1} ALL.chr*.vcf.gz
gunzip -c ${1} | wc -l > ${PRJDIR}/countVcf.txt
rsync ${TMPDIR}/${1} ${PRJDIR}/
rm -rf ${TMPDIR}



vcfToSnp.sh

#!/bin/sh
module load bcftools/1.5
if [ $# != 1 ]; then
  echo "$0 <vcf file>"
  exit
fi

if [ -e /bin/mktemp ]; then
  TMPDIR=`/bin/mktemp -d -p /scratch/` || exit
elif [ -e /usr/bin/mktemp ]; then
  TMPDIR=`/usr/bin/mktemp -d –p /scratch/` || exit
else
  echo “Error. Cannot find mktemp to create tmp directory”
  exit
fi

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE[\t%IUPACGT]\n' ${1} | vcfToSnp > ${TMPDIR}/${1}.snp
rsync ${TMPDIR}/${1}.snp .
rm -rf ${TMPDIR}

module unload bcftools/1.5


/mnt/22qFSS/data2/GenPro/1000g


/home/mxie/Projects/WGS/analysis1

/mnt/22qFSS/data2/GenPro/1000g/phase1_all


1/2/2018
in /home/mxie/Projects/WGS/analysis1 directory: submit 

qsub -cwd -pe smp 24 vcfToSnp.sh


Hi Maohua,

I’d like you to help map and merge WGS fastq files for the following directories:

/home/twingo/NYGC/PBNAS/twingo/upload/nygc

I’ve copied David Cutler (who wrote PEMapper/PECaller) and who may have a new version for us to use soon.

It’s not a high priority but needs to be done this month. We’ll use PEMapper/PECaller - https://github.com/wingolab-org/pecaller

Best wishes,
Thomas



vcfToSnp1.sh

#!/bin/sh
module load bcftools/1.5

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE[\t%IUPACGT]\n' MSBB.vcf.gz | vcfToSnp > MSBB.vcf.gz.snp

module unload bcftools/1.5



qsub -cwd -pe smp 24 vcfToSnp1.sh

## This one is working, you have to specify the vcfToSnp file
vcfToSnp2.sh

#!/bin/sh
module load bcftools/1.5

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE[\t%IUPACGT]\n' MSBB.vcf.gz | /home/mxie/Projects/WGS/analysis1/vcfToSnp > MSBB.vcf.gz.snp

module unload bcftools/1.5

qsub -cwd -pe smp 24 vcfToSnp2.sh

1/4/2018
nano vcfToSnp1.sh
#!/bin/sh
module load bcftools/1.5
for file in *.vcf.gz;
do
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE[\t%IUPACGT]\n' $file |/home/mxie/Projects/WGS/data/MSBB/vcfToSnp > $file.snp;
done

module unload module load bcftools/1.5

qsub -cwd -pe smp 24 vcfToSnp1.sh

1/5/2018
zcat
bcftools query -l MSBB.vcf.gz 

makepdbstep1.sh
#!/bin/sh
GenPro_make_perprotdb1.pl -v -s MSBB.vcf.gz.snp -l /mnt/22qFSS/data2/GenPro/hg19/idx -n hg19 -o /home/mxie/Projects/WGS/analysis1/db  -w /home/mxie/Projects/WGS/analysis1/ids.MSBB.batch_1.txt

qsub -cwd -pe smp 24 makepdbstep1.sh

runall_make


qsub -v USER -v PATH -cwd -t 1-26 runnall_make_refprot.sh \
MSBB.vcf.gz.snp \
  <binary genome index> \
  <output dir>


sed -i 's/[*]//g'  msbb.head.snp   # not working

sed 's/[[0-9]*]//g' msbb.head.snp > msbb.snp

sed 's/[[0-9]*]//g' MSBB.vcf.gz.snp  > MSBB.vcf.gz.m.snp 

1/6/2018
makepdbstep1.sh
#!/bin/sh
GenPro_make_perprotdb1.pl -v -s MSBB.vcf.gz.snp -l /mnt/22qFSS/data2/GenPro/hg19/idx -n hg19 -o /home/mxie/Projects/WGS/analysis1/db  -w /home/mxie/Projects/WGS/analysis1/ids.MSBB.batch_1.txt

qsub -cwd -pe smp 24 makepdbstep1.sh

1/7/2018
sed 's/:GT//g' MSBB.vcf.gz.m.snp > MSBB.vcf.gz.m1.snp

makepdbstepm1.sh
#!/bin/sh
GenPro_make_perprotdb1.pl -v -s MSBB.vcf.gz.snp -l /mnt/22qFSS/data2/GenPro/hg19/idx -n hg19 -o /home/mxie/Projects/WGS/analysis1/db  -w /home/mxie/Projects/WGS/analysis1/ids.MSBB.batch_1.txt

qsub -cwd -pe smp 24 makepdbstepm1.sh

sed 's/^\([0-9].*\)/char\1/g' (not try)

sed 's/:GT//g' msbb.head.10.snp | awk '{print "chr"$0}' > MSBB.snp

sed 's/:GT//g' MSBB.vcf.gz.m.snp | awk '{print "chr"$0}' > MSBB.snp
sed -i 's/chrFragment/Fragment/' MSBB.snp
awk '{if(NR==1){print $0} else{{print "chr"$0}}}' 



makepdbstep1m2.sh
#!/bin/sh
GenPro_make_perprotdb1.pl -v -s MSBB.snp -l /mnt/22qFSS/data2/GenPro/hg19/idx -n hg19 -o /home/mxie/Projects/WGS/analysis1/db2  -w /home/mxie/Projects/WGS/analysis1/ids.MSBB.batch_1.txt

qsub -cwd -pe smp 24 makepdbstep1m2.sh


node07:[GenPro] % cat makePdbStep2.sh
#!/bin/sh

if [ $# != 1 ]; then
  echo "$0 <id>"
  exit
fi

PRJDIR="/mnt/22qFSS/data2/GenPro"
REFDIR="/mnt/22qFSS/data2/GenPro/hg19/refProt"
PPDDIR="${PRJDIR}/hg19_phase1.snp"
OUTDIR="${PRJDIR}/hg19_phase1_${1}"

mkdir -p ${OUTDIR}

if [ -e /bin/mktemp ]; then
  TMPDIR=`/bin/mktemp -d -p /scratch/` || exit
elif [ -e /usr/bin/mktemp ]; then
  TMPDIR=`/usr/bin/mktemp -d â€“p /scratch/` || exit
else
  echo â€œError. Cannot find mktemp to create tmp directoryâ€
  exit
fi

GenPro_make_perprotdb2.pl -i ${1} -d ${PPDDIR} -r ${REFDIR} -o ${TMPDIR}
rsync ${TMPDIR}/* ${OUTDIR}/
rm -rf ${TMPDIR}

node07:[GenPro] % 

1/8/2018
makepdbstep2.sh
#!/bin/sh
GenPro_make_perprotdb2.pl -i 68411 -d db1 -r /mnt/22qFSS/data2/GenPro/hg19/refProt -o Json1

qsub -cwd -pe smp 24 makepdbstep2.sh

pigz -d -c Mayo.vcf.gz | cut -f7 | grep -P "^\.$" | wc -l


makepdbstep2m.sh
#!/bin/sh
for file in 68416 68419 68420 68422 68423; 
do
GenPro_make_perprotdb2.pl -i $file -d db1 -r /mnt/22qFSS/data2/GenPro/hg19/refProt -o Json1;
done

qsub -cwd -pe smp 24 makepdbstep2m.sh

1/9/2018
makepdbstep2m.sh
#!/bin/sh
for file in 68424 68425 68427 68428 68430;
do
GenPro_make_perprotdb2.pl -i $file -d db1 -r /mnt/22qFSS/data2/GenPro/hg19/refProt -o Json1;
done

qsub -cwd -pe smp 24 makepdbstep2m.sh

1/11/2018

/home/twingo/PBNAS/upload/nygc

makepdbstep3m.sh
#!/bin/sh
for file in $(cat ./ids.MSBB.batch_1_m.txt);
do
GenPro_make_perprotdb2.pl -i $file -d db1 -r /mnt/22qFSS/data2/GenPro/hg19/refProt -o Json1;
done

qsub -cwd -pe smp 24 makepdbstep3m.sh

1/16/2018

Here’s an example of the script running correctly. NOTE: I had to rename the file to get it to work (i.e., “.R?” to “_R?”).

hgcc:twingo@node00:[~/Projects/nygenomeComanyRemap/data] [Thu Jan 11 15:43:22]
==> map_directory_array2.pl HDA-01369 ${Hg38}

 Matching CGND-HDA-01369_TAATGCGC_HGKLJALXX_L002_001_R1 with CGND-HDA-01369_TAATGCGC_HGKLJALXX_L002_001_R2 

 Matching CGND-HDA-01369_TAATGCGC_HGKLJALXX_L003_001_R1 with CGND-HDA-01369_TAATGCGC_HGKLJALXX_L003_001_R2 

 Matching CGND-HDA-01369_TAATGCGC_HGKLJALXX_L004_001_R1 with CGND-HDA-01369_TAATGCGC_HGKLJALXX_L004_001_R2 

 Matching CGND-HDA-01369_TAATGCGC_HGKLJALXX_L005_001_R1 with CGND-HDA-01369_TAATGCGC_HGKLJALXX_L005_001_R2 

 Matching CGND-HDA-01369_TAATGCGC_HGKLJALXX_L006_001_R1 with CGND-HDA-01369_TAATGCGC_HGKLJALXX_L006_001_R2 

 Matching CGND-HDA-01369_TAATGCGC_HGKLJALXX_L007_001_R1 with CGND-HDA-01369_TAATGCGC_HGKLJALXX_L007_001_R2 

 Matching CGND-HDA-01369_TAATGCGC_HGKLJALXX_L008_001_R1 with CGND-HDA-01369_TAATGCGC_HGKLJALXX_L008_001_R2 
Your job 244530 ("HDA-01369.sh") has been submitted
hgcc:twingo@node00:[~/Projects/nygenomeComanyRemap/data] [Thu Jan 11 15:43:23]
==> qstat
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
 244530 0.21350 HDA-01369. twingo       r     01/11/2018 15:43:28 b.q@node02.local                  26        


Best wishes,

Thomas


From: "Wingo, Thomas" <thomas.wingo@emory.edu>
Date: Thursday, January 11, 2018 at 3:21 PM
To: "Xie, Maohua" <mxie6@emory.edu>
Cc: Aliza Wingo <aliza.wingo@emory.edu>
Subject: mapping genomes

Use `map_directory_array2.pl` to map all files in a single directory. Use Hg38 genome, which I set as an environmental variable: `export Hg38=/home/dcutler/hg38/hg38.sdx`. 
 
https://github.com/wingolab-org/pecaller
 
Best wishes,
Thomas
 
Thomas S. Wingo, MD
Physician Scientist, Atlanta Veterans Affairs Medical Center 
Assistant Professor, Departments of Neurology and Human Genetics
Emory University Center for Neurodegenerative Disease
office: 404-727-4905
website: http://wingolab.org
 

Hi Maohua,

See MSBB_WES_covariates.csv file in the /home/twingo/Projects/AG056533/Amp-AD/MSBB/Metadata directory to link the WGS IDs to the sample Ids for the proteomics. The proteomics directory is here: /home/twingo/Projects/AG056533/MSBB_Proteomics

Best wishes,

Thomas

1/21/2018

makepdbstep1b2.sh
#!/bin/sh
GenPro_make_perprotdb1.pl -v -s MSBB.snp -l /mnt/22qFSS/data2/GenPro/hg19/idx -n hg19 -o /home/mxie/Projects/WGS/analysis1/db2  -w /home/mxie/Projects/WGS/analysis1/ids.MSBB.batch_2.txt

qsub -cwd -pe smp 24 makepdbstep1b2.sh

makepdbstep2b2.sh
#!/bin/sh
for file in $(cat ./ids.MSBB.batch_2.txt);
do
GenPro_make_perprotdb2.pl -i $file -d db2 -r /mnt/22qFSS/data2/GenPro/hg19/refProt -o Json2;
done

qsub -cwd -pe smp 24 makepdbstep2b2.sh

1/23/2018
### create 5 tables and load csv files into tables

DROP TABLE IF EXISTS mssbwgs;
create table mssbwgs ( individualIdentifier varchar(100) NOT NULL, WGS_Barcode int NOT NULL, PRIMARY KEY (WGS_Barcode));

LOAD DATA LOCAL INFILE '/home/maohuaxie/Emoryjob/metadata/AMP-AD_MSBB_WGS__sample_barcode_brainBankID.csv' 
INTO TABLE mssbwgs 
FIELDS TERMINATED BY ',' 
ENCLOSED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS;

create table barcode10 (Barcode int NOT NULL, Proteomics_BM_10_Barcode int NULL , PRIMARY KEY (Barcode));

LOAD DATA LOCAL INFILE '/home/maohuaxie/Emoryjob/metadata/barcode_for_vivek.txt' 
INTO TABLE barcode10
FIELDS TERMINATED BY '\t' 
ENCLOSED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS;

DROP TABLE IF EXISTS barcode10;
create table barcode10 (Barcode int NOT NULL, Proteomics_BM_10_Barcode int default NULL, PRIMARY KEY (Barcode));

create table clinical (individualIdentifier varchar(100) NOT NULL,	
PMI int(11) NOT NULL, 
RACE varchar(2) NOT NULL,	
AOD varchar(10) NOT NULL,	
CDR double(10,2), 	
SEX varchar(2),
NP_1 int(4),
PlaqueMean double(10,2),  	
bbscore int(5) DEFAULT NULL,	
Apo1 int(5) DEFAULT NULL, 
Apo2 int(5) DEFAULT NULL,
PRIMARY KEY (individualIdentifier, PMI));

LOAD DATA LOCAL INFILE '/home/maohuaxie/Emoryjob/metadata/MSBB_clinical.csv' 
INTO TABLE clinical
FIELDS TERMINATED BY ',' 
ENCLOSED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS;

DROP TABLE IF EXISTS clinical;

create table clinical (individualIdentifier varchar(100) NOT NULL,	
PMI int(11) NOT NULL, 
RACE varchar(2) NOT NULL,	
AOD varchar(10) NOT NULL,	
CDR double(10,2), 	
SEX varchar(2),
NP_1 int(4),
PlaqueMean double(10,2),  	
bbscore int(5) DEFAULT NULL,	
Apo1 int(5) DEFAULT NULL, 
Apo2 int(5) DEFAULT NULL,
PRIMARY KEY (individualIdentifier, PMI));

create table msswes (
sampleIdentifier varchar(50) NOT NULL,
fileName varchar(50) NOT NULL,
synapseId varchar(50) NOT NULL,
Barcode int default NULL,
individualIdentifier varchar(50) NOT NULL,
BrodmannArea varchar(20) NOT NULL,
PRIMARY KEY (individualIdentifier));

LOAD DATA LOCAL INFILE '/home/maohuaxie/Emoryjob/metadata/MSBB_WES_covariates.csv' 
INTO TABLE msswes
FIELDS TERMINATED BY ',' 
ENCLOSED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS;
	
create table mssbdemo(
CaseID varchar(50) NOT NULL,
ID2 int (10) NOT NULL,
PMI int (10) NOT NULL,
AOD varchar(10) NOT NULL,
DxByPath  int (10) NOT NULL,
CDR double(10,3), 
SEX varchar(2),
NP1 int(4),
PlaqueMean	double(10,3),
MFG double(10,3),
STG double(10,3),	
Batch int(4) NOT NULL,
Remove_Noam varchar(10) NOT NULL,
Combat_Group varchar(10) NOT NULL,
PRIMARY KEY (PMI));
	
LOAD DATA LOCAL INFILE '/home/maohuaxie/Emoryjob/metadata/msbb.demo.txt' 
INTO TABLE mssbdemo
FIELDS TERMINATED BY '\t' 
ENCLOSED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS;	
	
1/26/2018
install shrimp for miRNA annotation: 
http://compbio.cs.toronto.edu/shrimp/releases/SHRiMP_2_2_3.src.tar.gz


1/27/2018
78  tail -n +31 ids.MSBB.batch_2.txt > ids.MSBB.batch_2m.txt 
79  less ids.MSBB.batch_2m.txt 
81  cp makepdbstep2b2.sh makepdbstep2bm2.sh
82  nano makepdbstep2bm2.sh 
   
makepdbstep2bm2.sh

#!/bin/sh
for file in $(cat ./ids.MSBB.batch_2m.txt);
do
GenPro_make_perprotdb2.pl -i $file -d db2 -r /mnt/22qFSS/data2/GenPro/hg19/refProt -o Json2;
done

qsub -cwd -pe smp 24 makepdbstep2bm21.sh


miRbase.sh
#!/bin/sh
set -euo pipefail
wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz   
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.zip    
wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.zip
wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3 
wget ftp://mirbase.org/pub/mirbase/CURRENT/miFam.dat.zip

gunzip hairpin.fa.gz
unzip mature.fa.zip
grep sapiens mature.fa |wc  　# 2588 
grep sapiens hairpin.fa |wc  # 1881 
# Homo sapiens
perl -alne '{if(/^>/){if(/Homo/){$tmp=1}else{$tmp=0}};next if $tmp!=1;s/U/T/g if !/>/;print }' hairpin.fa >hairpin.human.fa
perl -alne '{if(/^>/){if(/Homo/){$tmp=1}else{$tmp=0}};next if $tmp!=1;s/U/T/g if !/>/;print }' mature.fa >mature.human.fa


/*wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz   ##　28645　reads
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.zip   ##   35828 reads 
wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.zip
wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3 ##
wget ftp://mirbase.org/pub/mirbase/CURRENT/miFam.dat.zip
gunzip hairpin.fa.gz
unzip mature.fa.zip
grep sapiens mature.fa |wc  　# 2588 
grep sapiens hairpin.fa |wc  # 1881 
# Homo sapiens
perl -alne '{if(/^>/){if(/Homo/){$tmp=1}else{$tmp=0}};next if $tmp!=1;s/U/T/g if !/>/;print }' hairpin.fa >hairpin.human.fa
perl -alne '{if(/^>/){if(/Homo/){$tmp=1}else{$tmp=0}};next if $tmp!=1;s/U/T/g if !/>/;print }' mature.fa >mature.human.fa */



1/29/2018

adapter.fa
>TruSeq3_smRNA_IndexAdapter
TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC
>TruSeq3_smRNA_Universal
GATCGTCGGACTGTAGAACTCTGAACGTGTAGA

trimmomatic SE -threads 4 -phred33 S1_S1_L001_R1_001.fastq.gz S1_S1_L001_R1_001_miRNAtrimmed.fastq.gz ILLUMINACLIP:/home/emorycpl/project/miRNA/adapter.fa:2:30:10

trimmomatic SE -threads 4 -phred33 S2_S2_L001_R1_001.fastq.gz S2_S2_L001_R1_001_miRNAtrimmed.fastq.gz ILLUMINACLIP:/home/emorycpl/project/miRNA/adapter.fa:2:30:10

trimmomatic SE -threads 4 -phred33 S3_S3_L001_R1_001.fastq.gz S3_S3_L001_R1_001_miRNAtrimmed.fastq.gz ILLUMINACLIP:/home/emorycpl/project/miRNA/adapter.fa:2:30:10

trimmomatic SE -threads 4 -phred33 S4_S4_L001_R1_001.fastq.gz S4_S4_L001_R1_001_miRNAtrimmed.fastq.gz ILLUMINACLIP:/home/emorycpl/project/miRNA/adapter.fa:2:30:10

find . -name \*.gz -exec cp {} test \;

gmapper-ls --strata -o 1 --qv-offset 33 -Q -L /illumina/local/downloads/miRBase/miRBase_mature_hsa-ls - -N 14 --mode mirna -w 170% -E --un 1_R-110_unmapped-miRNA.fastq | samtools view -Su - | samtools sort -m 12000000000 - 1_R-110_miRNAtrimmed_to_MB_mt_sorted


ls *fastq.gz |while read id
do
echo $id
trimmomatic SE -threads 4 -phred33 $id ${id%%.*}_miRNAtrimmed.fq.gz ILLUMINACLIP:/home/emorycpl/project/miRNA/adapter.fa:2:30:10;
done


1/30/2018

ls *fastq.gz |while read id
do
echo $id
java -Xmx2g -jar /sw/hgcc/Pkgs/Trimmomatic/0.36/trimmomatic-0.36.jar SE -threads 4 -phred33 $id /home/mxie/Projects/miRNA/analysis/data/trimmeddata/${id%%.*}_miRNAtrimmed.fq.gz ILLUMINACLIP:/home/mxie/Projects/WGS/analysis2/adaptor.fa:2:30:10;
done


java -Xmx2g -jar /sw/hgcc/Pkgs/Trimmomatic/0.36/trimmomatic-0.36.jar SE -threads 4 -phred33 S1_S1_L001_R1_001.fastq.gz S1_S1_L001_R1_001_miRNAtrimmed.fastq.gz ILLUMINACLIP:/home/mxie/Projects/WGS/analysis2/adaptor.fa:2:30:10

java -Xmx2g -jar /sw/hgcc/Pkgs/Trimmomatic/0.36/trimmomatic-0.36.jar SE -threads 4 -phred33 S2_S2_L001_R1_001.fastq.gz S2_S2_L001_R1_001_miRNAtrimmed.fastq.gz ILLUMINACLIP:/home/mxie/Projects/WGS/analysis2/adaptor.fa:2:30:10

java -Xmx2g -jar /sw/hgcc/Pkgs/Trimmomatic/0.36/trimmomatic-0.36.jar SE -threads 4 -phred33 S3_S3_L001_R1_001.fastq.gz S3_S3_L001_R1_001_miRNAtrimmed.fastq.gz ILLUMINACLIP:/home/mxie/Projects/WGS/analysis2/adaptor.fa:2:30:10

java -Xmx2g -jar /sw/hgcc/Pkgs/Trimmomatic/0.36/trimmomatic-0.36.jar SE -threads 4 -phred33 S4_S4_L001_R1_001.fastq.gz S4_S4_L001_R1_001_miRNAtrimmed.fastq.gz ILLUMINACLIP:/home/mxie/Projects/WGS/analysis2/adaptor.fa:2:30:10



export SHRIMP_FOLDER=$PWD
$SHRIMP_FOLDER/utils/project-db.py --seed 00111111001111111100,00111111110011111100,00111111111100111100,00111111111111001100,00111111111111110000  --h-flag --shrimp-mode ls /home/emorycpl/project/miRNA/miRbase/mature.human.fa

$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/emorycpl/project/miRNA/miRbase/mature.human-ls S1_S1_L001_R1_001_miRNAtrimmed.fastq.gz  --qv-offset 33   \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 8  >map.out 2>map.log

$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/emorycpl/project/miRNA/miRbase/mature.human-ls S2_S2_L001_R1_001_miRNAtrimmed.fastq.gz  --qv-offset 33   \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 8  >map2.out 2>map2.log

$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/emorycpl/project/miRNA/miRbase/mature.human-ls S3_S3_L001_R1_001_miRNAtrimmed.fastq.gz  --qv-offset 33   \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 8  >map3.out 2>map3.log

$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/emorycpl/project/miRNA/miRbase/mature.human-ls S4_S4_L001_R1_001_miRNAtrimmed.fastq.gz  --qv-offset 33   \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 8  >map4.out 2>map4.log

export SHRIMP_FOLDER=$/sw/hgcc/Pkgs/SHRiMP/2.2.3
$SHRIMP_FOLDER/utils/project-db.py --seed 00111111001111111100,00111111110011111100,00111111111100111100,00111111111111001100,00111111111111110000  --h-flag --shrimp-mode ls /home/mxie/Projects/WGS/analysis2/miRbase/mature.human.fa


ls /sw/hgcc/Pkgs/SHRiMP/2.2.3

export SHRIMP_FOLDER=$PWD
$SHRIMP_FOLDER/utils/project-db.py --seed 00111111001111111100,00111111110011111100,00111111111100111100,00111111111111001100,00111111111111110000  --h-flag --shrimp-mode ls /home/emorycpl/project/miRNA/miRbase/hairpin.human.fa

$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/emorycpl/project/miRNA/miRbase/hairpin.human-ls S1_S1_L001_R1_001_miRNAtrimmed.fastq.gz  --qv-offset 33   \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 8  >maph.out 2>maph.log

$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/emorycpl/project/miRNA/miRbase/hairpin.human-ls S2_S2_L001_R1_001_miRNAtrimmed.fastq.gz  --qv-offset 33   \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 8  >map2h.out 2>map2h.log

$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/emorycpl/project/miRNA/miRbase/hairpin.human-ls S3_S3_L001_R1_001_miRNAtrimmed.fastq.gz  --qv-offset 33   \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 8  >map3h.out 2>map3h.log

$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/emorycpl/project/miRNA/miRbase/hairpin.human-ls S4_S4_L001_R1_001_miRNAtrimmed.fastq.gz  --qv-offset 33   \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 8  >map4h.out 2>map4h.log


#gmapper-ls --strata -o 1 --qv-offset 33 -Q -L /home/emorycpl/project/miRNA/miRbase/mature.human-ls -N 14 --mode mirna -w 170% -E --un S1_S1_L001_R1_001_miRNAtrimmed.fastq.gz  >maphh.out 2>maphh.log

export SHRIMP_FOLDER=$PWD

$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/emorycpl/project/miRNA/miRbase/mature.human-ls S1_S1_L001_R1_001_miRNAtrimmed.fq.gz  --qv-offset 33   \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 14 --mode mirna -w 170%  >maphh.out 2>maphh.log

$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/emorycpl/project/miRNA/miRbase/mature.human-ls S1_S1_L001_R1_001_miRNAtrimmed_miRNAtrimmed.fq.gz  --qv-offset 33   \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 8  >maphhh.out 2>maphhh.log

1/31/2018

samtools view  -SF 4 map.out |perl -alne '{$h{$F[2]}++}END{print "$_\t$h{$_}" foreach sort keys %h }'  > map.mature.counts

export SHRIMP_FOLDER=$PWD
ls *fq.gz |while read id
do
echo $id
$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/mxie/Projects/miRNA/analysis/miRbase/mature.human-ls $id --qv-offset 33   \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 14 --mode mirna -w 170%  >${id%%.*}.out 2>${id%%.*}.log;
done

2/1/2017
qrsh
module avail
module load samtools

nano counting.sh
#!/bin/sh
set -euo pipefail
module avail
module load samtools
ls *miRNAtrimmed.out | while read id
do
echo $id
samtools view  -SF 4 $id |perl -alne '{$h{$F[2]}++}END{print "$_\t$h{$_}" foreach sort keys %h }'  > ${id%%.*}.mature.counts;
done

2/2/2017




trimmomatic SE -threads 4 -phred33 --trimlog logfile S1_S1_L001_R1_001.fastq.gz S1_S1_L001_R1_001_miRNAtrimmed.fq.gz ILLUMINACLIP:/home/emorycpl/project/miRNA/adapter.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:14


2/3/2018 miRNA pipeline


find . -name \*.gz -exec cp {} ../test \;

nano removeadapter.sh
#!/bin/sh
set -euo pipefail
ls *fastq.gz |while read id
do
echo $id
trimmomatic SE -threads 4 -phred33 $id ${id%%.*}_miRNAtrimmed.fq.gz ILLUMINACLIP:/home/emorycpl/project/miRNA/adapter.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:14;
done


ls *fastq.gz |while read id
do
echo $id
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC -a ATCTCGTATGCCGTCTTCTGCTTG -e 0.1 -O 5 -m 15 -o ${id%%.*}_adaprm.fastq.gz $id > ${id%%.*}_adaprm.log
done





ls *adaprm.fastq.gz |while read id
do
echo $id
fastqc $id
done


ls *adaprm.fastq.gz |while read id
do
echo $id
$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/emorycpl/project/miRNA/miRbase/mature.human-ls $id --qv-offset 33   \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 14 --mode mirna -w 170%  >${id%%.*}_adaprm_mapping.out 2>${id%%.*}_adaprm_maping.log;
done


ls *adaprm_mapping.out | while read id
do
echo $id
samtools view  -SF 4 $id |perl -alne '{$h{$F[2]}++}END{print "$_\t$h{$_}" foreach sort keys %h }'  > ${id%%.*}.mature.counts;
done




cutadapt -a ADAPTER input.fastq | cutadapt -u 4 -u -4 - > output.fastq

cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC -a ATCTCGTATGCCGTCTTCTGCTTG -e 0.1 -O 5 -m 15 -o S1_S1_L001_R1_001_adaprm.fastq.gz S1_S1_L001_R1_001.fastq.gz > S1_S1_L001_R1_001.log

cutadapt -u 4 -a NNNNADAPTER -o output.fastq input.fastq

cutadapt -u 4 -u -4 -o S1_S1_L001_R1_001_adaprm_cut.fastq.gz S1_S1_L001_R1_001_adaprm.fastq.gz


gunzip -c S1_S1_L001_R1_001_adaprm.fastq.gz| awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c

gunzip -c S1_S1_L001_R1_001_miRNAtrimmed_cut.fq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c


nano mapping.sh
#!/bin/sh
set -euo pipefail
cd /home/emorycpl/src/SHRiMP/SHRiMP_2_2_3
export SHRIMP_FOLDER=$PWD
cd /home/emorycpl/project/miRNA/test2

$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/emorycpl/project/miRNA/miRbase/mature.human-ls  S4_S4_L001_R1_001_adaprm_filtered1.fastq --qv-offset 33   \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 14 --mode mirna -w 170%  > S4_S4_L001_R1_001_adaprm_filtered1.out 2> S4_S4_L001_R1_001_adaprm_filtered1.log;

cat reads.fastq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > read_length.txt

reads<-read.csv(file="read_length.txt", sep="", header=FALSE)
plot (reads$V2,reads$V1,type="l",xlab="read length",ylab="occurences",col="blue")
awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' file.fastq


cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG -e 0.1 -O 5 -m 15 -o S1_S1_L001_R1_001_adaprm.fastq.gz S1_S1_L001_R1_001.fastq.gz > S1_S1_L001_R1_001_cutadaptor.log


cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG -e 0.1 -O 5 -m 15 -o S2_S2_L001_R1_001_adaprm.fastq.gz S2_S2_L001_R1_001.fastq.gz > S2_S2_L001_R1_001_cutadaptor.log


cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG -e 0.1 -O 5 -m 15 -o S3_S3_L001_R1_001_adaprm.fastq.gz S3_S3_L001_R1_001.fastq.gz > S3_S3_L001_R1_001_cutadaptor.log


cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG -e 0.1 -O 5 -m 15 -o S3_S3_L001_R1_001_adaprm.fastq.gz S3_S3_L001_R1_001.fastq.gz > S3_S3_L001_R1_001_cutadaptor.log



2/7/2018




awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 10 && length(seq) <= 30) {print header, seq, qheader, qseq}}' <S4_S4_L001_R1_001_adaprm.fastq >S4_S4_L001_R1_001_adaprm_filtered1.fastq

awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 10 && length(seq) <= 30) {print header, seq, qheader, qseq}}' <S3_S3_L001_R1_001_adaprm.fastq >S3_S3_L001_R1_001_adaprm_filtered1.fastq

awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 10 && length(seq) <= 30) {print header, seq, qheader, qseq}}' <S2_S2_L001_R1_001_adaprm.fastq >S2_S2_L001_R1_001_adaprm_filtered1.fastq

awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 10 && length(seq) <= 30) {print header, seq, qheader, qseq}}' <S1_S1_L001_R1_001_adaprm.fastq >S1_S1_L001_R1_001_adaprm_filtered1.fastq


ls *adaprm_filtered1.fastq |while read id
do
echo $id
$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/emorycpl/project/miRNA/miRbase/mature.human-ls $id --qv-offset 33   \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 14 --mode mirna -w 170%  >${id%%.*}_adaprm_mapping.out 2>${id%%.*}_adaprm_maping.log;
done

ls *.out | while read id
do
echo $id
samtools view  -SF 4 $id |perl -alne '{$h{$F[2]}++}END{print "$_\t$h{$_}" foreach sort keys %h }'  > ${id%%.*}.mature.counts;
done


2/8/2018

cat file*.fastq > bigfile.fastq
cat file*.fastq.gz > bigfile.fastq.gz


less 11_R-152_S50_L004_R1_001.fastq.gz | grep TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC

cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC -a ATCTCGTATGCCGTCTTCTGCTTG -e 0.1 -O 5 -m 15 -o 11_R-152_S50_L004_R1_001_adaprm.fastq.gz 11_R-152_S50_L004_R1_001.fastq.gz > 11_R-152_S50_L004_R1_001.log


awk ‘NR%4==2’ in.fastq|sort |uniq -c|awk ‘{print $1″\t”$2}’ > out




2/9/2018


nano qc.sh
#!/bin/sh
set -euo pipefail
module load FastQC/0.11.4
ls *adaprm.fastq.gz |while read id
do
echo $id
fastqc $id
done



nano mapping.sh

#!/bin/sh
set -euo pipefail
cd /sw/hgcc/Pkgs/SHRiMP/2.2.3
export SHRIMP_FOLDER=$PWD
cd /home/mxie/Projects/miRNA/analysis/cutadapt
ls *adaprm.fastq.gz |while read id
do
echo $id
$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/mxie/Projects/miRNA/analysis/miRbase/mature.human-ls $id --qv-offset 33   \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 14 --mode mirna -w 170%  >${id%%.*}_adaprm_mapping.out 2>${id%%.*}_adaprm_maping.log;
done



nano counting.sh
#!/bin/sh
set -euo pipefail
module load samtools/1.5
ls *adaprm_mapping.out | while read id
do
echo $id
samtools view  -SF 4 $id |perl -alne '{$h{$F[2]}++}END{print "$_\t$h{$_}" foreach sort keys %h }'  > ${id%%.*}.mature.counts;
done

module load Anaconda3/4.2.0
nano cutadapter.sh
#!/bin/sh
set -euo pipefail
ls *fastq.gz |while read id
do
echo $id
python /sw/hgcc/Pkgs/Anaconda3/4.2.0/pkgs/cutadapt-1.15-py35_0/bin/cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC -a ATCTCGTATGCCGTCTTCTGCTTG -e 0.1 -O 5 -m 15 -o ${id%%.*}_adaprm.fastq.gz $id > ${id%%.*}_adaprm.log
done


2/12/2018
nano makepdbstep1b3.sh
#!/bin/sh
set -euo pipefail
GenPro_make_perprotdb1.pl -v -s MSBB.snp -l /mnt/22qFSS/data2/GenPro/hg19/idx -n hg19 -o /home/mxie/Projects/WGS/analysis1/db3  -w /home/mxie/Projects/WGS/analysis1/ids.MSBB.batch_3.txt

qsub -cwd -pe smp 24 makepdbstep1b3.sh

nano makepdbstep2b3.sh

#!/bin/sh
set -euo pipefail
for file in $(cat ./ids.MSBB.batch_3.txt);
do
GenPro_make_perprotdb2.pl -i $file -d db3 -r /mnt/22qFSS/data2/GenPro/hg19/refProt -o Json3;
done

qsub -cwd -pe smp 24 makepdbstep2b3.sh


nano makepdbstep1b4.sh
#!/bin/sh
set -euo pipefail
GenPro_make_perprotdb1.pl -v -s MSBB.snp -l /mnt/22qFSS/data2/GenPro/hg19/idx -n hg19 -o /home/mxie/Projects/WGS/analysis1/db4  -w /home/mxie/Projects/WGS/analysis1/ids.MSBB.batch_4.txt

qsub -cwd -pe smp 24 makepdbstep1b4.sh


2/13/2018

for ((i=204;i<=209;i++)) ;do wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP017/SRP017311/SRR620$i/SRR620$i.sra;done
ls *sra |while read id; do fastq-dump --split-3 $id;done



2/14/2018


sourcedir=/home/emorycpl/project/miRNA/pastsmaple
destdir=/home/emorycpl/project/miRNA/combinedsample

for f in $sourcedir/*(L\d\d\d)_R1_001.fastq.gz;
do
file=$(basename $id )
sample=${file%%.*}
echo $sample
zcat $id | gzip >$destdir/${sample}.fastq.gz; 
done	
	

	
	
	
sourcedir=/home/emorycpl/project/miRNA/pastsmaple
destdir=/home/emorycpl/project/miRNA/combinedsample

ls $sourcedir/*L*_R1_001.fastq.gz |while read id
do
file=$(basename $id )
sample=${file%%.*}
echo $id
echo $file
echo $sample
#zcat $id | gzip >$destdir/${sample}.fastq.gz 
done	
	
	
$file_name =~ /^(\d+)_R*S_(L\d\d\d)_(R1|R2)_(\d\d\d).fastq.gz$/.	
	
	
for id in 10_R-150_S49 
	

	
	
zcat 10_R-150_S49_L004_R1_001.fastq.gz 10_R-150_S49_L005_R1_001.fastq.gz 10_R-150_S49_L006_R1_001.fastq.gz | gzip > 10_R-150_S49.fastq.gz
		
zcat 11_R-152_S50_L004_R1_001.fastq.gz 11_R-152_S50_L005_R1_001.fastq.gz 11_R-152_S50_L006_R1_001.fastq.gz | gzip > 11_R-152_S50.fastq.gz	
	
zcat 12_R-156_S51_L004_R1_001.fastq.gz 12_R-156_S51_L005_R1_001.fastq.gz 12_R-156_S51_L006_R1_001.fastq.gz | gzip > 12_R-156_S51.fastq.gz

zcat 15_R-160_S52_L004_R1_001.fastq.gz 15_R-160_S52_L005_R1_001.fastq.gz 15_R-160_S52_L006_R1_001.fastq.gz | gzip > 15_R-160_S52.fastq.gz

zcat 40_R-243_S62_L004_R1_001.fastq.gz 40_R-243_S62_L005_R1_001.fastq.gz 40_R-243_S62_L006_R1_001.fastq.gz | gzip > 40_R-243_S62.fastq.gz


for name in *.fastq.gz; do
    rnum=${name##*_}
    rnum=${rnum%%.*}

    sample=${name#*_}
    sample=${sample%%_*}
    echo $sample
    #cat "$name" >>"${sample}_$rnum.fastq.gz"
done

R-150
R-150
R-150


for name in *.fastq.gz; do
    rnum=${name##*_}
    rnum=${rnum%%.*}

    sample=${name##*_}
    sample=${sample%%_*}
    echo $sample
    #cat "$name" >>"${sample}_$rnum.fastq.gz"
done
001.fastq.gz
001.fastq.gz


for name in *.fastq.gz; do
    rnum=${name##*_}
    rnum=${rnum%%.*}

    sample=${name#*_}
    sample=${sample%_*}
    echo $sample
    #cat "$name" >>"${sample}_$rnum.fastq.gz"
done

R-150_S49_L004_R1
R-150_S49_L005_R1

for name in *.fastq.gz; do
    rnum=${name##*_}
    rnum=${rnum%%.*}

    sample=${name#*_}
    sample=${sample%%_*}
    echo $sample
    #cat "$name" >>"${sample}_$rnum.fastq.gz"
done

~% FILE="example.tar.gz"
~% echo "${FILE%%.*}"
example
~% echo "${FILE%.*}"
example.tar
~% echo "${FILE#*.}"
tar.gz
~% echo "${FILE##*.}"
gz


cd /home/emorycpl/src/SHRiMP/SHRiMP_2_2_3
export SHRIMP_FOLDER=$PWD
cd /home/emorycpl/project/miRNA/pastsmaple/mergedfile


2/15/2018
nano parseid.py
import os
import csv
os.chdir("/home/mxie/Projects/miRNA/analysis/Past86samples")
filelist=[]
for f in os.listdir():
    
    f_name,f_ext = os.path.splitext(f)
    f1 = f_name.split('_')[0]
    f2 = f_name.split('_')[1]
    f3 = f_name.split('_')[2]
    newname= "_".join([f1,f2,f3])
    filelist.append(newname)
cleanlist = []
[cleanlist.append(x) for x in filelist if x not in cleanlist]
with open("/home/mxie/Projects/miRNA/analysis/pastsample/cleanlist.csv",'w') as resultFile:
    wr = csv.writer(resultFile, dialect='excel')
    wr.writerow(cleanlist)
with open ("/home/mxie/Projects/miRNA/analysis/pastsample/cleanlist.txt","w")as file:
   file.write('\n'.join(cleanlist))    

nano mergefastqfile.sh
#!/bin/sh
set -euo pipefail   
for file in $(cat /home/mxie/Projects/miRNA/analysis/pastsample/cleanlist.txt);
do
zcat $file* | pigz > /home/mxie/Projects/miRNA/analysis/pastsample/$file.fastq.gz; 
done   
   
qsub -cwd -pe smp 24 mergefastqfile.sh
   
   
   
for file in $(cat /home/emorycpl/project/miRNA/pastsmaple/parseidfile/cleanlist.txt);
do
echo $file
zcat $file* | gzip > /home/emorycpl/project/miRNA/pastsmaple/parseidforsmaples/mergedfile/$file.fastq.gz; 
done
   
    
import os
import csv
os.chdir("/home/emorycpl/project/miRNA/pastsmaple/parseidforsmaples")
filelist=[]
for f in os.listdir():
    
    f_name,f_ext = os.path.splitext(f)
    f1 = f_name.split('_')[0]
    f2 = f_name.split('_')[1]
    f3 = f_name.split('_')[2]
    newname= "_".join([f1,f2,f3])
    filelist.append(newname)
cleanlist = []
[cleanlist.append(x) for x in filelist if x not in cleanlist]
with open("/home/emorycpl/project/miRNA/pastsmaple/parseidfile/cleanlist.csv",'w') as resultFile:
    wr = csv.writer(resultFile, dialect='excel')
    wr.writerow(cleanlist)
with open ("/home/emorycpl/project/miRNA/pastsmaple/parseidfile/cleanlist.txt","w")as file:
   file.write('\n'.join(cleanlist)   




2/16/2018

module avail 
module load Anaconda2/4.2.0

nano cutadaptor.sh
#!/bin/sh
set -euo pipefail

ls *fastq.gz |while read id
do
echo $id
/sw/hgcc/Pkgs/Anaconda2/4.2.0/bin/cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC -a ATCTCGTATGCCGTCTTCTGCTTG -e 0.1 -O 5 -m 15 -o /home/mxie/Projects/miRNA/analysis/Data2repeatmerged/${id%%.*}_adaprm.fastq.gz $id > /home/mxie/Projects/miRNA/analysis/Data2repeatmerged/${id%%.*}_adaprm.log
done   
  
qsub -cwd -pe smp 24 cutadaptor.sh


2/19/2018

nano fastqc.sh

#!/bin/sh
set -euo pipefail
ls *adaprm.fastq.gz |while read id
do
echo $id
/sw/hgcc/Pkgs/FastQC/0.11.4/fastqc $id
done

qsub -cwd -pe smp 24 fastqc.sh

nano mapping.sh
#!/bin/sh
set -euo pipefail
cd /sw/hgcc/Pkgs/SHRiMP/2.2.3
export SHRIMP_FOLDER=$PWD
cd 
/home/mxie/Projects/miRNA/analysis/Data2repeatmerged
ls *adaprm.fastq.gz |while read id
do
echo $id
$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/mxie/Projects/miRNA/analysis/miRbase/mature.human-ls $id  \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 14 --mode mirna -w 170%  >${id%%.*}.sam 2>${id%%.*}.log;
done

qsub -cwd -pe smp 24 mapping.sh


2/20/2018

module load samtools/1.5

nano counting.sh

#!/bin/sh
set -euo pipefail
ls  *.out| while read id
do
echo $id
/sw/hgcc/Pkgs/samtools/1.5/bin/samtools view  -SF 4 $id |perl -alne '{$h{$F[2]}++}END{print "$_\t$h{$_}" foreach sort keys %h }'  > ${id%%.*}.mature.counts;
done

qsub -cwd -pe smp 24 counting.sh


module avail 
module load Anaconda2/4.2.0


nano cutadaptor.sh

#!/bin/sh
set -euo pipefail

ls *fastq.gz |while read id
do
echo $id
/sw/hgcc/Pkgs/Anaconda2/4.2.0/bin/cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC -a ATCTCGTATGCCGTCTTCTGCTTG -e 0.1 -O 5 -m 15 -o /home/mxie/Projects/miRNA/analysis/pastcutadaptsamples/${id%%.*}_adaprm.fastq.gz $id > /home/mxie/Projects/miRNA/analysis/Data2repeatmerged/${id%%.*}_adaprm.log
done   
  
qsub -cwd -pe smp 24 cutadaptor.sh



2/21/2018

nano comapping.sh
#!/bin/sh
set -euo pipefail
cd /sw/hgcc/Pkgs/SHRiMP/2.2.3
export SHRIMP_FOLDER=$PWD
cd /home/mxie/Projects/miRNA/analysis/pastcutadaptsamples
ls *adaprm.fastq.gz |while read id
do
echo $id
$SHRIMP_FOLDER/bin/gmapper-ls -L  /home/mxie/Projects/miRNA/analysis/miRbase/mature.human-ls $id \
-o 1 -H -E -a -1 -q -30 -g -30 --qv-offset 33 --strata -N 14 --mode mirna -w 170%  >/home/mxie/Projects/miRNA/analysis/pastcutadaptsamples/mappingresult/${id%%.*}.sam 2>/home/mxie/Projects/miRNA/analysis/pastcutadaptsamples/mappingresult/${id%%.*}.log;
done

qsub -cwd -pe smp 24 comapping.sh


nano combinedsam.sh

#!/bin/sh
set -euo pipefail
cd /sw/hgcc/Pkgs/SHRiMP/2.2.3
export SHRIMP_FOLDER=$PWD
cd /home/mxie/Projects/miRNA/analysis/pastcutadaptsamples/mappingresult
 
for id in $(cat /home/mxie/Projects/miRNA/analysis/pastsample/cleanlist.txt);
do
$SHRIMP_FOLDER/utils/merge-hits-diff-qr-same-db.py $id*_adaprm.sam > ./mergedout/$id.out;
done   
   
qsub -cwd -pe smp 24 combinedsam.sh


2/22/2018

module load samtools/1.5

nano counting.sh

#!/bin/sh
set -euo pipefail
ls  *.out| while read id
do
echo $id
/sw/hgcc/Pkgs/samtools/1.5/bin/samtools view  -SF 4 $id |perl -alne '{$h{$F[2]}++}END{print "$_\t$h{$_}" foreach sort keys %h }'  > ${id%%.*}.mature.counts;
done

qsub -cwd -pe smp 24 counting.sh


nano lengthdistribution.sh

#!/bin/sh
set -euo pipefail
ls *_adaprm.fastq.gz |while read id
do
echo $id
gunzip -c $id | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c >${id%%.*}_length.txt;
done

qsub -cwd -pe smp 24 lengthdistribution.sh


nano mappinglenthdistribution.sh

#!/bin/sh
set -euo pipefail
ls  *_adaprm.sam| while read id
do
echo $id
/sw/hgcc/Pkgs/samtools/1.5/bin/samtools view $id | gawk '{mlarr[length($10)]++}END{for (len = 1; len <=101; len++) printf("%d\t%d\n",len,mlarr[len])}' | sort -k 1n,1n  > ${id%%.*}_mapped_length_distribution_data.txt
done

qsub -cwd -pe smp 24 mappinglenthdistribution.sh


nano makepdbstep2b3m.sh

#!/bin/sh
set -euo pipefail
for file in $(cat ./ids.MSBB.batch_3m.txt);
do
GenPro_make_perprotdb2.pl -i $file -d db3 -r /mnt/22qFSS/data2/GenPro/hg19/refProt -o Json3;
done

qsub -cwd -pe smp 24 makepdbstep2b3m.sh


2/27/2018

nano removeadaptor.sh

#!/bin/sh
set -euo pipefail
ls *fastq.gz |while read id
do
echo $id
java -Xmx2g -jar /sw/hgcc/Pkgs/Trimmomatic/0.36/trimmomatic-0.36.jar SE -phred33 $id  ./trimmed/${id%%.*}_miRNAtrimmed.fastq.gz ILLUMINACLIP:/home/mxie/Projects/miRNA/analysis/illumina_RGA_3p_adapter.fa:2:30:10;
done   

qsub -cwd -pe smp 24 removeadaptor.sh


nano mapping.sh

#!/bin/sh
set -euo pipefail
cd /sw/hgcc/Pkgs/SHRiMP/2.2.3
export SHRIMP_FOLDER=$PWD
cd /home/mxie/Projects/miRNA/analysis/pastsample/trimmed
ls *miRNAtrimmed.fastq.gz |while read id
do
echo $id
$SHRIMP_FOLDER/bin/gmapper-ls --strata -o 1 --qv-offset 33 -Q -L /home/mxie/Projects/miRNA/analysis/miRbase/mature.human-ls $id \
-N 14 --mode mirna -w 170% -E --un ${sid}_unmapped-miRNA.fastq \
 >${id%%.*}_mapping.sam 2>${id%%.*}_maping.log;
done

qsub -cwd -pe smp 24 mapping.sh


nano mapping.sh

#!/bin/sh
set -euo pipefail
cd /sw/hgcc/Pkgs/SHRiMP/2.2.3
export SHRIMP_FOLDER=$PWD
cd /home/mxie/Projects/miRNA/analysis/Past86samples/trimmed
for tfq in *_L001_R1_001_miRNAtrimmed.fastq.gz; do sid=${tfq%%_S*}; spx=${tfq%%_L00*}; 
if [ ! -e ./trimmed/${sid}_mapping.sam ]; 
then cat ${spx}_L00?_R1_001_miRNAtrimmed.fastq.gz | $SHRIMP_FOLDER/bin/gmapper-ls --strata -o 1 --qv-offset 33 -Q -L /home/mxie/Projects/miRNA/analysis/miRbase/mature.human-ls \
-N 14 --mode mirna -w 170% -E --un ./trimmed/${sid}_unmapped-miRNA.fastq \
./trimmed/${sid}_mapping.sam 2>./trimmed/${sid}_mapping.log; fi; done

qsub -cwd -pe smp 24 mapping.sh


2/28/2018


nano mapping.sh

#!/bin/sh
set -euo pipefail
cd /sw/hgcc/Pkgs/SHRiMP/2.2.3
export SHRIMP_FOLDER=$PWD
cd /home/mxie/Projects/miRNA/analysis/Past86samples/trimmed
for tfq in *_L00?_R1_001_miRNAtrimmed.fastq.gz; do sid=${tfq%%_S*}; spx=${tfq%%_L00*}; 
if [ ! -e ./trimmed/${sid}_miRNAtrimmed_to_MB_mt_sorted.bam ]; 
then cat ${spx}_L00?_R1_001_miRNAtrimmed.fastq.gz | $SHRIMP_FOLDER/bin/gmapper-ls --strata -o 1 --qv-offset 33 -Q -L /home/mxie/Projects/miRNA/analysis/miRbase/mature.human-ls \
-N 14 --mode mirna -w 170% -E --un ./trimmed/${sid}_unmapped-miRNA.fastq | /sw/hgcc/Pkgs/samtools/1.5/bin/samtools view -Su - | /sw/hgcc/Pkgs/samtools/1.5/bin/samtools sort -m 12000000000 - ./trimmed/${sid}_miRNAtrimmed_to_MB_mt_sorted; 
/sw/hgcc/Pkgs/samtools/1.5/bin/samtools index ./trimmed/${sid}_miRNAtrimmed_to_MB_mt_sorted.bam; fi; done

qsub -cwd -pe smp 24 mapping.sh

tfq=10_R-150_S49_L001_R1_001_miRNAtrimmed.fastq.gz
sid=${tfq%%_S*}
$sid
10_R-150
spx=${tfq%%_L00*}
$spx
10_R-150_S49

###
nano removeadaptor.sh

#!/bin/sh
set -euo pipefail
ls *fastq.gz |while read id
do
echo $id
if [ ! -e ./trimmed/${id%%.*}_miRNAtrimmed.fastq.gz ];
then java -Xmx2g -jar /sw/hgcc/Pkgs/Trimmomatic/0.36/trimmomatic-0.36.jar SE -phred33 $id  ./trimmed/${id%%.*}_miRNAtrimmed.fastq.gz ILLUMINACLIP:/home/mxie/Projects/miRNA/analysis/illumina_RGA_3p_adapter.fa:2:30:10; fi; done   

qsub -cwd -pe smp 24 removeadaptor.sh


nano makepdbstep2b3m.sh

#!/bin/sh
set -euo pipefail
for file in $(cat ./ids.MSBB.batch_3m.txt);
do
GenPro_make_perprotdb2.pl -i $file -d db3 -r /mnt/22qFSS/data2/GenPro/hg19/refProt -o Json3;
done

qsub -cwd -pe smp 24 makepdbstep2b3m.sh




3/1/2018

nano counting.sh

#!/bin/sh
set -euo pipefail
ls  *.sam| while read id
do
echo $id
/sw/hgcc/Pkgs/samtools/1.5/bin/samtools view  -SF 4 $id |perl -alne '{$h{$F[2]}++}END{print "$_\t$h{$_}" foreach sort keys %h }'  > ${id%%.*}.mature.counts;
done

qsub -cwd -pe smp 24 counting.sh


nano samtobamsort.sh

#!/bin/sh
set -euo pipefail
ls  *.sam| while read id
do
echo $id
/sw/hgcc/Pkgs/samtools/1.5/bin/samtools view -Su $id > ./trimmed/${id}.bam;
/sw/hgcc/Pkgs/samtools/1.5/bin/samtools sort ./trimmed/${id}.bam > ./trimmed/${id%%.*}_miRNAtrimmed_to_MB_mt_sorted.bam
/sw/hgcc/Pkgs/samtools/1.5/bin/samtools index ./trimmed/${id}_miRNAtrimmed_to_MB_mt_sorted.bam
done

qsub -cwd -pe smp 24 samtobam.sh


/sw/hgcc/Pkgs/samtools/1.5/bin/samtools flagstat 10_R-150_S49_miRNAtrimmed_mapping.sam.bam


ls sed 's/miRNAtrimmed_mapping.sam_miRNAtrimmed_to_MB_mt_sorted/miRNAtrimmed_to_MB_mt_sorted/g'

ls *miRNAtrimmed_mapping.sam_miRNAtrimmed_to_MB_mt_sorted.bam | while read FILE ; do
    newfile="$(echo ${FILE} |sed -e 's/miRNAtrimmed_mapping.sam_miRNAtrimmed_to_MB_mt_sorted/miRNAtrimmed_to_MB_mt_sorted/')" ;
    mv "${FILE}" "${newfile}" ;
done 







## trim the reads
for fqf in *_R1_001.fastq.gz; do fid=${fqf%%_R1_001.fastq.gz}; if [ ! -e ${fid}_R1_001_miRNAtrimmed.fastq.gz ]; then echo ${fid}; qsub -pe threaded 1 -q threaded.q -cwd -N trim_${fid} -b y java -jar /illumina/local/downloads/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 ${fqf} ${fid}_R1_001_miRNAtrimmed.fastq.gz ILLUMINACLIP:/illumina/local/share/illumina/illumina_RGA_3p_adapter.fa:2:30:10; fi; done

# echo -n "length" > trimmed_length_distribution_header.txt; cp length_labels.txt trimmed_length_distribution_data.txt; for tfq in *_R1_001_miRNAtrimmed.fastq.gz ../402_EU_Wingo_86miRNASeq_2/*_R1_001_miRNAtrimmed.fastq.gz; do sbn=`basename ${tfq}`; sid=${sbn%%_R1*}; echo -n "	${sid}" >> trimmed_length_distribution_header.txt; mv trimmed_length_distribution_data.txt tmp_trimmed_length_distribution_data.txt; zcat ${tfq} | gawk '(NR %4) == 2{lenarr[length($0)]++}END{for (len = 1;len <= 101; len++) printf("%d\t%d\n",len,lenarr[len])}' | sort -k 1n,1n | join -t "	" -a 1 tmp_trimmed_length_distribution_data.txt - > trimmed_length_distribution_data.txt; done; echo "" >> trimmed_length_distribution_header.txt; cat trimmed_length_distribution_header.txt trimmed_length_distribution_data.txt > trimmed_length_distribution_table.txt
#
## run all the samples
# for tfq in *_L001_R1_001_miRNAtrimmed.fastq.gz; do sid=${tfq%%_S*}; spx=${tfq%%_L00*}; if [ ! -e ${sid}_miRNAtrimmed_to_MB_mt_sorted.bam ]; then qsub -pe threaded 1 -q threaded.q -cwd -V -N SHRiMP_${sid}_MB_mt <<< "cat ${spx}_L00?_R1_001_miRNAtrimmed.fastq.gz | $SHRIMP_FOLDER/bin/gmapper-ls --strata -o 1 --qv-offset 33 -Q -L /illumina/local/downloads/miRBase/miRBase_mature_hsa-ls - -N 14 --mode mirna -w 170% -E --un ${sid}_unmapped-miRNA.fastq | samtools view -Su - | samtools sort -m 12000000000 - ${sid}_miRNAtrimmed_to_MB_mt_sorted; samtools index ${sid}_miRNAtrimmed_to_MB_mt_sorted.bam"; fi; done
#
## extract counts for all bam files separately
# for tbf in *_miRNAtrimmed_to_MB_mt_sorted.bam; do sid=${tbf%%_miRNAtrimmed_to_MB_mt_sorted.bam}; if [ ! -e ${sid}_mt_subcount.txt ]; then echo ${sid}_mt_subcount.txt; for mirna in `samtools view -h ${tbf} | gawk '$1 == "@SQ"{printf("%s\n",substr($2,4))}' | sort`; do cnt=`samtools view -c ${tbf} ${mirna}`; echo "${mirna}	${cnt}"; done > ${sid}_mt_subcount.txt; fi; done
#
# for scp in ../*/*_mt_subcount.txt; do sid=`basename ${scp} _mt_subcount.txt`; echo ${sid}; done | sort -g > sampleIDList.txt
# gawk '{print $1}' 1_R-110_mt_subcount.txt > subcount_table.txt; echo -n "ID" > subcount_table_header.txt; for sid in `cat sampleIDList.txt`; do echo -n "	${sid}" >> subcount_table_header.txt; mv subcount_table.txt tmp_subcount_table.txt; join -t "	" tmp_subcount_table.txt ../*/${sid}_mt_subcount.txt > subcount_table.txt; done; echo "" >> subcount_table_header.txt; mv subcount_table.txt tmp_subcount_table.txt; cat subcount_table_header.txt tmp_subcount_table.txt > subcount_table.txt
#
## get mapped length distribution
# echo -n "length" > mapped_length_distribution_header.txt; cp length_labels.txt mapped_length_distribution_data.txt; for sid in `cat sampleIDList.txt`; do echo -n "	${sid}" >> mapped_length_distribution_header.txt; mv mapped_length_distribution_data.txt tmp_mapped_length_distribution_data.txt; samtools view ../402_EU_Wingo_86miRNASeq_*/${sid}_miRNAtrimmed_to_MB_mt_sorted.bam | gawk '{mlarr[length($10)]++}END{for (len = 1; len <=101; len++) printf("%d\t%d\n",len,mlarr[len])}' | sort -k 1n,1n | join -t "	" -a 1 tmp_mapped_length_distribution_data.txt - > mapped_length_distribution_data.txt; done; echo "" >> mapped_length_distribution_header.txt; cat mapped_length_distribution_header.txt mapped_length_distribution_data.txt > mapped_length_distribution_table.txt
#
# 















/home/mxie/Projects/miRNA/analysis/Past86samples

nano removeadaptor.sh
#!/bin/sh
set -euo pipefail
for fqf in *_R1_001.fastq.gz; do fid=${fqf%%_R1_001.fastq.gz}; 
if [ ! -e ${fid}_R1_001_miRNAtrimmed.fastq.gz ]; 
then echo ${fid}; 
qsub -cwd -pe smp 24 -N trim_${fid} -b y java -Xmx2g -jar /sw/hgcc/Pkgs/Trimmomatic/0.36/trimmomatic-0.36.jar SE -phred33 ${fqf} ./trimmed/${fid}_R1_001_miRNAtrimmed.fastq.gz ILLUMINACLIP:/home/mxie/Projects/miRNA/analysis/illumina_RGA_3p_adapter.fa:2:30:10; fi; done



echo -n "length" > trimmed_length_distribution_header.txt; 
cp length_labels.txt trimmed_length_distribution_data.txt; 
for tfq in *_R1_001_miRNAtrimmed.fastq.gz ../402_EU_Wingo_86miRNASeq_2/*_R1_001_miRNAtrimmed.fastq.gz; 
do sbn=`basename ${tfq}`; sid=${sbn%%_R1*}; 
echo -n "	${sid}" >> trimmed_length_distribution_header.txt; 
mv trimmed_length_distribution_data.txt tmp_trimmed_length_distribution_data.txt; 
zcat ${tfq} | gawk '(NR %4) == 2{lenarr[length($0)]++}END{for (len = 1;len <= 101; len++) printf("%d\t%d\n",len,lenarr[len])}' | sort -k 1n,1n | join -t "	" -a 1 tmp_trimmed_length_distribution_data.txt - > trimmed_length_distribution_data.txt; 
done; echo "" >> trimmed_length_distribution_header.txt; 
cat trimmed_length_distribution_header.txt trimmed_length_distribution_data.txt > trimmed_length_distribution_table.txt


#!/bin/sh
set -euo pipefail
cd /sw/hgcc/Pkgs/SHRiMP/2.2.3
export SHRIMP_FOLDER=$PWD
cd /home/mxie/Projects/miRNA/analysis/Past86samples/trimmed
for tfq in *_L001_R1_001_miRNAtrimmed.fastq.gz; do sid=${tfq%%_S*}; spx=${tfq%%_L00*}; 
if [ ! -e ${sid}_miRNAtrimmed_to_MB_mt_sorted.bam ]; 
then qsub -cwd -pe smp 24 -V -N SHRiMP_${sid}_MB_mt <<< "cat ${spx}_L00?_R1_001_miRNAtrimmed.fastq.gz | $SHRIMP_FOLDER/bin/gmapper-ls --strata -o 1 --qv-offset 33 -Q -L /home/mxie/Projects/miRNA/analysis/miRbase/mature.human-ls - -N 14 --mode mirna -w 170% -E --un ${sid}_unmapped-miRNA.fastq 
| /sw/hgcc/Pkgs/samtools/1.5/bin/samtools view -Su - | /sw/hgcc/Pkgs/samtools/1.5/bin/samtools sort -m 12000000000 - ${sid}_miRNAtrimmed_to_MB_mt_sorted; /sw/hgcc/Pkgs/samtools/1.5/bin/samtools index ${sid}_miRNAtrimmed_to_MB_mt_sorted.bam"; fi; done

echo -n "length" > mapped_length_distribution_header.txt; 
cp length_labels.txt mapped_length_distribution_data.txt; 
for sid in `cat sampleIDList.txt`; 
do echo -n "	${sid}" >> mapped_length_distribution_header.txt; 
mv mapped_length_distribution_data.txt tmp_mapped_length_distribution_data.txt; 
samtools view ../402_EU_Wingo_86miRNASeq_*/${sid}_miRNAtrimmed_to_MB_mt_sorted.bam | gawk '{mlarr[length($10)]++}END{for (len = 1; len <=101; len++) printf("%d\t%d\n",len,mlarr[len])}' | sort -k 1n,1n | join -t "	" -a 1 tmp_mapped_length_distribution_data.txt - > mapped_length_distribution_data.txt; done; echo "" >> mapped_length_distribution_header.txt; cat mapped_length_distribution_header.txt mapped_length_distribution_data.txt > mapped_length_distribution_table.txt
