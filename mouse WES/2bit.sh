# hg19.fa 是指 hg19的fasta文件。
# 它是很重要的参考文件，像bedtools，bowtie，bwa等，都需要它。

# 一般来说，使用的都是UCSC的hg19文件作为权威的。地址是：
# http://hgdownload.cse.ucsc.edu/downloads.html#human

# 然后下载
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/
# 中的 hg19.2bit

# 要注意的是，下载之后，获得的还不是fasta文件，需要一个程序(twoBitToFa )转化。
# 这个程序仅在linux中可运行。下载地址：
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/

# 之后进入twoBitToFa所在文件夹，运行如下命令：
twoBitToFa input.2bit output.fa   
# 即 
twoBitToFa hg19.2bit hg19.fa

# 就得到了hg19的fasta文件。

# 另外附上TwoBitToFa的其他参数：
# options:
# -seq=name - restrict this to just one sequence
# - - start at given position in sequence (zero-based)
# -end=X - end at given position in sequence (non-inclusive)
# -seqList=file - file containing list of the desired sequence names 
# in the format seqSpec[:start-end], e.g. chr1 or chr1:0-189
# where coordinates are half-open zero-based, i.e. [start,end)
# -noMask - convert sequence to all upper case
# -bpt=index.bpt - use bpt index instead of built in one
# -bed=input.bed - grab sequences specified by input.bed. Will exclude introns

# 最后再送上染色体的大小文件地址，也来自于UCSC.
# http://genome.ucsc.edu/goldenPath/help/hg19.chrom.sizes

please note: 
$ chmod +x twoBitToFa
$ ./twoBitToFa
