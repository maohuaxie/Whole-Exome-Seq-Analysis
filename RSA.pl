#!/usr/local/bin/perl -w
#######################################################################
#        Yingyao Zhou, zhou@gnf.org
#        April 29, 2007 
#        Last modified 6/13/08
#######################################################################
BEGIN{ unshift @INC, ".", "/home/yzhou/cheminfo/yzhou/Motif/" }

use strict;
use Getopt::Std;
use HGPValue;
use Data::Table;
use Data::Dumper;

#======= check input arguments ========
my %options=();
getopts("l:u:f:o:g:w:s:hrR:b", \%options) or usage();
usage() if ($options{h});

# lower bound and upper bound
my ($LB, $UB)=(0, 1);
$LB = $options{l} if (defined $options{l});
$UB = $options{u} if (defined $options{u});
my $outfile = $options{o} if (defined $options{o});
my $fmt = 'UNIX';
$fmt = uc $options{f} if (defined $options{f});
if ($fmt eq 'UNIX') {
  $fmt = 0;
} elsif ($fmt eq 'PC') {
  $fmt = 1;
} else {
  $fmt = 2;
}

my $REVERSE = (defined $options{r})?1:0;
### ZHOU
my $RANDOMIZE = (defined $options{R})?1:0;
###
my $BONFERRONI=(defined $options{b})?1:0;

my $fileName = $ARGV[0];
unless (defined($fileName)) {
  error("Input file is not specified!");
}

sub usage {
  print<<END;
Usage:
  RSA.pl [options] fileName
    -l: lower_bound, defaults to 0
    -u: upper bound, defaults to 1
    -r: reverse hit picking, the higher the score the better
        if -r flag is off, the lower the score the better
    -f: input file format: PC, UNIX or MAC, UNIX by default
    -o: output file name, STDOUT if not specified
### ZHOU
    -R: randomize score
###
    -h: help, this message

    filename: input file must a in CSV format
    the spreadsheet must contain at least three columns:
      Gene_ID: the gene identifier for the well
      Well_ID: the well identifier
      Score: numerical value for hit picking
    You can supply your own column names using the following three options

    -g: column name for gene ID, default "Gene_ID"
    -w: column name for well ID, default "Well_ID"
    -s: column name for score used for sorting, default "Score"
      E.g., RSA.pl -g myGene -w myWell
      will instructure RSA.pl to look for "myGene" instead of "Gene_ID",
      "myWell" instead of "Well_ID",
      but still expecting a column named "Score" (since -s is not used).
    -b: turn on Bonferroni correction, conceptually useful
        when there are different number of siRNAs per gene.

    Notice:
      1) the order of the these three columns can be arbitrary
      2) wells share the same Gene_ID are consider independent siRNAs for the same gene
      2) wells are ignored, if Gene_ID or Score is not defined
  Examples:

    RSA.pl -l 0.2 -u 0.8 -f PC -o output.csv input.csv
      wells with lower scores are considered more active
      wells <=0.2 are guaranteed hits, wells >0.8 are guaranteed non-hits
      wells (0.2,0.8] are determined by OPI algorithm
      input CSV file is a Windows format
      output results in output.csv file

    RSA.pl -l 1.2 -u 2.0 -f PC -r -o output.csv input.csv
      wells with higher scores are considered more active (specified by -r flag)
      wells >=2.0 are guaranteed hits, wells <1.2 are guaranteed non-hits
      wells [1.2,2.0) are determined by OPI algorithm
      input CSV file is a Windows format
      output results in output.csv file
END
  exit;
}

my $hgp = new HGPValue;

#my $t = Data::Table::fromCSV($fileName, 1, undef, {OS=>$fmt});
my $t = Data::Table::fromFile($fileName);

my $geneColumn = 'Gene_ID';
$geneColumn = $options{g} if defined($options{g});
my $wellColumn = 'Well_ID';
$wellColumn = $options{w} if defined($options{w});
my $scoreColumn = 'Score';
$scoreColumn = $options{s} if defined($options{s});

my $idx_gene = $t->colIndex($geneColumn);
my $idx_well = $t->colIndex($wellColumn);
my $idx_score = $t->colIndex($scoreColumn);
error("Missing column named $geneColumn!") if $idx_gene<0;
error("Missing column named $wellColumn!") if $idx_well<0;
error("Missing column named $scoreColumn!") if $idx_score<0;

# filter out undefined wells
$t = $t->match_pattern('defined($_->['.$idx_gene.']) && $_->['.$idx_gene.'] ne "" && defined($_->['.$idx_score.'])');
# Table must contains at least the following three columns
#"Gene_ID","Well_ID","Score"
my $N = $t->nofRow;

my %c_gene = ();
my %c_rank = ();
my %c_score = ();
my $R_logP = array($N);
my $R_bestActivity = array($N);
my $I_hit = array($N,0);
my $I_totWell = array($N);
my $I_hitWell = array($N);

### ZHOU
if ($RANDOMIZE) {
  my @rndScore=();
  for (my $i=0; $i<$N; $i++) {
    $rndScore[$i]=rand;
  }
  my @originScore=$t->col($scoreColumn);
  my $t2=new Data::Table([\@rndScore, \@originScore], ['rndScore', 'Score'],1);
  $t2->sort('rndScore',0,0);
  for (my $i=0; $i<$N; $i++) {
    $t->setElm($i,$scoreColumn,$t2->elm($i,'Score'));
  }
}
###

# build Gene_ID, well rank mapping
$t->sort($scoreColumn, 0, $REVERSE);
for (my $j = 0; $j < $N; $j++) {
  $c_score{$t->elm($j, $scoreColumn)}=$j;
}
for (my $j = 0; $j < $N; $j++) {
  my $s_gene = $t->elm($j, $geneColumn);
  unless (defined($c_gene{$s_gene})) {
    $c_gene{$s_gene}=[];
    $c_rank{$s_gene}=[];
  }
  push @{$c_gene{$s_gene}}, $j;
  # modified to deal with ties in the score, when tie happens, use the max rank
  #print $j, "<>", $c_rank{$t->elm($j, $scoreColumn)}, "\n";
  push @{$c_rank{$s_gene}}, $c_score{$t->elm($j, $scoreColumn)};
}
#print ">>>>>>>>> $N $LB $UB\n";
#print Dumper(%c_gene);
# Running OPI hit picking algorithm
foreach my $s_gene (keys %c_gene) {
  my $I_rank = $c_gene{$s_gene};
  my ($i_max, $i_min) = (undef, undef);
  #print ">>> $s_gene $LB $UB\n";
  #print join(" ", @{$I_rank}), "\n";
  for (my $k = 0; $k < scalar @{$I_rank}; $k++) {
    #print "== ", $t->elm($I_rank->[$k], $scoreColumn), "\n";
    if ($REVERSE) {
      $i_max = $k if $t->elm($I_rank->[$k], $scoreColumn)>=$LB;
      $i_min = $k if $t->elm($I_rank->[$k], $scoreColumn)>=$UB;
      # added in case all values are below $LB
      $i_max = $k-1 if ($t->elm($I_rank->[$k], $scoreColumn)<$LB && !defined($i_max));
    } else {
      $i_max = $k if $t->elm($I_rank->[$k], $scoreColumn)<=$UB;
      $i_min = $k if $t->elm($I_rank->[$k], $scoreColumn)<=$LB;
      # added in case all values are above $UB
      $i_max = $k-1 if ($t->elm($I_rank->[$k], $scoreColumn)>$UB && !defined($i_max));
    }
  }
  #print "> > > $i_min $i_max\n";
  #print Dumper($s_gene, $I_rank, $i_min, $i_max);
  # use c_rank instead of I_rank, c_rank considers ties
  my $rslt = OPIScore($c_rank{$s_gene}, $N, $i_min, $i_max, $BONFERRONI);
  #print Dumper($rslt);
  my $logP = $rslt->[0];
  my $cutoff = $rslt->[1];
  for (my $k = 0; $k < scalar @{$I_rank}; $k++) {
    $R_logP->[$I_rank->[$k]] = $logP;
    $R_bestActivity->[$I_rank->[$k]] = $t->elm($I_rank->[0], $scoreColumn);
    $I_hitWell->[$I_rank->[$k]] = $cutoff;
    $I_totWell->[$I_rank->[$k]] = scalar @{$I_rank};
    $I_hit->[$I_rank->[$k]] = 1 if $k < $cutoff;
  }
}

$t->addCol($R_logP, "LogP");
$t->addCol($R_bestActivity, "BestActivity");
$t->addCol($I_hit, "OPI_Hit");
$t->addCol($I_hitWell, "#hitWell");
$t->addCol($I_totWell, "#totalWell");

# Running simple cutoff algorithm
my $I_OPI_rank = array($N, 999999);
my $I_X_rank = array($N, 999999);
my $I_Exp_rank = array($N, 999999);
$t->addCol($I_OPI_rank, "OPI_Rank");
$t->addCol($I_X_rank, "Cutoff_Rank");
#$t->addCol($I_Exp_rank, "EXP_Rank");

# OPI hit ranking
my $cnt=0;
$t->sort('LogP',0,0,$scoreColumn,0,$REVERSE);
for (my $j=0; $j<$N; $j++) {
  if ($t->elm($j, 'OPI_Hit')>0) {
    $cnt+=1; $t->setElm($j, 'OPI_Rank', $cnt); 
  }
}

# Cutoff hit ranking
$cnt=0;
$t->sort($scoreColumn,0,$REVERSE,'LogP',0,0);
for (my $j=0; $j<$N; $j++) {
  $cnt+=1;
  if ($REVERSE) {
    last if ($t->elm($j, $scoreColumn)<$LB);
  } else {
    last if ($t->elm($j, $scoreColumn)>$UB);
  }
  $t->setElm($j, 'Cutoff_Rank', $cnt);
  #$t->setElm($j, 'EXP_Rank',
  #  ($t->elm($j, 'OPI_Rank') < $t->elm($j, 'Cutoff_Rank'))?
  #    $t->elm($j, 'OPI_Rank') : $t->elm($j, 'Cutoff_Rank')
  #);
}
#$t->sort('EXP_Rank',0,0);
$t->sort('LogP',0,0,'BestActivity',0,$REVERSE,$scoreColumn,0,$REVERSE);
if ($outfile) {
  $t->csv(1, {file=>$outfile});
} else {
  print $t->csv(1);
}
my $idx_opi=$t->colIndex('OPI_Hit');

print ">Total Number of Wells: ".$t->nofRow."\n";
my %gene=();
for (my $i=0; $i<$t->nofRow; $i++) {
  $gene{$t->elm($i, $geneColumn)}=1;
}
print ">Total Number of Genes: ".(scalar keys %gene)."\n";

$t=$t->match_pattern('$_->['.$idx_opi.']==1');
print ">Number of Hit Wells: ".$t->nofRow."\n";
%gene=();
for (my $i=0; $i<$t->nofRow; $i++) {
  $gene{$t->elm($i, $geneColumn)}=1;
}
print ">Number of Hit Genes: ".(scalar keys %gene)."\n";
exit;

# create an array of $n element, with default value to be $default
sub array {
  my ($n, $default) = @_;
  my @arr = ();
  unless (defined($default)) {
    $arr[$n-1] = undef;
  } else {
    for (my $i=0; $i<$n; $i++) {
      $arr[$i] = $default;
    }
  }
  return \@arr;
}

# core of the OPI algorithm
sub OPIScore {
  my ($I_rank, $N, $i_min, $i_max, $BONFERRONI) = @_;
  my ($i, $cutoff, $p_min, $logP);
  $cutoff=0; $p_min=1.0;
  my $n = scalar @{$I_rank};
  $i_max=($n - 1) unless defined($i_max);
  $i_min=0 unless defined($i_min);
  for (my $i=$i_min; $i<=$i_max; $i++) {
    next if (($i<$i_max) && ($I_rank->[$i]==$I_rank->[$i+1]));
    $logP = $hgp->pvalue($N, $I_rank->[$i]+1, $n, $i+1);
    if ($BONFERRONI) {
        $logP*=($i_max-$i_min+1);
    }
    $logP = ($logP < 1e-100)?1e-100:$logP;
    $logP=log($logP)/log(10);
    if ($logP < $p_min) {
      $p_min=sprintf('%.3f', $logP); $cutoff=$i+1;
    }
  }
  return [$p_min, $cutoff];
}

sub warning {
  my $msg = shift;
  print "WARNING> $msg\n";
}

sub error {
  my $msg = shift;
  print "ERROR> $msg\n";
  exit;
}
