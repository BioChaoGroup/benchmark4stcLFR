#!/usr/bin/env perl
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       stat bam from bwa alignment
# Author:            Chao | fangchao@genomics.cn
# Version:           V0.0
# Last modified:    10 Aug 2019 (since 10 Aug 2019)
# ===================================================================
# see detail below
use strict;

my($pID);
while(<>){
  chomp;
  my @s = split;
  next if $s[0] eq $pID;
  my ($mLen,$fLen) = ($s[5] eq "*")?(0,0):&cigarLen($s[5]);
  $s[11] =~ /NM:i:(\d+)/;
  my $NM = ($s[11] =~ /NM:i:(\d+)/)?$1:0;
  my @d = split("_",$s[0]);
  #$d[3] =~ s/flag=//;
  $d[4] =~ s/multi=//; $d[5] =~ s/len=//;
  print "$d[0]\t$d[1]\_$d[2]\t$d[3]\t$d[4]\t$d[5]\t$d[6]\t$d[8]\t$s[2]\t$s[5]\t$mLen\t$fLen\t$NM\n";
  $pID = $s[0];
  if($. % 1000 == 0 ){
    my $msg = sprintf("Reading %d \r",$.);
    print STDERR $msg;
  }
}
print STDERR "\n";
close;
exit;
##sub
sub cigarLen{
    my $cigar =shift;
    my ($mLen,$fLen) = (0,0);
    while($cigar){
        $cigar =~ s/^(\d*)([MIDNSHP=X])//;
        my ($mode,$n) = ($2,$1);
        $n ||=1;
        if($mode =~/[MINP=X]/){$mLen += $n;}
        if($mode =~/[HSMINP=X]/){$fLen += $n;}
    }
    return($mLen,$fLen)
}
