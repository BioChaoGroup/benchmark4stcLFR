#!/usr/bin/env perl
#use strict
my %HASH;
open GFF,"< $ARGV[0]";
while(<GFF>){
  next if /^##/;
  chomp;
  my @FS =split;
  my $relative_i = 0;
  my $length = $FS[4] - $FS[3] + 1;
  for(my $i = $FS[3];$i<=$FS[4];$i++){
    $relative_i += 1;
    $HASH{$i}{ID} = $FS[2];
    $HASH{$i}{RI} = ($FS[6] eq "+")?$relative_i:($length-$relative_i+1);
  }
}
close GFF;

my (@HEAD,@iPOST,@iREF,@iALT,@iFMT);
open VCF,"<$ARGV[1]";
while(<VCF>){
  next if /^##/;
  chomp;
  if(/^#/){       # header line
    $_ =~ s/bwa\///g;s/ToGenome\.bam//g;
    @HEAD = split;
    @iPOS = grep {$HEAD[$_] eq "POS"} keys @HEAD;
    @iREF = grep {$HEAD[$_] eq "REF"} keys @HEAD;
    @iALT = grep {$HEAD[$_] eq "ALT"} keys @HEAD;
    @iFMT = grep {$HEAD[$_] eq "FORMAT"} keys @HEAD;
  }elsif($_ =~ /Escherichia_coli_chromosome/){
    my @FS = split;
    next unless exists $HASH{$FS[$iPOS[0]]};
    unless(/INDEL/){
      $HASH{$FS[$iPOS[0]]}{REF} = $FS[$iREF[0]];
      $HASH{$FS[$iPOS[0]]}{ALT} = $FS[$iALT[0]];
    }
    for(my $i = $iFMT[0] + 1;$i<@FS;$i++){
      $FS[$i] =~ /(\d+):(\d+):([0-9,]+)/;
      my ($PL,$DP,$AD) = ($1, $2, $3);
      my @ADs = split ",", $AD;
      my $sum_alt = 0;
      foreach my $k (1 .. $#ADs){
        $sum_alt += $ADs[$k];
      }
      if(/INDEL/){
        $HASH{$FS[$iPOS[0]]}{$HEAD[$i]}{INDEL} = $sum_alt;
      }else{
        $HASH{$FS[$iPOS[0]]}{$HEAD[$i]}{DP} = $DP;
        $HASH{$FS[$iPOS[0]]}{$HEAD[$i]}{AD} = $sum_alt;
      }
    }
  }
}
close VCF;

#OUTPUT
foreach my $p (sort {$a<=>$b} keys %HASH){
  for(my $i = $iFMT[0] + 1;$i<@HEAD;$i++){
    next unless defined $HASH{$p}{REF};
    print "$HEAD[$i]\t$HASH{$p}{ID}\t$p\t$HASH{$p}{RI}\t$HASH{$p}{REF}\t$HASH{$p}{ALT}";
    print "\t$HASH{$p}{$HEAD[$i]}{DP}\t$HASH{$p}{$HEAD[$i]}{AD}\n";
  }
}

#exit
