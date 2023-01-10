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
    $HASH{$i}{ID} = "$FS[0]|$FS[2]";
    $HASH{$i}{RI} = ($FS[6] eq "+")?$relative_i:($length-$relative_i+1);
  }
}
close GFF;

my (@HEAD,@iPOST,@iREF,@iALT,@iFMT);
open VCF,"<$ARGV[1]";
while(<VCF>){
  next if /^##/;
  chomp;
  if(/^#CHROM/){       # header line
    $_ =~ s/bwa\///g;s/ToGenome(|.filter)\.bam//g;s/.genome.sort.filter.bam//g;

    @HEAD = split;
    @iPOS = grep {$HEAD[$_] eq "POS"} keys @HEAD;
    @iREF = grep {$HEAD[$_] eq "REF"} keys @HEAD;
    @iALT = grep {$HEAD[$_] eq "ALT"} keys @HEAD;
    @iFMT = grep {$HEAD[$_] eq "FORMAT"} keys @HEAD;

  }else{
    my @FS = split;
    next unless exists $HASH{$FS[$iPOS[0]]};
    my @ATs = split ",", $FS[$iALT[0]];
    unshift (@ATs, $FS[$iREF[0]]);

    for(my $i = $iFMT[0] + 1;$i<@FS;$i++){
      $FS[$i] =~ /(\d+):(\d+):([0-9,]+)/;
      my ($PL,$DP,$AD) = ($1, $2, $3);
      my @ADs = split ",", $AD;
      my $sum_alt = 0;
      my $ref_count =0;
      my $ref_type = $ATs[0];
      foreach my $k (0 .. $#ADs){
        $sum_alt += $ADs[$k];
        $ref_type  = ($ADs[$k]>$ref_count)?$ATs[$k]:$ref_type;
        $ref_count = ($ADs[$k]>$ref_count)?$ADs[$k]:$ref_count;
      }
      $sum_alt -= $ref_count;

      $HASH{$FS[$iPOS[0]]}{REF} = $ref_type;
      $HASH{$FS[$iPOS[0]]}{ALT} = join(",",@ATs);

      if(/INDEL/){
        my @indels=(0,0,0);
        foreach my $k (0 .. $#ADs){
          my $indel_type = (length($ATs[$k])>length($ref_type))?2:(length($ATs[$k])<length($ref_type))?0:1;
          $indels[$indel_type] += $ADs[$k];
        }
        $HASH{$FS[$iPOS[0]]}{$HEAD[$i]}{INDEL}{REF} = $ref_type;
        $HASH{$FS[$iPOS[0]]}{$HEAD[$i]}{INDEL}{ALT} = join(",",@ATs);
        $HASH{$FS[$iPOS[0]]}{$HEAD[$i]}{INDEL}{INS} = $indels[2];
        $HASH{$FS[$iPOS[0]]}{$HEAD[$i]}{INDEL}{DEL} = $indels[0];
        $HASH{$FS[$iPOS[0]]}{$HEAD[$i]}{INDEL}{DP}  = $DP;
      }else{
        $HASH{$FS[$iPOS[0]]}{$HEAD[$i]}{SUB}{REF} = $ref_type;
        $HASH{$FS[$iPOS[0]]}{$HEAD[$i]}{SUB}{ALT} = join(",",@ATs);
        $HASH{$FS[$iPOS[0]]}{$HEAD[$i]}{SUB}{DP} = $DP;
        $HASH{$FS[$iPOS[0]]}{$HEAD[$i]}{SUB}{AD} = $sum_alt;
      }
    }
  }
}
close VCF;

#OUTPUT
foreach my $p (sort {$a<=>$b} keys %HASH){
  for(my $i = $iFMT[0] + 1;$i<@HEAD;$i++){
    next unless defined $HASH{$p}{REF};
    my $info = "$HEAD[$i]\t$HASH{$p}{ID}\t$p\t$HASH{$p}{RI}";
    if($HASH{$p}{$HEAD[$i]}{SUB}{AD}>0){
      print "$info\t$HASH{$p}{$HEAD[$i]}{SUB}{REF}\t$HASH{$p}{$HEAD[$i]}{SUB}{ALT}";
      print "\tSUB\t$HASH{$p}{$HEAD[$i]}{SUB}{DP}\t$HASH{$p}{$HEAD[$i]}{SUB}{AD}\n";
    }
    if($HASH{$p}{$HEAD[$i]}{INDEL}{INS} > 0){
      print "$info\t$HASH{$p}{$HEAD[$i]}{INDEL}{REF}\t$HASH{$p}{$HEAD[$i]}{INDEL}{ALT}";
      print "\tINS\t$HASH{$p}{$HEAD[$i]}{INDEL}{DP}\t$HASH{$p}{$HEAD[$i]}{INDEL}{INS}\n";
    }
    if($HASH{$p}{$HEAD[$i]}{INDEL}{DEL} > 0){
      print "$info\t$HASH{$p}{$HEAD[$i]}{INDEL}{REF}\t$HASH{$p}{$HEAD[$i]}{INDEL}{ALT}";
      print "\tDEL\t$HASH{$p}{$HEAD[$i]}{INDEL}{DP}\t$HASH{$p}{$HEAD[$i]}{INDEL}{DEL}\n";
    }
  }
}

#exit
