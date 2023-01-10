#!/usr/bin/env perl
#use strict

my %UID;

while(<>){
  next if /^#/;
  chomp;
  my @FS = split;
  next unless $FS[5] =~ /^(\d+)M$/;
  my $len = $1;
  my $end = $FS[3] + $len - 1;
  my $CIGAR = &read_cigar($FS[1]);
  my $std = ($$CIGAR{4} == 1)?"-":"+";
  $FS[0] =~ /(\S+)_(5S|16S|23S|ITS)/;
  my($sp,$unit) = ($1,$2);
  $UID{$sp}{$unit} ++;
  print "$FS[2]\tbam2gff\t$unit\_$UID{$sp}{$unit}\t$FS[3]\t$end\t0\t$std\t.\tName=$sp\_$unit\_$UID{$sp}{$unit}\n"
}

exit;

# Escherichia_coli_chromosome     manually        rRNA_16S_1        2428583 2430125 0       -       .       Name=16S_rRNA;product=16S ribosomal RNA


sub read_cigar {
  my $cigar = shift;
  my %HS;
  for(my $c=11;$c>=0;$c--){
    my $code = 2**$c;
    my $remainder = $cigar % $code;
    $HS{$c} = ($remainder == $cigar)?0:1;
    $cigar = $remainder;
  }
  return(\%HS);
}


=pod
1 0x1 template having multiple segments in sequencing
2 0x2 each segment properly aligned according to the aligner
4 0x4 segment unmapped
8 0x8 next segment in the template unmapped
16 0x10 SEQ being reverse complemented
32 0x20 SEQ of the next segment in the template being reverse complemented
64 0x40 the first segment in the template
128 0x80 the last segment in the template
256 0x100 secondary alignment
512 0x200 not passing filters, such as platform/vendor quality controls
1024 0x400 PCR or optical duplicate
2048 0x800 supplementary alignmen
=cut
