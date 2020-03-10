#!/usr/bin/env perl
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Clipping RCA adaptor(s) from contigs (in a stupid way)
# Author:            Chao | fangchao@genomics.cn
# Version:           V0.1
# Last modified:     01 Nov 2019 (since 01 Nov 2019)
# ===================================================================
# see detail below
use strict;
use Getopt::Long qw(:config no_ignore_case); # For case sensitive

sub usage {
  my $msg = shift;
print <<USAGE;
$msg
usage:
  $0 -s <sam file> -o <output>
    -s   Sequence Alignment/Map format file
    -o   output filename
		-f   [0-4095]if specific, run a SAM flag convert function
		-b   if -f is specified, an explaination will return.
    -v   verbose
    -h   show help info
USAGE
}

&usage("Making a bed format file from SAM:") && exit unless @ARGV;

my ($samf,$out,$flag,$bit,$cigar,$verbose,$help);
GetOptions(
  "s=s" => \$samf,
  "o=s" => \$out,
	"f=i" => \$flag,
	"b=i" => \$bit,
	"g=s" => \$cigar,
  "v" => \$verbose,
  "h|help|?" => \$help,
);
&usage && exit if $help;

if(defined $flag){
	my @FEXP = (
	"the read is paired in sequencing",
	"the read is mapped in a proper pair",
	"the query sequence itself is unmapped",
	"the mate is unmapped",
	"strand of the query is reverse",
	"strand of the mate is reverse",
	"the read is the first read in a pair",
	"the read is the second read in a pair",
	"the alignment is not primary",
	"QC failure",
	"optical or PCR duplicate",
	"supplementary alignment"
	);
	my @bits = &flagExp($flag);
	if(defined $bit){
		print "[$bits[$bit]] $FEXP[$bit]\n";
	}else{
		print "@bits\n";
		for(my $i=0;$i<12;$i++){ print "$FEXP[$i]\n" if $bits[$i] eq 1 ;}
	}
	exit;
}

if(defined $cigar){
	my $len = &cigarLen($cigar);
	print $len;
	exit;
}

my (%GENE);
open SAMF, ($samf)?"<$samf":"<-" or die $!;
open OUT,  ($out) ?">$out":"> -" or die $!;
while(<SAMF>){
	$_ =~ /NM:i:(\d+).*AS:i:(\d+)/; next if $1 > 9 or $2 < 500;
	chomp;
	my @s = split;
	my @flags = &flagExp($s[1]);
	if((not defined $GENE{$s[2]}{$s[3]}{'name'}) || $flags[8]==0 ){
		$GENE{$s[2]}{$s[3]}{'name'} = $s[0];
		$GENE{$s[2]}{$s[3]}{'end'}  = $s[3] + &cigarLen($s[5]);
		$GENE{$s[2]}{$s[3]}{'strand'} = ($flags[4]==1)?"-":"+";
		$GENE{$s[2]}{$s[3]}{'primary'} = ($flags[8]==1)?0:1;
	}
}
close SAMF;
# summary output
foreach my $genome ( sort keys %GENE){
	foreach my $pos (sort {$a<=>$b} keys %{$GENE{$genome}}){
		my $end = $GENE{$genome}{$pos}{'end'};
		my $name= $GENE{$genome}{$pos}{'name'};
		my $std = $GENE{$genome}{$pos}{'strand'};
		my $pri = $GENE{$genome}{$pos}{'primary'};
		print OUT "$genome\t$pos\t$end\t$name\t$pri\t$std\t$pos\t$end\t50,50,05\n";
	}
}
close OUT;
exit;
################################################################################
# sub function
################################################################################
sub flagExp{
	return(split("",reverse sprintf("%012b", $_[0])));
}
# check list
=pod
Each bit in the FLAG field is defined as:
bit	Chr	Flag	Description
0	p	0x0001	the read is paired in sequencing
1	P	0x0002	the read is mapped in a proper pair
2	u	0x0004	the query sequence itself is unmapped
3	U	0x0008	the mate is unmapped
4	r	0x0010	strand of the query (1 for reverse)
5	R	0x0020	strand of the mate
6	1	0x0040	the read is the first read in a pair
7	2	0x0080	the read is the second read in a pair
8	s	0x0100	the alignment is not primary
9	f	0x0200	QC failure
X	d	0x0400	optical or PCR duplicate
    0x0800	supplementary alignment
The Please check <http://samtools.sourceforge.net> for the format specification and the tools for post-processing the alignment.
=cut

sub cigarLen{
	my $cigar =shift;
	my $len =0;
	while($cigar){
		$cigar =~ s/^(\d*)([MIDNSHP=X])//;
		my ($mode,$n) = ($2,$1);
		$n ||=1;
		if($mode =~/[MDNP=X]/){$len += $n;}
	}
	return($len)
}
