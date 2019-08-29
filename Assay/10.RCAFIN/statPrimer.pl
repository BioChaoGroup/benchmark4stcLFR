#!/usr/bin/env perl
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       stat primers blasted from clean reads
# Author:            Chao | fangchao@genomics.cn
# Version:           V0.0
# Last modified:    10 Aug 2019 (since 10 Aug 2019)
# ===================================================================
# see detail below
use strict;
use Getopt::Long qw(:config no_ignore_case); # For case sensitive

my ($mode,$verbose,$help);
sub usage {
  my $msg = shift;
print <<USAGE;
$msg
usage:
  $0 action [options]
action:
	count
	link
USAGE
}

&usage("stat primers blasted from clean reads") && exit unless @ARGV;

# Main start
$mode = shift @ARGV unless $ARGV[0] =~/^-/;
&run_count if $mode eq "count";
&run_link  if $mode eq "link";
# Main end
exit;

##########################
# Functions
##########################

sub usage4count {
  my $msg = shift;
print <<USAGE;
$msg
usage:
  $0 count -i input -o output
    -i  stLFR blast m6 format file
    -o  output stat file
    -v  verbose
    -h  show help info
USAGE
}


sub run_count{
&usage4count("Show count usage:") && exit unless @ARGV;
my ($inf,$out,$verbose,$help);
GetOptions(
  "i=s" => \$inf,
  "o=s" => \$out,
  "v" => \$verbose,
  "h|help|?" => \$help,
);
&usage4count("stat primers blasted from clean reads") && exit if $help;

open INF, ($inf)?"<$inf":"< -" or die $!;
open OUT, ($out)?">$out":"> -" or die $!;

my (%BEADS);
&verbose("[log] Start ... ");
while(<INF>){
	my @L = split(/\/|\t/);
	next if $L[5] < 15;
	$BEADS{$L[1]}{$L[3]} ++;
}
close INF;
##output
foreach my $bid (sort keys %BEADS){
	foreach my $pri (sort keys %{$BEADS{$bid}}){
		print OUT "$bid\t$pri\t$BEADS{$bid}{$pri}\n"
	}
}
close OUT;
# Main end

&verbose("Count done!\n");
}

###### stat link ###############
sub usage4link {
  my $msg = shift;
print <<USAGE;
$msg
usage:
  $0 count -i input -o output
    -i  stLFR blast m6 format file
    -o  output stat file
    -v  verbose
    -h  show help info
USAGE
}


sub run_link{
&usage4link("Show count usage:") && exit unless @ARGV;
my ($inf,$out,$verbose,$help);
GetOptions(
  "i=s" => \$inf,
  "o=s" => \$out,
  "v" => \$verbose,
  "h|help|?" => \$help,
);
&usage4link("stat primers-link-patterns blasted from clean reads") && exit if $help;

open INF, ($inf)?"<$inf":"< -" or die $!;
open OUT, ($out)?">$out":"> -" or die $!;
open LOG, ($out)?">$out.log":"> STD.log" or die $!;

my (%READS,$pL0,$pL1,%SB,%SR,%ST,%SID,@SIDM);
&verbose("[link] start ... ");
while(<INF>){
    my @L = split(/\/|\t/);
    next if $L[5] < 15 || $L[1]=~/0000/;
	if($L[0] ne $pL0 && keys %READS >0 ){
		my $pat = &sumR(\%READS);
		$ST{$pat} = 1;
		$SR{$pat} ++;
		$SID{$pat} = $pL0;
		if($pat=~/RCA\+\?|\?\+RCA/){print LOG "$pat\t$pL0\n"}
		if($L[1] ne $pL1){
			foreach my $p (sort keys %ST){
				$SB{$p} ++;
			}
			%READS = (); %ST = (); $pL1 = $L[1];
		}
	}
    $READS{$L[8]}{$L[9]} = $L[3];
	$pL0 = $L[0];
}
#save last one
my $pat = &sumR(\%READS);
$ST{$pat} = 1;
$SR{$pat} ++;
foreach my $p (sort keys %ST){$SB{$p} ++;}
#
close INF;
##output
foreach my $pat (sort keys %SB){
		my $pid = ($SB{$pat}==1)?$SID{$pat}:"multi";
        print OUT "$pat\t$SB{$pat}\t$SR{$pat}\t$pid\n";
}
close OUT;
# Main end
&verbose("[link] done!\n");
}

sub sumR{
my $RD = shift;
my %TMP;
my ($count,$ps,$pe,$pp,$p,$PAT) = (0,0,0,"","","");
foreach my $s (sort {$a<=>$b} keys %$RD){
#	my $me = 0;
	my $e = (sort {$b<=>$a} keys %{$$RD{$s}})[0];
#	foreach my $e (sort {$a<=>$b} keys %{$$RD{$s}}){
#		$me = ($e>$me)?$e:$me;
#	}
	$p = $$RD{$s}{$e}; # determined current primer
	# compare with previous primer
	if($pp eq ""){
		if($s<15){$PAT = "|+$p"}else{$PAT = "?+$p"};
		$count = 1;
	}else{
		if(abs($s-$pe-1) < 2){
			$PAT .= "+$p"; $count ++;
			if($pp eq "RCA"){
			#	if(100-$e>=15){$PAT.="+?"}else{$PAT.="+|"}
				last;
			}
		}elsif($s-$ps>15){
			$TMP{"$PAT+x"} = $count;
			($count,$PAT) = (1,"x+$p");
		}
	}
	$pp=$p; $ps=$s; $pe=$e;
}
if(100-$pe>=15){$PAT.="+?"}else{$PAT.="+|"};
foreach my $p (sort keys %TMP){
	if($TMP{$p} > $count){
		$PAT = $p; $count = $TMP{$p};
	}
}
return($PAT);
}

#################################
sub verbose{
  my $msg = shift;
  print STDERR $msg if $verbose;
}
