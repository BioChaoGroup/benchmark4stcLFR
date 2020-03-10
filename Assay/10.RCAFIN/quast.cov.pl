#!/usr/bin/env perl
use strict;
my $pfx = shift @ARGV;
my(%COV,%REG);
while(<>){
	chomp;
	my @a=split;
	my $ref = $a[4];
	for(my $i=$a[0];$i<=$a[1];$i++){
		$COV{$ref}{$i} ++;
	}
	unless(defined $REG{$ref}){
		$REG{$ref}{$a[0]}{E} = $a[1];
		$REG{$ref}{$a[0]}{M} = $a[5];
	}else{
		my @region = ($a[0],$a[1]);
		my @ss=(sort {$a<=>$b} keys %{$REG{$ref}});
		my $find = 0;
		for(my $i=0;$i<@ss;$i++){
			my $s = $ss[$i];
			if($a[1]<$s){
				$REG{$ref}{$a[0]}{E} = $a[1];
				$REG{$ref}{$a[0]}{M} = $a[5];
				$find = 1;
				last;
			}elsif($a[0]>$REG{$ref}{$s}{E}){
				next;
			}else{
				my @st = sort {$a<=>$b} ($a[0],$a[1],$s,$REG{$ref}{$s}{E});
				if($a[0]<$s){
					$REG{$ref}{$a[0]}{E} = ($REG{$ref}{$s}{E}<$a[1])?$a[1]:$REG{$ref}{$s}{E};
					$REG{$ref}{$a[0]}{M} = $REG{$ref}{$s}{M}.",".$a[5];
					delete $REG{$ref}{$s};
				}else{
					$REG{$ref}{$s}{E} = ($REG{$ref}{$s}{E}<$a[1])?$a[1]:$REG{$ref}{$s}{E};
					$REG{$ref}{$s}{M} = $REG{$ref}{$s}{M}.",".$a[5];
				}
				if($a[1] eq $st[3]){
					$i ++;
					$s = $ss[$i];
					while($st[3]>$s && $i <@ss){
						$REG{$ref}{$st[0]}{E} = ($st[3]<$REG{$ref}{$s}{E})?$REG{$ref}{$s}{E}:$st[3];
						$REG{$ref}{$st[0]}{M} .= ",".$REG{$ref}{$s}{M};
						delete $REG{$ref}{$s};
						$i ++;
						$s = $ss[$i];
					}
				}
				$find = 1;
				last;
			}
		}
		if($find==0){
			$REG{$ref}{$a[0]}{E} = $a[1];
			$REG{$ref}{$a[0]}{M} = $a[5];
		}
	}
	#
	print STDERR sprintf("Reading %8d lines ...\n",$.) if $. % 1000 == 0;
}
#
open COV,">$pfx.cov";
open REG,">$pfx.reg";
foreach my $r (sort keys %REG){
	foreach my $s (sort {$a<=>$b} keys %{$REG{$r}}){

		#print 1st COV
		my ($pPos,$pCov,$maxCov) = ($s,$COV{$r}{$s},$COV{$r}{$s});
		print COV "$r\t$s\t$s\t$pCov\n";
		#print when COV are differ
		for(my $i=$s+1;$i<$REG{$r}{$s}{E};$i++){
			if($COV{$r}{$i} ne $pCov){
				print COV "$r\t$s\t".($i-1)."\t$pCov\n" if $i > $pPos + 1;
				print COV "$r\t$s\t$i\t$COV{$r}{$i}\n";
				($pPos,$pCov) = ($i,$COV{$r}{$i});
				$maxCov = ($maxCov<$pCov)?$pCov:$maxCov;
			}
		}
		#print last COV
		print COV "$r\t$s\t$REG{$r}{$s}{E}\t$COV{$r}{$REG{$r}{$s}{E}}\n";
		#sum region info
		print REG "$r\t$s\t".($REG{$r}{$s}{E}-$s+1)."\t$maxCov\t$REG{$r}{$s}{M}\n";
	}
}
