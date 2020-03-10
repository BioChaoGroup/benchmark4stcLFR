use strict;
my ($pi,$region,$bn,%HS) = ("",0,0,);
while(<>){
  my @s = split;
  my @i = split("/",$s[0]);
  if($i[1] ne $pi){
    # summary
    foreach my $r (sort {$a<=>$b} keys %{$HS{$pi}}){
      print "$pi\t$r\t$HS{$pi}{$r}{scf}\t$HS{$pi}{$r}{min}\t$HS{$pi}{$r}{max}\t$HS{$pi}{$r}{num}\n";

    }
    # init new
    $region = 1; $bn=0;
    $HS{$i[1]}{$region}{scf}=$s[2];
    $HS{$i[1]}{$region}{min}=$s[3];
    $HS{$i[1]}{$region}{max}=$s[3];
    $HS{$i[1]}{$region}{num}=0; 
  }else{
    my $new = 1;
    foreach my $r (keys %{$HS{$i[1]}}){
      my ($min,$max,$scf) = ($HS{$i[1]}{$r}{min}, $HS{$i[1]}{$r}{max}, $HS{$i[1]}{$r}{scf});
      if($s[2] eq $scf && $min - 5000 < $s[3] && $max + 5000 > $s[3]){
        $HS{$i[1]}{$r}{min} = ($s[3]<$HS{$i[1]}{$r}{min})?$s[3]:$HS{$i[1]}{$r}{min};
        $HS{$i[1]}{$r}{max} = ($s[3]>$HS{$i[1]}{$r}{max})?$s[3]:$HS{$i[1]}{$r}{max};
        $HS{$i[1]}{$r}{num} ++; $new=0; last;
      }
    }
    if($new){
      $region ++;
      $HS{$i[1]}{$region}{scf}=$s[2];
      $HS{$i[1]}{$region}{min}=$s[3];
      $HS{$i[1]}{$region}{max}=$s[3];
      $HS{$i[1]}{$region}{num}=0; 
    }
  }
  #
  $pi = $i[1];
}
foreach my $r (sort {$a<=>$b} keys %{$HS{$pi}}){
  print "$pi\t$r\t$HS{$pi}{$r}{scf}\t$HS{$pi}{$r}{min}\t$HS{$pi}{$r}{max}\t$HS{$pi}{$r}{num}\n";

}
