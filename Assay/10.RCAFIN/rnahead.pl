#!/usr/bin/env perl
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       stat barnnap detected ID names
# Author:            Chao | fangchao@genomics.cn
# Version:           V0.0
# Last modified:    10 Aug 2019 (since 10 Aug 2019)
# ===================================================================
# see detail below
use strict;

my($pID);
while(<>){
  chomp;
  my @s = split /[>:_]/;
  $s[13] =~ /^(\d+)-(\d+)\(+|-\)$/;
  my($k5,$k3,$fr) = ($1,$2,$3);
  print "$s[4]\t$s[5]\_$s[6]\t$s[12]\t$k5\t$k3\t$fr\t$s[1]\n";
  if($. % 1000 == 0 ){
    my $msg = sprintf("Reading %d \r",$.);
    print STDERR $msg;
  }
}
print STDERR "DONE\n";
close;
exit;
