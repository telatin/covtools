#!/usr/bin/perl

use 5.012;
use warnings;
use File::Basename;
my ($file1, $file2) = @ARGV;

open(my $I, '<', "$file1") || die "File1?\n";
open(my $J, '<', "$file2") || die "File2?\n";

my $c = 0;
my $ok = 0;
while (my $i = readline($I)) {
  my $j = readline($J);
  chomp($i);
  chomp($j);
  $c++;
  my @fields_1 = split /\t/, $i;
  my @fields_2 = split /\t/, $j;

  if ($fields_1[0] eq $fields_2[0] and $fields_1[1] eq $fields_2[1] and $fields_1[2] eq $fields_2[2]) {
    $ok++;
  } else {
    print "$c\t$i\t$j\n";
    die;
  }

}