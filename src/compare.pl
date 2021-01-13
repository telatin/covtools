#!/usr/bin/perl

use 5.012;
use warnings;
use File::Basename;

my $dir = dirname($0);
my $nimbin = $dir . "/covtobed3";


for my $input_file (@ARGV) {
  if ( -e "$input_file" ) {
    my @covtobed = `covtobed -d "$input_file"`;
    my @covtonim = `$nimbin "$input_file"`;

    my $equal = 0;
    my $diff  = 0;
    for my $cov_line (@covtobed) {
      my $nim_line = shift @covtonim;
      chomp($cov_line);
      chomp($nim_line);
      my @cov = split /\t/, $cov_line;
      my @nim = split /\t/, $nim_line;

      if ($cov[0] eq $nim[0] and $cov[1]==$nim[1] and $cov[2]==$nim[2]) {
        $equal++;
      } else {
        $diff++;
        if ($diff == 1) {
          print "$equal\tCOV:$cov_line\tNIM:$nim_line\n";
        }
      }
    }
    say STDERR basename($input_file), "\t$equal OK; $diff Wrong";
  }

}