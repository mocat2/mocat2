#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

my ($carma, $coverage_header, $coverage_file, $temp_dir, $bin_dir);
my $OS = "_linux";
chomp(my $systemType = `uname -s`);
if ($systemType =~ m/Darwin/) {
  $OS = "_osx";
}

GetOptions(
	   'carma_output=s'     => \$carma,
	   'coverage_file=s'    => \$coverage_file,
	   'coverage_header=s'  => \$coverage_header,
	   'temp_folder=s'      => \$temp_dir,
	   'bin_dir=s'          => \$bin_dir
	  );

system "paste $coverage_header $coverage_file | awk 'NR > 2'> $temp_dir/cov; paste $carma $temp_dir/cov > $temp_dir/mod";

system "../gTP.pl -o $temp_dir/output -g /g/bork/kultima/projects/mock/temp/scaftigs/carma $temp_dir/mod";

#my @header = split "\t", `head -n 1 $temp_dir/cov`;
#my $columns = join (",", @header);
#system "$bin_dir/sqlite3$OS $temp_dir/db \"DROP TABLE IF EXISTS cov;\"";
#system "$bin_dir/sqlite3$OS $temp_dir/db \"CREATE TABLE IF NOT EXISTS cov ($columns)\";";
#system "$bin_dir/sqlite3$OS -separator '\t' $temp_dir/db \".import $temp_dir/cov cov\"";




