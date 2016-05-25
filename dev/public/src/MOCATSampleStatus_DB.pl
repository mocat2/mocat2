#!/usr/bin/env perl
use warnings;
use strict;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

#USAGE:parseSampleStatus2DB.pl <sqlite3.db> <samplefile> <bin_dir>
#Example: if you ran MOCAT.pl -sf <samplefile> -ss, then run: parseSampleStatus2DB.pl Cancer.sqlite samples.60

my $OS = "_linux";
chomp( my $systemType = `uname -s` );
if ( $systemType =~ m/Darwin/ ) {
	$OS = "_osx";
}

my $db      = $ARGV[0];
my $sf      = $ARGV[1];
my $bin_dir = $ARGV[2];
my @summary = split "\n", `ls $sf.*.summary`;

foreach my $i (@summary) {
	$i =~ m/$sf\.(.+)\.summary/;
	my $table = $1;
	$table =~ s/\./_/g;
	my @header = split "\t", `head -n 1 $i`;
	my $columns = join( ",", @header );
	chomp $columns;
	system "$bin_dir/sqlite3$OS $db \"DROP TABLE IF EXISTS $table;\"";
	system "$bin_dir/sqlite3$OS $db \"CREATE TABLE IF NOT EXISTS $table ($columns)\";";
	my $com2 = "$bin_dir/sqlite3$OS -separator '\t' $db \".import $i $table\"";
	system "$com2";
	my $com3 = "$bin_dir/sqlite3$OS $db \"DELETE from $table where sample like \\\"%ample\\\";\"";
	system "$com3";
}

exit 0;
