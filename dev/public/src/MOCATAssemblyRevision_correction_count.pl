#!/usr/bin/env perl
use warnings;
use strict;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

#####################################	VARIABLES	#######################################
die "$0 <log.list> <out> \n" unless ( @ARGV == 2 );
my ( $list, $out ) = @ARGV;
my ( $reduce_count, $add_count, $replace_count, $large_indel, $largest_add, $largest_delete );
#####################################	main	###########################################
my $date;

open LOG, "$list" or die "can not open $!\n";
open O,   ">$out" or die "can not open $!\n";
print O "name\treplace_count\tadd_count\tadd_base_num\treduce_count\treduce_base_num\tlarge_indel\tlargest_add\tlargest_delete\n";
$largest_add    = 0;
$largest_delete = 0;
my ( $add_num, $delete_num );
while (<LOG>) {
	my @t    = split /\//;
	my $name = $t[-3];
	open I, "$_" or die "can not open $!\n";
	( $reduce_count, $add_count, $replace_count, $large_indel, $largest_add, $largest_delete, $add_num, $delete_num ) = ( 0, 0, 0, 0, 0, 0, 0, 0 );

	while (<I>) {
		my @temp = split /\s+/;
		next if ( $temp[7] <= 10 );
		if ( @temp == 10 ) {
			$replace_count++;
		}
		elsif ( $temp[10] > $temp[11] ) {
			if ( $temp[8] =~ /\-/ ) {
				$reduce_count++;
				$delete_num += length( $temp[8] ) - 1;
				$large_indel++ if ( length( $temp[8] ) > 30 );
				$largest_delete = length( $temp[8] ) - 1 if ( $largest_delete < length( $temp[8] ) );
			}
			else {
				$add_count++;
				$add_num += length( $temp[8] ) - 1;
				$large_indel++ if ( length( $temp[8] ) > 30 );
				$largest_add = length( $temp[8] ) - 1 if ( $largest_add < length( $temp[8] ) );
			}
		}
		else {
			if ( $temp[9] =~ /\-/ ) {
				$reduce_count++;
				$delete_num += length( $temp[9] ) - 1;
				$large_indel++ if ( length( $temp[9] ) > 30 );
				$largest_delete = length( $temp[9] ) - 1 if ( $largest_delete < length( $temp[9] ) );
			}
			else {
				$add_count++;
				$add_num += length( $temp[9] ) - 1;
				$large_indel++ if ( length( $temp[9] ) > 30 );
				$largest_add = length( $temp[9] ) - 1 if ( $largest_add < length( $temp[9] ) );
			}
		}
	}
	print O "$name\t$replace_count\t$add_count\t$add_num\t$reduce_count\t$delete_num\t$large_indel\t$largest_add\t$largest_delete\n";

	#	$reduce_count=0;
	#	$add_count=0;
	#	$replace_count=0;
	#    $large_indel=0;
	#    $largest_add = 0;
	#    $largest_delete = 0;
	#    ($add_num, $delete_num) = (0,0);
	close I;
}
close LOG;
close O;
print STDERR "revision count program finished!\t", $date = localtime, "\n";
