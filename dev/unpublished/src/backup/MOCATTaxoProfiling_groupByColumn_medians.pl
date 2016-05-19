#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Statistics::Basic qw(:all);
use Number::Format;
my $format = new Number::Format( -thousands_sep => '' );



my ( $input_file, $group_column, $first_data_column, $header_rows, $number_of_columns );
my $usage = "
NOTE: THIS FUNCTION RETURNS MEDIAN VALUES, WHEN DISREGARDING 0s.
USAGE: groupByColumn.pl [-i|input_file] [-g|group_column] [-f|first_data_column] [-h|header_rows]
-i|input_file           file to parse
-g|group_column         non-unique values in this column will be grouped and rows in data rows summed up
-f|first_data_column    first column with data that should be summed up for each distinct value in group_column
-h|header_rows          number of rows that should be ignored for grouping
\n";

GetOptions(
	'i=s' => \$input_file,
	'g=i' => \$group_column,
	'f=i' => \$first_data_column,
	'h=i' => \$header_rows,
);

unless ( $input_file =~ /^\S+/ && $group_column >= 1 && $first_data_column >= 2 && $header_rows >= 0 ) { die $usage, "\nEXIT: check your input options\n\n" }
$first_data_column -= 1;
$group_column      -= 1;
open( IN, $input_file );

my $href = {};
while (<IN>) {
	unless ( $. <= $header_rows ) {
		my @line = split "\t", $_;
		chomp @line;
		my @val = @line[ $first_data_column .. $#line ];
		$number_of_columns = scalar @val;
		for ( my $i = 0 ; $i < scalar @val ; $i++ ) {
			if ( $val[$i] > 0 ) {
				push @{ $href->{ $line[$group_column] }{$i} }, $val[$i];
			}
		}
	}
	else {
		my @header = split "\t", $_;
		chomp @header;
		my @newheader;
		push( @newheader, $header[$group_column] );
		push( @newheader, @header[ $first_data_column .. $#header ] );
		print join( "\t", @newheader ), "\n";
	}
}

foreach my $key ( keys %{$href} ) {
	print "$key";
	for ( my $i = 0 ; $i < $number_of_columns ; $i++ ) {
		my $median;
		if ( $href->{$key}{$i} ) {
			$median = median( @{ $href->{$key}{$i} } );
			$median = $format->format_number($median);
		}
		else {
			$median = "0";
		}
		print "\t$median";
	}
	print "\n";
}
