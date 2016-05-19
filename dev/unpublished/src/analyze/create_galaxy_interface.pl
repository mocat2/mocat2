#!/usr/bin/perl
# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.
use strict;
use warnings;

# TOOLNAME
# TOOLDESC
# TOOLID

unless ( $ARGV[0] ) {
	die "Missing intput file with .R extension as input.\n";
}

my $file = $ARGV[0];

print "Processing $file...";

my $galaxy_location   = "/g/bork2/kultima/bin/galaxy/galaxy-dist/";
my $output            = "/g/bork2/kultima/bin/galaxy/galaxy-dist/tools/mocat_analyses/$file.xml";


my $input             = "$file.R";


my $METADATA_REQUIRED = 0;
chomp( my $DESC     = `grep '^# DESC' $input | sed 's/.*DESC//' | sed 's/^ *//'` );
chomp( my @REQUIRE  = `grep -h '^# REQUIRE' ../MOCATAnalyze_R.R $file.R | sed 's/.*REQUIRE//' | sed 's/^ *//'` );
chomp( my @SETTINGS = `grep -h '^# SETTINGS' ../MOCATAnalyze_R.R $file.R | sed 's/.*SETTINGS//' | sed 's/^ *//'` );
chomp( my @PRODUCES = `grep '^# PRODUCES' $file.R | sed 's/.*PRODUCES//' | sed 's/^ *//'` );
chomp( my @SUPPRESS = `grep '^# SUPPRESS' $file.R | sed 's/.*SUPPRESS//' | sed 's/^ *//'` );

my %REQUIRE;
foreach my $require (@REQUIRE) {
	$REQUIRE{$require} = 1;
}
my %SUPPRESS;
foreach my $require (@SUPPRESS) {
	$require =~ m/\s*(.*)\s*/;
	$SUPPRESS{$1} = 1;
}

system "cat part.1 > $output";

# Print parameters to MOCAT
foreach my $require (@REQUIRE) {
	if ( $require eq 'featureA' ) {
		open OUT, ">>$output";
		print OUT "    -analyze_condA '\$featureA'\n    -analyze_condB '\$featureB'\n";
		close OUT;
	}
	if ( $require =~ m/groups\s*\|\s*1/ ) {
		$METADATA_REQUIRED = 1;
		open OUT, ">>$output";
		print OUT "    -analyze_group '\$metadata_source.group'\n";
		print OUT "    -analyze_project_name '\$metadata_source.metadata'\n";
		close OUT;
	}
	if ( $require eq 'condA' ) {
		open OUT, ">>$output";
		print OUT "    -analyze_condA '\$condA'\n    -analyze_condB '\$condB'\n";
		close OUT;
	}
	if ( $require eq 'metadata' && !$REQUIRE{'groups|1'} ) {
		$METADATA_REQUIRED = 1;
		open OUT, ">>$output";
		print OUT "    -analyze_project_name '\$metadata'\n";
		close OUT;
	}
}

# Additional parameters to MOCAT printed
open OUT, ">>$output";
if ( defined( $PRODUCES[0] ) ) {
	print OUT "    -analyze_output ";
	foreach my $produce (@PRODUCES) {
		print OUT " '$produce' '\$output_$produce'";
	}
	print OUT "\n";
}
if ( defined( $SETTINGS[0] ) ) {
	print OUT "    -analyze_settings ";
	foreach my $setting (@SETTINGS) {
		$setting =~ m/\s*:\s*(\w*)\s*\[(.*)\]\s*\{(.*)\}\s*\{(.*)\}\s*\((.*)\)/;
		my $name        = $1;
		my $fields      = $2;
		my $check       = $3;
		my $default     = $4;
		my $description = $5;
		if ( !$SUPPRESS{$name} ) {
			print OUT " '$name' '\$output_$name' ";
		} elsif ($SUPPRESS{$name}) {
			print OUT " '$name' '$default' ";
		}
	}
	print OUT "\n";
}
close OUT;

# Beginning of input section
system "cat part.2 >> $output";

# Prints sections needed in galaxy
foreach my $require (@REQUIRE) {
	if ( $require eq 'featureA' ) {
		system "cat part.2feature >> $output";
	}
	if ( $require =~ m/groups\s*\|\s*1/ ) {
		system "cat part.group1 >> $output";
	}
	if ( $require eq 'condA' ) {
		system "cat part.condA >> $output";
	}
	if ( $require eq 'metadata' && !$REQUIRE{'groups|1'} ) {
		system "cat part.metadata >> $output";
	}
}

# Input section for SETTINGS
open OUT, ">>$output";
if ( defined( $SETTINGS[0] ) ) {
	foreach my $setting (@SETTINGS) {
		$setting =~ m/\s*:\s*(\w*)\s*\[(.*)\]\s*\{(.*)\}\s*\{(.*)\}\s*\((.*)\)/;
		my $name        = $1;
		my $fields      = $2;
		my $check       = $3;
		my $default     = $4;
		my $description = $5;
		unless ( $SUPPRESS{$name} ) {
			if ( $check eq 'CHECK' ) {
				print OUT "  <param name=\"output_$name\" type=\"select\" label=\"$description\">\n";
				my @fields = split ",", $fields;
				foreach my $field (@fields) {
					$field =~ s/^\s*//;
					$field =~ s/\s*$//;
					if ( $field eq $default ) {
						print OUT "     <option value=\"$field\" selected=\"true\" >$field</option>\n";
					}
					else {
						print OUT "     <option value=\"$field\" >$field</option>\n";
					}
				}
				print OUT "  </param>\n";
			}
			else {
				print OUT "  <param name=\"output_$name\" type=\"text\" value=\"$default\" size=\"40\" label=\"$description\" />\n";
			}
		}
	}
}
close OUT;

# OUTPUT section
open OUT, ">>$output";
print OUT "  \n</inputs>\n    <outputs>\n";
foreach my $produce (@PRODUCES) {
	print OUT "<data name=\"output_$produce\" format=\"$produce\" metadata_source=\"input_file\" input_dataset=\"input_file\" />\n";
}
print OUT "  </outputs>\n";
close OUT;

# If metadata is required, print which metadata is available
if ($METADATA_REQUIRED) {
	system "cat part.metadata.help >> $output";
}

system "cat part.end >> $output";

system "sed -i 's/TOOLNAME/$file/' $output";
system "sed -i 's/TOOLID/$file/' $output";
system "sed -i 's/TOOLDESC/$DESC/' $output";

print " OK!\n";

exit(0);

