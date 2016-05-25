#!/usr/bin/env perl
use strict;
use warnings;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

$|++;
my $prev_insert = "0";
my %fw;
my %rev;
my %fwc;
my %revc;
my $total = 0;
my $keep  = 0;
print STDERR getLoggingTime() . " MOCATFilter - filterPE : started\n";

while (<STDIN>) {
	my @line = split "\t", $_;
	$line[0] =~ m/(.+)\/([12])/ or die "ERROR & EXIT: Read ID not matching correct format (xxx/1 or xxx/2)";
	my $insert    = $1;
	my $direction = $2;
	my $ref_id    = $line[2] or die "ERROR & EXIT: Cannot set reference ID";
	$total++;
	if ( $insert ne $prev_insert ) {
		my %common = ();

		#put ref_ids that have been seen by both fw and rev into %common
		foreach my $k ( keys %fw ) {
			$common{$k}++ if exists $rev{$k};
		}
		if ( scalar( keys %common ) > 0 ) {
			foreach my $key ( keys %common ) {
				print STDOUT $fw{$key};
				$keep += $fwc{$key};
			}
			foreach my $key ( keys %common ) {
				print STDOUT $rev{$key};
				$keep += $revc{$key};
			}
		}
		else {
			foreach my $key ( keys %fw ) {
				print STDOUT $fw{$key};
				$keep += $fwc{$key};
			}
			foreach my $key ( keys %rev ) {
				print STDOUT $rev{$key};
				$keep += $revc{$key};
			}
		}
		%fw   = ();
		%rev  = ();
		%fwc  = ();
		%revc = ();
	}
	if ( $direction == 1 ) {    #a fw read has seen this ref_id
		$fw{$ref_id} .= $_;
		$fwc{$ref_id} += 1;
	}
	elsif ( $direction == 2 ) {    #a rev read has seen this ref_id
		$rev{$ref_id} .= $_;
		$revc{$ref_id} += 1;
	}
	$prev_insert = $insert;
}

my %common = ();
foreach my $k ( keys %fw ) {
	$common{$k}++ if exists $rev{$k};
}
if ( scalar( keys %common ) > 0 ) {
	foreach my $key ( keys %common ) {
		print STDOUT $fw{$key};
		$keep += $fwc{$key};
	}
	foreach my $key ( keys %common ) {
		print STDOUT $rev{$key};
		$keep += $revc{$key};
	}
}
elsif ( scalar( keys %common ) == 0 ) {
	foreach my $key ( keys %fw ) {
		print STDOUT $fw{$key};
		$keep += $fwc{$key};
	}
	foreach my $key ( keys %rev ) {
		print STDOUT $rev{$key};
		$keep += $revc{$key};
	}
}

print STDERR getLoggingTime() . " MOCATFilter - filterPE : [STATS] total_reads_in=$total | total_reads_out=$keep\n";
print STDERR getLoggingTime() . " MOCATFilter - filterPE : finished\n";

exit 0;

sub getLoggingTime {
	my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime(time);
	my $nice_timestamp = sprintf( "%04d%02d%02d %02d%02d%02d", $year + 1900, $mon + 1, $mday, $hour, $min, $sec );
	return $nice_timestamp;
}
