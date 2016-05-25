#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

# Takes a number of input files and references files as input and merges
# the input files in an ordered fashion, in comparion to the reference files
# Limitations:
# 1. reference files MUST be in the same order as they were given to SOAPaligner for mapping
# 2. when ran SOAPaligner, the files have to be in the order which is produced by 'ls *pair*gz *single*gz'
# 3. all input and reference files must be in gz format
# 4. input files MUST be in SOAP format
# 5. reference files MUST be FastQ format
# 6. input reference files must have files names like this: *pair.1.fq* *pair.2.fq* *single*

# NOTE! This version of the script is only part 2, which requires numbers at the STDIN

my ( @input, @input_open, @input_line, @input_line1, $temp_input_file, @input_direction );
my $printedCounter = 0;
my $refCounter     = 0;
my $lineInCounter  = 0;
GetOptions( 'temp_input_file:s' => \$temp_input_file, );

unless ($temp_input_file) {
	die "ERROR & EXIT: Missing -temp_input_file";
}

print STDERR getLoggingTime() . " MOCATFilter - part 2 : waiting for signal\n";

# Process
while (<STDIN>) {
	chomp( my $fileCount = $_ );
	last if $fileCount == -2;
	if ( $fileCount == -1 ) {
		die "ERROR & EXIT: Recevied signal to die from part 1";
	}
	print STDERR getLoggingTime() . " MOCATFilter - part 2 : received signal = $fileCount & processing secondary input files $temp_input_file.(1|2).$fileCount\n";
	for my $i ( 0 .. 1 ) {
		( -e "$temp_input_file.$i.$fileCount" ) or die "ERROR & EXIT: Cannot open $temp_input_file.$i.$fileCount: $!";
		open $input[$i], "gunzip -c $temp_input_file.$i.$fileCount | " or die "ERROR & EXIT: Cannot open input pipe 'gunzip -c $temp_input_file.$i.$fileCount | ': $!";
		$input_open[$i] = 1;
	}
	( -e "$temp_input_file.r.$fileCount" ) or die "ERROR & EXIT: Cannot open $temp_input_file.r.$fileCount: $!";
	open my $REF, "gunzip -c $temp_input_file.r.$fileCount | " or die "ERROR & EXIT: Cannot open input pipe 'gunzip -c $temp_input_file.r.$fileCount | ': $!";

	while (<$REF>) {
		$refCounter++;
		chomp( my $ref = $_ );
		$ref =~ s/.$//;

		foreach my $i ( 0 .. $#input ) {
			if ( $input_open[$i] ) {    # only process if there are lines left in this input file
				my $loop = 1;     # set that we will by default read more lines from file
				if ( $input_line1[$i] ) {    # we already have something stored from previous ref
					if ( $input_line1[$i] eq $ref ) {    # this new ref matches the old line ref
						$input_line1[$i] .= $input_direction[$i];
						print STDOUT $input_line[$i];    # print to /1 file
						delete $input_line1[$i];         # delete the ref
						$printedCounter++;
					}
					else {                                     # we still haven't gotten to a ref that exists, wait a bit more
						$loop = 0;                       # don't read any new liens from this input file
					}
				}
				while ($loop) {
					if ( my $line = readline( $input[$i] ) ) {    # this returns undef if the last line
						$lineInCounter++;
						my @line = split "\t", $line;       # split the line
						$line[0] =~ s/(.)$//;
						my $direction = $1;
						if ( $line[0] eq $ref ) {           # this line matches the currewnt ref, process it
							print STDOUT $line;       # print to /1 file
							$printedCounter++;
						}
						else {
							$input_direction[$i] = $direction;
							$input_line[$i]      = $line;        # this line reference don't match the current ref, we still have to read more ref before we get to this ref in the input file
							$input_line1[$i]     = $line[0];     # just so we don't have to split the line again next time we want to check if the ref is the same
							$loop                = 0;
						}
					}
					else {                                                   # last line, close this file by setting the checker to undef
						$input_open[$i] = undef;                       # set checker to undef
						$loop = 0;
						print STDERR getLoggingTime() . " MOCATFilter - part 2 :  end of input file " . ( $i + 1 ) . "\n";
					}
				}
			}
		}
	}
	foreach my $i ( 0 .. $#input ) {
		if ( $input_line1[$i] ) {                                                              # we already have something stored from previous ref
			$input_line1[$i] .= $input_direction[$i];
			print STDOUT $input_line[$i];                                                # print to /1 file
			$printedCounter++;
			delete $input_line1[$i];                                                     # delete the ref
		}
	}
	foreach my $i ( 0 .. $#input ) {
		if ( $input_open[$i] ) {
			my $loop = 1;
			while ($loop) {
				if ( my $line = readline( $input[$i] ) ) {                         # this returns undef if the last line
					$lineInCounter++;
					print STDOUT $line;                                      # print to /1 file
					$printedCounter++;
				}
				else {                                                             # last line, close this file by setting the checker to undef
					$input_open[$i] = undef;                                 # set checker to undef
					$loop = 0;
					print STDERR getLoggingTime() . " MOCATFilter - part 2 :  end of input file " . ( $i + 1 ) . "\n";
				}
			}
		}
	}
	unlink("$temp_input_file.r.$fileCount") or warn "WARNING: Could not remove $temp_input_file.r.$fileCount";
	unlink("$temp_input_file.1.$fileCount") or warn "WARNING: Could not remove $temp_input_file.1.$fileCount";
	unlink("$temp_input_file.0.$fileCount") or warn "WARNING: Could not remove $temp_input_file.0.$fileCount";
	print STDERR getLoggingTime() . " MOCATFilter - part 2 : finished, deleted files & waiting for signal\n";
}

print STDERR getLoggingTime() . " MOCATFilter - part 2 : [STATS] refrences_in=$refCounter | total_reads_in=$lineInCounter | total_reads_out=$printedCounter\n";
print STDERR getLoggingTime() . " MOCATFilter - part 2 : finished\n";
exit 0;

sub getLoggingTime {
	my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime(time);
	my $nice_timestamp = sprintf( "%04d%02d%02d %02d%02d%02d", $year + 1900, $mon + 1, $mday, $hour, $min, $sec );
	return $nice_timestamp;
}

