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

# NOTE! This is part 1, part 2 needs to be executed via a pipe after part 1

my ( @input_files, $refOut, @totalCounter, @allFiles, @reference_files, @input, $output2, $output, $output1, @input_line, @input_line1, @input_open, $temp_output_file, $direction, @input_direction );
my $zip            = 'gzip -2 -p 2';
my $fileCount      = 0;
my $printedCounter = 0;
GetOptions(
	'temp_output_file:s' => \$temp_output_file,
	'input_files=s{,}'   => \@input_files,
	'references=s{,}'    => \@reference_files,
	'zip=s'              => \$zip,
);

chomp( my $host = `hostname` );

# Open files
print STDERR getLoggingTime() . " MOCATFilter - part 1 : host=$host\n";
print STDERR getLoggingTime() . " MOCATFilter - part 1 : processing initial input & reference files\n";
for my $i ( 0 .. $#input_files ) {
	( -e "$input_files[$i]" ) or die "ERROR & EXIT: Cannot open $input_files[$i]: $!";
	open $input[$i], "gunzip -c $input_files[$i] | " or die "ERROR & EXIT: Cannot open input pipe 'gunzip -c $input_files[$i] | ': $!";
	$input_open[$i]   = 1;
	$totalCounter[$i] = 0;
}
unless ( $reference_files[0] ) {
	die "ERROR & EXIT: No reference files specified using -references";
}
foreach my $reference (@reference_files) {
	( -e $reference ) or die "ERROR & EXIT: Cannot open reference file $reference: $!";
}

# Process
foreach my $reference (@reference_files) {
	my ( $current, $match );
	if ( $reference =~ m/\.pair\.1\./ || $reference =~ m/\.single\./ ) {    # every time we've processed both either A. both apried files, or B. one single file, we move on to the next

		$fileCount++;
		print STDERR getLoggingTime() . " MOCATFilter - part 1 : writing to $temp_output_file.(0|1|r).$fileCount\n";
		close $output1 if $output1;
		close $output2 if $output2;
		close $refOut  if $refOut;

		open $output1, " | $zip -c > $temp_output_file.0.$fileCount" or die "ERROR & EXIT: Cannot open output pipe ' | $zip -c > $temp_output_file.0.$fileCount': $!";
		open $output2, " | $zip -c > $temp_output_file.1.$fileCount" or die "ERROR & EXIT: Cannot open output pipe ' | $zip -c > $temp_output_file.1.$fileCount': $!";
		open $refOut,  " | $zip -c > $temp_output_file.r.$fileCount" or die "ERROR & EXIT: Cannot open output pipe ' | $zip -c > $temp_output_file.r.$fileCount': $!";
		push @allFiles, $fileCount;

		if ( $fileCount > 1 ) {
			print STDERR getLoggingTime() . " MOCATFilter - part 1 : sending singal = " . ( $fileCount - 1 ) . "\n";
			print STDOUT ( $fileCount - 1 ) . "\n";
		}
	}
	if ( $reference =~ m/\.pair\.1\./ || $reference =~ m/\.single\./ ) {    # this flag is used for printing to the reference file only if we match pair.1 or single (pair.2 is redundant)
		$current = 1;
	}
	else {
		$current = undef;
	}

	print STDERR getLoggingTime() . " MOCATFilter - part 1 :  $reference\n";
	open my $REF, "gunzip -c $reference | awk '{if(NR%4==1){print}}' | " or die "ERROR & EXIT: Cannot open reference pipe 'gunzip -c $reference | awk '{if(NR%4==1){print}}' | ': $!";
	while (<$REF>) {
		my $ref = $_;
		$ref =~ s/^@// or die "ERROR & EXIT: Read ID should start with '\@' (line $.)";
		if ($current) {
			print {$refOut} $ref;
		}

		chomp($ref);
		foreach my $i ( 0 .. $#input ) {
			if ( $input_open[$i] ) {    # only process if there are lines left in this input file
				my $loop = 1;     # set that we will by default read more lines fmor file
				if ( $input_line1[$i] ) {    # we already have something stored from previous ref
					if ( $input_line1[$i] eq $ref ) {    # this new ref matches the old line ref
						if ( $ref =~ /1$/ ) {      # if run==2 we print to only one file, might as well be this one
							print {$output1} $input_line[$i];    # print to /1 file
							$printedCounter++;
						}
						elsif ( $ref =~ /2$/ ) {
							print {$output2} $input_line[$i];    # print to /2 file
							$printedCounter++;
						}
						else {
							die "ERROR & EXIT: Read ID does not end with /1 or /2";
						}
						delete $input_line1[$i];                       # delete the ref
					}
					else {                                                   # we still haven't gotten to a ref that exists, wait a bit more
						$loop = 0;                                     # don't read any enw liens from this input file
					}
				}
				while ($loop) {
					if ( my $line = readline( $input[$i] ) ) {               # this returns undef if the last line
						$totalCounter[$i]++;
						my @line = split "\t", $line;                  # split the line
						if ( $line[0] eq $ref ) {                      # this line matches the currewnt ref, process it
							if ( $ref =~ /1$/ ) {                # if run==2 we print to only one file, might as well be this one
								print {$output1} $line;    # print to /1 file
								$printedCounter++;
							}
							elsif ( $ref =~ /2$/ ) {
								print {$output2} $line;    # print to /2 file
								$printedCounter++;
							}
							else {
								die "ERROR & EXIT: Read ID does not end with /1 or /2";
							}
						}
						else {
							$input_line[$i]  = $line;            # this line reference don't match the current ref, we still have to read more ref before we get to this ref in the input file
							$input_line1[$i] = $line[0];         # just so we don't have to split the line again next time we want to check if the ref is the same
							$loop            = 0;
						}
					}
					else {                                                   # last line, close this file by setting the checker to undef
						$input_open[$i] = undef;                       # set checker to undef
						$loop = 0;
						print STDERR getLoggingTime() . " MOCATFilter - part 1 :  end of input file " . ( $i + 1 ) . "\n";
					}
				}
			}
		}
	}
}

# This section is needed to clear the last stored input from previous reference file (perhaps)
foreach my $i ( 0 .. $#input ) {
	if ( $input_line1[$i] ) {    # we already have something stored from previous ref
		if ( $input_line1[$i] =~ /1$/ ) {    # if run==2 we print to only one file, might as well be this one
			print {$output1} $input_line[$i];    # print to /1 file
			#print STDERR "printed end1.1 \{$output1\} $input_line[$i]\n";    ###REM###
			$printedCounter++;
		}
		elsif ( $input_line1[$i] =~ /2$/ ) {
			print {$output2} $input_line[$i];                                # print to /2 file
			#print STDERR "printed end1.2 \{$output2\} $input_line[$i]\n";    ###REM###
			$printedCounter++;
		}
		else {
			die "ERROR & EXIT: Read ID does not end with /1 or /2";
		}
		delete $input_line1[$i];                                                   # delete the ref
	}
}

foreach my $i ( 0 .. $#input ) {
	if ( $input_open[$i] ) {
		my $loop = 1;
		while ($loop) {
			if ( my $line = readline( $input[$i] ) ) {                       # this returns undef if the last line
				$totalCounter[$i]++;
				my @line = split "\t", $line;                          # split the line
				if ( $line[0] =~ /1$/ ) {                              # if run==2 we print to only one file, might as well be this one
					print {$output1} $line;                      # print to /1 file
					#print STDERR "printed end2.1 \{$output1\} $line\n";    ###REM###
					$printedCounter++;
				}
				elsif ( $line[0] =~ /2$/ ) {
					print {$output2} $line;                                # print to /2 file
					#print STDERR "printed end2.2 \{$output2\} $line\n";    ###REM###
					$printedCounter++;
				}
				else {
					die "ERROR & EXIT: Read ID does not end with /1 or /2";
				}
			}
			else {                                                                     # last line, close this file by setting the checker to undef
				$input_open[$i] = undef;                                         # set checker to undef
				$loop = 0;
				print STDERR getLoggingTime() . " MOCATFilter - part 1 :  end of input file " . ( $i + 1 ) . "\n";
			}
		}
	}
}

close $output1 if ($output1);
close $output2 if ($output2);
close $refOut  if ($refOut);
foreach my $file (@input) {
	close $file if $file;
}

print STDERR getLoggingTime() . " MOCATFilter - part 1 : sending last singal = $fileCount\n";
print STDOUT "$fileCount\n";
my $inputSum = 0;
for my $i ( 0 .. $#input_files ) {
	$inputSum += $totalCounter[$i];
}
print STDERR getLoggingTime() . " MOCATFilter - part 1 : [STATS] total_reads_in=$inputSum (" . join( ", ", @totalCounter ) . ") | total_reads_out=$printedCounter\n";
unless ( $printedCounter == $inputSum ) {
	print STDOUT "-1\n";
	die "ERROR & EXIT: input lines and printed lines arenot the same. :( Something is wrong, perhaps the order of input files wasn't correct? Notte that if you do ls *pair* *single* the single files will appear between the paired files, but this is not the order in which MOCAT pastes them, MOCAT pastes them like cat *pair* *single*, which adds all single files at the end";
}

print STDOUT "-2\n";
print STDERR getLoggingTime() . " MOCATFilter - part 1 : finished\n";
exit 0;

sub getLoggingTime {
	my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime(time);
	my $nice_timestamp = sprintf( "%04d%02d%02d %02d%02d%02d", $year + 1900, $mon + 1, $mday, $hour, $min, $sec );
	return $nice_timestamp;
}

