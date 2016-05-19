package MOCATImport;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

sub run {
	### DEFINE VARIABLES AND OPEN OUTPUT FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	
	my @org_files = <$import_path/*fq $import_path/*fastq $import_path/*fq.gz $import_path/*fastq.gz>;

	my $prevFile    = "";
	my $prevOrgFile = "";
	my $sample;
	my $lane;
	my $counter = -1;
	my @processed;
	for my $i ( 0 .. scalar @org_files - 1 ) {
		$processed[$i] = 0;
	}
	my @files = @org_files;
	foreach my $file (@files) {
		$counter++;
		my $org_file = $file;
		$file =~ s/$import_path\///;
		$file =~ s/^\///;

		if ( $file =~ m/_1_/ || $file =~ m/\.1\./ ) {
			$prevFile    = $file;
			$prevOrgFile = $org_file;
		}

		# Matches pair, enters here when apir 2 file comes
		if ( $file =~ m/_2_/ || $file =~ m/\.2\./ ) {
			if ( $prevFile ne "" ) {

				# True paired reads
				if ( $prevFile =~ m/_1_/ || $prevFile =~ m/.1./ ) {
					print "\tImported pair: $file & $prevFile\n";
					if ( length $file != length $prevFile ) {
						die "ERROR & EXIT: Found pair: $file & $prevFile\nBut they seem to have different file name lengths...\n";
					}
					$sample = "";
					$lane   = "";
					for my $i ( 0 .. length $file ) {
						my $fs  = substr $file,     0, $i;
						my $fps = substr $prevFile, 0, $i;
						if ( $fs eq $fps ) {
							$sample = $fs;
						}
						else {
							my $fn  = substr $fs,  $i - 1, 1;
							my $fpn = substr $fps, $i - 1, 1;
							if ( $fpn eq "1" && $fn eq "2" ) {
								my $lane2 = substr $file,     $i + 1, length $file;
								my $lane1 = substr $prevFile, $i + 1, length $prevFile;
								if ( $lane2 ne $lane1 ) {
									die "ERROR & EXIT: Lane names for $prevFile and $file different.\n";
								}
								$lane = $lane1;
							}
						}
					}
					if ( $sample eq "" ) {
						die "ERROR & EXIT: Sample name for paired files $file and $prevFile to short.\n";
					}
					$lane =~ s/^\.//;
					$lane =~ s/\.$//;
					$lane =~ s/_$//;
					$lane =~ s/^_//;
					my $lane1 = $lane;
					my $lane2 = $lane;
					if ( $lane =~ m/fastq$/ ) {
						$lane1 =~ s/fastq$/1.fq/;
						$lane2 =~ s/fastq$/2.fq/;
					}
					elsif ( $lane =~ m/fq$/ ) {
						$lane1 =~ s/fq$/1.fq/;
						$lane2 =~ s/fq$/2.fq/;
					}
					elsif ( $lane =~ m/fq.gz$/ ) {
						$lane1 =~ s/fq.gz$/1.fq.gz/;
						$lane2 =~ s/fq.gz$/2.fq.gz/;
					}
					elsif ( $lane =~ m/fastq.gz$/ ) {
						$lane1 =~ s/fastq.gz$/1.fq.gz/;
						$lane2 =~ s/fastq.gz$/2.fq.gz/;
					}
					$sample = substr $sample, 0, $import_max_sample_name;
					$sample =~ s/^\.//;
					$sample =~ s/\.$//;
					$sample =~ s/_$//;
					$sample =~ s/^_//;
					system "mkdir -p $cwd/$sample";
					system "ln -sf $org_file $cwd/$sample/$lane2";
					system "ln -sf $prevOrgFile $cwd/$sample/$lane1";
					$prevFile                  = "";
					$processed[$counter]       = 1;
					$processed[ $counter - 1 ] = 1;
					system "echo '$sample' >> samples.imported.paired";
				}
			}
		}    # end ture paired reads
	}

	# Process single files, those not processed above
	@files = @org_files;
	for my $i ( 0 .. scalar @files - 1 ) {
		if ( $processed[$i] == 1 ) {
			next;
		}
		my $file     = $files[$i];
		my $org_file = $file;
		$file =~ s/$import_path\///;
		$file =~ s/^\///;
		print "\tImported single: $file\n";
		$sample = substr $file, 0, $import_max_sample_name;
		$sample =~ s/^\.//;
		$sample =~ s/\.$//;
		$sample =~ s/_$//;
		$sample =~ s/^_//;
		my $lane = $file;
		$lane =~ s/$sample//;
		$lane =~ s/^\.//;
		$lane =~ s/\.$//;
		$lane =~ s/_$//;
		$lane =~ s/^_//;

		if ( $lane =~ m/fastq$/ ) {
			$lane =~ s/fastq$/fq/;
		}
		elsif ( $lane =~ m/fq$/ ) {
			$lane =~ s/fq$/fq/;
		}
		elsif ( $lane =~ m/fq.gz$/ ) {
			$lane =~ s/fq.gz$/fq.gz/;
		}
		elsif ( $lane =~ m/fastq.gz$/ ) {
			$lane =~ s/fastq.gz$/fq.gz/;
		}
		system "mkdir -p $cwd/$sample";
		system "ln -sf $org_file $cwd/$sample/$lane";
		system "echo '$sample' >> samples.imported.single";
	}
	system "sort -u samples.imported.paired > t; mv t samples.imported.paired";
	system "sort -u samples.imported.single > t; mv t samples.imported.single";
}

1;
