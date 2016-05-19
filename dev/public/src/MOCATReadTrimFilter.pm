package MOCATReadTrimFilter;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

sub pre_check_files {
	if ( $conf{readtrimfilter_use_precalc_5prime_trimming} eq 'yes' ) {
		if ( -s "$cwd/MOCAT.cutoff5prime" ) {
		}
		else {
			die "ERROR: Expected 5' cut off file in file in $cwd/MOCAT.cutoff5prime\n";
		}
	}

	# New: support pair.12 and single files
	if ( $conf{MOCAT_paired_end} eq "yes" ) {
		foreach my $sample (@samples) {
			chomp( my @fqs = `ls -1 $cwd/$sample/*fq $cwd/$sample/*fq.gz $cwd/$sample/*.fq.bz2 2>/dev/null | grep -v 'trimmed.filtered' | grep -P '\.single\.fq|\.pair\.1\.fq|\.pair\.2\.fq'` );
			if ( !( scalar @fqs % 3 == 0 ) && scalar @fqs > 0 ) {
				die "ERROR & EXIT: $sample seems to have .pair. and .single. files, but the number of files is not exactly 3 (number is " . ( scalar @fqs ) . "), please check that no files are mising";
			}
			if ( scalar @fqs % 3 == 0  && scalar @fqs > 0  ) {
				chomp( my @beginning = `ls -1 $cwd/$sample/*.pair.1.*fq.gz $cwd/$sample/*.pair.1.*.fq $cwd/$sample/*.pair.1.*fq.bz2 2>/dev/null | sed 's/.pair.1.*//'` );
				chomp( my $end       = `ls -1 $cwd/$sample/*.pair.1.*fq.gz $cwd/$sample/*.pair.1.*.fq $cwd/$sample/*.pair.1.*fq.bz2 2>/dev/null| sed 's/.*.pair.1.//' | sort -u` );
				foreach my $beginning (@beginning) {
					unless ( -e "$beginning.single.$end" ) {
						die "ERROR & EXIT: Expected $beginning.single.$end to exist. Do the files have different ending for the sample $sample?";
					}
					unless ( -e "$beginning.pair.1.$end" ) {
						die "ERROR & EXIT: Expected $beginning.single.$end to exist. Do the files have different ending for the sample $sample?";
					}
					unless ( -e "$beginning.pair.2.$end" ) {
						die "ERROR & EXIT: Expected $beginning.single.$end to exist. Do the files have different ending for the sample $sample?";
					}
				}
				$use3files{$sample} = 1;
			}
		}
	}

	# Check input files structure
	foreach my $sample (@samples) {
		unless ( $use3files{$sample} ) {
			chomp( my @fqs = `ls -1 $cwd/$sample/*fq $cwd/$sample/*fq.gz $cwd/$sample/*.fq.bz2 2>/dev/null | grep -v 'trimmed.filtered' | grep -v '.single.' | grep -v '.pair.'` );
			if ( $conf{MOCAT_paired_end} eq "yes" ) {
				if ( scalar @fqs % 2 == 1 || scalar @fqs == 0) {
					die <<"EOF";
ERROR & EXIT: Configuration 'MOCAT_paired_end' set to 'yes', but an odd number of files were found in sample $sample.
Make sure sample name files are in the format:
    - LANE_NAME1.1.fq(.gz)
    - LANE_NAME1.2.fq(.gz)
    - LANE_NAME2.1.fq(.gz)
    - LANE_NAME2.1.fq(.gz)
    
It is important with the .1. and .2. to identify the two paired end reads files.
EOF
				}
				foreach my $fq (@fqs) {
					$fq =~ m/(\S+)\.(1|2)\.(fq.gz|fq.bz2|fq)$/;
					unless ( defined $1 && defined $2 && defined $3 ) {
						die <<"EOF";
ERROR & EXIT: Files located in sample folder $sample seem to not be in format:
    - LANE_NAME1.1.fq(.gz)
    - LANE_NAME1.2.fq(.gz)
    - LANE_NAME2.1.fq(.gz)
    - LANE_NAME2.1.fq(.gz)

It is important with the .1. and .2. to identify the two paired end reads files.
EOF

					}
					if ( $2 eq '1' ) {
						unless ( -e "$1.1.$3" ) {
							die "ERROR & EXIT: Could not read file  $1.1.$3";
						}
						unless ( -e "$1.2.$3" ) {
							die "ERROR & EXIT: Found file:   $fq, also expected to be able to read file: $1.2.$3";
						}
					}

					if ( $2 eq '2' ) {
						unless ( -e "$1.2.$3" ) {
							die "ERROR & EXIT: Could not read file  $1.2.$3";
						}
						unless ( -e "$1.1.$3" ) {
							die "ERROR & EXIT: Found file:   $fq, also expected to be able to read file: $1.1.$3";
						}
					}
				}
			}
		}
	}
}

sub detect_format {

	# format = detect_format(sample)
	#
	#
	# Returns one of 's' (sanger), 'x' (solexa), or 'i' (solill)
	my $lane = shift;
	if ( $lane =~ /\.gz$/ ) {
		open( LANE, "gunzip -c $lane|" )
		  or die("Cannot open gunzip pipe to read '$lane': $!");
	}
	elsif ( $lane =~ /\.bz2$/ ) {
		open( LANE, "bunzip2 -c $lane|" )
		  or die("Cannot open bunzip2 pipe to read '$lane': $!");
	}
	else {
		open( LANE, '<', $lane )
		  or die("Cannot open '$lane': $!");
	}
	my $NR        = 0;
	my $is_solexa = 0;
	my $is_solill = 0;
	while (<LANE>) {
		++$NR;
		chomp;
		if ( ( $NR % 4 ) == 0 ) {
			if (/[!-:]/) { return 's' }

			#elsif (/[#-J]/) { return 'l' } # This doesn't have to be returned, because the v5+ script version should handle s format just fine even though new score scheme
			elsif (/[;-?]/) { $is_solexa = 1; }
			elsif (/[@-h]/) { $is_solill = 1; }

		}
		if ( $NR == 1000 ) { last }
	}
	close(LANE);
	if ($is_solexa) { return 'x'; }
	if ($is_solill) { return 'i'; }
	if ( $NR == 0 ) { return 's' }

	die "\nERROR & EXIT: Unknown format of files in sample $lane";
}

sub create_job {
	### DEFINE VARIABLES AND OPEN OUTPUT FILE ###
	my $job        = shift;
	my $processors = shift;
	my $jobfile    = "$jobdir/MOCATJob_$job\_$date";
	$ZIP =~ s/pigz.*/pigz -p $processors/;

	open JOB, '>', $jobfile or die "ERROR & EXIT: Cannot write $jobfile. Do you have permission to write in $jobdir?";

	my $now = localtime();
	print "$now: Creating $job jobs...";
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );

	### JOB SPECIFIC ###
	my $format;

	foreach my $sample (@samples) {
		&MOCATCore::mkdir_or_die("$temp_dir/$sample/temp");

		my @fqs = ();
		if ( $use3files{$sample} ) {
			chomp( @fqs = `ls -1 $cwd/$sample/*fq $cwd/$sample/*fq.gz $cwd/$sample/*.fq.bz2 2>/dev/null | grep -v 'trimmed.filtered' | grep -P '\.single\.|\.pair\.'` );
		}
		else {
			chomp( @fqs = `ls -1 $cwd/$sample/*fq $cwd/$sample/*fq.gz $cwd/$sample/*.fq.bz2 2>/dev/null | grep -v 'trimmed.filtered' | grep -v '.single.fq' | grep -v '.pair.'` );
			$use3files{$sample} = 0;
		}

		my @formats;
		if ( $conf{readtrimfilter_use_sanger_scale} eq 'yes' ) {
			foreach my $lane (@fqs) {
				$format = "s";
				push @formats, $format;
			}
		}
		elsif ( $conf{readtrimfilter_use_sanger_scale} eq 'auto' ) {
			foreach my $lane (@fqs) {
				$format = &detect_format($lane);
				push @formats, $format;
			}
		}
		elsif ( $conf{readtrimfilter_use_sanger_scale} eq 'no' ) {
			foreach my $lane (@fqs) {
				$format = "i";
				push @formats, $format;
			}
		}
		my $trim5prime = "yes";
		if ( $conf{readtrimfilter_trim_5_prime} eq "no" ) {
			$trim5prime = "no";
		}
		my $trim_paired_end = substr( $conf{MOCAT_paired_end}, 0, 1 );
		my $formats = join( '', @formats );

		if ($preprocessSOLiD) {
			print JOB "$scr_dir/filterConvertedColorspaceReads.pl " . " -indir=$cwd/$sample -outdir=$cwd/$sample -originaldir=$cwd/$sample/reads.raw.SOLiD -scoreThreshold=\"\$\" -lengthThreshold=40 " . " 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log ";
		}

		print JOB "exec 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log && ";
		print JOB "mkdir -p $temp_dir/$sample/temp && ";
		print JOB "$scr_dir/MOCATReadTrimFilter_aux.pl "
		  . "-sample $sample "
		  . "-trim_5prime_end $trim5prime "
		  . "-src_dir $scr_dir "
		  . "-paired_end_data $conf{MOCAT_paired_end} "
		  . "-file_formats_array $formats "
		  . "-length_cutoff $conf{readtrimfilter_length_cutoff} "
		  . "-qual_cutoff $conf{readtrimfilter_qual_cutoff} "
		  . "-solexaqa_or_fastx $conf{MOCAT_data_type} "
		  . "-bin_dir $bin_dir "
		  . "-cwd $cwd "
		  . "-use3files $use3files{$sample} "
		  . "-use_5prime_file $conf{readtrimfilter_use_precalc_5prime_trimming} "
		  . "-temp_dir $temp_dir "
		  . "-zcat \"$ZCAT\" "
		  . "2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log\n";
	}
	print " OK!\n";
	print localtime() . ": Temp directory is $temp_dir/SAMPLE/temp\n";
}

sub post_check_files {
	print localtime() . ": Checking files... ";
	my $p = 1;
	for my $sample (@samples) {
		if ( $conf{MOCAT_paired_end} eq "yes" ) {
			my @fqs = <$cwd/$sample/*.1.fq $cwd/$sample/*.1.fq.gz>;
			foreach my $fq (@fqs) {
				my @fq = split /\//, $fq;
				$fq = $fq[-1];
				$fq =~ s/.1.fq$//;
				$fq =~ s/.1.fq.gz$//;
				if ( $use3files{$sample} == 1 ) {
					$fq =~ s/\.pair//;
				}
				unless ( -s "$cwd/$sample/reads.processed.$conf{MOCAT_data_type}/$fq.pair.1.fq.gz" ) {
					print "\nERROR: Missing $cwd/$sample/reads.processed.$conf{MOCAT_data_type}/$fq.pair.1.fq.gz";
					$p = 0;
				}
				unless ( -s "$cwd/$sample/reads.processed.$conf{MOCAT_data_type}/$fq.pair.2.fq.gz" ) {
					print "\nERROR: Missing $cwd/$sample/reads.processed.$conf{MOCAT_data_type}/$fq.pair.2.fq.gz";
					$p = 0;
				}
				unless ( -s "$cwd/$sample/reads.processed.$conf{MOCAT_data_type}/$fq.single.fq.gz" ) {
					print "\nERROR: Missing $cwd/$sample/reads.processed.$conf{MOCAT_data_type}/$fq.single.fq.gz";
					$p = 0;
				}
			}
		}

		if ( $conf{MOCAT_paired_end} eq "no" ) {
			my @fqs = <$cwd/$sample/*fq $cwd/$sample/*fq.gz>;
			foreach my $fq (@fqs) {
				my @fq = split /\//, $fq;
				$fq = $fq[-1];
				$fq =~ s/.fq$//;
				$fq =~ s/.fq.gz$//;
				unless ( -s "$cwd/$sample/reads.processed.$conf{MOCAT_data_type}/$fq.pair.1.fq.gz" ) {
					print "\nERROR: Missing $cwd/$sample/reads.processed.$conf{MOCAT_data_type}/$fq.pair.1.fq.gz";
					$p = 0;
				}
				unless ( -s "$cwd/$sample/reads.processed.$conf{MOCAT_data_type}/$fq.pair.2.fq.gz" ) {
					print "\nERROR: Missing $cwd/$sample/reads.processed.$conf{MOCAT_data_type}/$fq.pair.2.fq.gz";
					$p = 0;
				}
				unless ( -s "$cwd/$sample/reads.processed.$conf{MOCAT_data_type}/$fq.single.fq.gz" ) {
					print "\nERROR: Missing $cwd/$sample/reads.processed.$conf{MOCAT_data_type}/$fq.single.fq.gz";
					$p = 0;
				}
			}
		}

	}
	if ( $p == 1 ) {
		print " OK!\n";
	}
	else {
		die "\nERROR & EXIT: Missing one or more files!";
	}
}

1;
