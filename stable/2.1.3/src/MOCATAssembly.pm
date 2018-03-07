package MOCATAssembly;
use strict;
use warnings;
use MOCATCore;
use Math::Round qw(nearest_floor);
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

sub pre_check_files {
	my @do_assembly = @reads;
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}
	print localtime() . ": Checking files... ";
	if ( $conf{MOCAT_paired_end} eq 'yes' ) {
		if ( $conf{assembly_calculate_insert_size_using} eq 'assembly' ) {
			unless ( -s "$cwd/MOCAT.preliminary_insert_sizes" ) {
				die "\nERROR & EXIT: Missing manually entered preliminary insert sizes in file '$cwd/MOCAT.preliminary_insert_sizes'.\nEither specify assembly_calculate_insert_sizes to be 'mapping' in the config file, or, for each sample specify <SAMPLE><TAB><LANE><TAB><insert size> in the file '$cwd/MOCAT.preliminary_insert_sizes'.";
			}
		}
	}
	if ( $do_assembly[0] eq '' ) {
		die "\nERROR & EXIT: Missing option -r|reads. Please set this option to either of: 'reads.processed', 'DATABASE NAME' or 'fastafile'.";
	}
	my $F;
	foreach my $sample (@samples) {
		my @F;
		if ( $do_assembly[0] eq 'reads.processed' ) {
			###NO ALL###      $F = "$cwd/$sample/reads.processed.$conf{MOCAT_data_type}/$sample.all.fq.gz"
			@F = `ls -1 $cwd/$sample/reads.processed.$conf{MOCAT_data_type}/*pair*fq.gz $cwd/$sample/reads.processed.$conf{MOCAT_data_type}/*single*fq.gz 2>/dev/null `;
		}
		else {
			foreach my $reads (@do_assembly) {
				###NO ALL###      $F = "$cwd/$sample/reads.screened.$reads.$conf{MOCAT_data_type}/$sample.screened.$reads.all.fq.gz"
				@F = `ls -1 $cwd/$sample/reads.$read_type.$reads.$conf{MOCAT_data_type}/*pair*fq.gz $cwd/$sample/reads.$read_type.$reads.$conf{MOCAT_data_type}/*single*fq.gz 2>/dev/null`;
			}
			unless ( ( scalar @F % 3 ) == 0 && scalar @F >= 3 ) {
				my @statsfiles = `ls -1 $cwd/$sample 2>/dev/null`;
				my $line       = "";
				foreach my $sf (@statsfiles) {
					if ( $sf =~ m/^reads.screened/ || $sf =~ m/^reads.processed/ ) {
						chomp($sf);
						$sf =~ s/^reads.screened.//;
						$sf =~ s/^reads.extracted.//;
						$sf =~ s/\.solexaqa$//;
						$sf =~ s/\.fastx$//;
						$line = $line . "- $sf\n";
					}

				}
				die "\nERROR & EXIT: Missing pair and/or single reads files for $sample\n-r|reads '$reads' is not any of the following:\na) 'reads.processed' (and the reads have been trimmed and filtered)\nb) 'DATABASE NAME' (and the reads have been screened using this database)\nc) 'FASTA FILE' (and the reads have been screened using a fasta file).\nPerhaps you meant either of these for the -r option:\n$line";
			}
		}
	}
	print " OK!\n";
}

sub create_job {
	### DEFINE VARIABLES AND OPEN OUTPUT FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	my @do_assembly = @reads;
	
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}

	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?";
	print localtime() . ": Creating $job jobs...";
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );

	### JOB SPECIFIC ###
	my $db = $databases[0];

	foreach my $sample (@samples) {
		
		print JOB "exec 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log && mkdir -p $temp_dir/$sample/temp; ";

		my ($max, $average, $kmer, $kmer_avg);
		my $max_max = 0;
		if (scalar @do_assembly == 1) {
			$reads = $do_assembly[0];
			( $max, $average, $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-a" );
		} else {
			foreach $reads (@do_assembly) {
				( $max, $average, $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-a" );
				$kmer_avg += $kmer;
				if ($max > $max_max) {
					$max_max = $max;
				}
			}
			$max = $max_max;
			$kmer = nearest_floor(1,$kmer_avg / scalar @do_assembly);
		}

		my $denovo;
		my $soapSpecific;

		if ( $conf{assembly_soap_version} eq "2.04" ) {
			my $OSX = '';
			if ( $systemType =~ m/Darwin/ ) {
				$OSX = '_OSX';
			}
			if ( $kmer <= 63 ) {
				$denovo = "SOAPdenovo2/SOAPdenovo-63mer$OSX";
			}
			elsif ( $kmer <= 127 ) {
				$denovo = "SOAPdenovo2/SOAPdenovo-127mer$OSX";
			}
			else {
				die "ERROR: Kmer for sample $sample is > 127. Kmer = $kmer.";
			}
		}
		elsif ( $conf{assembly_soap_version} eq "1.05" ) {
			$soapSpecific = "-L $max";
			if ( $kmer <= 31 ) {
				$denovo = "soap1.05/SOAPdenovo31mer";
			}
			elsif ( $kmer <= 63 ) {
				$denovo = "soap1.05/SOAPdenovo63mer";
			}
			elsif ( $kmer <= 127 ) {
				$denovo = "soap1.05/SOAPdenovo127mer";
			}
			else {
				die "ERROR: Kmer for sample $sample is > 127. Kmer = $kmer.";
			}
		}
		elsif ( $conf{assembly_soap_version} eq "1.06" ) {
			if ( $systemType =~ m/Darwin/ ) {

				# For some reason the original 1.06 doesn't work, so we'll use the nobamaio
				# version, which is fine for us. Only diff is that the nobamaio doesn't produce
				# specific output files we don't use anyway.

				$soapSpecific = "-L $average -F";
				if ( $kmer <= 31 ) {
					$denovo = "soap1.06OSXnobamaio/SOAPdenovo-31mer";
				}
				elsif ( $kmer <= 63 ) {
					$denovo = "soap1.06OSXnobamaio/SOAPdenovo-63mer";
				}
				elsif ( $kmer <= 127 ) {
					$denovo = "soap1.06OSXnobamaio/SOAPdenovo-127mer";
				}
				else {
					die "ERROR: Kmer for sample $sample is > 127. Kmer = $kmer.";
				}
			}
			else {
				$soapSpecific = "-L $average -F";
				if ( $kmer <= 31 ) {
					$denovo = "soap1.06/SOAPdenovo-31mer-static";
				}
				elsif ( $kmer <= 63 ) {
					$denovo = "soap1.06/SOAPdenovo-63mer-static";
				}
				elsif ( $kmer <= 127 ) {
					$denovo = "soap1.06/SOAPdenovo-127mer-static";
				}
				else {
					die "ERROR: Kmer for sample $sample is > 127. Kmer = $kmer.";
				}
			}
		}
		else {
			die "ERROR & EXIT: Incorrect SOAP version specified. Should be one of 'soap1.05' or 'soap1.06'.";
		}

		my $reads2 = join( "_AND_", @do_assembly );
		my $F;
		my $reads_pre  = "";
		my $reads_full = "";
		if ( $reads eq 'reads.processed' ) {
			$reads_pre  = 'assembly.processed';
			$reads_full = 'reads.processed';
			$F          = "$cwd/$sample/reads.processed.$conf{MOCAT_data_type}/";
		}
		else {
			$reads_pre  = "assembly.$reads2";
			$reads_full = "reads.$read_type.$reads2";
			foreach my $reads (@do_assembly) {
				$F .= "$cwd/$sample/reads.$read_type.$reads.$conf{MOCAT_data_type}/,";
			}
			$F =~ s/,$//;
			$F = "'$F'";
		}

		my $O   = "$cwd/$sample/assembly.$reads2.$conf{MOCAT_data_type}.K$kmer";
		my $A   = "assembly.$reads2.$conf{MOCAT_data_type}.K$kmer";
		my $LOG = " 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log ";

		my $insert_calc_way = $conf{assembly_calculate_insert_size_using};

		print JOB "perl $scr_dir/MOCATAssembly_aux.pl $F $data_dir/$db.index $max $bin_dir $kmer $cwd $sample_file $temp_dir $conf{MOCAT_data_type} $reads_full $reads2 $insert_calc_way $scr_dir $bin_dir $bin_dir/$denovo $sample $processors $conf{assembly_scaftig_min_length} \"$ZCAT\" $conf{MOCAT_paired_end} \"$LOG\" $LOG && ";

		if ( $conf{assembly_soap_version} eq "1.05" ) {
			print JOB "$bin_dir/$denovo all -s $O/$sample.config -K $kmer -o $temp_dir/$sample/temp/$sample -D 1 -M 3 -u $soapSpecific -p $processors $LOG && ";
		}
		if ( $conf{assembly_soap_version} eq "1.06beta" ) {
			print JOB "$bin_dir/$denovo pregraph -s $O/$sample.config -K $kmer -o $temp_dir/$sample/temp/$sample -p $processors $LOG && $bin_dir/$denovo contig -g $temp_dir/$sample/temp/$sample -D 1 -M 3 $LOG && $bin_dir/$denovo map -s $O/$sample.config -p $processors -g $temp_dir/$sample/temp/$sample $LOG && $bin_dir/$denovo scaff -L $average -p $processors -F -u -g $temp_dir/$sample/temp/$sample $LOG && ";
		}
		if ( $conf{assembly_soap_version} eq "1.06" ) {
			print JOB "$bin_dir/$denovo all -s $O/$sample.config -K $kmer -o $temp_dir/$sample/temp/$sample -p $processors -D 1 -M 3 -L $average -F -u $LOG && ";
		}
		if ( $conf{assembly_soap_version} eq "1.06nobamaio" ) {
			print JOB "$bin_dir/$denovo all -s $O/$sample.config -K $kmer -o $temp_dir/$sample/temp/$sample -p $processors -D 1 -M 3 -L $average -F -u $LOG && ";
		}
		if ( $conf{assembly_soap_version} eq "2.04" ) {
			print JOB "$bin_dir/$denovo all -s $O/$sample.config -K $kmer -o $temp_dir/$sample/temp/$sample -p $processors -D 1 -M 3 -L $average -F -u $LOG && ";
		}
		
		print JOB
		  " perl $scr_dir/MOCATAssembly_getScaf.pl $conf{assembly_scaftig_min_length} $temp_dir/$sample/temp/$sample.scafSeq $sample > $O/$sample.scaftig && perl $scr_dir/MOCATAssembly_stats.pl $temp_dir/$sample/temp/$sample.contig $temp_dir/$sample/temp/$sample.scafSeq $O/$sample.scaftig $O/$sample.stats $scr_dir && rsync -av --remove-sent-files $temp_dir/$sample/temp/*scafSeq $temp_dir/$sample/temp/*scafStatistics $temp_dir/$sample/temp/*contig $O $LOG && $scr_dir/MOCATAssembly_rmAssembly.pl $temp_dir/$sample/temp $LOG && mv $O/$sample.config $O/$sample.$A.config && mv $O/$sample.contig $O/$sample.$A.contig && mv $O/$sample.scafStatistics $cwd/$sample/stats/$sample.$A.scaftig.stats && mv $O/$sample.scafSeq $O/$sample.$A.scafSeq && mv $O/$sample.scaftig $O/$sample.$A.scaftig && mv $O/$sample.stats $cwd/$sample/stats/$sample.$A.assembly.stats $LOG && $ZIP -f -$ziplevel $O/$sample.$A.scaftig && $ZIP -f -$ziplevel $O/$sample.$A.contig && $ZIP -f -$ziplevel $O/$sample.$A.scafSeq\n";

	}
	print " OK!\n";
	print localtime() . ": Temp directory is $temp_dir/SAMPLE/temp\n";
}

sub post_check_files {
	print localtime() . ": Checking files... ";
	my @do_assembly = @reads;
	my $c = 1;
	my $t = "\nERROR: Missing file: ";
	my $lf;
	my $m;
	my @c;
	my $F;
	my $A;

	foreach my $sample (@samples) {

		my ($max, $average, $kmer, $kmer_avg);
		my $max_max = 0;
		if (scalar @do_assembly == 1) {
			$reads = $do_assembly[0];
			( $max, $average, $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-a" );
		} else {
			foreach $reads (@do_assembly) {
				( $max, $average, $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-a" );
				$kmer_avg += $kmer;
				if ($max > $max_max) {
					$max_max = $max;
				}
			}
			$max = $max_max;
			$kmer = nearest_floor(1,$kmer_avg / scalar @do_assembly);
		}
		my $reads2 = join( "_AND_", @do_assembly );
		
		$F  = "$cwd/$sample/assembly.$reads2.$conf{MOCAT_data_type}.K$kmer/$sample.";
		$A  = "assembly.$reads2.$conf{MOCAT_data_type}.K$kmer";
		$lf = "$A.scaftig.gz";
		$m  = "$t$F$lf";
		@c  = `ls $F$lf 2>/dev/null`;
		if ( scalar @c == 0 ) { print $m; $c = 0; }
		$lf = "$A.scafSeq.gz";
		$m  = "$t$F$lf";
		@c  = `ls $F$lf 2>/dev/null`;
		if ( scalar @c == 0 ) { print $m; $c = 0; }
		$lf = "$A.contig.gz";
		$m  = "$t$F$lf";
		@c  = `ls $F$lf 2>/dev/null`;
		if ( scalar @c == 0 ) { print $m; $c = 0; }
		$lf = "$A.config";
		$m  = "$t$F$lf";
		@c  = `ls $F$lf 2>/dev/null`;
		if ( scalar @c == 0 ) { print $m; $c = 0; }
		$F  = "$cwd/$sample/stats/$sample.";
		$lf = "$A.assembly.stats";
		$m  = "$t$F$lf";
		@c  = `ls $F$lf 2>/dev/null`;
		if ( scalar @c == 0 ) { print $m; $c = 0; }
		$lf = "$A.scaftig.stats";
		$m  = "$t$F$lf\n";
		@c  = `ls $F$lf 2>/dev/null`;
		if ( scalar @c == 0 ) { print $m; $c = 0; }
	}
	if ( $c == 0 ) {
		die "\nERROR & EXIT: One or more files are missing in the assembly folders.";
	}
	print " OK!\n";
}

1;
