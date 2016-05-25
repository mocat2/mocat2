package MOCATAssemblyRevision;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;
use Math::Round qw(nearest_floor);

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

sub pre_check {
	if ( $reads eq '' ) {
		die "ERROR & EXIT: Option -r|reads 'READ ORIGIN 'required. Please set it to either 'reads.processed', 'DATABASE NAME' or 'FASTA FILE'";
	}
}

sub create_job {
	### DEFINE VARIABLES AND OPEN OUTPUT FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	my @databases;
	if ( exists $_[2] ) {
		@databases = @{ $_[2] };
	}
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}
	
		
	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?";
	print localtime() . ": Creating $job jobs...";
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );

	### JOB SPECIFIC ###
	foreach my $sample (@samples) {
		
		print JOB "exec 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log && mkdir -p $temp_dir/$sample/temp; ";
		
		my $LOG = " 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log ";

		my ($max, $average, $kmer, $kmer_avg);
		my $max_max = 0;
		if (scalar @reads == 1) {
			$reads = $reads[0];
			( $max, $average, $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-ar" );
		} else {
			foreach $reads (@reads) {
				( $max, $average, $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-ar" );
				$kmer_avg += $kmer;
				if ($max > $max_max) {
					$max_max = $max;
				}
			}
			$max = $max_max;
			$kmer = nearest_floor(1,$kmer_avg / scalar @reads);
		}
		my $reads = join( "_AND_", @reads );

		# Check files
		unless ( -s "$cwd/$sample/assembly.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.assembly.$reads.$conf{MOCAT_data_type}.K$kmer.scaftig.gz" ) {
			die "\nERROR & EXIT: Missing assembly files for $sample.\nExpected this file to exist: $cwd/$sample/assembly.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.assembly.$reads.$conf{MOCAT_data_type}.K$kmer.scaftig.gz\nDid you set -r to one of 'reads.processed', 'FASTA FILE' or 'DATABASE NAME'?";
		}
		unless ( -s "$cwd/$sample/stats/$sample.assembly.$reads.$conf{MOCAT_data_type}.K$kmer.inserts.stats" ) {
			die "\nERROR & EXIT: Missing insert size file for $sample: $cwd/$sample/stats/$sample.assembly.$reads.$conf{MOCAT_data_type}.K$kmer.inserts.stats\nDid you set -r to one of 'reads.processed', 'FASTA FILE' or 'DATABASE NAME'?";
		}
		
		
		# Prepare files
		my @lanes;
		if ( $reads eq 'reads.processed' ) {
			@lanes = <$cwd/$sample/reads.processed.$conf{MOCAT_data_type}/*pair.*.fq.gz $cwd/$sample/reads.processed.$conf{MOCAT_data_type}/*single.fq.gz>;
		}
		else {
			foreach my $r (@reads) {
				my @l = <$cwd/$sample/reads.$read_type.$r.$conf{MOCAT_data_type}/*pair.*.fq.gz $cwd/$sample/reads.$read_type.$r.$conf{MOCAT_data_type}/*single.fq.gz>;
				push @lanes, @l;
			}
		}
		chomp(@lanes);
		
		system ("mkdir -p $cwd/$sample/temp");
		open OUT2, '>', "$cwd/$sample/temp/AC.lanes";
		print OUT2 join( "\n", @lanes );
		close OUT2;
		open OUT2, '>', "$cwd/$sample/temp/AC.lengths";
		print OUT2 "$sample\t$max\t$average\n";
		close OUT2;
		
		open IN2,  '<', "$cwd/$sample/stats/$sample.assembly.$reads.$conf{MOCAT_data_type}.K$kmer.inserts.stats";
		open OUT2, '>', "$cwd/$sample/temp/AC.inserts";
		while (<IN2>) {
			chomp($_);
			print OUT2 "$sample\_$_\n";
		}
		close IN2;
		close OUT2;
		
		
		# Print job
		print JOB " $ZIP -dc $cwd/$sample/assembly.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.assembly.$reads.$conf{MOCAT_data_type}.K$kmer.scaftig.gz > $cwd/$sample/assembly.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.assembly.$reads.$conf{MOCAT_data_type}.K$kmer.scaftig && ";
		print JOB " rm -rf $temp_dir/$sample/temp/assembly.revised && mkdir -p $temp_dir/$sample/temp/assembly.revised && ";
		print JOB " $scr_dir/MOCATAssemblyRevision.pl '$cwd/$sample/assembly.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.assembly.$reads.$conf{MOCAT_data_type}.K$kmer.scaftig' '$cwd/$sample/temp/AC.lanes' '$cwd/$sample/temp/AC.inserts' '$cwd/$sample/temp/AC.lengths' '$temp_dir/$sample/temp/assembly.revised' '$bin_dir' '$scr_dir' '$processors' '$conf{assembly_revision_scaftig_min_length}' '$sample' '$cwd/$sample/stats/$sample.assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer.scaftig.stats' $LOG && ";
		print JOB " cp $temp_dir/$sample/temp/assembly.revised/contig_correct.stat $cwd/$sample/stats/$sample.assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer.baseErrorAndIndelError.stats $LOG && mkdir -p $cwd/$sample/assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer/ && if [ -e $temp_dir/$sample/temp/assembly.revised/break_contig/$sample.scaftig.more$conf{assembly_revision_scaftig_min_length}.revised.break.fa ]; then echo \"Assembly revision OK.\"; else echo \"Assembly revision failed. No scaftigs longer that specified cutoff were created. Rerun and allow shorter scaftigs to be created. Possibly your sample it too small to have assembly revision run properly for this sample.\"; exit 111; fi && rsync -av --remove-sent-files $temp_dir/$sample/temp/assembly.revised/break_contig/$sample.scaftig.more$conf{assembly_revision_scaftig_min_length}.revised.break.fa $cwd/$sample/assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer.scaftig $LOG && ";
		print JOB " rm -fr $temp_dir/$sample/temp/assembly.revised && $ZIP -$ziplevel -f $cwd/$sample/assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer.scaftig\n";
	}
	print " OK!\n";
	print localtime() . ": Temp directory is $temp_dir/SAMPLE/temp\n";
}

sub post_check_files {
	
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}
	

	foreach my $sample (@samples) {

		# Get Kmer
		my ($max, $average, $kmer, $kmer_avg);
		my $max_max = 0;
		if (scalar @reads == 1) {
			$reads = $reads[0];
			( $max, $average, $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-ar" );
		} else {
			foreach $reads (@reads) {
				( $max, $average, $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-ar" );
				$kmer_avg += $kmer;
				if ($max > $max_max) {
					$max_max = $max;
				}
			}
			$max = $max_max;
			$kmer = nearest_floor(1,$kmer_avg / scalar @reads);
		}
		my $reads = join( "_AND_", @reads );

		unless ( -s "$cwd/$sample/assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer.scaftig.gz" ) {
			die "ERROR & EXIT: Missing revised assembly scaftig file: $cwd/$sample/assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer.scaftig.gz";
		}
		chomp(my $entries = `zgrep -c '>' $cwd/$sample/assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer.scaftig.gz`);
		if ($entries eq '0') {
			die "ERROR & EXIT: Revised assembly scaftig file $cwd/$sample/assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer.scaftig.gz has 0 scaftigs.\nMaybe this sample is too small for the assembly to be revised? It is normal that small sample assemblies cannot be revised.\nIt could also be that the sample had no revised scaftigs longer than $conf{assembly_revision_scaftig_min_length} (assembly_revision_scaftig_min_length setting in config file)";
		}
		unless ( -s "$cwd/$sample/stats/$sample.assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer.baseErrorAndIndelError.stats" ) {
			die "ERROR & EXIT: Missing revised assembly statistics file: $cwd/$sample/stats/$sample.assembly.revised.$reads.$conf{MOCAT_data_type}.K$kmer.baseErrorAndIndelError.stats";
		}
	}
}

1;
