package MOCATGenePredictionMetaGeneMark;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;
use Math::Round qw(nearest_floor);

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

sub pre_check {
	unless ( -s "$conf{MOCAT_dir}/ext/metagenemark/gmhmmp" ) {
		die "ERROR & EXIT: Missing external MetaGeneMark executable.\nPlease download MetaGeneMark and copy the content of the package file into $conf{MOCAT_dir}/ext/metagenemark.";
	}
	unless ( -s $ENV{"HOME"} . "/.gm_key" ) {
		print "Missing MetaGeneMark licence key file in " . $ENV{"HOME"} . "/.gm_key. Trying to copy file...";
		system "cp $conf{MOCAT_dir}/ext/metagenemark/gm_key " . $ENV{"HOME"} . "/.gm_key 2>/dev/null >/dev/null";
		if ( -s $ENV{"HOME"} . "/.gm_key" ) {
			print " Success!\n";
		}
		else {
			print " Fail!\n";
			die "ERROR & EXIT: Missing gm_key file in $conf{MOCAT_dir}/ext/metagenemark/gm_key.";
		}
	}
	if ( $assembly eq 'assembly' ) {
	}
	elsif ( $assembly eq 'assembly.revised' ) {
	}
	else {
		die "ERROR & EXIT: Option --gp|gene_prediction 'ASSEMBLY ORIGIN' must be either of 'assembly' or 'assembly.revised'";
	}
	if ( $reads eq '' ) {
		die "ERROR & EXIT: Option --r|reads 'READS ORIGIN' must be either 'reads.processed', 'fastafile' or 'DATABASE NAME'";
	}
}

sub create_job {
	### DEFINE VARIABLES AND OPEN OUTPUT FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?";
	
	print localtime() . ": Creating $job jobs...";
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
	my $temp = $assembly;
	$temp =~ s/assembly.revised/assembly_revision/;
	my $GP_min_len = "$temp\_scaftig_min_length";
	my $GP         = "MetaGeneMark.$conf{$GP_min_len}";


	### JOB SPECIFIC ###
	foreach my $sample (@samples) {
		
		print JOB "exec 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log && ";
		print JOB "mkdir -p $temp_dir/$sample/temp; ";
		
		my $LOG = " 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log ";
		
		# Support for multiple DBs
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

		# Create job
		my $lst     = "$cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$GP/$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.MetaGeneMark.lst";
		my $nt      = "$cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$GP/$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.MetaGeneMark.fna";
		my $aa      = "$cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$GP/$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.MetaGeneMark.faa";
		my $scaftig = "$cwd/$sample/$assembly.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$conf{gene_prediction_input}";

		unless ( -e "$scaftig.gz" ) {
			die "\nERROR & EXIT: Missing scaftig file $scaftig.gz\nDid you specify correct parameters for --gene_prediction and --reads?\nIn the config file you specified to predict genes on the '$conf{gene_prediction_input}' file.";
		}

		print JOB "$ZIP -dc $scaftig.gz > $scaftig &&";
		print JOB " mkdir -p $cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$GP && $conf{MOCAT_dir}/ext/metagenemark/gmhmmp -a -d -f L -m $conf{MOCAT_dir}/ext/metagenemark/MetaGeneMark_v1.mod -o $lst $scaftig $LOG && ";
		print JOB " $scr_dir/MOCATGenePredictionMetaGeneMark_aux.pl $sample $conf{MOCAT_data_type} $assembly $reads $cwd $conf{gene_prediction_input} $kmer $GP $LOG && rm -f $scaftig\n";
	}
	print " OK!\n";
}

sub post_check_files {
	print localtime() . ": Checking files... ";
	my $temp = $assembly;
	$temp =~ s/assembly.revised/assembly_revision/;
	my $GP_min_len = "$temp\_scaftig_min_length";
	my $GP         = "MetaGeneMark.$conf{$GP_min_len}";

	foreach my $sample (@samples) {

		# Support for multiple DBs
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


		# Create job
		my $lst     = "$cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$GP/$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.MetaGeneMark.lst";
		my $nt      = "$cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$GP/$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.MetaGeneMark.fna";
		my $aa      = "$cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$GP/$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.MetaGeneMark.faa";
		my $tab     = "$cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$GP/$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.MetaGeneMark.tab";
		
		unless ( -s $nt ) {
			die "\nERROR & EXIT: Missing predicted genes nucleotide sequence file $nt";
		}
		unless ( -s $aa ) {
			die "\nERROR & EXIT: Missing predicted proteins amino acid sequence file $aa";
		}
		unless ( -s $tab ) {
			die "\nERROR & EXIT: Missing predicted sequence file $tab";
		}
		unless ( -s $lst ) {
			die "\nERROR & EXIT: Missing predicted genes lst file $lst";
		}

	}
	print " OK!\n";
}

1;
