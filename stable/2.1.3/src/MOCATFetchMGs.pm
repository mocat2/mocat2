package MOCATFetchMGs;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;
use Math::Round qw(nearest_floor);

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

sub pre_check {
	if ( $assembly eq 'assembly' ) {
	}
	elsif ( $assembly eq 'assembly.revised' ) {
	}
	else {
		die "ERROR & EXIT: Option --fmg|fetch_mg 'ASSEMBLY ORIGIN' must be either of 'assembly' or 'assembly.revised'";
	}
	if ( $reads eq '' ) {
		die "ERROR & EXIT: Option --r|reads 'READS ORIGIN' must be either 'reads.processed', 'fastafile' or 'DATABASE NAME'";
	}
}

sub create_job {
	### DEFINE VARIABLES AND OPEN input FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?";
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
	
	print localtime() . ": Temp directory is $temp_dir\n";
	print localtime() . ": Creating $job jobs...\n";
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}
	my $temp = $assembly;
	$temp =~ s/assembly.revised/assembly_revision/;
	my $GP_min_len = "$temp\_scaftig_min_length";
	my $GP         = "$conf{gene_prediction_software}.$conf{$GP_min_len}";

	# Prechecks
	print localtime() . ": Checking HMM models...";
	my @hmms = (
		"COG0012.hmm", "COG0016.hmm", "COG0018.hmm", "COG0048.hmm", "COG0049.hmm", "COG0052.hmm", "COG0080.hmm", "COG0081.hmm", "COG0085.hmm", "COG0087.hmm", "COG0088.hmm", "COG0090.hmm", "COG0091.hmm", "COG0092.hmm", "COG0093.hmm", "COG0094.hmm", "COG0096.hmm", "COG0097.hmm", "COG0098.hmm", "COG0099.hmm",
		"COG0100.hmm", "COG0102.hmm", "COG0103.hmm", "COG0124.hmm", "COG0172.hmm", "COG0184.hmm", "COG0185.hmm", "COG0186.hmm", "COG0197.hmm", "COG0200.hmm", "COG0201.hmm", "COG0202.hmm", "COG0215.hmm", "COG0256.hmm", "COG0495.hmm", "COG0522.hmm", "COG0525.hmm", "COG0533.hmm", "COG0541.hmm", "COG0552.hmm"
	);
	foreach my $hmm (@hmms) {
		unless ( -e "$conf{MOCAT_dir}/lib/fetchMG/$hmm" ) {
			die "\nERROR & EXIT: Missing $conf{MOCAT_dir}/lib/fetchMG/$hmm";
		}
	}
	unless ( -e "$conf{MOCAT_dir}/lib/fetchMG/MG_BitScoreCutoffs.allhits.txt" ) {
		die "\nERROR & EXIT: Missing $conf{MOCAT_dir}/lib/fetchMG/MG_BitScoreCutoffs.allhits.txt";
	}
	print " OK!\n";
	
	### JOB SPECIFIC ###
	print localtime() . ": Continuing creating jobs...";
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
		my $tab           = "$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$conf{gene_prediction_software}.tab";
		my $faa           = "$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$conf{gene_prediction_software}.faa";
		my $fna           = "$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$conf{gene_prediction_software}.fna";
		my $input_folder  = "$cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$GP";
		my $output_folder = "$cwd/$sample/genes.fetched.MGs.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$GP/";

		unless ( -e "$input_folder/$faa" ) {
			die "\nERROR & EXIT: Missing $input_folder/$faa\nDid you perform gene predicion on these samples?";
		}
		unless ( -e "$input_folder/$fna" ) {
			die "\nERROR & EXIT: Missing $input_folder/$fna\nDid you perform gene predicion on these samples?";
		}
		unless ( -e "$input_folder/$tab" ) {
			die "\nERROR & EXIT: Missing $input_folder$tab\nDid you perform gene predicion on these samples?";
		}
		
		
		
		print JOB "mkdir -p $output_folder &&";
		print JOB "grep -P '\\tcomplete\$' $input_folder/$tab | cut -f 1 > $temp_dir/$sample/temp/complete.genes.list && ";
		print JOB "$bin_dir/faSomeRecords $input_folder/$faa $temp_dir/$sample/temp/complete.genes.list $temp_dir/$sample/temp/$faa && ";
		print JOB "$bin_dir/faSomeRecords $input_folder/$fna $temp_dir/$sample/temp/complete.genes.list $temp_dir/$sample/temp/$fna && ";		
		print JOB "$scr_dir/MOCATFetchMGs03.pl -m extraction -c all -o $output_folder -l $conf{MOCAT_dir}/lib/fetchMG -b $conf{MOCAT_dir}/lib/fetchMG/MG_BitScoreCutoffs.allhits.txt -t $processors -x $bin_dir $temp_dir/$sample/temp/$faa $LOG && ";
		print JOB "rm -fr $output_folder/temp $output_folder/hmmResults $temp_dir/$sample/temp/complete.gene.list $temp_dir/$sample/temp/$faa* $temp_dir/$sample/temp/$fna*\n"

	}
	print " OK!\n";
}

1;
