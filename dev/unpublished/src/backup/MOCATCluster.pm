package MOCATCluster;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

sub create_job {
	### DEFINE VARIABLES AND OPEN OUTPUT FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?\n";
	print localtime() . ": Creating $job jobs...";

	### JOB SPECIFIC ###
	my $db = $databases[0];
	foreach my $sample (@samples) {
		my $LOG = " 2>> $cwd/logs/MOCATJob_$job.$sample.log >> $cwd/logs/MOCATJob_$job.$sample.log ";

		( my $max, my $avg, my $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-r" );
		my $temp = "$assembly\_scaftig_min_length";
		$temp =~ s/assembly.revised/assembly_revision/;
		my $input = "$cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$conf{gene_prediction_software}.$conf{$temp}/$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$conf{gene_prediction_software}.fna";
		unless ( -e $input ) {
			die "\nERROR & EXIT: Missing expected input file $input\n";
		}
		my $output = "$cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$conf{gene_prediction_software}.$conf{$temp}/$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$conf{gene_prediction_software}.nr.fna";
		print JOB "$conf{cluster_bin} $conf{cluster_cmd} -T $processors -i $input -o $output $LOG\n";

	}
	print " OK!\n";

	#print localtime().": Temp directory is $temp_dir/SAMPLE/temp\n";
}

sub post_check_files {
	print localtime() . ": Checking files...";
	foreach my $sample (@samples) {
		( my $max, my $avg, my $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-r" );
		my $temp = "$assembly\_scaftig_min_length";
		$temp =~ s/assembly.revised/assembly_revision/;
		my $file = "$cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$conf{gene_prediction_software}.$conf{$temp}/$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$conf{gene_prediction_software}.nr.fna";
		unless ( -e $file ) {
			die "\nERROR & EXIT: Missing expected output file $file\n";
		}
	}
	print " OK!\n";
}

1;
