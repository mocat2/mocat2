package MOCATAnnotateUsingCarmaFile;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

sub create_job {
	
	### DEFINE VARIABLES AND OPEN input FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	my @databases;
	if ( exists $_[2] ) {
		@databases = @{ $_[2] };
	}
	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?\n";
	print localtime() . ": Creating $job jobs...";
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}

	### JOB SPECIFIC ###
	my @temp = split /\//, $carma_file;
	$carma_file = "$cwd/annotate.$temp[-1]/$temp[-1].carma.sorted";
	unless ( -e $carma_file ) {
			die "\nERROR & EXIT: Missing CARMA input file: $carma_file\n";
		}
		
	foreach my $sample (@samples) {
		my $LOG = " 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.log >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.log ";
		foreach my $start ( "base.coverage", "insert.coverage" ) {
			my $db = join( "_AND_", @do_annotate_using_carma_file );
			my $folder    = "$cwd/$sample/$start.$db.$conf{MOCAT_data_type}";
			my $input     = "$sample.filtered.$read_type.$reads.on.$db.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$start";
			unless ( -e "$folder/$input.norm" ) {
				die "\nERROR & EXIT: Missing coverage calculation results file: $folder/$input.norm\n";
			}
			unless ( -e "$folder/$input.count" ) {
				die "\nERROR & EXIT: Missing coverage calculation results file: $folder/$input.count\n";
			}
			foreach my $end ( "norm", "count" ) {
				my $inputfolder = "$cwd/$sample/taxonomic.compositions/$input/$start.$end";
				print JOB "mkdir -p $inputfolder && ";
				print JOB "paste $data_dir/$db.rownames $folder/$input.$end > $inputfolder/tax.tmp.$end.1 && ";
				print JOB "awk 'NR>2' $inputfolder/tax.tmp.$end.1 > $inputfolder/tax.tmp.$end.2 && cut -f 1 $carma_file > $inputfolder/tax.tmp.$end.3 && ";
				print JOB "paste $carma_file $inputfolder/tax.tmp.$end.2 > $inputfolder/tax.tmp.$end.5 && rm -f $inputfolder/$sample.$end.taxonomic.composition.tab && ";
				print JOB "$conf{MOCAT_DEV_DIR}/MOCATAnnotateUsingCarmaFile_aux.pl -o $inputfolder/$sample.$end.taxonomic.composition.tab -g $inputfolder -p $bin_dir/gnuplot $inputfolder/tax.tmp.$end.5 $LOG && ";
				print JOB "rm -f $inputfolder/tax.tmp.* && ";
			}
		}
		print JOB "echo COMPLETE $LOG\n";
	}
	close JOB;
	print " OK!\n";
}

1;
