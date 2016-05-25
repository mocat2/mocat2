package MOCATStatsFastqc;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

sub create_job {
	### DEFINE VARIABLES AND OPEN OUTPUT FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?";
	
	print localtime() . ": Creating $job jobs...";
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );

	### JOB SPECIFIC ###
	my $db = $databases[0];
	foreach my $sample (@samples) {
		
		print JOB "mkdir -p $temp_dir/$sample/temp; ";
		
		my $LOG = " 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log ";

		chomp( my @fqs = `ls -1 $cwd/$sample/*fq $cwd/$sample/*fq.gz 2>/dev/null | grep -v 'trimmed.filtered' | grep -v '.single.' | grep -v '.pair.'` );
		foreach my $lane (@fqs) {
			print JOB "$ext_dir/fastqc/fastqc -o $cwd/$sample/ -f fastq -t 1 $lane $LOG\n";
		}
	}
	print " OK!\n";

	#print localtime().": Temp directory is $temp_dir/SAMPLE/temp\n";
}

sub post_check_files {
	print localtime() . ": Checking files...";
	foreach my $sample (@samples) {
		
		chomp( my @fqs = `ls -1 $cwd/$sample/*fq $cwd/$sample/*fq.gz 2>/dev/null | grep -v 'trimmed.filtered' | grep -v '.single.' | grep -v '.pair.'` );
		foreach my $lane (@fqs) {
			$lane =~ s/.fq(.gz)//;
			unless ( -s "$lane\_fastqc.html" ) {
				die "\nERROR & EXIT: Missing FastQC summary file: $lane\_fastqc.html\nDid FastQC run correctly?";
			}
		}
	}
	print " OK!\n";
}

1;
