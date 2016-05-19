package MOCATRExample;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

sub pre_check_files {
}

sub create_job {
	### DEFINE VARIABLES AND OPEN input FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?\n";

	print localtime() . ": Creating $job jobs...\n";
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}

	### JOB SPECIFIC ###
	my $db = $databases[0];
	foreach my $sample (@samples) {

	}
	print " OK!\n";
	print localtime() . ": Temp directory is $temp_dir/SAMPLE/temp\n";
}

sub post_check_files {
	print localtime() . ": Checking files... ";
	print " OK!\n";
}

1;
