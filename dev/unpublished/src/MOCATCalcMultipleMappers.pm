package MOCATCalcMultipleMappers;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

sub create_job {
	### DEFINE VARIABLES AND OPEN OUTPUT FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}

	unless (     $multiple_mappers_summarize_level eq 'none'
		|| $multiple_mappers_summarize_level eq 'taxid'
		|| $multiple_mappers_summarize_level eq 'kingdom'
		|| $multiple_mappers_summarize_level eq 'phylum'
		|| $multiple_mappers_summarize_level eq 'class'
		|| $multiple_mappers_summarize_level eq 'order'
		|| $multiple_mappers_summarize_level eq 'family'
		|| $multiple_mappers_summarize_level eq 'genus'
		|| $multiple_mappers_summarize_level eq 'species'
		|| $multiple_mappers_summarize_level eq 'specI_clusters' )
	{
		die "\nERROR & EXIT: -summarize_level has to be of: 'none', 'taxid', 'kingdom', 'phylum', 'class' 'order', 'family', 'genus', 'species', 'specI_clusters'";
	}
	
	my $databases = join( "_AND_", @databases );
	my $taxo_profiling_map         = "$data_dir/$databases.refmg.map";
	my $taxo_profiling_map_tot_len = "$data_dir/$databases.refmg.map.total_length";
	my $taxo_profiling_motu_map    = "$data_dir/$databases.motu.map";
	if ( $multiple_mappers_summarize_level ne 'none' ) {
		unless ( -e "$taxo_profiling_map" ) {
			die "\nERROR & EXIT: Missing map file $taxo_profiling_map";
		}
	}
	else {
		$taxo_profiling_map = "not_used";
	}
	
		
	
	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?";

	print localtime() . ": Creating $job jobs...";
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );

	### JOB SPECIFIC ###
	my $db = $databases[0];
	foreach my $sample (@samples) {
		
		print JOB "mkdir -p $temp_dir/$sample/temp; ";
		
		my $LOG       = " 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.log >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.log ";
		my $folder    = "$cwd/$sample/reads.filtered.$databases.$conf{MOCAT_data_type}";
		my $output    = "$sample.filtered.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
		my $input;
		my $format;
		if ( $conf{MOCAT_mapping_mode} eq "allbest" ) {
			$format = 'BAM';
			$input  = "$folder/$output.bam";
		}
		elsif ( $conf{MOCAT_mapping_mode} eq "unique" || $conf{filter_mode} eq "random" ) {
			$format = 'SOAP';
			$input  = "$folder/$output.soap.gz";
		}
		unless ( -e $input ) {
			die "\nERROR & EXIT: Missing input file $input";
		}
		print JOB "$conf{MOCAT_DEV_DIR}/MOCATCalcMultipleMappers.pl -in $input -format $format -sample $sample -sum $multiple_mappers_summarize_level -map $taxo_profiling_map $LOG\n";
	}
	close JOB;
	print " OK!\n";
}

sub post_check_files {

	my $job = $_[0];

	print localtime() . ": Checking files... ";
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}
	my $fail = 0;
	my ( @read, @insert, @readstats, @insertstats );
	my $databases = join( "_AND_", @databases );
	foreach my $sample (@samples) {
		my $folder = "$cwd/$sample/reads.filtered.$databases.$conf{MOCAT_data_type}";
		my $output = "$sample.filtered.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
		my $input;
		my $format;
		if ( $conf{MOCAT_mapping_mode} eq "allbest" ) {
			$format = 'BAM';
			$input  = "$folder/$output.bam";
		}
		elsif ( $conf{MOCAT_mapping_mode} eq "unique" || $conf{filter_mode} eq "random" ) {
			$format = 'SOAP';
			$input  = "$folder/$output.soap.gz";
		}
		unless ( -e "$input.multiplemappers.read.allbyall.$multiple_mappers_summarize_level" ) {
			$fail = 1;
			print "MISSING $input.multiplemappers.read.allbyall.$multiple_mappers_summarize_level\n";
		}
		unless ( -e "$input.multiplemappers.insert.allbyall.$multiple_mappers_summarize_level" ) {
			$fail = 1;
			print "MISSING $input.multiplemappers.read.allbyall.$multiple_mappers_summarize_level\n";
		}
		unless ( -e "$input.multiplemappers.insert.stats" ) {
			$fail = 1;
			print "MISSING $input.multiplemappers.insert.stats\n";
		}
		unless ( -e "$input.multiplemappers.read.stats" ) {
			$fail = 1;
			print "MISSING $input.multiplemappers.read.stats\n";
		}
		push @read,        "$input.multiplemappers.read.allbyall.$multiple_mappers_summarize_level";
		push @insert,      "$input.multiplemappers.insert.allbyall.$multiple_mappers_summarize_level";
		push @readstats,   "$input.multiplemappers.read.stats";
		push @insertstats, "$input.multiplemappers.insert.stats";
	}
	if ($fail) {
		die "\nERROR & EXIT: One or more file(s) are missing";
	}
	print " OK!\n";

	print localtime() . ": Generating input files... ";
	my $folder = "$cwd/multiple.mappers/$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
	system "rm -fr $folder/$sample_file_basename.$multiple_mappers_summarize_level $folder/$sample_file_basename.$multiple_mappers_summarize_level.tar.gz && mkdir -p $folder/$sample_file_basename.$multiple_mappers_summarize_level";
	open B,  ">$folder/$sample_file_basename.$multiple_mappers_summarize_level/base"       or die "ERROR & EXIT: Cannot write $folder/$sample_file_basename/base";
	open O,  ">$folder/$sample_file_basename.$multiple_mappers_summarize_level/$sample_file_basename.$multiple_mappers_summarize_level.output.files"       or die "ERROR & EXIT: Cannot write $folder/$sample_file_basename.$multiple_mappers_summarize_level/$sample_file_basename.$multiple_mappers_summarize_level.output.files";
	open I,  ">$folder/$sample_file_basename.$multiple_mappers_summarize_level/$sample_file_basename.$multiple_mappers_summarize_level.insert.files"       or die "ERROR & EXIT: Cannot write $folder/$sample_file_basename.$multiple_mappers_summarize_level/$sample_file_basename.$multiple_mappers_summarize_level.insert.files";
	open R,  ">$folder/$sample_file_basename.$multiple_mappers_summarize_level/$sample_file_basename.$multiple_mappers_summarize_level.read.files"         or die "ERROR & EXIT: Cannot write $folder/$sample_file_basename.$multiple_mappers_summarize_level/$sample_file_basename.$multiple_mappers_summarize_level.read.files";
	open IS, ">$folder/$sample_file_basename.$multiple_mappers_summarize_level/$sample_file_basename.$multiple_mappers_summarize_level.insert.stats.files" or die "ERROR & EXIT: Cannot write $folder/$sample_file_basename.$multiple_mappers_summarize_level/$sample_file_basename.$multiple_mappers_summarize_level.insert.stats.files";
	open RS, ">$folder/$sample_file_basename.$multiple_mappers_summarize_level/$sample_file_basename.$multiple_mappers_summarize_level.read.stats.files"   or die "ERROR & EXIT: Cannot write $folder/$sample_file_basename.$multiple_mappers_summarize_level/$sample_file_basename.$multiple_mappers_summarize_level.read.stats.files";
	open S,  ">$folder/$sample_file_basename.$multiple_mappers_summarize_level/$sample_file_basename.$multiple_mappers_summarize_level.samples"            or die "ERROR & EXIT: Cannot write $folder/$sample_file_basename.$multiple_mappers_summarize_level/$sample_file_basename.$multiple_mappers_summarize_level.samples";
	foreach my $n (@read) {
		print R "$n\n";
	}
	foreach my $n (@insert) {
		print I "$n\n";
	}
	foreach my $n (@readstats) {
		print RS "$n\n";
	}
	foreach my $n (@insertstats) {
		print IS "$n\n";
	}
	foreach my $n (@samples) {
		print S "$n\n";
	}
	print O "$folder/$sample_file_basename.$multiple_mappers_summarize_level/$sample_file_basename.$multiple_mappers_summarize_level.results\n";
	print B "$sample_file_basename.$multiple_mappers_summarize_level\n";
	close I;
	close R;
	close IS;
	close RS;
	close S;
	close O;
	close B;
	print " OK!\n";

	print localtime() . ": Running R session...\n";
	my $LOG = " >> $cwd/logs/$job/samples/MOCATJob_$job.log 2>> $cwd/logs/$job/samples/MOCATJob_$job.log ";
	system "cd $folder/$sample_file_basename.$multiple_mappers_summarize_level && R --no-save --args '$sample_file_basename.$multiple_mappers_summarize_level' < '$conf{MOCAT_DEV_DIR}/MOCATCalcMultipleMappers.R' $LOG";
	print localtime() . ": R session completed.\n";
	unless ( -e "$folder/$sample_file_basename.$multiple_mappers_summarize_level/$sample_file_basename.$multiple_mappers_summarize_level.results.RData" ) {
		die "ERROR & EXIT: Seems like final file $folder/$sample_file_basename.$multiple_mappers_summarize_level/$sample_file_basename.$multiple_mappers_summarize_level.results.RData was not created correctly during the R session. Check out the log file $cwd/logs/$job/samples/MOCATJob_$job.log for details";
	}

	open OUT, ">>$cwd/MOCAT.results";
	print OUT "multiple mapper\t$job\_$date\t$sample_file_basename\t$sample_file\t$username\t$hostname\texecute\tcd $folder && R --interactive --no-save < $folder/$sample_file_basename.$multiple_mappers_summarize_level.execute.R\n";
	close OUT;

	# Generate shiny folder and tar.gz file
	print localtime() . ": Generating shiny files...";
	system "ln -s $conf{MOCAT_DEV_DIR}/shiny/multipleMapper/ui.R $folder/$sample_file_basename.$multiple_mappers_summarize_level 2>>/dev/null";
	system "ln -s $conf{MOCAT_DEV_DIR}/shiny/multipleMapper/server.R $folder/$sample_file_basename.$multiple_mappers_summarize_level 2>>/dev/null";
	system "cd $folder/$sample_file_basename.$multiple_mappers_summarize_level/ && tar -czvf $folder/$sample_file_basename.$multiple_mappers_summarize_level.tar.gz * 2>>/dev/null >>/dev/null";
	system "cp $conf{MOCAT_DEV_DIR}/shiny/multipleMapper/execute.R $folder/$sample_file_basename.$multiple_mappers_summarize_level.execute.R";
	open OUT, ">>$folder/$sample_file_basename.$multiple_mappers_summarize_level.execute.R";
	print OUT "setwd('$folder')\n";
	print OUT "runApp('$sample_file_basename.$multiple_mappers_summarize_level')\n";
	close OUT;
	print " OK!\n";

	print localtime() . ": RESULTS ARE SAVED IN $folder/$sample_file_basename/$sample_file_basename.$multiple_mappers_summarize_level.results.RData\n";
	print localtime() . ": TO VIEW RESULTS INTERACTIVELY, EXECUTE: 'R --interactive --no-save < $folder/$sample_file_basename.$multiple_mappers_summarize_level.execute.R'\n";
}

sub post_check_files2 {

	my $job       = $_[0];
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}
	my $fail = 0;
	my ( @read, @insert, @readstats, @insertstats );
	my $databases = join( "and", @databases );
	my $folder = "$cwd/multiple.mappers/$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
}
1;
