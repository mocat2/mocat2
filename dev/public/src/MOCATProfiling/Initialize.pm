package MOCATProfiling::Initialize;
use strict;
use warnings;
use MOCATProfiling::Variables;
use IO::Handle;
use File::Basename;
use Cwd;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

sub initialize {

	#############################################################################################################################################
	# PRINT INPUT LINES
	#############################################################################################################################################
	print STDERR "\nINPUT :";
	print STDERR " " . join( " ", @argv ) . "\n\n";
	#############################################################################################################################################

	#############################################################################################################################################
	# CHECK INPUT FILES
	#############################################################################################################################################
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Version is $version\n";
	if ($VERBOSE) {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Verbose mode is ON\n";
	} else {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Verbose mode is OFF\n";
	}
	$cwd = getcwd;
	chomp( my $host = `hostname` );
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Host = $host\n";
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Checking files and parameters\n";
	foreach my $input_file (@input_files) {
		( -e "$input_file" ) or die "ERROR & EXIT: Missing input file $input_file";
		if ($sam) {
			$input_file_format = "SAM";
			print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Format forced to SAM with option -sam\n";
		}
		elsif ( $input_file =~ m/\.bam$/ ) {
			$input_file_format = "BAM";
			print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Format set to BAM\n";
		}
		elsif ( $input_file =~ m/\.sam(.gz)$/ ) {
			$input_file_format = "SAM";
			print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Format set to SAM\n";
		}
		elsif ( $input_file =~ m/\.sam$/ ) {
			$input_file_format = "SAM";
			print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Format set to SAM\n";
		}
		elsif ( $input_file =~ m/\.soap.gz$/ ) {
			$input_file_format = "SOAP";
			print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Format set to SOAP\n";
		}
		else {
			die "ERROR & EXIT: Unknown file ending, cannot identify file format. File ending should be one of: .bam .sam .sam.gz";
		}
	}
	unless ($output_file) {
		die "ERROR & EXIT: Missing output file, please specify -output_file";
	}
	foreach my $coord_file (@coord_files) {
		( -e "$coord_file" || -e "$coord_file.gz" ) or die "ERROR & EXIT: Missing coord file $coord_file (or $coord_file.gz)";
		if ( -e "$coord_file.gz" ) {
			print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Found gz coord file ($coord_file.gz)\n";
		}
	}

	# used for checking if we have previously created HASH for these coord files
	my $name = MOCATProfiling::Misc::checkAndReturnDB( \@coord_files );
	$name = "$name.MOCAT_DATA";
	if ( -e "$name.HASH" && -e "$name.avg" ) {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Found MOCAT DATA file ($name.HASH)\n";
	}

	if ( $mode ne 'gene' ) {
		( -e "$map_file" || -e "$map_file.gz" ) or die "ERROR & EXIT: Missing $mode map file $map_file (or $map_file.gz)";
		if ( -e "$map_file.gz" ) {
			print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Found gz map file ($map_file.gz)\n";
		}
	}
	if ($rownames_file) {
		open my $OUT, ">$rownames_file" or die "ERROR & EXIT: Cannot write to $rownames_file";
		close $OUT or die "ERROR & EXIT: Closing pipe failed, this means something is wrong. Corrput files?";
	}

	# Parameter checks
	$PE_filter           or die "ERROR & EXIT: -PE_filter not specified";
	$sample_name         or die "ERROR & EXIT: -sample_name not specified";
	$zip_execute         or die "ERROR & EXIT: -zip_execute not specified";
	$THREADS             or die "ERROR & EXIT: -threads not specified";
	$BLOCKSIZE           or die "ERROR & EXIT: -blocksize not specified";
	$samtools_executable or die "ERROR & EXIT: -samtools_executable not specified";

	# Check and try to write stats files
	unless ($out_stats_file) {
		die "ERROR & EXIT: '' is an invalid name for the stats file, please specify -out_stats_file";
	}
	open my $out, ">", $out_stats_file or die "ERROR & EXIT: Cannot write to $out_stats_file: $!";
	close $out;
	unless ($out_PE_stats_file) {
		die "ERROR & EXIT: '' is an invalid name for the stats file, please specify -out_PE_stats_file";
	}
	open $out, ">", $out_PE_stats_file or die "ERROR & EXIT: Cannot write to $out_PE_stats_file: $!";
	close $out;
	#############################################################################################################################################

	#############################################################################################################################################
	# CHECK FOLDERS
	#############################################################################################################################################
	# Check folders
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Checking folders\n";

		unless ( $output_file =~ m/^\// ) {
			print STDERR "===== WARNING ===== : -output_file=$output_file set to $cwd/$output_file\n";
			$output_file = "$cwd/$output_file";
		}
		unless ( $rownames_file =~ m/^\// ) {
			print STDERR "===== WARNING ===== : -rownames_file=$rownames_file set to $cwd/$rownames_file\n";
			$rownames_file = "$cwd/$rownames_file";
		}
	$output_file_base = fileparse($output_file);
	$rownames_file_base = fileparse($rownames_file);
	if ($temp_folder) {
		MOCATProfiling::Misc::mkdir_or_die($temp_folder);
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Created & using main temporary folder $temp_folder\n";
	}
	MOCATProfiling::Misc::mkdir_or_die("$temp_file.output");
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Created & using specific temporary folder $temp_file.output\n";
	
	#############################################################################################################################################

	#############################################################################################################################################
	# DEFINE levels, INITALIZE TEMP FILES
	#############################################################################################################################################
	# note that for functional, the multiple mapper level is on the genes, and not on the cog, ko, mod and pathway level
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Checking mode = $mode. Valid options are NCBI, mOTU, functional & gene\n";
	if ( $mode eq 'NCBI' ) {
		@levels = ( 'gene', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'specI_cluster', 'taxaid' );
	}
	elsif ( $mode eq 'mOTU' ) {
		@levels = ( 'gene', 'mOTU' );
	}
	elsif ( $mode eq 'functional' ) {
		@levels = ('gene');    # this will be changed later to cog,ko,mod,path
	}
	elsif ( $mode eq 'gene' ) {
		@levels = ('gene');
	}
	else {
		die "ERROR & EXIT: Incorrect mode $mode specified";
	}

	# set values and initialize temp files
	if ($temp_file) {
		foreach my $level (@levels) {
			unless ($zip_execute) {
				die "ERROR & EXIT: if -temp_file is specified, -zip must also be specified";
			}
			for my $i ( 1 .. $THREADS ) {
				open my $OUT, ">", "$temp_file.MM.$level.$i.data.gz" or die "ERROR & EXIT: Cannot write to temp file $temp_file.MM.$level.$i.data.gz";
				close $OUT or die "ERROR & EXIT: Closing pipe failed, this means something is wrong. Corrput files?";
				open my $FH, "| $zip_execute -c > $temp_file.MM.$level.$i.data.gz" or die "ERROR & EXIT: Cannot write to OUTPUT PIPE | $zip_execute -c > $temp_file.MM.$level.$i.data";
				$open_temp_files{"$level.$i"} = $FH;
			}
		}
		for my $i ( 1 .. $THREADS ) {
			open my $OUT, ">", "$temp_file.HASH.$i.data.gz" or die "ERROR & EXIT: Cannot write to temp file $temp_file.HASH.$i.data.gz";
			close $OUT or die "ERROR & EXIT: Closing pipe failed, this means something is wrong. Corrput files?";
			open my $FH, "| $zip_execute -c > $temp_file.HASH.$i.data.gz" or die "ERROR & EXIT: Cannot write to OUTPUT PIPE | $zip_execute -c > $temp_file.HASH.$i.data";
			$open_temp_files{"HASH.$i"} = $FH;
		}
	}

	# Set paired end filtering variables
	if ( $PE_filter eq 'yes' ) {
		foreach my $i (@levels) {
			$PEaffectedInserts{$i} = 0;
		}
	}
	else {
		foreach my $i (@levels) {
			$PEaffectedInserts{$i} = 'PE_filter_set_to_off';
		}
	}
	#############################################################################################################################################
}

1;
