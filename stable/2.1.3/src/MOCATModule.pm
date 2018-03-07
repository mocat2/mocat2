package MOCATModule;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

sub create_job {
	### DEFINE VARIABLES AND OPEN input FILE ###
	my $job        = $_[0];
	my $single_job = $_[1];
	my $processors = $_[2];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?";

	print localtime() . ": Creating $job jobs...";
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}
	my $mc_e = 'no';
	if ($use_extracted_reads) {
		my $mc_e = 'yes';
	}

	my $databases = MOCATCore::checkAndReturnDB( \@MODULE_DB );
	
	unless ($databases) {
		$databases = "";
	}
	
	if ( $databases eq 's' || $databases eq 'c' || $databases eq 'f' || $databases eq 'r' ) {
		die "ERROR & EXIT: Sorry, currently these options are not supported for modules";
	}
	
	my $mc_m;
	my $mc_tax;
	
	unless ($mod_requests{type}) {
		$mod_requests{type} = "";
	}
	
	if ($mod_requests{type} eq 'NCBI') {
		$mc_m='taxonomic';
		$mc_tax = "NCBI";
	}
elsif ($mod_requests{type} eq 'mOTU') {
		$mc_m='taxonomic';
		$mc_tax = "mOTU";
	}
elsif ($mod_requests{type} eq 'functional') {
		$mc_m='functional';
		$mc_tax = "NA";
	} elsif ($mod_overrides{no_type}) {
		$mc_m='NA';
		$mc_tax='NA';
	}
	else {
		die "\nERROR & EXIT: -type '$mod_requests{type}' was specified, but expecte one of functional, NCBI or mOTU"
	}

	### JOB SPECIFIC ###
	if ( $single_job eq 'TRUE' ) {

		print JOB "export MC_E=$mc_e; ";
		print JOB "export MC_SF=$sample_file; ";
		print JOB "export MC_OUT=$cwd/$mod_overrides{outfolder}/$job; ";
		print JOB "export MC_TMP=$cwd/$mod_overrides{outfolder}/$job; ";
		print JOB "export MC_WD=$cwd; ";
		print JOB "export MC_EXT=$ext_dir; ";
		print JOB "export MC_R=$reads; ";
		print JOB "export MC_DB=$databases; ";
		print JOB "export MC_ID=95; ";
		print JOB "export MC_M=$mc_m; ";
		print JOB "export MC_TAX=$mc_tax; ";
		print JOB "export MC_MODE=$conf{MOCAT_data_type}; ";
		print JOB "export MC_MAPPING=$conf{MOCAT_mapping_mode}; ";
		print JOB "export MC_LEN=$conf{filter_length_cutoff}; ";
		print JOB "export MC_ID=$conf{filter_percent_cutoff}; ";
		print JOB "export MC_BIN=$bin_dir; ";
		print JOB "export MC_SRC=$src_dir; ";
		print JOB "export MC_MOD=$mod_dir/$job/; ";
		print JOB "export MC_MODNAME=$job; ";
		print JOB "export MC_NAME='$modcomment'; ";
		print JOB "export MC_MODCFG=$MODCFG; ";
		print JOB "export MC_DATE=$date; ";
		print JOB "export MC_DATA=$data_dir; ";
		print JOB "export MC_CPU=$cpu_module; ";
		
		foreach my $k (sort keys %mod_optionals) {
			print JOB "export $k='$mod_optionals{$k}'; ";
		}
		print JOB "export MC_FUNCT_MAP=$data_dir/$databases.functional.map; ";
		
		unless($mod_requests{metadata}) {
			$mod_requests{metadata}= "NA";
		}
		unless($mod_requests{groups}) {
			$mod_requests{groups} = "NA"
		}
		unless ($mod_requests{grouping}) {
			$mod_requests{grouping} = "NA";
		}
		unless ($mod_requests{assembly_type}) {
			$mod_requests{assembly_type} = "NA"
		}
                unless ($mod_requests{blasttype}) {
		    $mod_requests{blasttype} = "NA"
                }
		print JOB "export MC_ASS=$mod_requests{assembly_type}; ";
		print JOB "export MC_OPT=\"BLASTTYPE=$mod_requests{blasttype} METADATA=$mod_requests{metadata} GROUPS=$mod_requests{groups} GROUPING=$mod_requests{grouping}\"; ";
		#print JOB "mkdir -p $cwd/$mod_overrides{outfolder}/temp/$job; ";
		print JOB "bash $mod_dir/$job/$job.sh >$cwd/logs/$job/samples/MOCATJob_${job}.$date.log 2>$cwd/logs/$job/samples/MOCATJob_${job}.$date.log\n";
		@samples = ($date);
	}
	else {
		die "SINGLE_JOB=FALSE currently not suppported. Sorry.";
		foreach my $sample (@samples) {

			print JOB "mkdir -p $temp_dir/$sample/temp; ";

		}
	}

	print " OK!\n";
	
	$MODLOG = "$cwd/logs/$job/samples/MOCATJob_${job}.$date.log";
}

1;
