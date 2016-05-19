package MOCATAnalyze;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;
use Term::ANSIColor;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

sub run {

	print localtime() . ": MOCATAnalyze : starting\n";

	### DEFINE VARIABLES AND OPEN input FILE ###

	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	my @databases;
	my $downsample;
	my $desc = "";
	my $execute;
	my @line;
	my $GROUP = "";
	my %analyze_settings;
	our $load_metadata = "no";
	my $iTOL_used  = 0;
	my $iTOL_label = "Label";

	# metadata
	unless ( defined( $analyze_metadata[0] ) ) {
		if ( $conf{analyze_metadata} ) {
			@analyze_metadata = split( ",", $conf{analyze_metadata} );
		}
	}

	# feature annotations
	unless ( defined( $analyze_feature_annotations[0] ) ) {
		if ( $conf{analyze_feature_annotations} ) {
			@analyze_feature_annotations = split( ",", $conf{analyze_feature_annotations} );
		}
	}

	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}
	my $databases = join( "_AND_", @do_analyze );
	if ($analyze_file) {
		unless ( -e $analyze_file ) {
			die "ERROR & EXIT: Missing file $analyze_file";
		}
		chomp( $databases = `basename $analyze_file` );
		@samples = ();
		
		chomp( our $systemType = `uname -s` );
		chomp( @samples = `head -1 $analyze_file | sed 's/\\t/\\n/g' | awk 'NR>1'` );
	}
	chomp( my $analyze_date = `date +\%Y\%b\%d_\%H\%M\%S` );

	# Prepare data
	print localtime() . ": MOCATAnalyze : creating settings\n";
	my @i    = split ",", $conf{analyze_base_and_insert};
	my @end  = split ",", $conf{analyze_normalization};
	my @name = split ",", $conf{analyze_taxonomic_level};

	# Choose function
	my @functions = `ls -1 $conf{MOCAT_DEV_DIR}/analyze/*.R`;
	my @basenames;
	foreach my $f (@functions) {
		chomp $f;
		chomp( my $bn   = `basename $f` );
		chomp( my $desc = `grep '^# DESC' $f | sed 's/.*DESC//' | sed 's/^ *//'` );
		$bn =~ s/.R$//;
		chomp($bn);
		push @basenames, $bn;
	}
	unless ($analyze_function) {
		print "$sep\nAVAILABLE FUNCTIONS TO EXECUTE ARE:\n";
		foreach my $bn (@basenames) {
			chomp( my $desc = `grep '^# DESC' $conf{MOCAT_DEV_DIR}/analyze/$bn.R | sed 's/.*DESC//' | sed 's/^ *//'` );
			print " - ";
			print color 'bold';
			print "$bn";
			print color 'reset';
			print "\t($desc)\n";
		}
		print "\nADDITIONAL OPTIONS: ";
		print color 'bold';
		print "<FUNCTION>,<TYPE>,<DESCRIPTION>,<DOWNSAMPLING>\n";
		print color 'reset';
		print "   - TYPE: eg. 135.\n";
		print "           D1: 1=base, 2=insert\n";
		print "           D2: 1=raw, 2=norm, 3=scaled, 4=alignment fraction norm\n";
		print "           D3: 1=kingd, 2=phyl, 3=class, 4=order, 5=fam, 6=genus, 7=sp, 8=cur.sp, 9=gene, a=COG, b=KO, c=module, d=pathway, e=mOTU\n";
		print "           (limitations: gene, COG, KO, module, pathway cannot use 'scaled')\n";
		print "   - DESCRIPTION: Short description of this particular analysis\n";
		print "   - DOWNSAMPLING: rrarefy.min OR rrarefy.top90percent OR rrarefy OR forslund (additionally specify size to downsample to later)\n";
		print "   - Example: use function 'sum' and downsample to min: 'sum,,,rrarefy.min'\n\n";
		print color 'bold';
		print "Run function: ";
		print color 'reset';
		chomp( $execute = <> );

		if ( $execute eq '' ) {
			die "ERROR & EXIT: No funciton specified.";
		}
		@line = split ",", $execute;
		$line[0] =~ s/^\s+//;
		$line[0] =~ s/\s+$//;
		$execute    = $line[0];
		$downsample = "none";
	}
	else {
		$execute = $analyze_function;
		if ($analyze_downsample) {
			$downsample = $analyze_downsample;
		}
		else {
			$downsample = "none";
		}
	}
	my $exists = 0;
	foreach my $fun (@basenames) {
		if ( $fun eq $execute ) {
			$exists = 1;
		}
	}
	if ( $exists == 0 ) {
		die "ERROR & EXIT: Function $execute is not defined.";
	}

	if ($analyze_desc) {
		$analyze_desc =~ s/\s+/_/g;
		$desc = ".$analyze_desc";
	}
	else {
		$desc = "";
	}

	if ( $analyze_type && $analyze_norm && $analyze_taxa ) {
		@i    = ();
		@end  = ();
		@name = ();
		push @i,    $analyze_type;
		push @end,  $analyze_norm;
		push @name, $analyze_taxa;
	}
	if ($analyze_file) {
		if ( scalar @i > 1 || scalar @end > 1 || scalar @name > 1 ) {
			@i    = ();
			@end  = ();
			@name = ();
			push @i,    'unknown';
			push @end,  'unknown';
			push @name, 'unknown';
		}
	}
	else {
		if ( $line[1] ) {
			$line[1] =~ s/^\s+//;
			$line[1] =~ s/\s+$//;
			unless ( $line[1] =~ m/^\d\d\d|\d\d(a|b|c|d|e)$/ ) {
				die "ERROR & EXIT: Type is not valid.";
			}
			unless ( length( $line[1] ) == 3 ) {
				die "ERROR & EXIT: Type length is exactly 3.";
			}
			my $counter = 0;
			my @options = split /|/, $line[1];
			@i    = ();
			@end  = ();
			@name = ();
			while ( scalar @options > 0 ) {
				$counter++;
				my $value = shift @options;
				if ( $counter == 1 ) {
					if ( $value eq '1' ) {
						push @i, "base";
					}
					if ( $value eq '2' ) {
						push @i, "insert";
					}
				}
				if ( $counter == 2 ) {
					if ( $value eq '1' ) {
						push @end, "raw";
					}
					if ( $value eq '2' ) {
						push @end, "norm";
					}
					if ( $value eq '3' ) {
						push @end, "scaled";
					}
					if ( $value eq '4' ) {
						push @end, "algfrac.norm";
					}
				}
				if ( $counter == 3 ) {
					if ( $value eq '9' ) {
						push @name, "gene";
					}
					if ( $value eq '1' ) {
						push @name, "kingdom";
					}
					if ( $value eq '2' ) {
						push @name, "phylum";
					}
					if ( $value eq '3' ) {
						push @name, "class";
					}
					if ( $value eq '4' ) {
						push @name, "order";
					}
					if ( $value eq '5' ) {
						push @name, "family";
					}
					if ( $value eq '6' ) {
						push @name, "genus";
					}
					if ( $value eq '7' ) {
						push @name, "species";
					}
					if ( $value eq '8' ) {
						push @name, "specI_clusters";
					}
					if ( $value eq 'a' ) {
						push @name, "cog";
					}
					if ( $value eq 'b' ) {
						push @name, "ko";
					}
					if ( $value eq 'c' ) {
						push @name, "module";
					}
					if ( $value eq 'd' ) {
						push @name, "pathway";
					}
					if ( $value eq 'e' ) {
						push @name, "mOTU";
					}

				}
				if ( $counter == 4 ) {
					$counter = 0;
				}
			}
		}
		if ( $line[2] ) {
			$line[2] =~ s/^\s+//;
			$line[2] =~ s/\s+$//;
			$desc = $line[2];
			$desc =~ s/\s+/_/g;
			$desc = ".$desc";
		}
		if ( $line[3] ) {
			$line[3] =~ s/^\s+//;
			$line[3] =~ s/\s+$//;
			unless ( $line[3] eq 'rrarefy.min' || $line[3] eq 'rrarefy.top90percent' || $line[3] eq 'forslund' || $line[3] eq 'rrarefy' ) {
				die "ERROR & EXIT: Not a valid downsample method.";
			}
			$downsample = $line[3];
		}
	}

	print "$sep\n";
	print localtime() . ": MOCATAnalyze : creating folders\n";

	#chomp(my $usr = `whoami`);
	my $usr                = "";
	my $main_output_folder = "$cwd/analysis/$usr/$databases/$sample_file_basename/$analyze_date.$execute$desc.".join('.', @name).".".join('.', @end).".".join('.', @i);
	my $calculations       = "$main_output_folder/settings";
	my $tables             = "$main_output_folder/tables";
	my $figures            = "$main_output_folder/figures";
	my $data               = "$main_output_folder/data";
	my $excel              = "$main_output_folder/excel";
	system "mkdir -p $data $calculations $tables $figures $excel";

	open SETTINGS, ">$calculations/R.settings" or die "ERROR & EXIT: Cannot write to $calculations/R.settings";
	print SETTINGS "type\tsetting\n";

	# Add file names to settings file
	print localtime() . ": MOCATAnalyze : printing settings\n";
	open FILES, ">$calculations/R.files" or die "ERROR & EXIT: Cannot write to $calculations/R.files";
	print FILES "index\tfile\tinsert_or_base\tnorm_type\ttax_level\tsettings\tall\tzipped\n";
	my $counter = 0;
	for my $i (@i) {
		for my $end (@end) {
			for my $name (@name) {
				$counter++;
				my $file;
				if ($analyze_file) {
					$file = $analyze_file;
				}
				else {
					if ( $name eq 'gene' ) {
						$file = "$cwd/abundance.tables/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}/$sample_file.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$i.coverage.$end";
					}
					elsif ( $name eq 'cog' || $name eq 'ko' || $name eq 'module' || $name eq 'pathway' ) {
						my $output_folder = "$cwd/functional.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}/";
						$file = "$output_folder/$sample_file_basename.functional.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$i.$end.$name";
						unless ( -e "$data_dir/$name" . "2name.txt" ) {
							close SETTINGS;
							system "rm -fr $main_output_folder";
							die "ERROR & EXIT: If you want to use $name in your analysis, the following file should exist: $data_dir/$name" . "2name.txt";
						}
						print SETTINGS "$name" . "2name\t$data_dir/$name" . "2name.txt\n";
					}
					elsif ($name eq 'mOTU') {
						$file = "$cwd/taxonomic.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}/$sample_file.mOTU.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$i.$end.$name";
					}
					else {
						$file = "$cwd/taxonomic.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}/$sample_file.taxonomic.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$i.$end.$name";
					}
				}
				my $zipped = 'no';
				if ( -e "$file.gz" ) {
					$file = "$file.gz";
					$zipped = 'yes';
				}
				unless ( -e $file ) {
					close FILES;
					close SETTINGS;
					system "rm -f $calculations/R.settings; rm -fr $main_output_folder";
					die "ERROR & EXIT: Missing file $file\nDid you run coverage calculations and then taxonomic profiling pipeline using this specific $sample_file sample file?\nAnd if you specified gene, you also need to summarize the coverage calculations using -pcf";
				}
				print FILES "$counter\t$file\t$i\t$end\t$name\t$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}\t$i.$end.$name.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}\t$zipped\n";
			}
		}
	}
	close FILES;

	# Check if metadata is required
	chomp( my @require  = `grep -h '^# REQUIRE' $conf{MOCAT_DEV_DIR}/MOCATAnalyze_R.R $conf{MOCAT_DEV_DIR}/analyze/$execute.R | sed 's/.*REQUIRE//' | sed 's/^ *//'` );
	chomp( my @settings = `grep -h '^# SETTINGS' $conf{MOCAT_DEV_DIR}/MOCATAnalyze_R.R $conf{MOCAT_DEV_DIR}/analyze/$execute.R | sed 's/^# *//'` );
	chomp( my @SUPPRESS = `grep -h '^# SUPPRESS' $conf{MOCAT_DEV_DIR}/MOCATAnalyze_R.R $conf{MOCAT_DEV_DIR}/analyze/$execute.R | sed 's/.*SUPPRESS//' | sed 's/^ *//'` );
	my %SUPPRESS;
	foreach my $require (@SUPPRESS) {
		$require =~ m/\s*(.*)\s*/;
		$SUPPRESS{$1} = 1;
	}

	# Write some settings, keep file open for more writing later
	print SETTINGS "working_dir\t$calculations\n";
	print SETTINGS "tables_dir\t$tables\n";
	print SETTINGS "data_dir\t$data\n";
	print SETTINGS "figures_dir\t$figures\n";
	print SETTINGS "downsample\t$downsample\n";
	print SETTINGS "load_functions\t$conf{MOCAT_DEV_DIR}/MOCATAnalyze_functions.R\n";

	# If downsample is an integer, select method 
	if ( $analyze_downsample_size && ( $downsample eq 'rarefy' || $downsample eq 'forslund' ) ) {
		print SETTINGS "downsample_size\t$analyze_downsample_size\n";
	}
	elsif ( $downsample eq 'rarefy' || $downsample eq 'forslund' ) {
		print "$sep\nPlease specify a number to downsample to (you specified method $downsample).\nNOTE, method 'forslund' requires the R script to load a METADATA file containing the column 'average_read_length':\n";
		chomp( my $downsample_size = <> );
		if ( $downsample_size =~ m/^\d+$/ ) {
			print SETTINGS "downsample_size\t$downsample_size\n";
		}
		else {
			close SETTINGS;
			system "rm -fr $main_output_folder";
			die "ERROR & EXIT: Incorrect downsampling size.";
		}
	}

	# Set some variables that may change, these will be written to the config file
	my $load_groups              = "no";
	my $load_feature_annotations = 'no';

	if ( $analyze_settings[0] ) {
		for ( my $i = 0 ; $i <= scalar @analyze_settings - 1 ; $i += 2 ) {
			unless ( defined( $analyze_settings[ ( $i + 1 ) ] ) ) {
				close SETTINGS;
				system "rm -fr $main_output_folder";
				die "ERROR & EXIT: -analyze_settings does not have an even number of entries";
			}
			$analyze_settings{ $analyze_settings[$i] } = $analyze_settings[ $i + 1 ];
		}
	}
	foreach my $setting (@settings) {
		my @info = split /\[/, $setting;
		$info[0] =~ m/SETTINGS\s*:\s*(\w*)/;
		my $sett = $1;
		if ( !$SUPPRESS{$sett} ) {
			if ( defined( $analyze_settings{$sett} ) ) {
				print localtime() . ": MOCATAnalyze : SPECIFIED $sett = $analyze_settings{$sett}\n";
				print SETTINGS "$sett\t$analyze_settings{$sett}\n";
			}
			elsif ($analyze_running_galaxy) {
				close SETTINGS;
				system "rm -fr $main_output_folder";
				die "ERROR & EXIT: Missing settings value $sett";
			}
			else {
				print localtime() . ": MOCATAnalyze : REQUIRE $info[0]\n";
				$setting =~ m/(.*)\s*\[(.*)\]\s*\{(.*)\}\s*\{(.*)\}\s*\((.*)\)/;
				print "$sep\n$5 ($2) (default: '$4'):\n";
				chomp( my $input = <> );
				if ( $input eq '' ) {
					print SETTINGS "$sett\t$4\n";
				}
				else {
					if ( $3 eq 'CHECK' ) {
						my @opts = split ",", $2;
						my $checkOK = 0;
						foreach my $opt (@opts) {
							$opt =~ s/^\s*//;
							$opt =~ s/\s*$//;
							if ( $opt eq $input ) {
								$checkOK = 1;
							}
						}
						unless ($checkOK) {
							close SETTINGS;
							system "rm -fr $main_output_folder";
							die "ERROR & EXIT: $input is not a valid option for $sett";
						}
					}
					print SETTINGS "$sett\t$input\n";
				}
				print "$sep\n";
			}
		}
	}

	foreach my $req (@require) {
		my @split_req = split /\|/, $req;
		print localtime() . ": MOCATAnalyze : checking required $split_req[0]\n";

		# cond A or B
		if ( $split_req[0] eq 'condA' || $split_req[0] eq 'condB' || $split_req[0] eq 'featureA' || $split_req[0] eq 'featureB' ) {
			$load_metadata = load_metadata( $load_metadata, $main_output_folder );
			if ( $split_req[0] eq 'condA' || $split_req[0] eq 'featureA' ) {
				if ($analyze_condA) {
				}
				elsif ($analyze_running_galaxy) {
					close SETTINGS;
					system "rm -fr $main_output_folder";
					die "ERROR & EXIT: Missing $split_req[0] without default";
				}
				else {
					if ( $GROUP ne "" ) {
						chomp( my @fields = `cut -f \` head -1 $load_metadata | sed 's/\\t/\\n/g' | grep -n '$GROUP' | cut -f 1 -d':' \` $load_metadata | sort -u | grep -v $GROUP` );
						print "$sep\nAvailable conditions: " . join( ", ", @fields ) . "\n";
					}
					else {
						print "$sep\n";
					}
					print "Please enter condidtion/feature A:\n";
					chomp( $analyze_condA = <> );
					print "$sep\n";
				}
				print SETTINGS "condA\t$analyze_condA\n";
			}
			if ( $split_req[0] eq 'condB' || $split_req[0] eq 'featureB' ) {
				if ($analyze_condB) {
				}
				elsif ($analyze_running_galaxy) {
					close SETTINGS;
					system "rm -fr $main_output_folder";
					die "ERROR & EXIT: Missing $split_req[0] without default";
				}
				else {
					if ( $GROUP ne "" ) {
						chomp( my @fields = `cut -f \` head -1 $load_metadata | sed 's/\\t/\\n/g' | grep -n '$GROUP' | cut -f 1 -d':' \` $load_metadata | sort -u | grep -v $GROUP` );
						print "$sep\nAvailable conditions: " . join( ", ", @fields ) . "\n";
					}
					else {
						print "$sep\n";
					}
					print "Please enter condidtion/feature B:\n";
					chomp( $analyze_condB = <> );
					print "$sep\n";
				}
				print SETTINGS "condB\t$analyze_condB\n";
			}
		}

		# end cond A or B

		# iTOL support
		if ( $split_req[0] eq 'iTOL' ) {
			print SETTINGS "iTOL\tyes\n";
			$iTOL_used = 1;
		}

		# iTOL support

		# metadata
		if ( $split_req[0] eq 'metadata' ) {
			$load_metadata = load_metadata($load_metadata, $main_output_folder);
		}

		# end metadata

		# group
		if ( $split_req[0] eq 'groups' ) {
			$load_groups = "yes";
			if ($analyze_group) {
				open GROUPS, ">$calculations/R.groups" or die "ERROR & EXIT: Cannot write to $calculations/R.groups";
				my @grps = split ",", $analyze_group;
				$GROUP = $grps[0];
				foreach my $grp (@grps) {
					print GROUPS "$grp\n";
				}
				close GROUPS;
			}
			elsif ($analyze_running_galaxy) {
				close SETTINGS;
				system "rm -fr $main_output_folder";
				die "ERROR & EXIT: Missing groups without default";
			}
			else {
				$load_metadata = load_metadata($load_metadata, $main_output_folder);
				unless ( $split_req[1] ) {
					close SETTINGS;
					system "rm -fr $main_output_folder";
					die "ERROR & EXIT: Error in R script. Expected numeric after REQUIRE group in R script.";
				}
				my $groups = $split_req[1];
				$groups =~ s/^\s+//;
				$groups =~ s/\s+$//;
				if ( $groups > 1 ) {
					print "$sep\nThe $execute function requires $groups different groups to be entered.\n";
					print "Please choose groups from the metadata. You can choose between:\n";
				}
				else {
					print "$sep\nThe $execute function requires you to choose a metadata field to seperate\n";
					print "the data into two (or more) groups. Please choose one of the following fields to seperate by:\n";
				}
				my @metadata = `head -1 $load_metadata | sed 's/\\t/\\n/g'`;
				for my $md ( 1 .. scalar @metadata - 1 ) {
					print " - $metadata[$md]";
				}
				print "\n";
				open GROUPS, ">$calculations/R.groups" or die "ERROR & EXIT: Cannot write to $calculations/R.groups";
				for my $number ( 1 .. $groups ) {
					if ( $groups > 1 ) {
						print "Group $number: ";
					}
					elsif ( $groups == 1 ) {
						print "Metadata field: ";
					}
					else {
						close SETTINGS;
						system "rm -fr $main_output_folder";
						die "ERROR & EXIT: $groups doesn't seem to be a positive integer? Check R script.";
					}
					chomp( my $grp = <> );
					print GROUPS "$grp\n";
					$GROUP = $grp;
					my $exists = 0;
					foreach my $md (@metadata) {
						chomp $md;
						if ( $md eq $grp ) {
							$exists = 1;
						}
					}
					unless ($exists) {
						close SETTINGS;
						system "rm -fr $main_output_folder";
						die "ERROR & EXIT: $grp is no a valid option for metadata group.";
					}
				}
				close GROUPS;
			}
		}

		# end group
	}

	# feature annotations do not need to be specified, but if they're provided, we'll load them
	if ( defined( $analyze_feature_annotations[0] ) ) {
		unless ( $analyze_feature_annotations[0] eq '' ) {
			$load_feature_annotations = 'yes';
			open FEATURES, ">$calculations/R.feature_annotations" or die "ERROR & EXIT: Cannot write to $calculations/R.feature_annotations";
			for ( my $i = 0 ; $i <= scalar @analyze_feature_annotations - 1 ; $i += 2 ) {
				unless ( defined( $analyze_feature_annotations[ ( $i + 1 ) ] ) ) {
					close SETTINGS;
					system "rm -fr $main_output_folder";
					die "ERROR & EXIT: -analyze_feature_annotations does not have an even number of entries";
				}
				my $file = $analyze_feature_annotations[ $i + 1 ];
				unless ( -e $file ) {
					close SETTINGS;
					system "rm -fr $main_output_folder";
					die "ERROR & EXIT: Missing feature data file $file";
				}

				chomp( my $bn = `basename $file` );
				print localtime() . ": MOCATAnalyze : counting rows in $bn feature data\n";
				chomp( my $rows = `grep -c . $file` );
				print localtime() . ": MOCATAnalyze : counting columns in $bn feature data\n";
				chomp( my $columns = `awk -F\"\\t\" \'\{if (NF>max)\{max=NF\}\}; END\{print max\}\' $file` );
				print FEATURES "$analyze_feature_annotations[$i]\t$file\t$columns\t$rows\n";
			}
			close FEATURES;
		}

	}

	# Print some more settings
	print SETTINGS "load_groups\t$load_groups\n";
	print SETTINGS "load_metadata\t$load_metadata\n";
	print SETTINGS "load_feature_annotations\t$load_feature_annotations\n";
	close SETTINGS;

	open FUNCTIONS, ">$calculations/R.functions" or die "ERROR & EXIT: Cannot write to $calculations/R.functions";
	print FUNCTIONS "name\tlocation\n";

	# Running all is disabled.
	#if ( $execute eq 'all' ) {
	#	foreach my $f (@functions) {
	#		chomp $f;
	#		chomp( my $bn = `basename $f` );
	#		$bn =~ s/.R$//;
	#		print FUNCTIONS "$bn\t$f\n";
	#	}
	#}
	#else {
	#}

	unless ( -e "$conf{MOCAT_DEV_DIR}/analyze/$execute.R" ) {
		close SETTINGS;
		system "rm -fr $main_output_folder";
		die "ERROR & EXIT: Function $execute does not exist in the library.";
	}
	print FUNCTIONS "$execute\t$conf{MOCAT_DEV_DIR}/analyze/$execute.R\n";

	close FUNCTIONS;

	print "$sep\n" . localtime() . ": MOCATAnalyze : DATA and RESULTS are located in $main_output_folder/\n";
	print "$sep\n" . localtime() . ": MOCATAnalyze : ADDITIONAL ERROR MESSAGES ARE LOCATED IN $calculations/R.errors\n$sep\n";
	my $to_execute = "cd $calculations; R --silent --no-save -f $conf{MOCAT_DEV_DIR}/MOCATAnalyze_R.R 2>$calculations/R.errors | grep -ve '^>' -ve '^\+' -ve '^txt>' -ve '^txt+'";

	@samples = ();
	push @samples, "analysis";

	create_job( "analyze.$execute", 1, "$to_execute" );

	print localtime() . ": MOCATAnalyze : executing R-script\n$sep\n";

	MOCATCore::execute_job( "analyze.$execute", 1, 1, "15gb" );

	chomp( my $error = `grep -c 'Execution halted' $calculations/R.errors` );
	if ( $error > 0 ) {
		open IN, "<$calculations/R.errors";
		print "$sep\n";
		print "ERROR MESSAGES FORM R SCRIPT:\n";
		while (<IN>) {
			print;
		}
		die "$sep\nEXECUTION FAILED!\n$sep\n";
	}
	print "$sep\n";
	print localtime() . ": MOCATAnalyze : creating xls file\n";

	my @files    = `ls -1 $tables/*.table 2>/dev/null`;
	my $workbook = Spreadsheet::WriteExcel->new("$excel/results.xls");
	my $count    = 0;
	foreach my $file (@files) {
		$count++;
		chomp( my $name = `basename $file` );
		$name =~ s/\.table//;
		open IN, "<$file";
		my $i = -1;
		my $worksheet;
		if ( length($name) > 30 ) {
			$worksheet = $workbook->add_worksheet($count);
			$worksheet->write( 0, 0, $name );
			$i = 1;
		}
		else {
			$worksheet = $workbook->add_worksheet($name);
		}
		while (<IN>) {
			$i++;
			chomp;
			my @line = split ",";
			for my $j ( 0 .. scalar @line - 1 ) {
				$worksheet->write( $i, $j, $line[$j] );
			}
		}
		close IN;
	}
	$workbook->close();

	# iPATH support
	if ( $name[0] eq 'cog' || $name[0] eq 'ko' || $name[0] eq 'module' || $name[0] eq 'pathway' ) {
		print localtime() . ": MOCATAnalyze : initializing iPTAH\n";

		chomp( my @datasets = `ls -1 $data/*.multibar.data.iTOL $data/abundances.data.iTOL 2>/dev/null` );
		my $counter = 0;
		for my $dataset (@datasets) {
			my @name = split /\//, $dataset;
			my $name = $name[-1];
			$name =~ s/.data.iTOL//;
			print localtime() . ": MOCATAnalyze : generating $name iPTAH map\n";
			chomp( my @minmax = ` awk -F\",\" 'NR == 2 {max=\$2 ; min=\$2} \$2 >= max {max = \$2} \$2 <= min {min = \$2} END { print min\"\\n\"max }' $dataset` );
			my $min = $minmax[0];
			my $max = $minmax[1];
			open IN,  "awk 'NR>1' $dataset |" or die "ERROR & EXIT: Cannot open $dataset";
			open OUT, ">$dataset.iPATH"       or die "ERROR & EXIT: Cannot open $dataset.iPATH";
			my $header = <IN>;

			while (<IN>) {
				chomp;
				my @line = split /,/;
				$line[0] =~ s/ko//;
				if ( $line[1] < 0 ) {
					my $number = -$line[1] / -$min * 17 + 3;
					$number =~ s/\..*//;
					my $color = "#009ECE";
					print OUT "$line[0]\tW$number\t$color\n";
				}
				if ( $line[1] > 0 ) {
					my $number = $line[1] / $max * 17 + 3;
					$number =~ s/\..*//;
					my $color = "#9CCF31";
					print OUT "$line[0]\tW$number\t$color\n";

				}
			}
			system "echo 'GENERATE $figures/$analyze_date.$execute$desc.$name.iPATH.png' >> $calculations/iPATH.log";
			system "echo '$conf{MOCAT_DEV_DIR}/iPATH.pl --infile $dataset.iPATH --outfile $figures/$analyze_date.$execute$desc.$name.iPATH.png' >> $calculations/iPATH.log";
			system "$conf{MOCAT_DEV_DIR}/iPATH.pl --infile $dataset.iPATH --outfile $figures/$analyze_date.$execute$desc.$name.iPATH.png >> $calculations/iPATH.log";

		}

	}

	# iPATH support

	# iTOL support
	if ( $name[0] eq 'species' || $name[0] eq 'genus' || $name[0] eq 'family' || $name[0] eq 'order' || $name[0] eq 'class' || $name[0] eq 'phylum' || $name[0] eq 'kingdom' ) {
		if ($iTOL_used) {

			# Convert taxa names to taxa IDs
			print localtime() . ": MOCATAnalyze : initializing iTOL tree\n";
			print localtime() . ": MOCATAnalyze : converting iTOL taxa to IDs\n";
			system "cut -f 1 -d\" \" $data/taxa.iTOL | sort -u > $data/taxa.uniq.iTOL";
			
			chomp( our $systemType = `uname -s` );
			if ( $systemType =~ m/Darwin/ ) {
				system "sed -i '' -e 's/^/\\t\\|\\t/' -e 's/\$/\\t\\|/' $data/taxa.iTOL";
			}
			else {
				system "sed -i -e 's/^/\\t\\|\\t/' -e 's/\$/\\t\\|/' $data/taxa.iTOL";
			}
			
			system "fgrep -f $data/taxa.iTOL $conf{analyze_itol_ncbi_names_dmp} | cut -f 1 > $data/id.iTOL";
			system "fgrep -f $data/taxa.iTOL $conf{analyze_itol_ncbi_names_dmp} | cut -f 1,3 > $data/id2taxa.iTOL";

			print localtime() . ": MOCATAnalyze : generating color file\n";
			open IN, "<$data/taxa.uniq.iTOL";
			my @color_scale = ( "#FFB6C1", "#FFA07A", "#FFFACD", "#FFE4B5", "#E6E6FA", "#87CEEB", "#AFEEEE", "#B0E0E6", "#66CDAA", "#FFF8DC", "#DEB887", "#DAA520", "#F0FFF0", "#FFE4E1", "#DCDCDC", "#778899" );
			my %colors;
			my $c = 0;
			while (<IN>) {
				chomp;
				$colors{$_} = $color_scale[$c];
				$c++;
				if ( $c > scalar @color_scale - 1 ) {
					$c = 0;
				}
			}
			close IN;
			open IN,  "<$data/id2taxa.iTOL";
			open OUT, ">$data/colors.iTOL";
			while (<IN>) {
				chomp;
				my @line = split "\t";
				my @l2 = split /\s+/, $line[1];
				my $col;
				if ( $colors{ $l2[0] } ) {
					$col = $colors{ $l2[0] };
				}
				else {
					$col = "#FFFFFF";
				}
				print OUT "$line[0]\trange\t$col\trandom color\n";
			}
			close IN;
			close OUT;

			open OUT, ">$calculations/iTOL.config" or die "ERROR & EXIT: Cannot create $calculations/iTOL.config";
			print OUT "treeName = MyTree\n";
			print OUT "assignTaxonomy = 1\n";
			print OUT "ncbiFormat = idsFull\n";
			print OUT "ncbiFile = $data/id.iTOL\n\n";

			chomp( my @datasets = `ls -1 $data/*.data.iTOL` );
			my $counter = 0;

			my @colors  = ( "#0000FF",   "#FFCC00",   "#003366",   "#FFBAD2",   "#34ACAF",   "#9D538E",   "#000000",   "#009966" );
			my @colors2 = ( "\\#0000FF", "\\#FFCC00", "\\#003366", "\\#FFBAD2", "\\#34ACAF", "\\#9D538E", "\\#000000", "\\#009966" );
			my @datasetlist;
			foreach my $dataset (@datasets) {
				chomp( my $iTOL_label = `basename $dataset` );
				$iTOL_label =~ s/.data.iTOL//;
				print localtime() . ": MOCATAnalyze : processing $iTOL_label\n";
				system "awk 'NR>1' $dataset | cut -f 1 -d',' | sed -e 's/^/\\t\\|\\t/' -e 's/\$/\\t\\|/' | fgrep -f - $conf{analyze_itol_ncbi_names_dmp} | sed 's/\\t|\\t/\\t/g' | cut -f 1,2 > $dataset.map";
				open MAP, "<$dataset.map";
				my %map;
				while (<MAP>) {
					chomp;
					@line = split /\t/;
					$map{ $line[1] } = $line[0];
				}
				close MAP;
				open DATA,    "<$dataset"    or die;
				open DATAOUT, ">$dataset.id" or die;
				while (<DATA>) {
					chomp;
					my $line = $_;
					my @line = split /,/, $line;
					unless ( $line[0] eq "" ) {
						if ( $map{ $line[0] } ) {
							$line =~ s/$line[0]/$map{$line[0]}/;
							print DATAOUT "$line\n";
						}
						else {
							print "MISSING MAPPING NAME FOR $line[0]. THIS TAXON WILL NOT BE DISPLAYED IN iTOL TREE!\n";
						}
					}
				}
				close DATA;
				close DATAOUT;

				$counter++;
				if ( $counter > 10 ) {
					die "ERROR & EXIT: Only 10 datasets for each iTOL tree is supported. Please change the number of output tables in the iTOL variable in your R script.";
				}
				print OUT "dataset$counter" . "File = $dataset.id\n";
				print OUT "dataset$counter" . "Label = $iTOL_label\n";
				print OUT "dataset$counter" . "Type = simplebar\n";
				print OUT "dataset$counter" . "Separator = comma\n";
				print OUT "dataset$counter" . "Color = $colors[($counter-1)]\n";
				print OUT "dataset$counter" . "BarSizeMax = 1000\n";
				print OUT "dataset$counter" . "BranchColoringType = both\n\n";
				push @datasetlist, "dataset$counter";
			}
			close OUT;

			my @colors_to_use;
			print localtime() . ": MOCATAnalyze : generating combined data file\n";
			my %main;
			my $main_sum = 0;
			open IN, "<$data/abundances.data.iTOL.id";
			while (<IN>) {
				chomp;
				my @line = split /,/;
				if ( $line[1] > $main_sum ) {
					$main_sum = $line[1];
				}
			}
			close IN;

			open IN, "<$data/abundances.data.iTOL.id";
			while (<IN>) {
				chomp;
				my @line = split /,/;
				push @{ $main{ $line[0] } }, $line[1] / $main_sum * 100;

			}
			close IN;

			# multibar data begin
			chomp( @datasets = `ls -1 $data/*.multibar.data.iTOL.id` );
			my @labels;
			$counter = -1;
			foreach my $dataset (@datasets) {
				my $temp_max = 0.00000000001;
				my $temp_min = 0.00000000001;
				$counter++;
				push @colors_to_use, ( "#FFFFFF", $colors[$counter] );
				open IN, "<$dataset";
				my %temp;
				chomp( my $iTOL_label = `basename $dataset` );
				$iTOL_label =~ s/.multibar.data.iTOL.id//;
				push @labels, ( "empty", $iTOL_label );

				while (<IN>) {
					chomp;
					my @line = split /,/;
					$temp{ $line[0] } = $line[1];
					if ( $line[1] =~ m/-(\d+)/ ) {
						if ( $1 > $temp_min ) {
							$temp_min = $1;
						}
					}
					elsif ( $line[1] =~ m/(\d+)/ ) {
						if ( $line[1] > $temp_max ) {
							$temp_max = $line[1];
						}
					}
					else {
						die "ERROR & EXIT: Unexpected value: $line[1] in $dataset";
					}

				}
				foreach my $key ( keys %main ) {
					if ( $temp{$key} ) {
						if ( $temp{$key} =~ m/-(\d+)/ ) {
							push @{ $main{$key} }, ($temp_min) / ( $temp_max + $temp_min ) * 100, "-" . $1 / $temp_min * ( ($temp_min) / ( $temp_max + $temp_min ) * 100 );
						}
						else {
							push @{ $main{$key} }, ($temp_min) / ( $temp_max + $temp_min ) * 100, $temp{$key} / $temp_max * ( ($temp_max) / ( $temp_max + $temp_min ) * 100 );
						}
					}
					else {
						push @{ $main{$key} }, ($temp_min) / ( $temp_max + $temp_min ) * 100, "0";
					}
				}
				close IN;
			}
			open OUT, ">$calculations/iTOL.config.v2" or die "ERROR & EXIT: Cannot create $calculations/iTOL.config.v2";
			print OUT "treeName = MyTree\n";
			print OUT "assignTaxonomy = 1\n";
			print OUT "ncbiFormat = idsFull\n";
			print OUT "ncbiFile = $data/id.iTOL\n";
			print OUT "colorDefinitionFile = $data/colors.iTOL\n\n";
			print OUT "dataset1File = $data/data.iTOL.v2\n";
			print OUT "dataset1Label = Data\n";
			print OUT "dataset1Type = multibar\n";
			print OUT "dataset1MultiBarAlign = 1\n";
			print OUT "dataset1Separator = comma\n";
			print OUT "dataset1BarSizeMax = 1000\n";

			chomp( @datasets = `ls -1 $data/*.colorscale.data.iTOL.id 2>/dev/null` );
			$counter = 1;
			foreach my $dataset (@datasets) {
				chomp( my $lines = `grep -c . $dataset` );
				if ( $lines > 0 ) {
					$counter++;
					chomp( my $iTOL_label = `basename $dataset` );
					system "cut -f 1,2 -d\",\" $dataset > $dataset.tmp; mv $dataset.tmp $dataset";
					chomp( my $color = `head -1 $dataset | cut -f 2 -d','` );
					$iTOL_label =~ s/.colorscale.data.iTOL.id//;
					print OUT "\ndataset$counter" . "File = $dataset\n";
					print OUT "dataset$counter" . "Label = $iTOL_label\n";
					print OUT "dataset$counter" . "Type = colorstrip\n";
					print OUT "dataset$counter" . "Separator = comma\n";
					print OUT "dataset$counter" . "Color = \\$color\n";
				}
			}

			close OUT;
			open DATA, ">$data/data.iTOL.v2" or die "ERROR & EXIT: Cannot create $data/data.iTOL.v2";
			print DATA "LABELS,abundances," . join( ",", @labels ) . "\n";
			print DATA "COLORS,#FF0000," . join( ",", @colors_to_use ) . "\n";
			foreach my $key ( keys %main ) {
				print DATA "$key," . join( ",", @{ $main{$key} } ) . "\n";
			}
			close DATA;

			# multibar data end

			print localtime() . ": MOCATAnalyze : generating and uploading iTOL tree\n";
			my $executeiTOL = "$conf{MOCAT_DEV_DIR}/iTOL_uploader.pl --config $calculations/iTOL.config.v2";
			system "echo ' == $executeiTOL ==' > $calculations/iTOL.log";
			system "$executeiTOL >> $calculations/iTOL.log";

			open IN, "<$calculations/iTOL.log" or die "ERROR & EXIT: MISSING $calculations/iTOL.log";
			my $ID = 0;
			while (<IN>) {
				chomp;
				if (m/^(\d+)$/) {
					$ID = $1;
				}
			}
			if ( $ID == 0 ) {
				die "ERROR & EXIT: Seems like iTOL tree was not uploaded correctly. Could not fetch ID from $calculations/iTOL.log";
			}
			print localtime() . ": MOCATAnalyze : downloading and saving iTOL tree\n";
			print localtime() . ": MOCATAnalyze : iTOL tree ID is $ID\n";

			$executeiTOL = "$conf{MOCAT_DEV_DIR}/iTOL_downloader.pl --tree $ID --outputFile $figures/$analyze_date.$execute$desc.iTOL.tree.pdf --format pdf --displayMode circular --colorBranches --rotation 200 --omitDashedLines 1 --fontSize 24 --arc 355 --datasetList " . join( ",", @datasetlist );
			system "echo ' == $executeiTOL ==' >> $calculations/iTOL.log";
			system "$executeiTOL >> $calculations/iTOL.log";

			$executeiTOL = "$conf{MOCAT_DEV_DIR}/iTOL_downloader.pl --tree $ID --outputFile $figures/$analyze_date.$execute$desc.iTOL.tree.svg --format svg --displayMode circular --colorBranches --rotation 200 --omitDashedLines 1 --fontSize 24 --arc 355 --datasetList " . join( ",", @datasetlist );
                        system "echo ' == $executeiTOL ==' >> $calculations/iTOL.log";
                        system "$executeiTOL >> $calculations/iTOL.log";
		}
	}

	# iTOL support

	if ( $analyze_output[0] ) {
		for ( my $i = 0 ; $i <= scalar @analyze_output - 1 ; $i += 2 ) {
			unless ( defined( $analyze_output[ ( $i + 1 ) ] ) ) {
				die "ERROR & EXIT: -analyze_output does not have an even number of entries";
			}
			my $path;
			if ( $analyze_output[$i] eq 'xls' ) {
				$path = "$excel/results.xls";
			}
			elsif ( $analyze_output[$i] eq 'pdf' ) {
				$path = "$figures/*.pdf";
			}
			elsif ( $analyze_output[$i] eq 'jpg' ) {
				$path = "$figures/*.jpg";
			}
			elsif ( $analyze_output[$i] eq 'csv' ) {
				$path = "$tables/*.table";
			}
			elsif ( $analyze_output[$i] eq 'RData' ) {
				$path = "$data/session.RData";
			}
			else {
				die "ERROR & EXIT: No export function defined for file format '$analyze_output[$i]'. This need to be coded into MOCATAnalyze.pm script.";
			}
			my @n = `ls -1 $path`;
			if ( scalar @n > 1 ) {
				print localtime() . ": MOCATAnalyze : CANNOT EXPORT $analyze_output[$i] to $path because there are multiple files of the type\n";
			}
			else {
				print localtime() . ": MOCATAnalyze : exporting $analyze_output[$i] ($path -> " . $analyze_output[ $i + 1 ] . ")\n";
				system "cp $path " . $analyze_output[ $i + 1 ];
			}
		}
	}

	print "$sep\n";
	print localtime() . ": MOCATAnalyze : done\n";
}

sub load_metadata {
	my $load_metadata      = $_[0];
	my $main_output_folder = $_[1];
	my $metadata_ok        = 0;
	if ( scalar @analyze_metadata == 1 && $analyze_project_name eq '' ) {
		$load_metadata = $analyze_metadata[0];
		$metadata_ok   = 1;
	}
	elsif ( scalar @analyze_metadata % 2 == 0 && $analyze_project_name eq '' ) {
		close SETTINGS;
		system "rm -fr $main_output_folder";
		die "ERROR & EXIT: Metadata seems to be defined in config or command line, but not -analyze_project_name to use to match metadata";
	}
	elsif ( $analyze_project_name ne '' ) {
		for ( my $i = 0 ; $i <= scalar @analyze_metadata - 1 ; $i += 2 ) {
			if ( $analyze_metadata[$i] eq $analyze_project_name ) {
				$metadata_ok   = 1;
				$load_metadata = $analyze_metadata[ $i + 1 ];
			}
		}
	}
	else {
		close SETTINGS;
		system "rm -fr $main_output_folder";
		die "ERROR & EXIT: Do either: 1. Define a set of metadata in config file (separated by comma in config, but whitespace on command line) or as analyze_metadata with\n[project name] [metadata file] [project name] [metadata file]...\nAND define a project name using -analyze_project_name NAME, OR 2. set the analyze_metadata to only [metadata file].";
	}

	if ( $metadata_ok == 0 ) {
		close SETTINGS;
		system "rm -fr $main_output_folder";
		die "ERROR & EXIT: Metadata wasn't loaded correctly. Did you specify it in either config or on commandline using -analyze_metadata field?";
	}
	unless ( -e $load_metadata ) {
		close SETTINGS;
		system "rm -fr $main_output_folder";
		die "ERROR & EXIT: Metadata file $load_metadata does not exist.";
	}
	my $list = join @samples, "|";
	my @samples_in_meta = `cut -f 1 $load_metadata | awk 'NR>1' | grep -E \"$list\" - | sort -u`;
	my %samples_in_meta;
	foreach my $sample (@samples_in_meta) {
		chomp $sample;
		$samples_in_meta{$sample} = 1;
	}
	my $all_exists = 1;
	my %samples_not_in_meta;
	foreach my $sample (@samples) {
		unless ( $samples_in_meta{$sample} ) {
			$samples_not_in_meta{$sample} = 1;
			$all_exists = 0;
		}
	}
	if ( $all_exists == 0 ) {
		print "ERROR & EXIT: The folowing samples do not have metadata in $load_metadata:";
		foreach my $sample ( sort keys %samples_not_in_meta ) {
			print " - $sample\n";
		}
		close SETTINGS;
		system "rm -fr $main_output_folder";
		die "Please add metadata for these samples.\n";
	}
	return ($load_metadata);
}

sub create_job {
	### DEFINE VARIABLES AND OPEN OUTPUT FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	my $command = $_[2];

	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?";
	print localtime() . ": Creating $job jobs...";

	### JOB SPECIFIC ###

	my $LOG = "";
	if ( $conf{MOCAT_qsub_system} eq 'SGE' ) {
		$LOG = " >> $cwd/logs/$job/samples/MOCATJob_$job.$samples[0].log ";
	}

	print JOB "$command $LOG\n";

	print " OK!\n";

}

1;
