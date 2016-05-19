package MOCATProfiling::Print;
use strict;
use warnings;
use MOCATProfiling::Variables;

sub printFiles {
	#############################################################################################################################################
	# LOAD VARIABLES
	#############################################################################################################################################
	my $HASH_pointer        = shift;
	my $level               = shift;
	my $horizontal_position = shift;    # this is set to 1=INIT_NUM+1=11 for genes, and 4 for others4=INIT_NUM+6=16
	my $position            = 2;        # the position variable here keeps track of where in the HASH_pointer, and ultimately the @array we should find the count values to use
	                                    # the first two values 0,1 are start and stop, the third value is the first base raw thus we should set position to 2
	#############################################################################################################################################

	#############################################################################################################################################
	# open files for output and save their file handles in a hash
	#############################################################################################################################################
	my %open_files = ();
	for my $distType ( '', '.only.unique', '.mm.dist.among.unique' ) {
		for my $countType ( 'base', 'insert' ) {
			for my $file ( 'raw', 'norm', 'scaled' ) {
				open my $FH, '>', "$temp_file.output/$output_file_base.$countType$distType.$file.$level" or die "ERROR & EXIT: Cannot write to file $temp_file.output/$output_file_base.$countType$distType.$file.$level $!";
				$open_files{"$countType$distType.$file"} = $FH;
				print $FH "$sample_name\n";
				if ( ( $mode ne 'mOTU' && $mode ne 'functional' || $level eq 'gene' ) || ( $mode eq 'functional' && $level eq 'gene' ) ) {    # funcitonal profiles have a different header, see below
					if ( $file eq 'raw' || $file eq 'scaled' ) {                                                                        # currently, raw and scaled is the only one that has the -1 fraction printed                                                                                                                                                       # print the -1 fraction for raw files, and NA for other files, except for NCBI mode, print norm as well
						print $FH $TOTAL{$countType} - $TOTAL_MAPPED{$countType}{$level} . "\n";
					}

					# a note here, the norm is only defined if we have first mapped against mOTU DB and then RefMG db, in this second step we'd use the NCBI mode
					# we can't control that the user will always do it in this order, but at least we can only calculate the norm -1 if mode is NCBI
					elsif ( $mode eq 'NCBI' && $file eq 'norm' && $level ne 'gene' ) {
						print $FH ( $TOTAL{$countType} - $TOTAL_MAPPED{$countType}{$level} ) / $ALL_TAXA_MEAN_LENGTH . "\n";
					}
					else {
						print $FH "NA\n";
					}
				}
			}
		}
	}
	foreach my $end ( '', '.uniq' ) {
		for my $countType ( 'base', 'insert' ) {
			if ($print_rownames_file) {    # when running MOCAT, only the first job in a main job array should have a rownames file specified, this way this file is only printed once
				                     # create links for rowname files
				open my $FH, '>', "$rownames_file.$level$end.$date" or die "ERROR & EXIT: Cannot write to file $rownames_file$end.$date";
				$open_files{"rownames$end"} = $FH;
				if ( $firstColumnName ne "" ) {
					print $FH "$firstColumnName\n";    # support for custom names for scaftigs, contigs, ...
				}
				else {
					print $FH "$level\n";
				}

				#these are the two first header rows
				if ( ( $mode ne 'mOTU' && $mode ne 'functional' || $level eq 'gene' ) || ( $mode eq 'functional' && $level eq 'gene' ) ) {
					print $FH "-1\n";
				}
			}
			system_("ln -fs $rownames_file.$level$end $temp_file.output/$output_file_base.$level.rownames$end");
		}
	}
	if ($CALCULATE_HORIZONTAL_COVERAGE) {
		open my $FH, '>', "$temp_file.output/$output_file_base.horizontal.$level" or die "ERROR & EXIT: Cannot write to file $temp_file.output/$output_file_base.horizontal.$level $!";
		$open_files{horizon} = $FH;
		print $FH "$sample_name\n";
		print $FH "NA\n";
		if ( $mode eq 'functional' && $level ne 'gene' ) {
			print $FH "NA\n";
		}
	}
#############################################################################################################################################

#############################################################################################################################################
	# add header lines for functional profiles
#############################################################################################################################################
	if ( $mode eq 'functional' && $level ne 'gene' ) {
		for my $distType ( '', '.only.unique', '.mm.dist.among.unique' ) {
			for my $countType ( 'base', 'insert' ) {
				for my $normType ( 'raw', 'norm', 'scaled' ) {
					my $value;    # this is used for the total mapped below

					# this is a bit tricky, for mm dist we have to sum the unique and mm dist of the genes
					if ( $distType eq '.mm.dist.among.unique' ) {
						if ( $normType eq 'norm' ) {
							$value = $TOTAL_MAPPED_LENGTH_NORM{"$countType.only.unique"}{'gene'} + $TOTAL_MAPPED_LENGTH_NORM{"$countType.mm.dist.among.unique"}{'gene'};    # this should always be gene
						}
						elsif ( $normType eq 'raw' || $normType eq 'scaled' ) {
							$value = $TOTAL_MAPPED{"$countType.only.unique"}{'gene'} + $TOTAL_MAPPED{"$countType.mm.dist.among.unique"}{'gene'};                            # this should always be gene
						}

					}
					else {                                                                                                                                                              # do this for '' and only.unique
						if ( $normType eq 'norm' ) {
							$value = $TOTAL_MAPPED_LENGTH_NORM{"$countType$distType"}{'gene'};                                                                              # this should always be gene
						}
						elsif ( $normType eq 'raw' || $normType eq 'scaled' ) {
							$value = $TOTAL_MAPPED{"$countType$distType"}{'gene'};                                                                                          # this should always be gene
						}
					}

					unless ( $TOTAL_MAPPED_ANNOTATED{$level}{"$countType$distType.$normType"} ) {
						$TOTAL_MAPPED_ANNOTATED{$level}{"$countType$distType.$normType"} = 0;
					}
					unless ( $TOTAL_MAPPED_NOT_ANNOTATED{$level}{"$countType$distType.$normType"} ) {
						$TOTAL_MAPPED_NOT_ANNOTATED{$level}{"$countType$distType.$normType"} = 0;
					}

					# scaled files are same as raw files, but the file name is different, that's why they're outside they loop above
					my $normType2 = $normType;
					if ( $normType eq 'scaled' ) {
						$normType2 = "raw";
					}

					# In the new headers we only want to print the mapped and unassigned
					if ($USE_OLD_HEADERS) {
						print { $open_files{"$countType$distType.$normType"} } ( $TOTAL{"$countType$distType"} . "\n" );
					}

					print { $open_files{"$countType$distType.$normType"} } ( $value . "\n" );

					if ($USE_OLD_HEADERS) {
						print { $open_files{"$countType$distType.$normType"} } ( $TOTAL_MAPPED_ANNOTATED{$level}{"$countType$distType.$normType2"} + $TOTAL_MAPPED_NOT_ANNOTATED{$level}{"$countType$distType.$normType2"} . "\n" );
					}

					print { $open_files{"$countType$distType.$normType"} } ( $TOTAL_MAPPED_NOT_ANNOTATED{$level}{"$countType$distType.$normType2"} . "\n" );

					if ($USE_OLD_HEADERS) {
						print { $open_files{"$countType$distType.$normType"} } ( $TOTAL_MAPPED_ANNOTATED{$level}{"$countType$distType.$normType2"} . "\n" );
					}
				}
			}
		}
		if ($print_rownames_file) {
			foreach my $end ( '', '.uniq' ) {
				if ($USE_OLD_HEADERS) {
					print { $open_files{"rownames$end"} } ("total\n");
				}
				print { $open_files{"rownames$end"} } ("mapped\n");
				if ($USE_OLD_HEADERS) {
					print { $open_files{"rownames$end"} } ("assigned_and_unassigned\n");
				}
				print { $open_files{"rownames$end"} } ("unassigned\n");
				if ($USE_OLD_HEADERS) {
					print { $open_files{"rownames$end"} } ("assigned\n");
				}
			}
		}
	}
#############################################################################################################################################

#############################################################################################################################################
	# calculate weights
#############################################################################################################################################
	my %weighted_mean;    # here we store the weighted means for abses and inserts that are used firther down
	for my $distType ( '', '.only.unique', '.mm.dist.among.unique' ) {
		for my $countType ( 'base', 'insert' ) {
			unless ( $TOTAL_MAPPED{"$countType$distType"}{$level} ) {
				$TOTAL_MAPPED{"$countType$distType"}{$level}             = 0;
				$TOTAL_MAPPED_LENGTH_NORM{"$countType$distType"}{$level} = 0;
				$weighted_mean{"$countType$distType"}                    = 0;
			}

			if ( $distType eq '.mm.dist.among.unique' ) {    # mm dist is the sum of mm dist and only unique
				if ( ( $TOTAL_MAPPED_LENGTH_NORM{"$countType$distType"}{$level} + $TOTAL_MAPPED_LENGTH_NORM{"$countType.only.unique"}{$level} ) > 0 ) {
					$weighted_mean{"$countType$distType"} = ( $TOTAL_MAPPED{"$countType$distType"}{$level} + $TOTAL_MAPPED{"$countType.only.unique"}{$level} ) / ( $TOTAL_MAPPED_LENGTH_NORM{"$countType$distType"}{$level} + $TOTAL_MAPPED_LENGTH_NORM{"$countType.only.unique"}{$level} );
				}
				else {
					$weighted_mean{"$countType$distType"} = 0;
				}
				if ($VERBOSE) {
					print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : level=$level : weighted_mean{$countType$distType} = " . ( $TOTAL_MAPPED{"$countType$distType"}{$level} + $TOTAL_MAPPED{"$countType.only.unique"}{$level} ) . " / " . ( $TOTAL_MAPPED_LENGTH_NORM{"$countType$distType"}{$level} + $TOTAL_MAPPED_LENGTH_NORM{"$countType.only.unique"}{$level} ) . " = " . $weighted_mean{"$countType$distType"} . "\n";    #TEXT#
				}
			}
			else {
				if ( $TOTAL_MAPPED_LENGTH_NORM{"$countType$distType"}{$level} > 0 ) {
					$weighted_mean{"$countType$distType"} = $TOTAL_MAPPED{"$countType$distType"}{$level} / $TOTAL_MAPPED_LENGTH_NORM{"$countType$distType"}{$level};
				}
				else {
					$weighted_mean{"$countType$distType"} = 0;
				}
				if ($VERBOSE) {
					print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : level=$level : weighted_mean{$countType$distType} = " . $TOTAL_MAPPED{"$countType$distType"}{$level} . " / " . $TOTAL_MAPPED_LENGTH_NORM{"$countType$distType"}{$level} . " = " . $weighted_mean{"$countType$distType"} . "\n";                                                                                                                            #TEXT#
				}
			}
		}
	}
#############################################################################################################################################

#############################################################################################################################################
	# Actually print the files
#############################################################################################################################################
	if ($VERBOSE) {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : printing files for level=$level to $temp_file.output\n";
	}

	# loop over the entries in the HASH and print to each file
	foreach my $ref_id_index ( sort keys %{$HASH_pointer} ) {
		my @array = @{ $HASH_pointer->{$ref_id_index} };

		### EXTRACT REF ID ###
		if ( $level eq 'gene' ) {
			$ref_id_index =~ /^(.*)\.[0-9]+$/;
			$ref_id = $1;
		}
		else {
			$ref_id = $ref_id_index;
		}
		### EXTRACT REF ID ###

		my $length;
		if ( $level eq 'gene' ) {
			$length = $array[1] - $array[0] + 1;
		}
		my %values;
		for my $countType ( 'base', 'insert' ) {
			my $position2;    # we have to get another position index becaus ethe insert values are one position further than the base values
			if ( $countType eq 'base' ) {
				$position2 = $position;
			}
			elsif ( $countType eq 'insert' ) {
				$position2 = $position + 1;
			}

			for my $fill ( 0 .. 12 ) {    # we have to do this because when we initialize the large hash we don't store 0s in it :) 0-6 is the normal non norm, then there are 2 mm dist norm and then 4 more
				unless ( $array[ $position2 + $fill ] ) {
					$array[ $position2 + $fill ] = 0;
				}
			}

			# good luck! have a look at the hash layout at the top of the script, note the +2 for only.unique and the +4 and +6 for mm.dist
			print { $open_files{"$countType.raw"} }                      ( $array[$position2] . "\n" );
			print { $open_files{"$countType.only.unique.raw"} }          ( $array[ $position2 + 2 ] . "\n" );
			print { $open_files{"$countType.mm.dist.among.unique.raw"} } ( ( $array[ $position2 + 2 ] ) + ( $array[ $position2 + 4 ] ) . "\n" );

			# unfortunately there are two scenarios, for genes we can calculate the norm values, but for other ones that are merged form gene norm values they are stored in HASH_VALUES position 10-13
			if ( $level eq 'gene' ) {
				print { $open_files{"$countType.norm"} }                      ( $array[$position2] / $length . "\n" );
				print { $open_files{"$countType.only.unique.norm"} }          ( $array[ $position2 + 2 ] / $length . "\n" );
				print { $open_files{"$countType.mm.dist.among.unique.norm"} } ( ( $array[ $position2 + 2 ] ) / $length + ( $array[ $position2 + 6 ] ) . "\n" );

				print { $open_files{"$countType.scaled"} }                      ( $array[$position2] / $length * $weighted_mean{"$countType"} . "\n" );
				print { $open_files{"$countType.only.unique.scaled"} }          ( $array[ $position2 + 2 ] / $length * $weighted_mean{"$countType.only.unique"} . "\n" );
				print { $open_files{"$countType.mm.dist.among.unique.scaled"} } ( ( ( $array[ $position2 + 2 ] ) / $length + ( $array[ $position2 + 6 ] ) ) * $weighted_mean{"$countType.mm.dist.among.unique"} . "\n" );
			}
			else {
				print { $open_files{"$countType.norm"} }                      ( $array[ $position2 + $INIT_NUM - 2 ] . "\n" );
				print { $open_files{"$countType.only.unique.norm"} }          ( $array[ $position2 + 2 + $INIT_NUM - 2 ] . "\n" );
				print { $open_files{"$countType.mm.dist.among.unique.norm"} } ( ( $array[ $position2 + 2 + $INIT_NUM - 2 ] ) + ( $array[ $position2 + 6 ] ) . "\n" );    # the mm dist norm values are already stored

				print { $open_files{"$countType.scaled"} }                      ( $array[ $position2 + $INIT_NUM - 2 ] * $weighted_mean{"$countType"} . "\n" );
				print { $open_files{"$countType.only.unique.scaled"} }          ( $array[ $position2 + 2 + $INIT_NUM - 2 ] * $weighted_mean{"$countType.only.unique"} . "\n" );
				print { $open_files{"$countType.mm.dist.among.unique.scaled"} } ( ( ( $array[ $position2 + 2 + $INIT_NUM - 2 ] ) + ( $array[ $position2 + 6 ] ) ) * $weighted_mean{"$countType.mm.dist.among.unique"} . "\n" );
			}
		}
		if ( $open_files{"rownames"} ) {
			print { $open_files{"rownames"} } ("$ref_id\n");
			if ( $level eq 'gene' ) {
				print { $open_files{"rownames.uniq"} } ("$ref_id\_$array[0]\_$array[1]\n");
			}
			else {

				# we dont print these files. They are created, but we leave them empty, easiest fix, maybe someone will ask why they dont exist...
			}
		}
		if ($CALCULATE_HORIZONTAL_COVERAGE) {
			if ( $array[ $INIT_NUM + $horizontal_position ] ) {
				if ( $horizontal_position == 1 ) {
					print { $open_files{horizon} } ( ( $array[ $INIT_NUM + $horizontal_position ] =~ tr/1// ) . "\n" );
				}
				else {
					print { $open_files{horizon} } ( ( $array[ $INIT_NUM + $horizontal_position ] / $array[ $INIT_NUM + $horizontal_position + 1 ] ) . "\n" );
				}
			}
			else {
				print { $open_files{horizon} } ("0\n");
			}
		}
	}
#############################################################################################################################################

#############################################################################################################################################
	# close files
#############################################################################################################################################
	for my $countType ( 'base', 'insert' ) {
		for my $file ( 'raw', 'norm', 'scaled' ) {
			close $open_files{"$countType.$file"} or die "ERROR & EXIT: Closing pipe failed, this means something is wrong. Corrput files?";
		}
	}
	if ($CALCULATE_HORIZONTAL_COVERAGE) {
		close $open_files{horizon} or die "ERROR & EXIT: Closing pipe failed, this means something is wrong. Corrput files?";
	}
#############################################################################################################################################

#############################################################################################################################################
	# move header files
#############################################################################################################################################
	if ($print_rownames_file) {
		foreach my $end ( '', '.uniq' ) {
			system_("mv $rownames_file.$level$end.$date $rownames_file.$level$end");
		}
	}
#############################################################################################################################################

}

sub printStats {
	#############################################################################################################################################
	# PRINT MAIN STATS FILE
	#############################################################################################################################################
	open my $out, ">", $out_stats_file or die "ERROR & EXIT: Cannot write to $out_stats_file: $!";
	print $out "total_inserts\tmapped_inserts\tfraction_mapped_inserts\ttotal_bases\tmapped_bases\tfraction_mapped_bases\tdb_average_entry_length\n";
	my $frI = 0;
	my $frB = 0;
	if ( $TOTAL{insert} ) {
		$frI = $TOTAL_MAPPED{insert}{gene} / $TOTAL{insert};
	}
	if ( $TOTAL{base} ) {
		$frB = $TOTAL_MAPPED{base}{gene} / $TOTAL{base};
	}
	print $out "$TOTAL{insert}\t$TOTAL_MAPPED{insert}{gene}\t$frI\t$TOTAL{base}\t$TOTAL_MAPPED{base}{gene}\t$frB\t$db_average_entry_length\n";
	close $out;

	#############################################################################################################################################
	# PRINT PE STATS FILE
	#############################################################################################################################################
	for my $i ( 6 .. 8 ) {    # support for old stats files that didn't have these columns
		unless ( $STATS_FILE_DATA[$i] ) {
			$STATS_FILE_DATA[$i] = "NA";
		}
	}
	if ($VERBOSE) {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : printed $out_stats_file\n";
	}
	open $out, ">", $out_PE_stats_file or die "ERROR & EXIT: Cannot write to $out_PE_stats_file: $!";
	print $out "Reads\tBases\tMax\tAvg\tKmer\tInserts Min % identity\tMin length\tSOAP max mismatches\n";
	print $out "NA\t$TOTAL_MAPPED{base}{gene}\t$STATS_FILE_DATA[2]\t$STATS_FILE_DATA[3]\t$STATS_FILE_DATA[4]\t$TOTAL_MAPPED{insert}{gene}\t$STATS_FILE_DATA[6]\t$STATS_FILE_DATA[7]\t$STATS_FILE_DATA[8]\n";
	close $out;
	if ($VERBOSE) {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : printed $out_PE_stats_file\n";
	}
	#############################################################################################################################################
}

sub system_ {
	my $cmd = shift;
	system($cmd) == 0 or die("\nERROR & EXIT: system($cmd) failed: $!");
}

1;
