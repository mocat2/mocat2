package MOCATTaxoProfiling;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

sub create_job {

	### DEFINE VARIABLES AND OPEN input FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?";
	print localtime() . ": Creating $job jobs...\n";
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}
	
	$taxo_profiling_mode = $profiling_mode;

	# Reset arrays used to save file names
	@taxo_profiling_base   = ();
	@taxo_profiling_insert = ();
	my $databases = join( "_AND_", @do_taxo_profiling );

	# This was stored in the config previously, now we require the files to be DB specific
	
	my $old = '';
	if ($use_old_version) {
		$old = '.old';
	}
	my $taxo_profiling_map         = "$data_dir/$databases.refmg.map$old";
	my $taxo_profiling_map_tot_len = "$data_dir/$databases.refmg.map$old.total_length";
	my $taxo_profiling_motu_map    = "$data_dir/$databases.motu.map$old";

	### CHECKS ###
	print localtime() . ": Checking mode...";
	if ( $taxo_profiling_mode eq 'mOTU' ) {
		$taxo_profiling_map = $taxo_profiling_motu_map;
	}
	elsif ( $taxo_profiling_mode eq 'RefMG' ) {
		unless ( -e "$taxo_profiling_map_tot_len" ) {
			die "\nERROR & EXIT: Missing length file $taxo_profiling_map_tot_len";
		}
	}
	else {
		die "\nERROR & EXIT: -mode must be either 'mOTU' or 'RefMG'";
	}
	print " OK!\n";
	print localtime() . ": Checking map file...";
	unless ( -e "$taxo_profiling_map" ) {
		die "\nERROR & EXIT: Missing map file $taxo_profiling_map";
	}
	my @tmp = `awk -F\"\\t\" '{print NF}' $taxo_profiling_map | sort -u`;
	my %map;
	foreach my $tmp (@tmp) {
		chomp($tmp);
		if ( $taxo_profiling_mode eq 'RefMG' ) {
			unless ( $tmp eq '9' ) {
				die "\nERROR & EXIT: Mapping file ($taxo_profiling_map) does not have exactly 9 columns on all lines.";
			}
			open IN, "<$taxo_profiling_map" or die "ERROR & EXIT: Missing file $taxo_profiling_map";
			while (<IN>) {
				chomp;
				my @line = split "\t";
				$map{ $line[0] } = 1;
			}
			close IN;
		}
		if ( $taxo_profiling_mode eq 'mOTU' ) {
			unless ( $tmp eq '4' ) {
				die "\nERROR & EXIT: Mapping file ($taxo_profiling_map) does not have exactly 4 columns on all lines.";
			}
			open IN, "<$taxo_profiling_map" or die "ERROR & EXIT: Missing file $taxo_profiling_map";
			while (<IN>) {
				chomp;
				my @line = split "\t";
				$map{ $line[0] } = 1;
			}
			close IN;
		}
		print " OK!\n";
	}
	if ( $taxo_profiling_mode eq 'RefMG' ) {
		print localtime() . ": Checking map total length file...";
		unless ( -e "$taxo_profiling_map_tot_len" ) {

			print "\n" . localtime() . ": $taxo_profiling_map_tot_len was missing, create it from the .coord file(s)...\n";
			my %h;
			foreach my $db (@do_taxo_profiling) {
				unless ( -e "$data_dir/$db.coord" ) {
					die "\nERROR & EXIT: Missing .coord file $data_dir/$db.coord";
				}
				open IN,  "<$data_dir/$db.coord";				
				while (<IN>) {
					chomp;
					if (m/^(\d+)\.\S+\t(\d+)\t(\d+)$/) {
						$h{$1} += $3 - $2 + 1;
					}
					else {
						die "ERROR & EXIT: $data_dir/$db.coord has the incorrect format. It should be: <taxid>.XXXXXX<TAB><start><TAB><stop>";
					}
				}
				close IN;
			}
			open OUT, ">$taxo_profiling_map_tot_len" or die "ERROR & EXIT: Cannot write to $taxo_profiling_map_tot_len";
			foreach my $key ( keys %h ) {
				print OUT "$key\t$h{$key}\n";
			}
			close OUT;
			print localtime() . ": $taxo_profiling_map_tot_len was genereated...";
		}
		my @tmp = `awk -F\"\\t\" '{print NF}' $taxo_profiling_map_tot_len | sort -u`;
		foreach my $tmp (@tmp) {
			chomp($tmp);
			unless ( $tmp eq '2' ) {
				die "\nERROR & EXIT: Length file doesn't have exactly 2 columns on all lines.";
			}
		}
		print " OK!\n";
	}
	print localtime() . ": Checking database...";
	unless ( -e "$data_dir/$databases.rownames" ) {
		die "\nERROR & EXIT: Missing database rownames file $data_dir/$databases.rownames";
	}
	print " OK!\n";

	print localtime() . ": Checking taxa id mapping in database...";
	open IN, "<$data_dir/$databases.rownames" or die "ERROR & EXIT: Missing file $data_dir/$databases.rownames";
	<IN>;
	<IN>;
	while (<IN>) {
		chomp;
		my @line = split "\t";
		if ( $taxo_profiling_mode eq 'RefMG' ) {
			@line = split /\./, $line[0];
		}
		unless ( $map{ $line[0] } ) {
			die "\nERROR & EXIT: From $data_dir/$databases.rownames\n$line[0] is not a valid taxa id. \nThis id does not exist in the $taxo_profiling_map mapping file.\nDid you map to a database ($databases) that only has genes, which has identifiers on the form TAXAID.xxx?";
		}
	}
	print " OK!\n";

	print localtime() . ": Initial checks OK, continuing creating jobs...\n";

	### JOB SPECIFIC ###
	foreach my $sample (@samples) {
		
		print JOB "mkdir -p $temp_dir/$sample/temp; ";
		
		my $LOG = " 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.log >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.log ";

		# Check insert and base coverages

		my $insert = "$cwd/$sample/insert.coverage.$databases.$conf{MOCAT_data_type}/$sample.filtered.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.insert.coverage.count";
		my $base   = "$cwd/$sample/base.coverage.$databases.$conf{MOCAT_data_type}/$sample.filtered.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.base.coverage.count";
		unless ( -e $insert ) {
			die "ERROR & EXIT: Missing $insert";
		}
		unless ( -e $base ) {
			die "ERROR & EXIT: Missing $base";
		}

		# Paste rownames
		print JOB "paste $data_dir/$databases.rownames $insert > $insert.tmp && paste $data_dir/$databases.rownames $base > $base.tmp && ";

		if ( $taxo_profiling_mode eq 'RefMG' ) {

			# Actual jobs
			print JOB "echo 'RUNNING INSERT CALCULATIONS' >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.log && ";
			print JOB "$conf{MOCAT_DEV_DIR}/MOCATTaxoProfiling_generateSummary.pl -m $taxo_profiling_mode -c $insert.tmp -t $taxo_profiling_map -l $taxo_profiling_map_tot_len $LOG && ";
			print JOB "echo 'RUNNING BASE CALCULATIONS' >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.log && ";
			print JOB "$conf{MOCAT_DEV_DIR}/MOCATTaxoProfiling_generateSummary.pl -m $taxo_profiling_mode -c $base.tmp -t $taxo_profiling_map -l $taxo_profiling_map_tot_len $LOG && ";

			# Create links to rownames file
			print JOB "ln -fs $taxo_profiling_map.rownames $insert.by.taxonomy.rownames && ln -fs $taxo_profiling_map.rownames $base.by.taxonomy.rownames &&";

		}

		if ( $taxo_profiling_mode eq 'mOTU' ) {

			# Actual jobs
			print JOB "echo 'RUNNING INSERT CALCULATIONS' >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.log && ";
			print JOB "$conf{MOCAT_DEV_DIR}/MOCATTaxoProfiling_generateSummary.pl -m $taxo_profiling_mode -c $insert.tmp -g $taxo_profiling_map $LOG && ";
			print JOB "echo 'RUNNING BASE CALCULATIONS' >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.log && ";
			print JOB "$conf{MOCAT_DEV_DIR}/MOCATTaxoProfiling_generateSummary.pl -m $taxo_profiling_mode -c $base.tmp -g $taxo_profiling_map $LOG && ";

			# Create links to rownames file
			print JOB "ln -fs $taxo_profiling_map.rownames $insert.by.mOTU.rownames && ln -fs $taxo_profiling_map.rownames $base.by.mOTU.rownames &&";

		}

		# Remove tmp files
		print JOB " rm -f $base.tmp $insert.tmp";

		# Push filenames into array, used in second step
		push @taxo_profiling_insert, $insert;
		push @taxo_profiling_base,   $base;

		# End line
		print JOB "\n";

	}
	close JOB;
	print localtime() . ": Jobs created\n";
}

sub summarize {
	
	my $job        = $_[0];
	my $LOG = " 2>> $cwd/logs/$job/samples/MOCATJob_$job.ALL_SAMPLES.log";

	# Variables
	my $databases = join( "_AND_", @do_taxo_profiling );
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}

	chomp( my $host = `hostname` );
	my $groupByColumnExecutable = "$bin_dir/groupByColumn.$host";
	unless ( -e "$groupByColumnExecutable" ) {
		die "ERROR & EXIT: Missing host specific build of groupByColumn executable '$groupByColumnExecutable'. Please run on a different server or make this version.";
	}

	my @mg = ( 'COG0012', 'COG0049', 'COG0052', 'COG0048', 'COG0016', 'COG0018', 'COG0080', 'COG0088', 'COG0081', 'COG0087', 'COG0090', 'COG0085', 'COG0091', 'COG0092', 'COG0093', 'COG0094', 'COG0096', 'COG0097', 'COG0098', 'COG0099', 'COG0100', 'COG0102', 'COG0103', 'COG0124', 'COG0172', 'COG0184', 'COG0185', 'COG0186', 'COG0197', 'COG0200', 'COG0201', 'COG0202', 'COG0215', 'COG0256', 'COG0522', 'COG0495', 'COG0533', 'COG0525', 'COG0552', 'COG0541' );

	my ( @ends, @names, $start );
	if ( $taxo_profiling_mode eq 'mOTU' ) {
		@ends = ( "by.mOTU.norm", "by.mOTU.raw", "by.mOTU.scaled" );
		@names = ( "NA", "NA", "NA", "mOTU" );
		$start = 3;
	}
	elsif ( $taxo_profiling_mode eq 'RefMG' ) {
		@ends = ( "by.taxonomy.norm", "by.taxonomy.raw", "by.taxonomy.scaled" );
		@names = ( "NA", "NA", "kingdom", "phylum", "class", "order", "family", "genus", "species", "specI_clusters" );
		$start = 2;
	}

	# Check files
	print localtime() . ": Check mapping rownames...";
	my $taxo_profiling_map         = "$data_dir/$databases.refmg.map";
	my $taxo_profiling_motu_map    = "$data_dir/$databases.motu.map";
	if ( $taxo_profiling_mode eq 'mOTU' ) {
		$taxo_profiling_map = $taxo_profiling_motu_map;
	}
	elsif ( $taxo_profiling_mode eq 'RefMG' ) {
	}
	else {
		die "\nERROR & EXIT: -mode must be either 'mOTU' or 'RefMG'";
	}
	unless ( -e "$taxo_profiling_map.rownames" ) {
		die "ERROR & EXIT: Missing rownames file of mapping file.\nExpected this file to exist: $taxo_profiling_map.rownames\nPerhaps create it from the $taxo_profiling_map file?\nAdd 2 rows at the top, but make sure the sort order is correct with respect to what is in the data files.\nThis file should have been created in the previous step, though. Somethign wrong?";
	}
	print " OK!\n";
	print localtime() . ": Taxonomy files created. Checking files...";
	foreach my $file ( @taxo_profiling_base, @taxo_profiling_insert ) {
		foreach my $end (@ends) {
			unless ( -e "$file.$end" ) {
				die "ERROR & EXIT: Missing $file.$end";
			}
		}
	}
	print " OK!\n";

	# Create structure
	print localtime() . ": Creating folders...";
	if ( $taxo_profiling_mode eq 'mOTU' ) {
		system "mkdir -p $cwd/taxonomic.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}/COGs";
	}
	elsif ( $taxo_profiling_mode eq 'RefMG' ) {
		system "mkdir -p $cwd/taxonomic.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
	}
	print " OK!\n";

	# Paste
	for my $i ( "insert", "base" ) {
		my @boi;
		if ( $i eq "insert" ) {
			@boi = @taxo_profiling_insert;
		}
		if ( $i eq "base" ) {
			@boi = @taxo_profiling_base;
		}
		for my $end ( "raw", "norm", "scaled" ) {
			my @temp;
			foreach my $t (@boi) {
				if ( $taxo_profiling_mode eq 'mOTU' ) {
					push @temp, "$t.by.mOTU.$end";
				}
				elsif ( $taxo_profiling_mode eq 'RefMG' ) {
					push @temp, "$t.by.taxonomy.$end";
				}
			}
			print localtime() . ": Pasting $i:$end files into all file...";
			chomp( my $basename = `basename $sample_file` );
			my ( $output, $output2 );

			if ( $taxo_profiling_mode eq 'mOTU' ) {
				$output  = "$cwd/taxonomic.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}/$basename.mOTU.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$i.$end";
				$output2 = "$cwd/taxonomic.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}/COGs/$basename.mOTU.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$i.$end";
			}
			elsif ( $taxo_profiling_mode eq 'RefMG' ) {
				$output = "$cwd/taxonomic.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}/$basename.taxonomic.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$i.$end";
			}

			#system "paste $taxo_profiling_map.rownames " . join( " ", @temp ) . " > $output.all";
			#produces: "Can't exec "/bin/sh": Argument list too long"
			#
			#fixed here:
			my %hash;
			my $counter1 = 0;
			my $counter2 = 0;
			for my $k (@temp) {
				$counter1++;
				unless ( $hash{$counter2} ) {
					$hash{$counter2} = "";
				}
				$hash{$counter2} = $hash{$counter2} . " $k ";
				if ( $counter1 == 100 ) {
					$counter1 = 0;
					$counter2++;
				}
			}

			system "cp $taxo_profiling_map.rownames $output.all";
			foreach my $f ( sort { $a <=> $b } keys %hash ) {
				system "paste $output.all $hash{$f} > $output.all.tmp && mv $output.all.tmp $output.all";
			}
			system "rm -f $output.all.tmp";

			print " OK!\n";
			for my $k ( $start .. ( scalar @names - 1 ) ) {
				print localtime() . ": Grouping $i by $names[$k]...";
				if ( $taxo_profiling_mode eq 'mOTU' ) {
					system "$groupByColumnExecutable -g $k -f 4 -h 1 $output.all > $output.$names[$k] $LOG";
				}
				elsif ( $taxo_profiling_mode eq 'RefMG' ) {
					system "$groupByColumnExecutable -g $k -f 10 -h 2 $output.all > $output.$names[$k] $LOG";
				}

				# It doesn't make sense to calculate fractions for mOTUs
				if ( $taxo_profiling_mode ne 'mOTU' ) {
					open IN,  "<$output.$names[$k]"          or die "ERROR & EXIT: Missing $output.$names[$k]";
					open OUT, ">$output.$names[$k].fraction" or die "ERROR & EXIT: Cannot write to $output.$names[$k].fraction";
					my $line = <IN>;
					print OUT $line;
					my @sums;
					while (<IN>) {
						my $line = $_;
						chomp $line;
						my @line = split /\t/, $line;
						for my $col ( 1 .. scalar @line - 1 ) {
							$sums[$col] += $line[$col];
						}
					}
					close IN;
					open IN, "<$output.$names[$k]" or die "ERROR & EXIT: Missing $output.$names[$k]";
					$line = <IN>;
					while (<IN>) {
						my $line = $_;
						chomp $line;
						my @line = split /\t/, $line;
						for my $col ( 0 .. scalar @line - 1 ) {
							if ( $col == 0 ) {
								print OUT $line[$col];
							}
							else {
								print OUT $line[$col] / $sums[$col];
							}
							unless ( scalar @line - 1 == $col ) {
								print OUT "\t";
							}
						}
						print OUT "\n";
					}
					close IN;
					close OUT;
					system "$ZIP -$ziplevel -f $output.$names[$k]";
					system "$ZIP -$ziplevel -f $output.$names[$k].fraction";
				}
				print " OK!\n";

				if ( $taxo_profiling_mode eq 'mOTU' ) {
					print localtime() . ": Grouping $i by COGs...";

					#Old for when there was a file specifying the COGs, right now let's hard code it above instead.
					#open( MG, "<$taxo_profiling_ids" ) or die "ERROR & EXIT: Cannot open $taxo_profiling_ids for input";
					#my @mg = <MG>;
					#chomp @mg;
					#close MG;

					foreach my $i (@mg) {
						open IN,  "<$output.$names[$k]"     or die "ERROR & EXIT: Missing $output.$names[$k]";
						open OUT, ">$output2.$names[$k].$i" or die "ERROR & EXIT: Cannot write to $output2.$names[$k].$i\n";
						while (<IN>) {
							if ( $. == 1 || $_ =~ /^$i\./ ) {
								print OUT $_;
							}
						}
						close IN;
						close OUT;
						system "$ZIP -$ziplevel -f $output2.$names[$k].$i";
					}
					print " OK!\n";
					system "$ZIP -$ziplevel -f $output.$names[$k]";
				}
			}
			system "$ZIP -$ziplevel -f $output.all";
		}
	}
	print localtime() . ": All output files created.\n";
	print localtime() . ": OUTPUT SAVED IN $cwd/taxonomic.profiles/$databases\n";
}

1;
