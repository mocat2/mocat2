package MOCATProfiling::Functional;
use strict;
use warnings;
use MOCATProfiling::Variables;
use threads::shared;

sub summarizeFunctionalLevels {

	# for funcitonal stuff, not all genes have a cog, but instead of first saving the maps for each gene in a hash
	# we can just loop over the actual map file, and while doing that we get the data from that gene, this saves us from
	# even storing the map into memory, pretty cool, eh :)
	my @list  = ( '', '', 'base.raw', 'insert.raw', 'base.only.unique.raw', 'insert.only.unique.raw', 'base.mm.dist.among.unique.raw', 'insert.mm.dist.among.unique.raw', 'base.mm.dist.among.unique.norm', 'insert.mm.dist.among.unique.norm', 'base.norm', 'insert.norm', 'base.only.unique.norm', 'insert.only.unique.norm' );    # this will be used further down, note the use of 'raw', which ahs to do with final file names in printFile and the order for the last 6 units
	my @list2 = ( '', '', 'base',     'insert',     'base.only.unique',     'insert.only.unique',     'base.mm.dist.among.unique',     'insert.mm.dist.among.unique',     'base.mm.dist.among.unique.norm', 'insert.mm.dist.among.unique.norm', 'base.norm', 'insert.norm', 'base.only.unique.norm', 'insert.only.unique.norm' );    # this will be used further down, note the use of 'raw', which ahs to do with final file names in printFile and the order for the last 6 units
	my $MAP;
	my $index;                                                                                                                                                                                                                                                                                                                       # each gene gets stored with a gene.index = gene.0, gene.1, ...
	my $modified_name_without_raw;

	if ( -e "$map_file.gz" ) {
		open $MAP, "gunzip -c $map_file.gz | " or die;
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : summarizing functional levels based on maps in $map_file.gz\n";
	}
	else {
		open $MAP, "<", "$map_file" or die;
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : summarizing functional levels based on maps in $map_file\n";
	}
	chomp( my $header = <$MAP> );
	my @header = split /\s+/, $header;
	my @ORDER;
	foreach my $lvl (@LEVELS) {
		foreach my $h ( 0 .. scalar @header - 1 ) {

			#print "$lvl = $header[$h]?\n";
			if ( $header[$h] eq $lvl ) {
				push @ORDER, $h;
			}
		}
	}

	#print join( " ", @ORDER ) . " : " . join( " ", @header ) . " : " . join( " ", @LEVELS ) . " X\n";

	while (my $line = <$MAP>) {
		chomp;
		my @line = split /\t/, $line;

		$index = 0;
		while ( exists $HASH{"$line[0].$index"} ) {
			my @array  = @{ $HASH{"$line[0].$index"} };    # get the one line array, map data is stored in INIT_NUM and to the end
			my $length = $array[1] - $array[0] + 1;        # this is the length of the current gene
			for my $i (@ORDER) {                           # here we won't run over the map file in order, but rather in the order of the @levels, it's a bit tricky
				my $level = $header[$i];
				unless ( $line[$i] ) {
					$line[$i] = '';
				}

				#print join(" ", @line) . " - $level - $i - $line[$i]\n";

				my @cogs = split /[,|]/, $line[$i];
				my $cog_numbers = scalar @cogs;
				foreach my $cog (@cogs) {
					if ( !($cog eq '' || $cog eq 'NA') ) {
						@{ $HASH_LEVELS{$level}{$cog} }[1] += $length;    # saves the length in the uniq rownames files
						for my $j ( 2 .. $INIT_NUM - 1 ) {                # note that later we have to get the norm data, and this we can get by storing the concatenated length in position 1
							if ( $array[$j] ) {                     # this prevents us from storing 0s
								@{ $HASH_LEVELS{$level}{$cog} }[$j] += $array[$j];    # this does the actual distr to the new levels, we do it by sum, another method would've been by mean or something but we always use sum

								# here it again becomes tricky, for the HASH storing genes, it was fine not to store the gene length norm values, but here we have to keep them
								# but in position 8,9 we already have them for base and insert norm mm dist
								if ( $j < $INIT_NUM - 2 ) {
									@{ $HASH_LEVELS{$level}{$cog} }[ $j + $INIT_NUM - 2 ] += $array[$j] / $length;
								}

								# this appreas four times, 2 here (*NOTE number 2* and 2 just further down when going through the genes that weren't annotated)
								# for genes and non funcitonal entities these tables are filled out in the main processing while going line by line
								# but here we have to populate them for the funcitonal levels
								# raw base and read counts are stored in array position 2,3
								# NOTE that just here I used j for counting, but the other three has t for counting
								# *NOTE number 2*: I've uncommented them below because I don't think what wasnt annotated should be included to calculate the scaling factors
								$TOTAL_MAPPED{ $list2[$j] }{$level}             += $array[$j];
								$TOTAL_MAPPED_LENGTH_NORM{ $list2[$j] }{$level} += $array[$j] / $length;
								$TOTAL_MAPPED_ANNOTATED{$level}{ $list[$j] }    += $array[$j];

								# the values that go into only.unique, also have ot go into the mm.dist
								# but it doesn't have to be done for the only unique norm
								if ( $j == 4 || $j == 5 || $j == 12 || $j == 13 ) {
									$TOTAL_MAPPED_ANNOTATED{$level}{ $list[ $j + 2 ] } += $array[$j];
								}
							}
						}

						for my $t ( 10 .. 13 ) {    # 2-7 are normal non normalized, 8-9 is mm dist norm, 10-13 is norm and only unique norm
							if ( $array[ $t - $INIT_NUM + 2 ] ) {
								$TOTAL_MAPPED_ANNOTATED{$level}{ $list[$t] } += $array[ $t - $INIT_NUM + 2 ] / $length;    # note that we get thevalues from elsewhere in the array

								# for genes and non funcitonal entities these tables are filled out in the main processing while going line by line
								# but here we have to populate them for the funcitonal levels
								# raw base and read counts are stored in array position 2,3
								$TOTAL_MAPPED_LENGTH_NORM{$level}{ $list2[$t] } += $array[ $t - $INIT_NUM + 2 ];           # these values are already length normalized
							}
						}

						# Here we add the horizontal coverage, note that we have to add up the percentages and then at some point divide by the number of genes in each COG
						if ($CALCULATE_HORIZONTAL_COVERAGE) {
							unless ( @{ $HASH_LEVELS{$level}{$cog} }[ $INIT_NUM + 6 ] ) {
								@{ $HASH_LEVELS{$level}{$cog} }[ $INIT_NUM + 6 ] = 0;
								@{ $HASH_LEVELS{$level}{$cog} }[ $INIT_NUM + 7 ] = 0;
							}
							if ( $array[ $INIT_NUM + 1 ] ) {
								@{ $HASH_LEVELS{$level}{$cog} }[ $INIT_NUM + 6 ] += ( $array[ $INIT_NUM + 1 ] =~ tr/1// );
							}
							@{ $HASH_LEVELS{$level}{$cog} }[ $INIT_NUM + 7 ]++;

							#print scalar $HASH_LEVELS{$level}{$cog} . " : HASH_LEVELS{$level}{$cog} : " . join (" - ", @{$HASH_LEVELS{$level}{$cog}}) . "\n";
						}
					}
				}

				# we need to add to the total, for stats later on
				for my $t ( 2 .. 7 ) {    # 2-7 are normal non normalized, 8-9 is mm dist norm, 10-13 is norm and only unique norm
					if ( $array[$t] ) {
						if ( $cog_numbers == 0 ) {

							# here, if we didn't see a cog/ko/whatever for this particular gene we add it's value to the not annotated sum
							# these values are level dependent as well
							$TOTAL_MAPPED_NOT_ANNOTATED{$level}{ $list[$t] } += $array[$t];

							# the values that go into only.unique, also have ot go into the mm.dist
							# but it doesn't have to be done for the only unique norm
							#if ( $list[$t] =~ m/only.unique/ ) {
							if ( $t == 4 || $t == 5 || $t == 12 || $t == 13 ) {
								$TOTAL_MAPPED_NOT_ANNOTATED{$level}{ $list[ $t + 2 ] } += $array[$t];
							}
						}
					}
				}
				for my $t ( 10 .. 13 ) {    # 2-7 are normal non normalized, 8-9 is mm dist norm, 10-13 is norm and only unique norm
					if ( $array[ $t - $INIT_NUM + 2 ] ) {
						if ( $cog_numbers == 0 ) {
							$TOTAL_MAPPED_NOT_ANNOTATED{$level}{ $list[$t] } += $array[ $t - $INIT_NUM + 2 ] / $length;    # note that we get thevalues from elsewhere in the array
						}
					}
				}
			}

			# here we delete the entry in the HASH, so that we can run through the ones we didn't delete later down to get the stats
			delete $HASH{"$line[0].$index"};
			$index++;
		}
	}
	close $MAP or die "ERROR & EXIT: Closing pipe failed, this means something is wrong. Corrput files?";
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : summarizing genes that were not in map file\n";

	# we have to loop over the remaining objects in the hash to get the mapped but not annotated stats
	foreach my $key ( keys %HASH ) {

		# we nned to add to the total, for stats later on
		my @array  = @{ $HASH{$key} };
		my $length = $array[1] - $array[0] + 1;
		for my $t ( 2 .. 7 ) {    # 2-7 are normal non normalized, 8-9 is mm dist norm, 10-13 is norm and only unique norm
			if ( $array[$t] ) {

				#print "saving to not annotated $key : $list[$t]=$array[$t]\n";    #TEXT
				foreach my $level (@levels) {
					$TOTAL_MAPPED_NOT_ANNOTATED{$level}{ $list[$t] } += $array[$t];

					# the values that go into only.unique, also have ot go into the mm.dist
					# but it doesn't have to be done for the only unique norm
					#if ( $list[$t] =~ m/only.unique/ ) {
					if ( $t == 4 || $t == 5 || $t == 12 || $t == 13 ) {
						$TOTAL_MAPPED_NOT_ANNOTATED{$level}{ $list[ $t + 2 ] } += $array[$t];
					}

					# *NOTE number 2*: I've uncommented this below because I don't think what wasnt annotated should be included to calculate the scaling factors
					# This is a difference between version 1.4.0, compared to version 1.3 and below
					#my $modified_name_without_raw = $list[$t];
					#$modified_name_without_raw =~ s/.raw//;
					#$TOTAL_MAPPED{$modified_name_without_raw}{$level} += $array[$t];
					#$TOTAL_MAPPED_LENGTH_NORM{$modified_name_without_raw}{$level} += $array[$t] / $length;
				}
			}
		}

		for my $t ( 10 .. 13 ) {    # 2-7 are normal non normalized, 8-9 is mm dist norm, 10-13 is norm and only unique norm
			foreach my $level (@levels) {
				if ( $array[ $t - $INIT_NUM + 2 ] ) {
					$TOTAL_MAPPED_NOT_ANNOTATED{$level}{ $list[$t] } += $array[ $t - $INIT_NUM + 2 ] / $length;    # note that we get thevalues from elsewhere in the array
				}

				# *NOTE number 2*: I've uncommented this below because I don't think what wasnt annotated should be included to calculate the scaling factors
				# This is a difference between version 1.4.0, compared to version 1.3 and below#my $modified_name_without_raw = $list[$t];
				#$modified_name_without_raw =~ s/.raw//;
				#$TOTAL_MAPPED_LENGTH_NORM{$modified_name_without_raw}{$level} += $array[ $t - $INIT_NUM + 2 ];    # these values are already length normalized
			}
		}

		# empty hash
		delete $HASH{"$key"};
	}

	# this is a dirty fix, but this is the last one that hasn't been added through out the loops above
	# for the others these are summarized at run time just above
	foreach my $level (@levels) {
		if ( $TOTAL_MAPPED_NOT_ANNOTATED{$level}{'base.only.unique.norm'} ) {
			$TOTAL_MAPPED_NOT_ANNOTATED{$level}{'base.mm.dist.among.unique.norm'} += $TOTAL_MAPPED_NOT_ANNOTATED{$level}{'base.only.unique.norm'};
		}
		if ( $TOTAL_MAPPED_NOT_ANNOTATED{$level}{'insert.only.unique.norm'} ) {
			$TOTAL_MAPPED_NOT_ANNOTATED{$level}{'insert.mm.dist.among.unique.norm'} += $TOTAL_MAPPED_NOT_ANNOTATED{$level}{'insert.only.unique.norm'};
		}
		if ( $TOTAL_MAPPED_ANNOTATED{$level}{'base.only.unique.norm'} ) {
			$TOTAL_MAPPED_ANNOTATED{$level}{'base.mm.dist.among.unique.norm'} += $TOTAL_MAPPED_ANNOTATED{$level}{'base.only.unique.norm'};
		}
		if ( $TOTAL_MAPPED_ANNOTATED{$level}{'insert.only.unique.norm'} ) {
			$TOTAL_MAPPED_ANNOTATED{$level}{'insert.mm.dist.among.unique.norm'} += $TOTAL_MAPPED_ANNOTATED{$level}{'insert.only.unique.norm'};
		}
	}

}

1;
