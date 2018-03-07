package MOCATProfiling::DistMM;
use strict;
use warnings;
use MOCATProfiling::Variables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2014
# This code is released under GNU GPL v3.

# Main threader module, this sub then calls the subs below
sub distributeMultipleMappersThreadWrapper {
	my $thread = shift;
	my $i      = $thread - 1;

	# time to distribute multiple mappers
	foreach my $level (@levels) {
		MOCATProfiling::DistMM::distributeMultipleMappersWrapper( $multipleMapper{$level}, $level, $thread );    # here we distribute the multiple mappers

		if ($VERBOSE) {
			print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] FINISHED distributing MMs for $level [total mapped base length norm=" . $TOTAL_MAPPED_LENGTH_NORM{'base'}{$level} . " | total mapped insert length norm=" . $TOTAL_MAPPED_LENGTH_NORM{'insert'}{$level} . "]\n";
		}

		# we need these lines because to print them below they have to at least have a value :)
		unless ( $TOTAL_MAPPED_LENGTH_NORM{'base.mm.dist.among.unique'}{$level} ) {
			$TOTAL_MAPPED_LENGTH_NORM{'base.mm.dist.among.unique'}{$level} = 0;
		}
		unless ( $TOTAL_MAPPED_LENGTH_NORM{'insert.mm.dist.among.unique'}{$level} ) {
			$TOTAL_MAPPED_LENGTH_NORM{'insert.mm.dist.among.unique'}{$level} = 0;
		}
	}
}

# This wrapper function had to be implemented because of the two ways we can store the data from the multiple mappers
# But as of the current version we don't support non temp files, we die in the else statement below
sub distributeMultipleMappersWrapper {
	my $multipleMapper = shift;
	my $level          = shift;
	my $thread         = shift;

	my $NCBInumber;
	if ( $level eq 'taxaid' ) {
		$NCBInumber = 0;
	}
	if ( $level eq 'kingdom' ) {
		$NCBInumber = 1;
	}
	if ( $level eq 'phylum' ) {
		$NCBInumber = 2;
	}
	if ( $level eq 'class' ) {
		$NCBInumber = 3;
	}
	if ( $level eq 'order' ) {
		$NCBInumber = 4;
	}
	if ( $level eq 'family' ) {
		$NCBInumber = 5;
	}
	if ( $level eq 'genus' ) {
		$NCBInumber = 6;
	}
	if ( $level eq 'species' ) {
		$NCBInumber = 7;
	}
	if ( $level eq 'specI_cluster' ) {
		$NCBInumber = 8;
	}

	if ($temp_file) {
		if ($VERBOSE) {
			print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] STARTING distributing multiple mappers for level=$level using TEMPORARY FILE=$temp_file.MM.$level.$thread.data.gz\n";
		}
		open my $IN, "gunzip -c $temp_file.MM.$level.$thread.data.gz |" or die "ERROR & EXIT: Cannot open IN PIPE gunzip -c $temp_file.MM.$level.$thread.data.gz |";
		while (<$IN>) {
			chomp;
			my @line = split "\t";

			my @insert;
			$insert[0] = $line[0];
			$insert[1] = $line[1];
			$insert[2] = $line[2];
			for my $i ( 3 .. scalar @line - 1 ) {
				my @array = split ",", $line[$i];
				$insert[$i] = \@array;
			}

			#MOCATProfiling::DistMM::distributeMultipleMappers( \@insert, $thread, $level );
			# new in v1.5.2 is that I moved this inside the loop here, which should save time
			#my $temp   = shift;      # we get an array reference that we dereference in two steps
			#my @insert = @{$temp};
			#my $thread = shift;
			#my $level  = shift;

			my $J;
			for my $j ( 0 .. scalar @levels - 1 ) {
				if ( $levels[$j] eq $level ) {
					$J = $j;
				}
			}

			for my $i ( 0 .. 2 ) {

				#this loops over base.1 base.2 and insert in that order
				if ( $insert[$i] > 0 ) {    # we only run it if we actually have something to give

					# we have to change from i which is 0,1 for base.1 base.2 and insert which is 2, to k=0 for i=0,1 and k=1 if i=2
					my $k;
					if ( $i == 0 || $i == 1 ) {
						$k = 0;
					}
					else {
						$k = 1;
					}

					# first we need to get the total sum of the taxa based on the unique abundances we have already allocated
					my @array_gene               = ();
					my @array                    = ();
					my %total_unique_sum_of_taxa = ();    # this is the total som of the taxa based by the unique counts
					my %taxaToGiveTo             = ();    # this hash contains IDs of the taxa we should give to
					my %j;                                # j determines where in the @array we should get the values from, 4=only unique, 2=normal original because no one has unique

					#foreach my $level (@levels) {
					$j{$level} = 4;                       # have to be set to 4 each time we get to a new level
					my %numberOfGeneMatches;

					##########################################################################
					### process gene first, to establish which genes we have actually seen ###
					# test 1, using unique abundances
					foreach my $current ( @{ $insert[ $i + 3 ] } ) {
						@array_gene = @{ $HASH{$current} };
						if ( $array_gene[ $k + 4 ] ) {
							if ( $array_gene[ $k + 4 ] > 0 ) {    # this means it is a unique taxa
								$total_unique_sum_of_taxa{gene}{'raw'} += $array_gene[ $k + 4 ];
								$total_unique_sum_of_taxa{gene}{'norm'} += $array_gene[ $k + 4 ] / ( $array_gene[1] - $array_gene[0] + 1 );
								$taxaToGiveTo{gene}{$current} = ( $array_gene[1] - $array_gene[0] + 1 );    # we are smart and store the length of the taxa here
							}
						}
					}

					# test 2, if 1 failed, we use abundances from raw
					if ( scalar keys %{ $taxaToGiveTo{gene} } == 0 ) {
						%{ $total_unique_sum_of_taxa{gene} } = ();                                                      # this is the total som of the taxa based by the unique counts
						%{ $taxaToGiveTo{gene} }             = ();                                                      # this hash contains IDs of the taxa we should give to
						$j{gene} = 2;
						foreach my $current ( @{ $insert[ $i + 3 ] } ) {
							@array_gene = @{ $HASH{$current} };
							$total_unique_sum_of_taxa{gene}{'raw'} += $array_gene[ $k + 2 ];
							$total_unique_sum_of_taxa{gene}{'norm'} += $array_gene[ $k + 2 ] / ( $array_gene[1] - $array_gene[0] + 1 );
							$taxaToGiveTo{gene}{$current} = ( $array_gene[1] - $array_gene[0] + 1 );              # we are smart and store the length of the taxa here
						}
					}
					##########################################################################

					##########################################################################
					# then process any additional taxa
					if ( $level ne 'gene' ) {

						#test 1
						foreach my $current ( @{ $insert[ $i + 3 ] } ) {
							if ( $taxaToGiveTo{gene}{$current} ) {                                                # only process if the gene has been seen from above
								my $cref;
								if ( $level eq 'mOTU' ) {
									$cref  = @{ $HASH{$current} }[ $INIT_NUM + $J - 1 ];
									@array = @{ $HASH_LEVELS{$level}{$cref} };
								}
								else {
									$current =~ m/^(\d+)\..*/;
									my @l = split "\t", $MAP{$1};
									$cref  = $l[$NCBInumber];
									@array = @{ $HASH_LEVELS{$level}{$cref} };
								}
								@array_gene = @{ $HASH{$current} };
								if ( $array[ $k + 4 ] ) {
									if ( $array[ $k + 4 ] > 0 ) {    # this means it is a unique taxa
										$numberOfGeneMatches{$cref}{$current} = 1;
										unless ( $taxaToGiveTo{$level}{$cref} ) {    # we only want to add the taxa once
											$total_unique_sum_of_taxa{$level}{'raw'}  += $array[ $k + 4 ];
											$total_unique_sum_of_taxa{$level}{'norm'} += $array[ $k + 4 + $INIT_NUM - 2 ];

											#print "$level : $current : total_unique_sum_of_taxa{$level}{'norm'} += array[ $k + 4 + $INIT_NUM - 2 ];" . ( $array[ $k + 4 + $INIT_NUM - 2 ] ) . "\n";
											#$taxaToGiveTo{$current} = $insert[6];    # we are smart and store the length of the taxa here
											#$taxaToGiveTo{$level}{$cref} = @{ $insert[ 5 + $J ] }[$K];    # we are smart and store the length of the taxa here
										}

										# new in v1.5.2 is that we use taxaid length here, which is correct rather than gene length
										if ( $mode eq 'NCBI' ) {
											my @split = split /\./, $current;
											$taxaToGiveTo{$level}{$cref}{$current} = ( $NCBI_TAXAID_LENGTH{ $split[0] } );    # we are smart and store the length of the taxa here # we need to add length for each taxaid, and this has to come after the check just above
										}
										elsif ( $level eq 'mOTU' ) {
											$taxaToGiveTo{$level}{$cref}{$current} = $array_gene[1] - $array_gene[0] + 1;     # we need to add length for each gene, and this has to come after the check just above
										}
										else {
											die "INTERNAL ERROR: Unexpected combination. Contact MOCAT support.";
										}

									}
								}
							}
						}

						# test 2
						if ( scalar keys %{ $taxaToGiveTo{$level} } == 0 ) {
							%{ $total_unique_sum_of_taxa{$level} } = ();    # this is the total som of the taxa based by the unique counts
							%{ $taxaToGiveTo{$level} }             = ();    # this hash contains IDs of the taxa we should give to
							$j{$level} = 2;

							# this is the same as before, but not getting the data from the original fields and not unique
							foreach my $current ( @{ $insert[ $i + 3 ] } ) {
								if ( $taxaToGiveTo{gene}{$current} ) {    # only process if the gene has been seen from above
									my $cref;
									if ( $level eq 'mOTU' ) {
										$cref  = @{ $HASH{$current} }[ $INIT_NUM + $J - 1 ];
										@array = @{ $HASH_LEVELS{$level}{$cref} };
									}
									else {
										$current =~ m/^(\d+)\..*/;
										my @l = split "\t", $MAP{$1};
										$cref  = $l[$NCBInumber];
										@array = @{ $HASH_LEVELS{$level}{$cref} };
									}
									@array_gene = @{ $HASH{$current} };
									$numberOfGeneMatches{$cref}{$current} = 1;

									unless ( $taxaToGiveTo{$level}{$cref} ) {
										$total_unique_sum_of_taxa{$level}{'raw'}  += $array[ $k + 2 ];
										$total_unique_sum_of_taxa{$level}{'norm'} += $array[ $k + 2 + $INIT_NUM - 2 ];

										#$taxaToGiveTo{$current} = $insert[6];    # we are smart and store the length of the taxa here
										#$taxaToGiveTo{$level}{$cref} = @{ $insert[ 5 + $J ] }[$K];    # we are smart and store the length of the taxa here
									}

									# new in v1.5.2 is that we use taxaid length here, which is correct rather than gene length
									if ( $mode eq 'NCBI' ) {
										my @split = split /\./, $current;
										$taxaToGiveTo{$level}{$cref}{$current} = ( $NCBI_TAXAID_LENGTH{ $split[0] } );    # we are smart and store the length of the taxa here # we need to add length for each taxaid, and this has to come after the check just above
									}
									elsif ( $level eq 'mOTU' ) {
										$taxaToGiveTo{$level}{$cref}{$current} = $array_gene[1] - $array_gene[0] + 1;     # we need to add length for each gene, and this has to come after the check just above
									}
									else {
										die "INTERNAL ERROR: Unexpected combination. Contact MOCAT support.";
									}
								}
							}
						}
					}
					##########################################################################

					# Here we do the actual giving
					my $giveToGene;
					foreach my $giveTo ( @{ $insert[ $i + 3 ] } ) {
						my $give = 0;
						if ( $level eq 'gene' ) { $give = 1; }
						if ( $level ne 'gene' ) {
							if ( $taxaToGiveTo{gene}{$giveTo} ) {
								$giveToGene = $giveTo;
								if ( $level eq 'mOTU' ) {
									$giveTo = @{ $HASH{$giveTo} }[ $INIT_NUM + $J - 1 ];
								}
								else {
									$giveTo =~ m/^(\d+)\..*/;
									my @l = split "\t", $MAP{$1};
									$giveTo = $l[$NCBInumber];
								}
								$give = 1;
							}
							else {
								$taxaToGiveTo{$level}{$giveTo} = undef;
							}
						}

						if ( $taxaToGiveTo{$level}{$giveTo} && $give ) {
							my $numberOfGeneMatches = scalar keys %{ $numberOfGeneMatches{$giveTo} };
							if ( $level eq 'gene' ) {
								@{ $HASH{$giveTo} }[ $k + 6 ] += @{ $HASH{$giveTo} }[ $k + $j{$level} ] / $total_unique_sum_of_taxa{$level}{'raw'} * $insert[$i];
								@{ $HASH{$giveTo} }[ $k + 8 ] += ( ( @{ $HASH{$giveTo} }[ $k + $j{$level} ] / $taxaToGiveTo{$level}{$giveTo} ) / $total_unique_sum_of_taxa{$level}{'norm'} ) * ( $insert[$i] / $taxaToGiveTo{$level}{$giveTo} );    # $taxaToGiveTo{$giveTo} = the length of the taxa
							}
							else {

								@{ $HASH_LEVELS{$level}{$giveTo} }[ $k + 6 ] += @{ $HASH_LEVELS{$level}{$giveTo} }[ $k + $j{$level} ] / $total_unique_sum_of_taxa{$level}{'raw'} / $numberOfGeneMatches * $insert[$i];
								unless ( $taxaToGiveTo{$level}{$giveTo}{$giveToGene} ) {
									die "INTERNAL ERROR: missing taxaToGiveTo{$level}{$giveTo}{$giveToGene}";
								}

								@{ $HASH_LEVELS{$level}{$giveTo} }[ $k + 8 ] += @{ $HASH_LEVELS{$level}{$giveTo} }[ $k + $j{$level} + $INIT_NUM - 2 ] / $total_unique_sum_of_taxa{$level}{'norm'} / $numberOfGeneMatches * $insert[$i] / $taxaToGiveTo{$level}{$giveTo}{$giveToGene};
							}

							# first here we start making the total sum for genes and mm dist norm
							# note, there's a diff between gene and non gene, gene uses HASH and none gene uses HASH_LEVELS
							if ( $level eq 'gene' ) {
								if ( $k == 0 ) {
									$TOTAL_MAPPED_LENGTH_NORM{'base.mm.dist.among.unique'}{$level} += ( ( @{ $HASH{$giveTo} }[ $k + $j{$level} ] / $taxaToGiveTo{$level}{$giveTo} ) / $total_unique_sum_of_taxa{$level}{'norm'} ) * ( $insert[$i] / $taxaToGiveTo{$level}{$giveTo} );
								}
								if ( $k == 1 ) {
									$TOTAL_MAPPED_LENGTH_NORM{'insert.mm.dist.among.unique'}{$level} += ( ( @{ $HASH{$giveTo} }[ $k + $j{$level} ] / $taxaToGiveTo{$level}{$giveTo} ) / $total_unique_sum_of_taxa{$level}{'norm'} ) * ( $insert[$i] / $taxaToGiveTo{$level}{$giveTo} );
								}
							}
							else {
								if ( $k == 0 ) {
									$TOTAL_MAPPED_LENGTH_NORM{'base.mm.dist.among.unique'}{$level} += @{ $HASH_LEVELS{$level}{$giveTo} }[ $k + $j{$level} + $INIT_NUM - 2 ] / $total_unique_sum_of_taxa{$level}{'norm'} / $numberOfGeneMatches * $insert[$i] / $taxaToGiveTo{$level}{$giveTo}{$giveToGene};
								}
								if ( $k == 1 ) {
									$TOTAL_MAPPED_LENGTH_NORM{'insert.mm.dist.among.unique'}{$level} += @{ $HASH_LEVELS{$level}{$giveTo} }[ $k + $j{$level} + $INIT_NUM - 2 ] / $total_unique_sum_of_taxa{$level}{'norm'} / $numberOfGeneMatches * $insert[$i] / $taxaToGiveTo{$level}{$giveTo}{$giveToGene};
									#OLD: $TOTAL_MAPPED_LENGTH_NORM{'insert.mm.dist.among.unique'}{$level} += ( ( @{ $HASH_LEVELS{$level}{$giveTo} }[ $k + $j{$level} ] / $taxaToGiveTo{$level}{$giveTo}{$giveToGene} ) / $total_unique_sum_of_taxa{$level}{'norm'} ) / $numberOfGeneMatches * ( $insert[$i] / $taxaToGiveTo{$level}{$giveTo}{$giveToGene} );
								}
							}
						}
					}
				}
			}

		}
		close $IN;
	}

	# THIS HAS TO BE TESTED!
	else {
		die "INTERNAL ERROR, this needs to be updated";

		#		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] STARTING distributing multiple mappers for level=$level using HASH IN MEMORY\n";
		#		foreach my $insert ( keys %{$multipleMapper} ) {    # each key has an array form each insert, in this array there are the values and also arrays of taxa to deal to
		#			my @insert = @{ $multipleMapper->{$insert} };
		#			MOCATProfiling::DistMM::distributeMultipleMappers( \@insert, $thread, $level );
		#		}
	}
}

# This is the sub that does the actual job
#sub distributeMultipleMappers {
#}

1;
