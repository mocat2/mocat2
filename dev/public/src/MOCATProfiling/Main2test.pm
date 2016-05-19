package MOCATProfiling::Main2test;
use strict;
use warnings;
use File::Slurp;
use MOCATProfiling::Variables;
use File::Sync qw(fsync sync);
use Storable;

#use Data::Dumper;

# this is the routine started for each thred
sub main {

	my $thread = threads->tid();
	my $i      = $thread - 1;
	my %NCBI;

	# Load map
	#MOCATProfiling::Load::loadCoordFile($thread);    # we just set -1 because of semantic for not printing text THREAD X later on, but also because we want to loop over all entries, see Loader script.

	# Main loop
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] Collecting inserts from MAIN thread\n";

#	while ( !$TERMINATE && defined( my $file = $request_q->dequeue() ) ) {
#
#		# not needed anymore because we are back to loading one map per process...
#		#		while ( $LOADED_MAP == 0 ) {    # wait until the maps have been loaded
#		#			sleep 1;
#		#		}
#
#		sync();
#		open my $fileIn, "<$file" or die "ERROR & EXIT: Missing $file";
#		fsync($fileIn);
#		close $fileIn;
#		open $fileIn, "<$file.hash" or die "ERROR & EXIT: Missing $file.hash";
#		fsync($fileIn);
#		close $fileIn;
#
#		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] RECEIVED JOB ($file)\n";
#		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] LOADING $file.hash(levels)\n";
#		%HASH = %{ retrieve("$file.hash") };
#		if ( $mode eq 'mOTU' || $mode eq 'NCBI' ) {
#			%HASH_LEVELS = %{ retrieve("$file.hashlevels") };
#		}
#		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] LOADED HASH, PROESSING $file\n";
#
#		open $fileIn, "gunzip -c $file |" or die "ERROR & EXIT: Cannot read from INPUT PIPE gunzip -c $file |";
#
#		#open $fileIn, "<$file" or die "ERROR & EXIT: Cannot read from INPUT PIPE gunzip -c $file |";

		# we reset the local variables for each new files
		$prev_insert = "0";
		$prev_read   = "0";

		#my @lines = read_file($file);
		my $line;
		my $i_ref;
		my $countThisTaxa;
		my $j;
		my %miniNCBImap;
		my %preCount;

		#				if ($temp_file) {    # if we use a temp file, clear the multipleMapper hash
		%multipleMapper = ();

		#				}

		#while ( defined( $line = pop @lines ) ) {

		while (<$READ>) {
			
			

			##### READLINE #####
			##### READLINE #####
			##### READLINE #####
			##### READLINE #####
			#chomp($line);    # not used
			chomp(my $line = $_);
			my @line = split "\t";
			$read = $line[0];
			if ( $input_file_format eq 'BAM' || $input_file_format eq 'SAM' ) {
				$line[0] =~ m/(.+)\/([12])$/;
				$insert     = $1;
				$direction  = $2;
				$ref_id     = $line[2];
				$first_base = $line[3];
				$length     = length $line[9];
				$last_base  = $first_base + $length - 1;    # read will be revcomp if direction is "-", so $first_base is always left-most base
				if ( $mode eq 'NCBI' ) {
					$ref_id =~ m/^(\d+)\..*/;
					$line .= "\t" . $MAP{$1} . "\n";
					@line = split "\t", $line;
					$NCBI{taxaid}        = $line[13];
					$NCBI{kingdom}       = $line[14];
					$NCBI{phylum}        = $line[15];
					$NCBI{class}         = $line[16];
					$NCBI{order}         = $line[17];
					$NCBI{family}        = $line[18];
					$NCBI{genus}         = $line[19];
					$NCBI{species}       = $line[20];
					$NCBI{specI_cluster} = $line[21];
					$NCBI{"length"}        = $line[22];
				}
			}
			elsif ( $input_file_format eq 'SOAP' ) {
				$line[0] =~ m/(.+)\/([12])$/;
				$insert     = $1;
				$direction  = $2;
				$ref_id     = $line[7];
				$first_base = $line[8];
				$length     = $line[5];
				$last_base  = $first_base + $length - 1;    # read will be revcomp if direction is "-", so $first_base is always left-most base
				if ( $mode eq 'NCBI' ) {
					$NCBI{taxaid}        = $line[12];
					$NCBI{kingdom}       = $line[13];
					$NCBI{phylum}        = $line[14];
					$NCBI{class}         = $line[15];
					$NCBI{order}         = $line[16];
					$NCBI{family}        = $line[17];
					$NCBI{genus}         = $line[18];
					$NCBI{species}       = $line[19];
					$NCBI{specI_cluster} = $line[20];
					$NCBI{"length"}        = $line[21];
				}
			}
			##### READLINE #####
			##### READLINE #####
			##### READLINE #####
			##### READLINE #####

			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			# If we see a new insert, do the actual processing
			if ( $insert ne $prev_insert && $prev_insert ne "" ) {

				##### PROCESS INSERT ##### COPY FROM HERE
				##### PROCESS INSERT ##### COPY FROM HERE
				##### PROCESS INSERT ##### COPY FROM HERE
				##### THIS LARGE SECTION APPEARS TWICE,
				##### ONCE FOR ALL LINES EXCEPT THE LAST,
				##### AND THEN THE LAST OUTSIDE THE LOOP

				#############################################################################################################################################
				### Paired end filtering --> ###
				#############################################################################################################################################
				my %intersection = ();
				my %actual_gene;
				my ( $gene_id1, $gene_id2 );
				my %preCount2;    # preCount is a bit tricky, it has to do with when we get only a aubset of taxa, and we need to then decrease the number of genes that are used for counting, and this value comes from adding values in a hash, these values are calculated below in the inner loop
				if ( $PE_filter eq 'yes' ) {
					my ( $b1, $b2, %i );    # %i=intermediate hashes to save some look up time
					foreach my $i (@levels) {
						$b1 = scalar keys %{ $seen_types{'base.1'}{$i} };
						$b2 = scalar keys %{ $seen_types{'base.2'}{$i} };
						if ( $b1 > 0 && $b2 > 0 ) {
							my %intersect;
							if ( $i eq 'gene' ) {
								while ( my $k = each %{ $seen_types{'base.1'}{$i} } ) {
									if ( exists $seen_types{'base.2'}{$i}{$k} ) {
										$intersect{$k} = 1;
									}
								}
							}
							else {
								while ( my $k = each %{ $seen_types{'base.1'}{$i} } ) {
									$gene_id1 = $seen_types{'base.1'}{$i}{$k};
									$gene_id2 = $seen_types{'base.2'}{$i}{$k};
									if ( $gene_id2 && exists $i{$gene_id1} && exists $i{$gene_id2} ) {
										$intersect{$k} = $seen_types{'base.2'}{$i}{$k};    # the value of the hash is the real gene name, this is needed below for the multipleMapper hash, which needs the gene length
										$preCount2{insert}{$i}   += scalar keys %{ $preCount{insert}{$i}{$k} };
										$preCount2{'base.2'}{$i} += scalar keys %{ $preCount{'base.2'}{$i}{$k} };
										$preCount2{'base.1'}{$i} += scalar keys %{ $preCount{'base.1'}{$i}{$k} };
									}
								}
							}

							%{ $intersection{'insert'}{$i} } = %intersect;                                   # the insert intersection is the intersect of the two bases

							# but the bases we have to check the intersection of each base with the insert
							# if the insert intersection is 0, we revert to using the original seen_types hashes
							if ( scalar keys %intersect == 0 ) {
								%{ $intersection{'base.1'}{$i} } = %{ $seen_types{'base.1'}{$i} };
								%{ $intersection{'base.2'}{$i} } = %{ $seen_types{'base.2'}{$i} };
								%{ $intersection{'insert'}{$i} } = %{ $seen_types{'insert'}{$i} };
								$preCount2{'base.1'}{$i} = 0;
								$preCount2{'base.2'}{$i} = 0;
								$preCount2{'insert'}{$i} = 0;
							}
							else {
								%{ $intersection{'base.1'}{$i} } = %{ $intersection{'insert'}{$i} };
								%{ $intersection{'base.2'}{$i} } = %{ $intersection{'insert'}{$i} };
							}
							if ( $i eq 'gene' ) {
								foreach my $H ( \%{ $seen_types{'base.2'}{gene} }, \%{ $seen_types{'base.1'}{gene} } ) {
									while ( my $k = each %$H ) {
										$i{$k} = 1;
									}
								}
							}
						}
						else {
							%intersection = %seen_types;
						}
					}
				}
				else {
					%intersection = %seen_types;
				}
				### <-- Paired end filtering ###
				#############################################################################################################################################

				#############################################################################################################################################
				# THIS IS PROBABLY THE HEART OF EVERYTHING. HERE WE GIVES ABUNDANCES TO EACH GENE, MOTU, COG, ...
				#############################################################################################################################################
				$INSERT_COUNTER++;    # this is just a cunter from creating new entries in the multipleMapper hash
				my %count;            # this hash saves the count from gene iteration, which needs ot be used for the other taxa

				# These were originally inside while ref id loop. but moved out to save time
				my ( $taxaToGiveTo, $count, $position, $length );

				# LOOP OVER EACH TYPE (base.1, base.2, insert) IN THE HASH, WHICH CONTAINS EACH GENE FOR EACH TYPE
				while ( my $type = each %hash ) {
					foreach my $i (@levels) {
						my %tempHash = %{ $hash{$type}{gene} };
						while ( my $ref_id_index = each %tempHash ) {    # hash{gene} contains gene IDs, the other hash{level} are empty

							#############################################################################################################################################
							# GET VALUE AND ID
							#############################################################################################################################################
							my $value = $tempHash{$ref_id_index};
							if ( $i eq 'gene' ) {
								$i_ref        = "$ref_id_index";
								$length       = @{ $HASH{$ref_id_index} }[1] - @{ $HASH{$ref_id_index} }[0] + 1;    # the stop-start+1
								$count        = scalar keys %{ $intersection{$type}{$i} };
								$count{$type} = $count;                                                             # we have to save the count for the genes for using with the other taxa in the next iteration
								$taxaToGiveTo = $count;

							}
							elsif ( $i eq 'mOTU' ) {
								$length       = @{ $HASH{$ref_id_index} }[1] - @{ $HASH{$ref_id_index} }[0] + 1;    # the stop-start+1
								$i_ref        = @{ $HASH{"$ref_id_index"} }[$INIT_NUM];
								$taxaToGiveTo = scalar keys %{ $intersection{$type}{$i} };
								if ( $preCount2{$type}{$i} ) {
									$count = $preCount2{$type}{$i};
								}
								else {
									$count = $count{$type};
								}

							}
							else {
								$length       = $miniNCBImap{"length"}{$ref_id_index};
								$i_ref        = $miniNCBImap{$i}{$ref_id_index};
								$taxaToGiveTo = scalar keys %{ $intersection{$type}{$i} };
								if ( $preCount2{$type}{$i} ) {
									$count = $preCount2{$type}{$i};
								}
								else {
									$count = $count{$type};
								}
							}
							#############################################################################################################################################

							#############################################################################################################################################
							# HERE IS THE ACTUAL GIVING
							#############################################################################################################################################
							if ( exists $intersection{$type}{$i}{$i_ref} ) {    #COUNT THIS TAXA

								#############################################################################################################################################
								# ADD TO THE TOTAL AND ALSO TO THE HASH/HASH_LEVELS
								#############################################################################################################################################
								# so these lines just below appear three times within the next 50ish lines
								if ( $type eq 'insert' ) {
									$TOTAL_MAPPED_LENGTH_NORM{'insert'}{$i} += ( 1 / $count ) * $value / $length;    # this value is used further down when calculating the scaled values
									$TOTAL_MAPPED{'insert'}{$i} += 1 / $count * $value;
									$position = 3;
								}
								else {
									$TOTAL_MAPPED_LENGTH_NORM{'base'}{$i} += ( 1 / $count ) * $value / $length;      # this value is used further down when calculating the scaled values
									$TOTAL_MAPPED{'base'}{$i} += ( 1 / $count ) * $value;
									$position = 2;
								}

								if ( $i eq 'gene' ) {
									@{ $HASH_local{$i_ref} }[$position] += 1 / $count * $value;
								}
								else {
									@{ $HASH_LEVELS_local{$i}{$i_ref} }[$position] += 1 / $count * $value;

									# we also have to add the gene length norm, for functional, these are added in the Functional subroutine,
									# but for mOTU these have to be done on runtime here
									@{ $HASH_LEVELS_local{$i}{$i_ref} }[ $position + $INIT_NUM - 2 ] += 1 / $count * $value / $length;
								}
								#############################################################################################################################################

								#############################################################################################################################################
								# IF THE GENE IS A SINGLE MAPPER ON THIS TAXONOMIC LEVEL
								#############################################################################################################################################
								if ( $taxaToGiveTo == 1 ) {

									#print "UNIQUE === $count | $taxaToGiveTo | LEVEL=$i || $prev_insert || $i_ref += 1 / $count * $value %" . ( 1 / $count * $value ) . "\n";

									if ( $type eq 'insert' ) {
										$TOTAL_MAPPED_LENGTH_NORM{'insert.only.unique'}{$i} += ( 1 / $count ) * $value / $length;    # this value is used further down when calculating the scaled values
										$TOTAL_MAPPED{'insert.only.unique'}{$i} += ( 1 / $count ) * $value;
									}
									else {
										$TOTAL_MAPPED_LENGTH_NORM{'base.only.unique'}{$i} += ( 1 / $count ) * $value / $length;      # this value is used further down when calculating the scaled values
										$TOTAL_MAPPED{'base.only.unique'}{$i} += ( 1 / $count ) * $value;
									}

									if ( $i eq 'gene' ) {
										@{ $HASH_local{$i_ref} }[ $position + 2 ] += 1 / $count * $value;
									}
									else {
										@{ $HASH_LEVELS_local{$i}{$i_ref} }[ $position + 2 ] += 1 / $count * $value;

										# we also have to add the gene length norm, for functional, these are added in the Functional subroutine, but for mOTU & NCBI these have to be done on runtime here
										@{ $HASH_LEVELS_local{$i}{$i_ref} }[ $position + $INIT_NUM - 2 + 2 ] += 1 / $count * $value / $length;

									}
								}
								#############################################################################################################################################

								#############################################################################################################################################
								# IF THE GENE IS A MULTIPLE MAPPER ON THIS TAXONOMIC LEVEL
								#############################################################################################################################################
								else {
									my $position2;
									if ( $type eq 'base.1' ) {
										$position2 = 0;
									}
									if ( $type eq 'base.2' ) {
										$position2 = 1;
									}
									if ( $type eq 'insert' ) {
										$position2 = 2;
									}

									# here we save the inserts in a hash (could it be done in an array), 0->2 has the values for base.1 base.2 and insert, 3->5 has the arrays of the taxa to dist. to
									my ( @array, @length, $ref_hash );
									if ( scalar keys %{ $intersection{$type}{$i} } == 0 ) {    # this means that the intersection hash is empty and we need to get all from the seen hash
										$ref_hash = \%seen_types;
									}
									else {
										$ref_hash = \%intersection;
									}
									my %tempHash = %{ $ref_hash->{$type}{$i} };
									if ( $i eq 'gene' ) {
										while ( my $key = each %tempHash ) {
											push @length, @{ $HASH{$key} }[1] - @{ $HASH{$key} }[0] + 1;
											push @array,  $key;
										}

									}
									else {
										while ( my $key = each %tempHash ) {
											push @array,  $tempHash{$key};
											push @length, @{ $HASH{ $tempHash{$key} } }[1] - @{ $HASH{ $tempHash{$key} } }[0] + 1;
										}
									}

									# if we have specified a temp file we save the MM datat in that and read it line by line later, rather than
									# saving the data in a larger and larger array, this is done further doen at the end of each insert
									unless ( $multipleMapper{$i}{"$INSERT_COUNTER"} ) {
										@{ $multipleMapper{$i}{"$INSERT_COUNTER"} } = ( 0, 0, 0 );
									}
									@{ $multipleMapper{$i}{"$INSERT_COUNTER"} }[$position2]       = $value;     # here the annotion is a bit tricky
									@{ $multipleMapper{$i}{"$INSERT_COUNTER"} }[ $position2 + 3 ] = \@array;    # here the annotion is a bit tricky
									@{ $multipleMapper{$i}{"$INSERT_COUNTER"} }[6]                = \@length;

									if ( $type eq 'insert' ) {

										# because mm dist norm comes form multiple sources it's a bit tricky and have to be summarized later, when distr. MMs
										# so here we only summarize the TOTAL_MAPPED, and not the TOTAL_MAPPED_LENGTH_NORM
										$TOTAL_MAPPED{'insert.mm.dist.among.unique'}{$i} += ( 1 / $count ) * $value;
									}
									else {

										# because mm dist norm comes form multiple sources it's a bit tricky and have to be summarized later, when distr. MMs
										# so here we only summarize the TOTAL_MAPPED, and not the TOTAL_MAPPED_LENGTH_NORM
										$TOTAL_MAPPED{'base.mm.dist.among.unique'}{$i} += ( 1 / $count ) * $value;
									}
								}
							}
							else {

								# the read matches a gene that isn't in the coord file
							}
						}
					}
				}    # end type

				%hash        = ();
				%seen_types  = ();
				%miniNCBImap = ();
				%preCount    = ();
				##### PROCESS INSERT ##### COPY TO HERE
				##### PROCESS INSERT ##### COPY TO HERE
				##### PROCESS INSERT ##### COPY TO HERE

			}
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####
			##### PROCESS INSERT #####

			#############################################################################################################################################
			# MAIN LOOP THAT LOOPS OVER EACH LINE
			#############################################################################################################################################
			##### INNER LOOP #####
			##### INNER LOOP #####
			##### INNER LOOP #####
			##### INNER LOOP #####
			if ( exists $HASH{"$ref_id.0"} ) {
				$reference_for_read_exists_in_coord_file++;
				for ( my $index = 0 ; exists $HASH{"$ref_id.$index"} ; ++$index ) {
					my @bounds = @{ $HASH{"$ref_id.$index"} }[ 0 .. 1 ];
					if ( !( $last_base < $bounds[0] || $first_base > $bounds[1] ) ) {
						my $subtract = 0;
						if ( $first_base < $bounds[0] ) {
							$subtract = $bounds[0] - $first_base;
						}
						if ( $last_base > $bounds[1] ) {
							$subtract = $subtract + $last_base - $bounds[1];
						}
						$i_ref                                           = "$ref_id.$index";
						$hash{'insert'}{gene}{"$ref_id.$index"}          = 1;
						$hash{"base.$direction"}{gene}{"$ref_id.$index"} = $length;
						foreach my $level (@levels) {

							### GET ID ###
							### GET ID ###
							### GET ID ###
							if ( $level eq 'mOTU' ) {
								$i_ref                                                         = @{ $HASH{"$ref_id.$index"} }[$INIT_NUM];
								$preCount{"base.$direction"}{$level}{$i_ref}{"$ref_id.$index"} = 1;
								$preCount{"insert"}{$level}{$i_ref}{"$ref_id.$index"}          = 1;
							}
							elsif ( $level ne 'gene' ) {
								$i_ref                                                         = $NCBI{$level};
								$miniNCBImap{$level}{"$ref_id.$index"}                         = $i_ref;
								$miniNCBImap{"length"}{"$ref_id.$index"}                       = $NCBI{"length"};
								$preCount{"base.$direction"}{$level}{$i_ref}{"$ref_id.$index"} = 1;
								$preCount{"insert"}{$level}{$i_ref}{"$ref_id.$index"}          = 1;
							}
							### GET ID ###
							### GET ID ###
							### GET ID ###

							$seen_types{"base.$direction"}{$level}{$i_ref} = "$ref_id.$index";
							$seen_types{'insert'}{$level}{$i_ref}          = "$ref_id.$index";
						}

						$seen_this_insert = 1;
						$TOTAL_LENGTH_MAPPED += $length;    # REMEMBER TO ADD THIS WHEN ADDING BETTER BASE SUPPORT# - $subtract ;
					}
					$TOTAL_LENGTH += $length;
				}
			}
			else {
				$reference_for_read_do_not_exists_in_coord_file++;
			}
			##### INNER LOOP #####
			##### INNER LOOP #####
			##### INNER LOOP #####
			##### INNER LOOP #####

			# After reading a line a processing, if needed, set this insert to prev_insert
			$prev_insert = $insert;
			#############################################################################################################################################

		}

		#############################################################################################################################################
		# End of main IN loop. Here we run the process insert one last time
		#############################################################################################################################################

		##### PROCESS INSERT ##### COPY FROM HERE
		##### PROCESS INSERT ##### COPY FROM HERE
		##### PROCESS INSERT ##### COPY FROM HERE
		##### THIS LARGE SECTION APPEARS TWICE,
		##### ONCE FOR ALL LINES EXCEPT THE LAST,
		##### AND THEN THE LAST OUTSIDE THE LOOP

		#############################################################################################################################################
		### Paired end filtering --> ###
		#############################################################################################################################################
		my %intersection = ();
		my %actual_gene;
		my ( $gene_id1, $gene_id2 );
		my %preCount2;    # preCount is a bit tricky, it has to do with when we get only a aubset of taxa, and we need to then decrease the number of genes that are used for counting, and this value comes from adding values in a hash, these values are calculated below in the inner loop
		if ( $PE_filter eq 'yes' ) {
			my ( $b1, $b2, %i );    # %i=intermediate hashes to save some look up time
			foreach my $i (@levels) {
				$b1 = scalar keys %{ $seen_types{'base.1'}{$i} };
				$b2 = scalar keys %{ $seen_types{'base.2'}{$i} };
				if ( $b1 > 0 && $b2 > 0 ) {
					my %intersect;
					if ( $i eq 'gene' ) {
						while ( my $k = each %{ $seen_types{'base.1'}{$i} } ) {
							if ( exists $seen_types{'base.2'}{$i}{$k} ) {
								$intersect{$k} = 1;
							}
						}
					}
					else {
						while ( my $k = each %{ $seen_types{'base.1'}{$i} } ) {
							$gene_id1 = $seen_types{'base.1'}{$i}{$k};
							$gene_id2 = $seen_types{'base.2'}{$i}{$k};
							if ( $gene_id2 && exists $i{$gene_id1} && exists $i{$gene_id2} ) {
								$intersect{$k} = $seen_types{'base.2'}{$i}{$k};    # the value of the hash is the real gene name, this is needed below for the multipleMapper hash, which needs the gene length
								$preCount2{insert}{$i}   += scalar keys %{ $preCount{insert}{$i}{$k} };
								$preCount2{'base.2'}{$i} += scalar keys %{ $preCount{'base.2'}{$i}{$k} };
								$preCount2{'base.1'}{$i} += scalar keys %{ $preCount{'base.1'}{$i}{$k} };
							}
						}
					}

					%{ $intersection{'insert'}{$i} } = %intersect;                                   # the insert intersection is the intersect of the two bases

					# but the bases we have to check the intersection of each base with the insert
					# if the insert intersection is 0, we revert to using the original seen_types hashes
					if ( scalar keys %intersect == 0 ) {
						%{ $intersection{'base.1'}{$i} } = %{ $seen_types{'base.1'}{$i} };
						%{ $intersection{'base.2'}{$i} } = %{ $seen_types{'base.2'}{$i} };
						%{ $intersection{'insert'}{$i} } = %{ $seen_types{'insert'}{$i} };
						$preCount2{'base.1'}{$i} = 0;
						$preCount2{'base.2'}{$i} = 0;
						$preCount2{'insert'}{$i} = 0;
					}
					else {
						%{ $intersection{'base.1'}{$i} } = %{ $intersection{'insert'}{$i} };
						%{ $intersection{'base.2'}{$i} } = %{ $intersection{'insert'}{$i} };
					}
					if ( $i eq 'gene' ) {
						foreach my $H ( \%{ $seen_types{'base.2'}{gene} }, \%{ $seen_types{'base.1'}{gene} } ) {
							while ( my $k = each %$H ) {
								$i{$k} = 1;
							}
						}
					}
				}
				else {
					%intersection = %seen_types;
				}
			}
		}
		else {
			%intersection = %seen_types;
		}
		### <-- Paired end filtering ###
		#############################################################################################################################################

		#############################################################################################################################################
		# THIS IS PROBABLY THE HEART OF EVERYTHING. HERE WE GIVES ABUNDANCES TO EACH GENE, MOTU, COG, ...
		#############################################################################################################################################
		$INSERT_COUNTER++;    # this is just a cunter from creating new entries in the multipleMapper hash
		my %count;            # this hash saves the count from gene iteration, which needs ot be used for the other taxa

		# These were originally inside while ref id loop. but moved out to save time
		my ( $taxaToGiveTo, $count, $position, $length );

		# LOOP OVER EACH TYPE (base.1, base.2, insert) IN THE HASH, WHICH CONTAINS EACH GENE FOR EACH TYPE
		while ( my $type = each %hash ) {
			foreach my $i (@levels) {
				my %tempHash = %{ $hash{$type}{gene} };
				while ( my $ref_id_index = each %tempHash ) {    # hash{gene} contains gene IDs, the other hash{level} are empty

					#############################################################################################################################################
					# GET VALUE AND ID
					#############################################################################################################################################
					my $value = $tempHash{$ref_id_index};
					if ( $i eq 'gene' ) {
						$i_ref        = "$ref_id_index";
						$length       = @{ $HASH{$ref_id_index} }[1] - @{ $HASH{$ref_id_index} }[0] + 1;    # the stop-start+1
						$count        = scalar keys %{ $intersection{$type}{$i} };
						$count{$type} = $count;                                                             # we have to save the count for the genes for using with the other taxa in the next iteration
						$taxaToGiveTo = $count;

					}
					elsif ( $i eq 'mOTU' ) {
						$length       = @{ $HASH{$ref_id_index} }[1] - @{ $HASH{$ref_id_index} }[0] + 1;    # the stop-start+1
						$i_ref        = @{ $HASH{"$ref_id_index"} }[$INIT_NUM];
						$taxaToGiveTo = scalar keys %{ $intersection{$type}{$i} };
						if ( $preCount2{$type}{$i} ) {
							$count = $preCount2{$type}{$i};
						}
						else {
							$count = $count{$type};
						}

					}
					else {
						$length       = $miniNCBImap{"length"}{$ref_id_index};
						$i_ref        = $miniNCBImap{$i}{$ref_id_index};
						$taxaToGiveTo = scalar keys %{ $intersection{$type}{$i} };
						if ( $preCount2{$type}{$i} ) {
							$count = $preCount2{$type}{$i};
						}
						else {
							$count = $count{$type};
						}
					}
					#############################################################################################################################################

					#############################################################################################################################################
					# HERE IS THE ACTUAL GIVING
					#############################################################################################################################################
					if ( exists $intersection{$type}{$i}{$i_ref} ) {    #COUNT THIS TAXA

						#############################################################################################################################################
						# ADD TO THE TOTAL AND ALSO TO THE HASH/HASH_LEVELS
						#############################################################################################################################################
						# so these lines just below appear three times within the next 50ish lines
						if ( $type eq 'insert' ) {
							$TOTAL_MAPPED_LENGTH_NORM{'insert'}{$i} += ( 1 / $count ) * $value / $length;    # this value is used further down when calculating the scaled values
							$TOTAL_MAPPED{'insert'}{$i} += 1 / $count * $value;
							$position = 3;
						}
						else {
							$TOTAL_MAPPED_LENGTH_NORM{'base'}{$i} += ( 1 / $count ) * $value / $length;      # this value is used further down when calculating the scaled values
							$TOTAL_MAPPED{'base'}{$i} += ( 1 / $count ) * $value;
							$position = 2;
						}

						if ( $i eq 'gene' ) {
							@{ $HASH_local{$i_ref} }[$position] += 1 / $count * $value;
						}
						else {
							@{ $HASH_LEVELS_local{$i}{$i_ref} }[$position] += 1 / $count * $value;

							# we also have to add the gene length norm, for functional, these are added in the Functional subroutine,
							# but for mOTU these have to be done on runtime here
							@{ $HASH_LEVELS_local{$i}{$i_ref} }[ $position + $INIT_NUM - 2 ] += 1 / $count * $value / $length;
						}
						#############################################################################################################################################

						#############################################################################################################################################
						# IF THE GENE IS A SINGLE MAPPER ON THIS TAXONOMIC LEVEL
						#############################################################################################################################################
						if ( $taxaToGiveTo == 1 ) {

							#print "UNIQUE === $count | $taxaToGiveTo | LEVEL=$i || $prev_insert || $i_ref += 1 / $count * $value %" . ( 1 / $count * $value ) . "\n";

							if ( $type eq 'insert' ) {
								$TOTAL_MAPPED_LENGTH_NORM{'insert.only.unique'}{$i} += ( 1 / $count ) * $value / $length;    # this value is used further down when calculating the scaled values
								$TOTAL_MAPPED{'insert.only.unique'}{$i} += ( 1 / $count ) * $value;
							}
							else {
								$TOTAL_MAPPED_LENGTH_NORM{'base.only.unique'}{$i} += ( 1 / $count ) * $value / $length;      # this value is used further down when calculating the scaled values
								$TOTAL_MAPPED{'base.only.unique'}{$i} += ( 1 / $count ) * $value;
							}

							if ( $i eq 'gene' ) {
								@{ $HASH_local{$i_ref} }[ $position + 2 ] += 1 / $count * $value;
							}
							else {
								@{ $HASH_LEVELS_local{$i}{$i_ref} }[ $position + 2 ] += 1 / $count * $value;

								# we also have to add the gene length norm, for functional, these are added in the Functional subroutine, but for mOTU & NCBI these have to be done on runtime here
								@{ $HASH_LEVELS_local{$i}{$i_ref} }[ $position + $INIT_NUM - 2 + 2 ] += 1 / $count * $value / $length;

							}
						}
						#############################################################################################################################################

						#############################################################################################################################################
						# IF THE GENE IS A MULTIPLE MAPPER ON THIS TAXONOMIC LEVEL
						#############################################################################################################################################
						else {
							my $position2;
							if ( $type eq 'base.1' ) {
								$position2 = 0;
							}
							if ( $type eq 'base.2' ) {
								$position2 = 1;
							}
							if ( $type eq 'insert' ) {
								$position2 = 2;
							}

							# here we save the inserts in a hash (could it be done in an array), 0->2 has the values for base.1 base.2 and insert, 3->5 has the arrays of the taxa to dist. to
							my ( @array, @length, $ref_hash );
							if ( scalar keys %{ $intersection{$type}{$i} } == 0 ) {    # this means that the intersection hash is empty and we need to get all from the seen hash
								$ref_hash = \%seen_types;
							}
							else {
								$ref_hash = \%intersection;
							}
							my %tempHash = %{ $ref_hash->{$type}{$i} };
							if ( $i eq 'gene' ) {
								while ( my $key = each %tempHash ) {
									push @length, @{ $HASH{$key} }[1] - @{ $HASH{$key} }[0] + 1;
									push @array,  $key;
								}

							}
							else {
								while ( my $key = each %tempHash ) {
									push @array,  $tempHash{$key};
									push @length, @{ $HASH{ $tempHash{$key} } }[1] - @{ $HASH{ $tempHash{$key} } }[0] + 1;
								}
							}

							# if we have specified a temp file we save the MM datat in that and read it line by line later, rather than
							# saving the data in a larger and larger array, this is done further doen at the end of each insert
							unless ( $multipleMapper{$i}{"$INSERT_COUNTER"} ) {
								@{ $multipleMapper{$i}{"$INSERT_COUNTER"} } = ( 0, 0, 0 );
							}
							@{ $multipleMapper{$i}{"$INSERT_COUNTER"} }[$position2]       = $value;     # here the annotion is a bit tricky
							@{ $multipleMapper{$i}{"$INSERT_COUNTER"} }[ $position2 + 3 ] = \@array;    # here the annotion is a bit tricky
							@{ $multipleMapper{$i}{"$INSERT_COUNTER"} }[6]                = \@length;

							if ( $type eq 'insert' ) {

								# because mm dist norm comes form multiple sources it's a bit tricky and have to be summarized later, when distr. MMs
								# so here we only summarize the TOTAL_MAPPED, and not the TOTAL_MAPPED_LENGTH_NORM
								$TOTAL_MAPPED{'insert.mm.dist.among.unique'}{$i} += ( 1 / $count ) * $value;
							}
							else {

								# because mm dist norm comes form multiple sources it's a bit tricky and have to be summarized later, when distr. MMs
								# so here we only summarize the TOTAL_MAPPED, and not the TOTAL_MAPPED_LENGTH_NORM
								$TOTAL_MAPPED{'base.mm.dist.among.unique'}{$i} += ( 1 / $count ) * $value;
							}
						}
					}
					else {

						# the read matches a gene that isn't in the coord file
					}
				}
			}
		}    # end type

		%hash        = ();
		%seen_types  = ();
		%miniNCBImap = ();
		%preCount    = ();
		##### PROCESS INSERT ##### COPY TO HERE
		##### PROCESS INSERT ##### COPY TO HERE
		##### PROCESS INSERT ##### COPY TO HERE

		#close $fileIn;
		my $file ="";
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] CLOSE $file for reading\n";

		# here we have processed base.1 base.2 and insert, if we have a temp file we print to it and also reset the hash
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] print to MM file\n";
#		if ($temp_file) {
#
#			#foreach my $level (@levels) {
#			foreach my $level (@levels) {
#				foreach my $INSERT_COUNTER ( keys %{ $multipleMapper{$level} } ) {
#					my @array = @{ $multipleMapper{$level}{$INSERT_COUNTER} };
#					print { $open_temp_files{"$level.$thread"} } "$array[0]\t$array[1]\t$array[2]";
#					for my $i ( 3 .. 5 ) {
#						if ( $array[$i] ) {
#							print { $open_temp_files{"$level.$thread"} } "\t" . join( ",", @{ $array[$i] } );
#						}
#						else {
#							print { $open_temp_files{"$level.$thread"} } "\t";
#						}
#					}
#					print { $open_temp_files{"$level.$thread"} } "\t" . join( ",", @{ $array[6] } );
#					print { $open_temp_files{"$level.$thread"} } "\n";
#				}
#			}
#		}
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] finished printing to MM file\n";

		# here we lock the main HASH and enter the values for the hashes
		# HASH_LEVELS
		#		print STDERR MOCATProfiling::Misc::getLoggingTime() . " : [THREAD $thread] LOCK HASH_LEVELS\n";
		#		#lock(%HASH_LEVELS);
		#		foreach my $key1 ( keys %HASH_LEVELS_local ) {
		#			foreach my $key2 ( keys %{ $HASH_LEVELS_local{$key1} } ) {
		#				my @array = @{ $HASH_LEVELS_local{$key1}{$key2} };
		#				if ( exists $array[2] ) {    # only add it if it has data
		#					                   #print join (" ", @array) . "\n";
		#					                   # the length is empty in the local HASH
		#					                   #@{ $HASH_LEVELS{$key1}{$key2} }[0] = $array[0];    # start is always the same
		#					                   #@{ $HASH_LEVELS{$key1}{$key2} }[1] = $array[1];    # stop is always the same
		#					for my $i ( 2 .. scalar @array - 1 ) {
		#						if ( $array[$i] ) {
		#							@{ $HASH_LEVELS{$key1}{$key2} }[$i] += $array[$i];
		#						}
		#					}
		#				}
		#			}
		#		}

		# HASH
		#print STDERR MOCATProfiling::Misc::getLoggingTime() . " : [THREAD $thread] SAVING LOCAL HASH to $file.returned\n";
		#store \%HASH_local, "$file.returned";
		# 		#lock(%HASH);
		#		my $toPrint;
		#		foreach my $key1 ( keys %HASH_local ) {
		#			my @array = @{ $HASH_local{$key1} };
		#			my @array2 = @{ $HASH{$key1} };
		#			if ( exists $array[2] ) {    # only add it if it has data
		#				                   # the length is empty in the local HASH
		#				                   #@{ $HASH{$key1} }[0] = $array[0];    # start is always the same
		#				                   #@{ $HASH{$key1} }[1] = $array[1];    # stop is always the same
		#				                   #print "$thread array " . join( " ", @array ) . "\n";
		#				$toPrint .= "$key1\t$array2[0]\t$array2[1]";
		#				for my $i ( 2 .. scalar @array - 1 ) {
		#					if ( $array[$i] ) {
		#						#@{ $HASH{$key1} }[$i] += $array[$i];
		#						$toPrint .= "\t$array[$i]";
		#					} else {
		#						$toPrint .= "\t0";
		#					}
		#				}
		#				$toPrint .= "\n";
		#			}
		#		}

		# these are the files we need to process later on
		#lock(@RETURNED);
		#push @RETURNED, "$file.returned";
		#print STDERR MOCATProfiling::Misc::getLoggingTime() . " : [THREAD $thread] SAVED LOCAL HASH\n";

		# delete files when we are done with them
		#unlink "$file"      or warn "ERROR: Could not unlink $file: $!\n";
		#unlink "$file.hash" or warn "ERROR: Could not unlink $file.hash: $!\n";
		if ( $mode eq 'mOTU' || $mode eq 'NCBI' ) {

			#unlink "$file.hashlevels" or warn "ERROR: Could not unlink $file.hashlevels: $!\n";
		}

		# here we save the HASH and HASH_LEVELS to files
		# by doing it here, we avoid that it grows bigger and bigger
		# and we can summarize in another queue in the main loop at runtime
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] SAVING $file.return.HASH(_LEVELS)\n";
		store \%HASH_local,        "$file.return.HASH";
		store \%HASH_LEVELS_local, "$file.return.HASH_LEVELS";
		%HASH_local        = ();    # becasue we save them, we reset them
		%HASH_LEVELS_local = ();    # becasue we save them, we reset them
		$return_q->enqueue("$file.return");

		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] SAVED $file.return.HASH(_LEVELS)\n";

	#}

	for my $ref ( \$TOTAL_LENGTH, \$TOTAL_LENGTH_MAPPED, \$reference_for_read_do_not_exists_in_coord_file, \$reference_for_read_exists_in_coord_file ) {
		$$ref = 0 unless $$ref;
	}

	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] RETURN [mapped=$TOTAL_LENGTH_MAPPED | total length=$TOTAL_LENGTH | ref in coord=$reference_for_read_exists_in_coord_file | ref not in coord=$reference_for_read_do_not_exists_in_coord_file]\n";
	lock $FINISHED_THREADS;
	$FINISHED_THREADS++;
	my @return = ( $TOTAL_LENGTH, $TOTAL_LENGTH_MAPPED, \%TOTAL_MAPPED, \%TOTAL_MAPPED_LENGTH_NORM, $reference_for_read_do_not_exists_in_coord_file, $reference_for_read_exists_in_coord_file );
	return ( \@return );
}

1;
