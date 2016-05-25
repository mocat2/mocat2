package MOCATProfiling::Main;
use strict;
use warnings;
use File::Slurp;
use threads;
use threads::shared;
use MOCATProfiling::Variables;
use File::Sync qw(fsync sync);
use Math::Round qw(round);
use List::Util qw[min max];
use Storable;

# this is the routine started for each thread
sub main
{
  #############################################################################################################################################
  # INITALIZE
  #############################################################################################################################################
  my $thread = threads->tid();
  my $i      = $thread - 1;
  my ( %NCBI, %horizon );
  #############################################################################################################################################

  #############################################################################################################################################
  # MAIN LOOP
  #############################################################################################################################################
  print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] Collecting inserts from MAIN thread\n";
  my $file;
  while ( defined( $file = $request_q->dequeue() ) )
  {
    sync();
    open my $fileIn, "<$file" or die "ERROR & EXIT: Missing $file";
    fsync($fileIn);
    close $fileIn;
    open $fileIn, "<$file.hash" or die "ERROR & EXIT: Missing $file.hash";
    fsync($fileIn);
    close $fileIn;
    if ($VERBOSE)
    {
      print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] LOADING JOB $file.hash(levels)\n";
    }
    my %HASH = %{ retrieve("$file.hash") };

    if ($VERBOSE)
    {
      print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] LOADED HASH, PROESSING $file\n";
    }
    open $fileIn, "gunzip -c $file | " or die "ERROR & EXIT: Cannot read from INPUT PIPE gunzip -c $file |";

    # we reset the local variables for each new file that is processed
    $prev_insert = "0";
    my $line;
    my $i_ref;
    my $countThisTaxa;
    my %miniNCBImap;
    my %preCount;
    %multipleMapper = ();
    while (<$fileIn>)
    {

      ##### READLINE #####
      ##### READLINE #####
      ##### READLINE #####
      ##### READLINE #####
      my @line = split "\t";
      $read = $line[0];
      if ( $input_file_format eq 'BAM' || $input_file_format eq 'SAM' )
      {
        $line[0] =~ m/(.+)\/([12])$/ or die "ERROR & EXIT: The BAM/SAM file does not have /1 or /2 at the end of the read ID line. Please corect this.";
        $insert    = $1;
        $direction = $2;
        $ref_id    = $line[2];

        # old length
        #$length     = length $line[9];
        # new length
        #$line[5] =~ /(\d+)M/;
        #$length = $1;

        $first_base = $line[3];

        ##################### PROCESS CIGAR #####################
        my $CIGAR             = $line[5];
        my @CIGAR             = ( $CIGAR =~ /([0-9]+)([A-Z]+)/gi );
        my $lengthBeforeMatch = 0;
        my $matchLength       = 0;
        my $match             = 0;
        my $addToMatch        = 0;
        for ( my $i = 0 ; $i < $#CIGAR ; $i += 2 )
        {

          if ( $CIGAR[ $i + 1 ] eq "M" || $CIGAR[ $i + 1 ] eq "=" )
          {
            $match++;
            $matchLength = $CIGAR[$i] + $addToMatch + $matchLength;
            $addToMatch  = 0;
          } else
          {
            if ( $match > 0 )
            {
              if ( ( $CIGAR[ $i + 1 ] eq "I" ) )
              {
                $addToMatch += $CIGAR[$i];
              } elsif ( $CIGAR[ $i + 1 ] eq "D" || $CIGAR[ $i + 1 ] eq "H" || $CIGAR[ $i + 1 ] eq "S" )
              {
                # no nothing
              } else
              {
                die "INTERNAL ERROR 1: No support for CIGAR letter $CIGAR[$i+1] ($i+1; $CIGAR[$i]$CIGAR[$i+1]). Please correct MOCAT2 accordingly. CIGAR=$CIGAR match=$match";
              }
            }    # end match
            else
            {    # this is before matching begins
              if ( !( $CIGAR[ $i + 1 ] eq "H" || $CIGAR[ $i + 1 ] eq "S" || $CIGAR[ $i + 1 ] eq "I" ) )
              {
                die "INTERNAL ERROR 2: No support for CIGAR letter $CIGAR[$i+1] ($i+1; $CIGAR[$i]$CIGAR[$i+1]). Please correct MOCAT2 accordingly. CIGAR=$CIGAR match=$match";
              }
              $lengthBeforeMatch += $CIGAR[$i];
            }    # end not matched yet
          }    # end not matching M or =
        }    # end loop over cigar
        $length     = $matchLength;
        $first_base = $first_base + $lengthBeforeMatch;
##################### PROCESS CIGAR #####################

        # continue
        $last_base = $first_base + $length - 1;    # read will be revcomp if direction is "-", so $first_base is always left-most base
        if ( $mode eq 'NCBI' )
        {
          $NCBI{taxaid}        = $line[13];
          $NCBI{kingdom}       = $line[14];
          $NCBI{phylum}        = $line[15];
          $NCBI{class}         = $line[16];
          $NCBI{order}         = $line[17];
          $NCBI{family}        = $line[18];
          $NCBI{genus}         = $line[19];
          $NCBI{species}       = $line[20];
          $NCBI{specI_cluster} = $line[21];
          $NCBI{"length"}      = $line[22];
        }
      } elsif ( $input_file_format eq 'SOAP' )
      {
        $line[0] =~ m/(.+)\/([12])$/ or die "ERROR & EXIT: The SOAP file does not have /1 or /2 at the end of the read ID line. Please corect this.";
        $insert     = $1;
        $direction  = $2;
        $ref_id     = $line[7];
        $first_base = $line[8];
        $length     = $line[5];
        $last_base  = $first_base + $length - 1;    # read will be revcomp if direction is "-", so $first_base is always left-most base
        if ( $mode eq 'NCBI' )
        {
          $NCBI{taxaid}        = $line[12];
          $NCBI{kingdom}       = $line[13];
          $NCBI{phylum}        = $line[14];
          $NCBI{class}         = $line[15];
          $NCBI{order}         = $line[16];
          $NCBI{family}        = $line[17];
          $NCBI{genus}         = $line[18];
          $NCBI{species}       = $line[19];
          $NCBI{specI_cluster} = $line[20];
          $NCBI{"length"}      = $line[21];
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
      if ( $insert ne $prev_insert && $prev_insert ne "" )
      {

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
        if ( $PE_filter eq 'yes' )
        {
          my ( $b1, $b2, %i );    # %i=intermediate hashes to save some look up time
          foreach my $i (@levels)
          {
            $b1 = scalar keys %{ $seen_types{'base.1'}{$i} };
            $b2 = scalar keys %{ $seen_types{'base.2'}{$i} };
            if ( $b1 > 0 && $b2 > 0 )
            {
              my %intersect;
              if ( $i eq 'gene' )
              {
                while ( my $k = each %{ $seen_types{'base.1'}{$i} } )
                {
                  if ( exists $seen_types{'base.2'}{$i}{$k} )
                  {
                    $intersect{$k} = 1;
                  }
                }
              } else
              {
                while ( my $k = each %{ $seen_types{'base.1'}{$i} } )
                {
                  $gene_id1 = $seen_types{'base.1'}{$i}{$k};
                  $gene_id2 = $seen_types{'base.2'}{$i}{$k};
                  if ( $gene_id2 && exists $i{$gene_id1} && exists $i{$gene_id2} )
                  {
                    $intersect{$k} = $seen_types{'base.2'}{$i}{$k};    # the value of the hash is the real gene name, this is needed below for the multipleMapper hash, which needs the gene length
                    $preCount2{insert}{$i}   += scalar keys %{ $preCount{insert}{$i}{$k} };
                    $preCount2{'base.2'}{$i} += scalar keys %{ $preCount{'base.2'}{$i}{$k} };
                    $preCount2{'base.1'}{$i} += scalar keys %{ $preCount{'base.1'}{$i}{$k} };
                  }
                }
              }

              %{ $intersection{'insert'}{$i} } = %intersect;           # the insert intersection is the intersect of the two bases

              # but the bases we have to check the intersection of each base with the insert
              # if the insert intersection is 0, we revert to using the original seen_types hashes
              if ( scalar keys %intersect == 0 )
              {
                %{ $intersection{'base.1'}{$i} } = %{ $seen_types{'base.1'}{$i} };
                %{ $intersection{'base.2'}{$i} } = %{ $seen_types{'base.2'}{$i} };
                %{ $intersection{'insert'}{$i} } = %{ $seen_types{'insert'}{$i} };
                $preCount2{'base.1'}{$i} = 0;
                $preCount2{'base.2'}{$i} = 0;
                $preCount2{'insert'}{$i} = 0;
              } else
              {
                %{ $intersection{'base.1'}{$i} } = %{ $intersection{'insert'}{$i} };
                %{ $intersection{'base.2'}{$i} } = %{ $intersection{'insert'}{$i} };
              }
              if ( $i eq 'gene' )
              {
                foreach my $H ( \%{ $seen_types{'base.2'}{gene} }, \%{ $seen_types{'base.1'}{gene} } )
                {
                  while ( my $k = each %$H )
                  {
                    $i{$k} = 1;
                  }
                }
              }
            } else
            {
              %intersection = %seen_types;
            }
          }
        } else
        {
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
        while ( my $type = each %hash )
        {
          foreach my $i (@levels)
          {
            my %tempHash = %{ $hash{$type}{gene} };

            while ( my $ref_id_index = each %tempHash )
            {    # hash{gene} contains gene IDs, the other hash{level} are empty

              #############################################################################################################################################
              # GET VALUE AND ID
              #############################################################################################################################################
              my $value = $tempHash{$ref_id_index};
              my $start = @{ $HASH{$ref_id_index} }[0];
              my $stop  = @{ $HASH{$ref_id_index} }[1];
              if ( $i eq 'gene' )
              {
                $i_ref        = "$ref_id_index";
                $length       = $stop - $start + 1;                          # the stop-start+1
                $count        = scalar keys %{ $intersection{$type}{$i} };
                $count{$type} = $count;                                      # we have to save the count for the genes for using with the other taxa in the next iteration
                $taxaToGiveTo = $count;

              } elsif ( $i eq 'mOTU' )
              {
                $length       = $stop - $start + 1;                          # the stop-start+1
                $i_ref        = @{ $HASH{"$ref_id_index"} }[$INIT_NUM];
                $taxaToGiveTo = scalar keys %{ $intersection{$type}{$i} };
                if ( $preCount2{$type}{$i} )
                {
                  $count = $preCount2{$type}{$i};
                } else
                {
                  $count = $count{$type};
                }

              } else
              {
                $length       = $miniNCBImap{"length"}{$ref_id_index};
                $i_ref        = $miniNCBImap{$i}{$ref_id_index};
                $taxaToGiveTo = scalar keys %{ $intersection{$type}{$i} };
                if ( $preCount2{$type}{$i} )
                {
                  $count = $preCount2{$type}{$i};
                } else
                {
                  $count = $count{$type};
                }
              }
              #############################################################################################################################################

              #############################################################################################################################################
              # HERE IS THE ACTUAL GIVING
              #############################################################################################################################################
              if ( exists $intersection{$type}{$i}{$i_ref} )
              {    #COUNT THIS TAXA

                #############################################################################################################################################
                # ADD TO THE TOTAL AND ALSO TO THE HASH/HASH_LEVELS
                #############################################################################################################################################
                # so these lines just below appear three times within the next 50ish lines
                if ( $type eq 'insert' )
                {
                  $TOTAL_MAPPED_LENGTH_NORM{'insert'}{$i} += ( 1 / $count ) * $value / $length;    # this value is used further down when calculating the scaled values
                  $TOTAL_MAPPED{'insert'}{$i} += 1 / $count * $value;
                  $position = 3;

                  #############################################################################################################################################
                  # HORIONTAL GENE COVERAGE
                  #############################################################################################################################################
                  if ($CALCULATE_HORIZONTAL_COVERAGE)
                  {
                    if ( $i eq 'gene' )
                    {
                      my $first = round( max( ( $horizon{start}{$i_ref} - $start ), 0 ) * 99 / $length );
                      my $last = round( min( ( $horizon{stop}{$i_ref} - $start ), $stop - $start ) * 99 / $length );
                      unless ( length @{ $HASH_local{$i_ref} }[ $INIT_NUM + 1 ] )
                      {
                        @{ $HASH_local{$i_ref} }[ $INIT_NUM + 1 ] = "0" x 100;
                      }
                      substr( @{ $HASH_local{$i_ref} }[ $INIT_NUM + 1 ], $first, $last - $first + 1 ) = "1" x ( $last - $first + 1 );
                    }
                  }
                  #############################################################################################################################################

                } else
                {
                  $TOTAL_MAPPED_LENGTH_NORM{'base'}{$i} += ( 1 / $count ) * $value / $length;    # this value is used further down when calculating the scaled values
                  $TOTAL_MAPPED{'base'}{$i} += ( 1 / $count ) * $value;
                  $position = 2;
                }

                if ( $i eq 'gene' )
                {
                  @{ $HASH_local{$i_ref} }[$position] += 1 / $count * $value;

                } else
                {
                  @{ $HASH_LEVELS_local{$i}{$i_ref} }[$position] += 1 / $count * $value;

                  # we also have to add the gene length norm, for functional, these are added in the Functional subroutine,
                  # but for mOTU these have to be done on runtime here
                  @{ $HASH_LEVELS_local{$i}{$i_ref} }[ $position + $INIT_NUM - 2 ] += 1 / $count * $value / $length;
                }
                #############################################################################################################################################

                #############################################################################################################################################
                # IF THE GENE IS A SINGLE MAPPER ON THIS TAXONOMIC LEVEL
                #############################################################################################################################################
                if ( $taxaToGiveTo == 1 )
                {
                  if ( $type eq 'insert' )
                  {
                    $TOTAL_MAPPED_LENGTH_NORM{'insert.only.unique'}{$i} += ( 1 / $count ) * $value / $length;    # this value is used further down when calculating the scaled values
                    $TOTAL_MAPPED{'insert.only.unique'}{$i} += ( 1 / $count ) * $value;
                  } else
                  {
                    $TOTAL_MAPPED_LENGTH_NORM{'base.only.unique'}{$i} += ( 1 / $count ) * $value / $length;      # this value is used further down when calculating the scaled values
                    $TOTAL_MAPPED{'base.only.unique'}{$i} += ( 1 / $count ) * $value;
                  }

                  if ( $i eq 'gene' )
                  {
                    @{ $HASH_local{$i_ref} }[ $position + 2 ] += 1 / $count * $value;
                  } else
                  {
                    @{ $HASH_LEVELS_local{$i}{$i_ref} }[ $position + 2 ] += 1 / $count * $value;

                    # we also have to add the gene length norm, for functional, these are added in the Functional subroutine, but for mOTU & NCBI these have to be done on runtime here
                    @{ $HASH_LEVELS_local{$i}{$i_ref} }[ $position + $INIT_NUM - 2 + 2 ] += 1 / $count * $value / $length;

                  }
                }
                #############################################################################################################################################

                #############################################################################################################################################
                # IF THE GENE IS A MULTIPLE MAPPER ON THIS TAXONOMIC LEVEL
                #############################################################################################################################################
                else
                {
                  my $position2;
                  if ( $type eq 'base.1' )
                  {
                    $position2 = 0;
                  }
                  if ( $type eq 'base.2' )
                  {
                    $position2 = 1;
                  }
                  if ( $type eq 'insert' )
                  {
                    $position2 = 2;
                  }

                  # here we save the inserts in a hash (could it be done in an array), 0->2 has the values for base.1 base.2 and insert, 3->5 has the arrays of the taxa to dist. to
                  my ( @array, @length, $ref_hash );
                  if ( scalar keys %{ $intersection{$type}{$i} } == 0 )
                  {    # this means that the intersection hash is empty and we need to get all from the seen hash
                    $ref_hash = \%seen_types;
                  } else
                  {
                    $ref_hash = \%intersection;
                  }
                  my %tempHash = %{ $ref_hash->{$type}{$i} };
                  if ( $i eq 'gene' )
                  {
                    while ( my $key = each %tempHash )
                    {
                      push @length, @{ $HASH{$key} }[1] - @{ $HASH{$key} }[0] + 1;
                      push @array,  $key;
                    }

                  } else
                  {
                    while ( my $key = each %tempHash )
                    {
                      push @array,  $tempHash{$key};
                      push @length, @{ $HASH{ $tempHash{$key} } }[1] - @{ $HASH{ $tempHash{$key} } }[0] + 1;
                    }
                  }

                  # if we have specified a temp file we save the MM datat in that and read it line by line later, rather than
                  # saving the data in a larger and larger array, this is done further doen at the end of each insert
                  unless ( $multipleMapper{$i}{"$INSERT_COUNTER"} )
                  {
                    @{ $multipleMapper{$i}{"$INSERT_COUNTER"} } = ( 0, 0, 0 );
                  }
                  @{ $multipleMapper{$i}{"$INSERT_COUNTER"} }[$position2] += ( 1 / $count ) * $value;    # here the annotion is a bit tricky
                  @{ $multipleMapper{$i}{"$INSERT_COUNTER"} }[ $position2 + 3 ] = \@array;               # here the annotion is a bit tricky
                  @{ $multipleMapper{$i}{"$INSERT_COUNTER"} }[6] = \@length;

                  if ( $type eq 'insert' )
                  {

                    # because mm dist norm comes form multiple sources it's a bit tricky and have to be summarized later, when distr. MMs
                    # so here we only summarize the TOTAL_MAPPED, and not the TOTAL_MAPPED_LENGTH_NORM
                    $TOTAL_MAPPED{'insert.mm.dist.among.unique'}{$i} += ( 1 / $count ) * $value;
                  } else
                  {

                    # because mm dist norm comes form multiple sources it's a bit tricky and have to be summarized later, when distr. MMs
                    # so here we only summarize the TOTAL_MAPPED, and not the TOTAL_MAPPED_LENGTH_NORM
                    $TOTAL_MAPPED{'base.mm.dist.among.unique'}{$i} += ( 1 / $count ) * $value;
                  }
                }
              } else
              {

                # the read matches a gene that isn't in the coord file
                #print STDERR "NON MATCHING READ : $type | $i | $i_ref\n";
              }
            }
          }
        }    # end type

        %horizon     = ();
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
      if ( exists $HASH{"$ref_id.0"} )
      {
        $reference_for_read_exists_in_coord_file++;
        for ( my $index = 0 ; exists $HASH{"$ref_id.$index"} ; ++$index )
        {
          my @bounds = @{ $HASH{"$ref_id.$index"} }[ 0 .. 1 ];
          if ( !( $last_base < $bounds[0] || $first_base > $bounds[1] ) )
          {
            my $subtract = 0;
            if ( $first_base < $bounds[0] )
            {
              $subtract = $bounds[0] - $first_base;
            }
            if ( $last_base > $bounds[1] )
            {
              $subtract = $subtract + $last_base - $bounds[1];
            }
            $i_ref = "$ref_id.$index";
            $hash{'insert'}{gene}{"$ref_id.$index"} = 1;
            if ( $hash{"base.$direction"}{gene}{"$ref_id.$index"} )
            {
              if ( ( $length - $subtract ) > $hash{"base.$direction"}{gene}{"$ref_id.$index"} )
              {
                $hash{"base.$direction"}{gene}{"$ref_id.$index"} = $length - $subtract;
              }
            } else
            {
              $hash{"base.$direction"}{gene}{"$ref_id.$index"} = $length - $subtract;
            }

            # These two lines are used to the horizontal gene coverage
            $horizon{start}{"$ref_id.$index"} = $first_base;
            $horizon{stop}{"$ref_id.$index"}  = $last_base;

            foreach my $level (@levels)
            {

              ### GET ID ###
              ### GET ID ###
              ### GET ID ###
              if ( $level eq 'mOTU' )
              {
                $i_ref                                                         = @{ $HASH{"$ref_id.$index"} }[$INIT_NUM];
                $preCount{"base.$direction"}{$level}{$i_ref}{"$ref_id.$index"} = 1;
                $preCount{"insert"}{$level}{$i_ref}{"$ref_id.$index"}          = 1;
              } elsif ( $level ne 'gene' )
              {
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

            if ($profiling_v2_counting)
            {
              $TOTAL_LENGTH_MAPPED += $length;    # REMEMBER TO ADD THIS WHEN ADDING BETTER BASE SUPPORT# - $subtract ;
            } else
            {
              $TOTAL_LENGTH_MAPPED += $length - $subtract;
            }
          } else
          {

            #print STDERR "===== WARNING ===== : This read matches outside of the COORD region: $read -> $ref_id ($bounds[0]->$bounds[1] | read was $first_base -> $last_base)\n";
          }
          $TOTAL_LENGTH += $length;
        }
      } else
      {
        $reference_for_read_do_not_exists_in_coord_file++;
        if ($VERBOSE)
        {
          print STDERR "===== WARNING ===== : This read does not have a reference: $read -> $ref_id (length=$length)\n";
        }

        #die "ERROR & EXIT: This read does not have a reference: $read -> $ref_id (length=$length)"
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
    if ( $PE_filter eq 'yes' )
    {
      my ( $b1, $b2, %i );    # %i=intermediate hashes to save some look up time
      foreach my $i (@levels)
      {
        $b1 = scalar keys %{ $seen_types{'base.1'}{$i} };
        $b2 = scalar keys %{ $seen_types{'base.2'}{$i} };
        if ( $b1 > 0 && $b2 > 0 )
        {
          my %intersect;
          if ( $i eq 'gene' )
          {
            while ( my $k = each %{ $seen_types{'base.1'}{$i} } )
            {
              if ( exists $seen_types{'base.2'}{$i}{$k} )
              {
                $intersect{$k} = 1;
              }
            }
          } else
          {
            while ( my $k = each %{ $seen_types{'base.1'}{$i} } )
            {
              $gene_id1 = $seen_types{'base.1'}{$i}{$k};
              $gene_id2 = $seen_types{'base.2'}{$i}{$k};
              if ( $gene_id2 && exists $i{$gene_id1} && exists $i{$gene_id2} )
              {
                $intersect{$k} = $seen_types{'base.2'}{$i}{$k};    # the value of the hash is the real gene name, this is needed below for the multipleMapper hash, which needs the gene length
                $preCount2{insert}{$i}   += scalar keys %{ $preCount{insert}{$i}{$k} };
                $preCount2{'base.2'}{$i} += scalar keys %{ $preCount{'base.2'}{$i}{$k} };
                $preCount2{'base.1'}{$i} += scalar keys %{ $preCount{'base.1'}{$i}{$k} };
              }
            }
          }

          %{ $intersection{'insert'}{$i} } = %intersect;           # the insert intersection is the intersect of the two bases

          # but the bases we have to check the intersection of each base with the insert
          # if the insert intersection is 0, we revert to using the original seen_types hashes
          if ( scalar keys %intersect == 0 )
          {
            %{ $intersection{'base.1'}{$i} } = %{ $seen_types{'base.1'}{$i} };
            %{ $intersection{'base.2'}{$i} } = %{ $seen_types{'base.2'}{$i} };
            %{ $intersection{'insert'}{$i} } = %{ $seen_types{'insert'}{$i} };
            $preCount2{'base.1'}{$i} = 0;
            $preCount2{'base.2'}{$i} = 0;
            $preCount2{'insert'}{$i} = 0;
          } else
          {
            %{ $intersection{'base.1'}{$i} } = %{ $intersection{'insert'}{$i} };
            %{ $intersection{'base.2'}{$i} } = %{ $intersection{'insert'}{$i} };
          }
          if ( $i eq 'gene' )
          {
            foreach my $H ( \%{ $seen_types{'base.2'}{gene} }, \%{ $seen_types{'base.1'}{gene} } )
            {
              while ( my $k = each %$H )
              {
                $i{$k} = 1;
              }
            }
          }
        } else
        {
          %intersection = %seen_types;
        }
      }
    } else
    {
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
    while ( my $type = each %hash )
    {
      foreach my $i (@levels)
      {
        my %tempHash = %{ $hash{$type}{gene} };

        while ( my $ref_id_index = each %tempHash )
        {    # hash{gene} contains gene IDs, the other hash{level} are empty

          #############################################################################################################################################
          # GET VALUE AND ID
          #############################################################################################################################################
          my $value = $tempHash{$ref_id_index};
          my $start = @{ $HASH{$ref_id_index} }[0];
          my $stop  = @{ $HASH{$ref_id_index} }[1];
          if ( $i eq 'gene' )
          {
            $i_ref        = "$ref_id_index";
            $length       = $stop - $start + 1;                          # the stop-start+1
            $count        = scalar keys %{ $intersection{$type}{$i} };
            $count{$type} = $count;                                      # we have to save the count for the genes for using with the other taxa in the next iteration
            $taxaToGiveTo = $count;

          } elsif ( $i eq 'mOTU' )
          {
            $length       = $stop - $start + 1;                          # the stop-start+1
            $i_ref        = @{ $HASH{"$ref_id_index"} }[$INIT_NUM];
            $taxaToGiveTo = scalar keys %{ $intersection{$type}{$i} };
            if ( $preCount2{$type}{$i} )
            {
              $count = $preCount2{$type}{$i};
            } else
            {
              $count = $count{$type};
            }

          } else
          {
            $length       = $miniNCBImap{"length"}{$ref_id_index};
            $i_ref        = $miniNCBImap{$i}{$ref_id_index};
            $taxaToGiveTo = scalar keys %{ $intersection{$type}{$i} };
            if ( $preCount2{$type}{$i} )
            {
              $count = $preCount2{$type}{$i};
            } else
            {
              $count = $count{$type};
            }
          }
          #############################################################################################################################################

          #############################################################################################################################################
          # HERE IS THE ACTUAL GIVING
          #############################################################################################################################################
          if ( exists $intersection{$type}{$i}{$i_ref} )
          {    #COUNT THIS TAXA

            #############################################################################################################################################
            # ADD TO THE TOTAL AND ALSO TO THE HASH/HASH_LEVELS
            #############################################################################################################################################
            # so these lines just below appear three times within the next 50ish lines
            if ( $type eq 'insert' )
            {
              $TOTAL_MAPPED_LENGTH_NORM{'insert'}{$i} += ( 1 / $count ) * $value / $length;    # this value is used further down when calculating the scaled values
              $TOTAL_MAPPED{'insert'}{$i} += 1 / $count * $value;
              $position = 3;

              #############################################################################################################################################
              # HORIONTAL GENE COVERAGE
              #############################################################################################################################################
              if ($CALCULATE_HORIZONTAL_COVERAGE)
              {
                if ( $i eq 'gene' )
                {
                  my $first = round( max( ( $horizon{start}{$i_ref} - $start ), 0 ) * 99 / $length );
                  my $last = round( min( ( $horizon{stop}{$i_ref} - $start ), $stop - $start ) * 99 / $length );
                  unless ( length @{ $HASH_local{$i_ref} }[ $INIT_NUM + 1 ] )
                  {
                    @{ $HASH_local{$i_ref} }[ $INIT_NUM + 1 ] = "0" x 100;
                  }
                  substr( @{ $HASH_local{$i_ref} }[ $INIT_NUM + 1 ], $first, $last - $first + 1 ) = "1" x ( $last - $first + 1 );
                }
              }
              #############################################################################################################################################

            } else
            {
              $TOTAL_MAPPED_LENGTH_NORM{'base'}{$i} += ( 1 / $count ) * $value / $length;    # this value is used further down when calculating the scaled values
              $TOTAL_MAPPED{'base'}{$i} += ( 1 / $count ) * $value;
              $position = 2;
            }

            if ( $i eq 'gene' )
            {
              @{ $HASH_local{$i_ref} }[$position] += 1 / $count * $value;

            } else
            {
              @{ $HASH_LEVELS_local{$i}{$i_ref} }[$position] += 1 / $count * $value;

              # we also have to add the gene length norm, for functional, these are added in the Functional subroutine,
              # but for mOTU these have to be done on runtime here
              @{ $HASH_LEVELS_local{$i}{$i_ref} }[ $position + $INIT_NUM - 2 ] += 1 / $count * $value / $length;
            }
            #############################################################################################################################################

            #############################################################################################################################################
            # IF THE GENE IS A SINGLE MAPPER ON THIS TAXONOMIC LEVEL
            #############################################################################################################################################
            if ( $taxaToGiveTo == 1 )
            {
              if ( $type eq 'insert' )
              {
                $TOTAL_MAPPED_LENGTH_NORM{'insert.only.unique'}{$i} += ( 1 / $count ) * $value / $length;    # this value is used further down when calculating the scaled values
                $TOTAL_MAPPED{'insert.only.unique'}{$i} += ( 1 / $count ) * $value;
              } else
              {
                $TOTAL_MAPPED_LENGTH_NORM{'base.only.unique'}{$i} += ( 1 / $count ) * $value / $length;      # this value is used further down when calculating the scaled values
                $TOTAL_MAPPED{'base.only.unique'}{$i} += ( 1 / $count ) * $value;
              }

              if ( $i eq 'gene' )
              {
                @{ $HASH_local{$i_ref} }[ $position + 2 ] += 1 / $count * $value;
              } else
              {
                @{ $HASH_LEVELS_local{$i}{$i_ref} }[ $position + 2 ] += 1 / $count * $value;

                # we also have to add the gene length norm, for functional, these are added in the Functional subroutine, but for mOTU & NCBI these have to be done on runtime here
                @{ $HASH_LEVELS_local{$i}{$i_ref} }[ $position + $INIT_NUM - 2 + 2 ] += 1 / $count * $value / $length;

              }
            }
            #############################################################################################################################################

            #############################################################################################################################################
            # IF THE GENE IS A MULTIPLE MAPPER ON THIS TAXONOMIC LEVEL
            #############################################################################################################################################
            else
            {
              my $position2;
              if ( $type eq 'base.1' )
              {
                $position2 = 0;
              }
              if ( $type eq 'base.2' )
              {
                $position2 = 1;
              }
              if ( $type eq 'insert' )
              {
                $position2 = 2;
              }

              # here we save the inserts in a hash (could it be done in an array), 0->2 has the values for base.1 base.2 and insert, 3->5 has the arrays of the taxa to dist. to
              # bugfix in v1.5.1+, in the multipleMapper hash below, we havent used array or length
              # but we need to divide by the taxiid length in DistMM, and in previous versions we accidentially
              # used the gene length. reason for not storing in hash is because we need to save it to and load form tet file
              my ( @array, @length, $ref_hash );
              if ( scalar keys %{ $intersection{$type}{$i} } == 0 )
              {    # this means that the intersection hash is empty and we need to get all from the seen hash
                $ref_hash = \%seen_types;
              } else
              {
                $ref_hash = \%intersection;
              }
              my %tempHash = %{ $ref_hash->{$type}{$i} };
              if ( $i eq 'gene' )
              {
                while ( my $key = each %tempHash )
                {
                  push @length, @{ $HASH{$key} }[1] - @{ $HASH{$key} }[0] + 1;
                  push @array,  $key;
                }

              } else
              {
                while ( my $key = each %tempHash )
                {
                  push @array,  $tempHash{$key};
                  push @length, @{ $HASH{ $tempHash{$key} } }[1] - @{ $HASH{ $tempHash{$key} } }[0] + 1;
                }
              }

              # if we have specified a temp file we save the MM datat in that and read it line by line later, rather than
              # saving the data in a larger and larger array, this is done further doen at the end of each insert
              unless ( $multipleMapper{$i}{"$INSERT_COUNTER"} )
              {
                @{ $multipleMapper{$i}{"$INSERT_COUNTER"} } = ( 0, 0, 0 );
              }
              @{ $multipleMapper{$i}{"$INSERT_COUNTER"} }[$position2] += ( 1 / $count ) * $value;    # here the annotion is a bit tricky
              @{ $multipleMapper{$i}{"$INSERT_COUNTER"} }[ $position2 + 3 ] = \@array;               # here the annotion is a bit tricky
              @{ $multipleMapper{$i}{"$INSERT_COUNTER"} }[6] = \@length;

              if ( $type eq 'insert' )
              {

                # because mm dist norm comes form multiple sources it's a bit tricky and have to be summarized later, when distr. MMs
                # so here we only summarize the TOTAL_MAPPED, and not the TOTAL_MAPPED_LENGTH_NORM
                $TOTAL_MAPPED{'insert.mm.dist.among.unique'}{$i} += ( 1 / $count ) * $value;
              } else
              {

                # because mm dist norm comes form multiple sources it's a bit tricky and have to be summarized later, when distr. MMs
                # so here we only summarize the TOTAL_MAPPED, and not the TOTAL_MAPPED_LENGTH_NORM
                $TOTAL_MAPPED{'base.mm.dist.among.unique'}{$i} += ( 1 / $count ) * $value;
              }
            }
          } else
          {

            # the read matches a gene that isn't in the coord file
            #print STDERR "NON MATCHING READ : $type | $i | $i_ref\n";
          }
        }
      }
    }    # end type

    %horizon     = ();
    %hash        = ();
    %seen_types  = ();
    %miniNCBImap = ();
    %preCount    = ();
    ##### PROCESS INSERT ##### COPY TO HERE
    ##### PROCESS INSERT ##### COPY TO HERE
    ##### PROCESS INSERT ##### COPY TO HERE

    #############################################################################################################################################
    #############################################################################################################################################
    # BEGINNING OF THE LAST PART. NOW WE HAVE PROCESSED THE WHOLE FILE
    #############################################################################################################################################
    close $fileIn;
    if ($VERBOSE)
    {
      print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] CLOSE $file for reading\n";
      print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] print to MM file\n";
    }

    # here we have processed base.1 base.2 and insert, if we have a temp file we print to it and also reset the hash (CURRENTLY WE ONLY SUPPORT TEMP FILES)
    if ($temp_file)
    {
      foreach my $level (@levels)
      {
        foreach my $INSERT_COUNTER ( keys %{ $multipleMapper{$level} } )
        {
          my @array = @{ $multipleMapper{$level}{$INSERT_COUNTER} };
          print { $open_temp_files{"$level.$thread"} } "$array[0]\t$array[1]\t$array[2]";
          for my $i ( 3 .. 5 )
          {
            if ( $array[$i] )
            {
              print { $open_temp_files{"$level.$thread"} } "\t" . join( ",", @{ $array[$i] } );
            } else
            {
              print { $open_temp_files{"$level.$thread"} } "\t";
            }
          }
          print { $open_temp_files{"$level.$thread"} } "\t" . join( ",", @{ $array[6] } );
          print { $open_temp_files{"$level.$thread"} } "\n";
        }
      }
    }
    if ($VERBOSE)
    {
      print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] finished printing to MM file\n";
    }

    # delete files when we are done with them
    chomp( my $file_size_all = `du -hsLc $temp_file* | tail -n 1 | cut -f 1` );
    chomp( my $file_size     = `du -hsLc $file* | tail -n 1 | cut -f 1` );
    unlink "$file"      or warn "ERROR: Could not unlink $file: $!\n";
    unlink "$file.hash" or warn "ERROR: Could not unlink $file.hash: $!\n";
    if ($VERBOSE)
    {
      print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] using $file_size_all ($file_size from this thread) in temp folder\n";
    }

    # here we save the HASH and HASH_LEVELS to files
    # by doing it here, we avoid that it grows bigger and bigger
    # and we can summarize in another queue in the main loop at runtime
    if ($VERBOSE)
    {
      print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] SAVING $file.return.HASH(_LEVELS)\n";
    }
    store \%HASH_local,        "$file.return.HASH";
    store \%HASH_LEVELS_local, "$file.return.HASH_LEVELS";
    %HASH_local        = ();    # becasue we save them, we reset them
    %HASH_LEVELS_local = ();    # becasue we save them, we reset them
    $return_q->enqueue("$file.return");
    my $waiting  = $request_q->pending();
    my $waiting2 = $return_q->pending();

    if ($VERBOSE)
    {
      print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] FINISHED $file.return.HASH(_LEVELS) [$waiting waiting for submission | $waiting2 waiting for merging]\n";
    }
  }    # END OF PROCESS, next queued file

  #############################################################################################################################################
  # RETURN STUFF
  #############################################################################################################################################

  #print "TM mm:$TOTAL_MAPPED{'base.mm.dist.among.unique'}{gene} to:$TOTAL_MAPPED{'base'}{gene} un:$TOTAL_MAPPED{'base.only.unique'}{gene}\n";

  for my $ref ( \$TOTAL_LENGTH, \$TOTAL_LENGTH_MAPPED, \$reference_for_read_do_not_exists_in_coord_file, \$reference_for_read_exists_in_coord_file )
  {
    $$ref = 0 unless $$ref;
  }
  print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD $thread] RETURN [mapped=$TOTAL_LENGTH_MAPPED | total length=$TOTAL_LENGTH | ref in coord=$reference_for_read_exists_in_coord_file | ref not in coord=$reference_for_read_do_not_exists_in_coord_file]\n";
  {
    lock $FINISHED_THREADS;
    $FINISHED_THREADS++;
  }
  my @return = ( $TOTAL_LENGTH, $TOTAL_LENGTH_MAPPED, \%TOTAL_MAPPED, \%TOTAL_MAPPED_LENGTH_NORM, $reference_for_read_do_not_exists_in_coord_file, $reference_for_read_exists_in_coord_file );
  return ( \@return );
  #############################################################################################################################################
}

1;
