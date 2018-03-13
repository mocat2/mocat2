

##################### PROCESS CIGAR #####################
  my $CIGAR             = $F[5];
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
  $length = $matchLength;
##################### PROCESS CIGAR #####################

  head cigars | perl -wlane '
chomp;
my $CIGAR = $_;

@CIGAR = ($CIGAR =~ /([0-9]+)([A-Z]+)/gi);

  my @CIGAR             = ( $CIGAR =~ /([0-9]+)([A-Z]+)/gi );
  my $lengthBeforeMatch = 0;
  my $matchLength       = 0;
  my $match             = 0;
  my $addToMatch        = 0;
  my $addToMatchEnd     = 0;
  my $atTheEnd = 0;
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
       if ($atTheEnd && !($CIGAR[ $i + 1 ] eq "H" || $CIGAR[ $i + 1 ] eq "S")) {
        die "INTERNAL ERROR 3: No support for CIGAR letter $CIGAR[$i+1] ($i+1; $CIGAR[$i]$CIGAR[$i+1]). Please correct MOCAT2 accordingly. CIGAR=$CIGAR match=$match";
       }
        if ( ( $CIGAR[ $i + 1 ] eq "I" ) )
        {
          $addToMatch += $CIGAR[$i];
        } elsif ( $CIGAR[ $i + 1 ] eq "D" )
        {
          # no nothing
        } elsif ( $CIGAR[ $i + 1 ] eq "H" || $CIGAR[ $i + 1 ] eq "S") {
         $atTheEnd = 1;
         $addToMatchEnd += $CIGAR[$i]; 
        }else
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
  $totalLength = $lengthBeforeMatch + $matchLength + $addToMatchEnd;
  
print join ("", @CIGAR) . " matchLength=$matchLength :: before=$lengthBeforeMatch :: length=$totalLength";

'