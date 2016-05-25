#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright Jens Roat Kultima, 2016
# This code is released under GNU GPL v3.

# $conf{screen_percent_cutoff} = $ARGV[0]
# $conf{screen_length_cutoff} = $ARGV[1]

my ( %inserts, @F, $mm, @CIGAR, $CIGAR, $lengthBeforeMatch, $matchLength, $match, $addToMatch, $as );

while ( my $line = <STDIN> )
{
  chomp $line;
  @F = split "\t", $line;
  unless ( $line =~ m/.*\s*NM:i:(\d+)\s.*/ )
  {
    die "ERROR & EXIT: Missing NM field";
  }

  $mm                = $1;

  $CIGAR             = $F[5];
  @CIGAR             = ( $CIGAR =~ /([0-9]+)([A-Z]+)/gi );
  $lengthBeforeMatch = 0;
  $matchLength       = 0;
  $match             = 0;
  $addToMatch        = 0;
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
  $as = 100 - ( $mm / $matchLength ) * 100;
  if ( $as >= $ARGV[0] && $matchLength >= $ARGV[1] )
  {
    $F[0] =~ /(.+)\/[12]$/;
    print STDOUT "$line\n";
    print STDERR "$1\n";
  }
}
