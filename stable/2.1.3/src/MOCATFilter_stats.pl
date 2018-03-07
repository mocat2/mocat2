#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

my $prev_insert = "";
my $prev_read   = "";
my $INSERTS     = 0;
my $BASES       = 0;
my $READS       = 0;
my $MAX         = 0;
my $IN          = 0;
my $insert;
my $format;
my $identity;
my $length;
my $screened;
my $extracted;
my $stats;
my $prev_length;
my $read;
my $db;

GetOptions(
            'format:s'    => \$format,
            'identity:s'  => \$identity,
            'length:s'    => \$length,
            'screened:s'  => \$screened,
            'extracted:s' => \$extracted,
            'stats:s'     => \$stats,
            'db:s'        => \$db
);

print STDERR getLoggingTime() . " MOCATFilter - stats : started\n";
while (<STDIN>)
{
  $IN++;
  chomp( my $line = $_ );
  my @line = split "\t", $_;
  my $length;

  if ( $format eq 'SAM' )
  {
    $line[0] =~ m/.*\/([12])$/;
    $read   = $line[0];
    $insert = $1;

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
    $length = $matchLength;
##################### PROCESS CIGAR #####################

    # continue
  } elsif ( $format eq 'SOAP' )
  {
    $line[0] =~ m/.*\/([12])$/;
    $read   = $line[0];
    $insert = $1;
    $length = $line[5];
  } else
  {
    die "INTERNAL ERROR UNKNOWN FORMAT";
  }
  unless ($insert)
  {
    die "INTERNAL ERROR CANNOT DETECT /1 OR /2";
  }
  if ( $insert ne $prev_insert && $prev_insert ne '' )
  {
    $INSERTS++;
  }
  if ( $read ne $prev_read && $prev_read ne '' )
  {
    $READS++;
    $BASES += $prev_length;
  }
  if ( $length > $MAX )
  {
    $MAX = $length;
  }
  $prev_read   = $read;
  $prev_insert = $insert;
  $prev_length = $length;
  print STDOUT "$line\n";
}
$INSERTS++;
$READS++;
$BASES += $prev_length;

print STDERR getLoggingTime() . " MOCATFilter - stats : done processing, printing stats\n";
open IN, "<$stats" or die "ERROR & EXIT: Cannot open file $stats (error was $!)";
chomp( my $S1 = <IN> );
chomp( my $S2 = <IN> );
close IN;
my @S2 = split "\t", $S2;
if ( !$S2[6] )
{
  $S2[6] = '';
  $S2[7] = '';
  $S2[8] = '';
}

my $ostats = $stats;
$ostats =~ s/.stats$//;
$ostats = "$ostats.filtered.on.$db.l$length.p$identity.stats";
open STATS, ">$ostats" or die("Could not open output stats file ($ostats): $!");
print STATS "Reads\tBases\tMax\tAvg\tKmer\tInserts\tMin % identity\tMin length\tSOAP max mismatches\n";
print STATS "$READS\t$BASES\t$MAX\t" . $BASES / $READS . "\t$S2[4]\t$INSERTS\t$S2[6]\t$S2[7]\t$S2[8]\n";
close STATS;
print STDERR getLoggingTime() . " MOCATFilter - stats : [STATS] total_reads_in=$IN | unique_reads=$READS | total_bases=$BASES\n";
print STDERR getLoggingTime() . " MOCATFilter - stats : finished\n";

exit 0;

sub getLoggingTime
{
  my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime(time);
  my $nice_timestamp = sprintf( "%04d%02d%02d %02d%02d%02d", $year + 1900, $mon + 1, $mday, $hour, $min, $sec );
  return $nice_timestamp;
}

