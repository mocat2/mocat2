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

    # old length
    #$length     = length $line[9];
    # new length
    #$line[5] =~ /.*(\d+)M.*/;
    #$length = $1;

    my $cigar = $line[5];
    while ( $cigar !~ /^$/ )
    {
      if ( $cigar =~ /^([0-9]+[MIDNSHP])/ )
      {
        my $cigar_part = $1;
        if ( $cigar_part =~ /(\d+)M/ )
        {
          $length = $1;
        }
        $cigar =~ s/^$cigar_part//;
      } else
      {
        die "ERROR & EXIT: Unexpected cigar = $cigar\n";
      }
    }

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

