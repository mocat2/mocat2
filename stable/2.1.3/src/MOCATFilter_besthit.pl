#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

$|++;

my $prev_read = "";
my ( %fw, %rev, @keep_reads, $max_score );
my $total = 0;
my $keep  = 0;

print STDERR getLoggingTime() . " MOCATFilter - besthit : started\n";
while (<STDIN>)
{
  my $line = $_;
  my @line = split "\t", $line;
  $total++;
  if ( $line[0] ne $prev_read )
  {
    $keep += scalar @keep_reads;
    foreach my $read (@keep_reads)
    {
      print STDOUT $read;
    }
    @keep_reads = ();
    $max_score  = -1;
  }
  my $AS;
  if ( $line[11] =~ m/AS:.*:(.*)/ )
  {
    $AS = $1;
  } elsif ( $line[12] =~ m/AS:.*:(.*)/ )
  {
    $AS = $1;
  } elsif ( $line[13] =~ m/AS:.*:(.*)/ )
  {
    $AS = $1;
  } elsif ( $line[14] =~ m/AS:.*:(.*)/ )
  {
    $AS = $1;
  } elsif ( $line[15] =~ m/AS:.*:(.*)/ )
  {
    $AS = $1;
  } else
  {
    die "ERROR & EXIT: Not valid SAM Format for this script to handle (may still be a valid SAM format that this script cannot handle) (line $. : $line)";
  }
  if ($AS < 0){
   $AS = -$AS;
  }
  if ( $AS > $max_score )
  {
    @keep_reads = ();
    push @keep_reads, $line;
    $max_score = $AS;
  } elsif ( $AS == $max_score )
  {
    push @keep_reads, $line;
  }
  $prev_read = $line[0];
}
$keep += scalar @keep_reads;
foreach my $read (@keep_reads)
{
  print STDOUT $read;

}
print STDERR getLoggingTime() . " MOCATFilter - besthit : [STATS] total_reads_in=$total | total_reads_out=$keep\n";
print STDERR getLoggingTime() . " MOCATFilter - besthit : finished\n";
exit 0;

sub getLoggingTime
{
  my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime(time);
  my $nice_timestamp = sprintf( "%04d%02d%02d %02d%02d%02d", $year + 1900, $mon + 1, $mday, $hour, $min, $sec );
  return $nice_timestamp;
}
