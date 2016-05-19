#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

my ( $insert, @db, $read, %position, $ids );
my $total = 0;
my $keep  = 0;
my $SAM   = 0;

print STDERR getLoggingTime() . " MOCATFilter - remove_in_padded : started, loading coord files\n";

GetOptions( 'db:s{,}' => \@db,
            'sam'     => \$SAM, );

# Load hash
foreach my $file (@db)
{
  open POS, "<$file.coord" or killme("Cannot read $file.coord");
  while (<POS>)
  {
    chomp;
    my @line = split /\s+/;
    $position{ $line[0] } = [] unless ( exists $position{ $line[0] } );
    push @{ $position{ $line[0] } }, $line[1], $line[2];
  }
}
print STDERR getLoggingTime() . " MOCATFilter - remove_in_padded : loaded coord files, waiting\n";

# Process STDIN
my $ref_id;
my @line;
my $line;
my ( $first_base, $last_base, $length );
while (<STDIN>)
{
  $line = $_;
  @line = split "\t", $line;
  $total++;
  if ($SAM)
  {
    $ref_id = $line[2];

    $first_base = $line[3];
    my $cigar  = $line[5];
    my $seen_M = 0;
    while ( $cigar !~ /^$/ )
    {
      if ( $cigar =~ /^([0-9]+[MIDNSHP])/ )
      {
        my $cigar_part = $1;
        if ( $cigar_part =~ /(\d+)M/ )
        {
          $length = $1;
          $seen_M = 1;
        } elsif ( $cigar_part =~ /(\d+)I/ || $cigar_part =~ /(\d+)D/ || $cigar_part =~ /(\d+)N/ || $cigar_part =~ /(\d+)S/ || $cigar_part =~ /(\d+)H/ || $cigar_part =~ /(\d+)P/ )
        {
          unless ($seen_M)
          {
            $first_base += $1;
          }
        }
        $cigar =~ s/^$cigar_part//;
      } else
      {
        die "ERROR & EXIT: Unexpected cigar = $cigar\n";
      }
    }

    # old length
    #$length     = length $line[9];
    # new length
    #$line[5] =~ /.*(\d+)M.*/;
    #$length = $1;

    # continue
    $last_base = $first_base + $length - 1;    # read will be revcomp if direction is "-", so $first_base is always left-most base

  } else
  {
    $ref_id     = $line[7];
    $first_base = $line[8];
    $length     = $line[5];
    $last_base  = $first_base + $length - 1;
  }

  my $bounds = $position{$ref_id};
  ( defined $bounds ) or killme("$ref_id was not found in coord file ($line)\n");
  my $n = 0;
  while ( $n < scalar @{$bounds} )
  {
    my $lowb = $bounds->[$n];
    my $upb  = $bounds->[ $n + 1 ];
    $n += 2;

    # Here the read matches, print it if it matches, otherwise not
    if ( !( $last_base < $lowb || $first_base > $upb ) )
    {
      $keep++;
      print STDOUT "$line";
    } else {
     #print STDOUT "THROW AWAY $first_base->$last_base  AND $lowb->$upb  $line";
    }
  }
}

print STDERR getLoggingTime() . " MOCATFilter - remove_in_padded : [STATS] total_reads_in=$total | total_reads_out=$keep\n";
print STDERR getLoggingTime() . " MOCATFilter - remove_in_padded : finished\n";

exit 0;

sub getLoggingTime
{
  my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime(time);
  my $nice_timestamp = sprintf( "%04d%02d%02d %02d%02d%02d", $year + 1900, $mon + 1, $mday, $hour, $min, $sec );
  return $nice_timestamp;
}

sub killme
{
  my $err = $_[0];
  my $gip = getpgrp();
  print STDERR "ERROR & EXIT: $err\n";
  kill -9, $gip;
}
