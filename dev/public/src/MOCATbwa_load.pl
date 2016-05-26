#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# This code is part of MOCAT2
# Code is (c) Copyright Jens Roat Kultima, 2016
# This code is released under GNU GPL v3.

my ( $line, $read );

foreach my $file ( @ARGV[ 1 .. $#ARGV ] )
{
  open my $IN, "gunzip -c $file | ";
  while (<$IN>)
  {
    chomp( $line = $_ );
    if ( $. % 4 == 1 )
    {
      if ( $line =~ m/^\S+\/([12])$/ )
      {
        print "$line/$1\n";
      } else
      {
        die "INTERNAL ERROR: FastQ format '$line' not supported. Please coorect MOCAT2 source code.";
      }
    } elsif ( $. % 4 == 2 && $ARGV[0] == 2 )
    {
      print "N\n";
    } elsif ( $. % 4 == 0 && $ARGV[0] == 2 )
    {
      print "E\n";
    } else
    {
      print "$line\n";
    }
  }
}

