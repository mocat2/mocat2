#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# This code is part of MOCAT2
# Code is (c) Copyright Jens Roat Kultima, 2016
# This code is released under GNU GPL v3.

my ( $line, $read );

foreach my $file (@ARGV)
{
  open my $IN, "gunzip -c $file | ";
  while ( chomp( $line = <$IN> ) )
  {
    if ( $. % 4 == 1 )
    {
      if ( $line =~ m/^\S+\/([12])$/ )
      {
        print "$line/$1\n";
      } else
      {
        die "INTERNAL ERROR: FastQ format '$line' not supported. Please coorect MOCAT2 source code.";
      }
    } else
    {
      print "$line\n";
    }
  }
}

