#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

my ( $minlength, $identity, $soapmaxmm, $toRemoveFile, $inputFile, $outputFile, $extractedFile, $ZIP, $ziplevel, %hash, $print, @inputLane, @outputLane, @extractedLane, $statsFile, $EstatsFile );
my $max      = 0;
my $length   = 0;
my $sum      = 0;
my $Emax     = 0;
my $Elength  = 0;
my $Esum     = 0;
my $inserts  = 0;
my $Einserts = 0;
my $bwa;
my $bindir;

print STDERR "MOCATScreen_filter.pl: STARTING.\n";

GetOptions(
            'toremove=s'  => \$toRemoveFile,
            'in=s{,}'     => \@inputLane,
            'out=s{,}'    => \@outputLane,
            'ex=s{,}'     => \@extractedLane,
            'stats=s'     => \$statsFile,
            'estats=s'    => \$EstatsFile,
            'zip=s'       => \$ZIP,
            'ziplevel=i'  => \$ziplevel,
            'soapmaxmm=s' => \$soapmaxmm,
            'identity=s'  => \$identity,
            'length=s'    => \$minlength,
            'bwa=s'       => \$bwa,
            'bindir=s'    => \$bindir
);

unless ($bwa)
{
 print STDERR "MOCATScreen_filter.pl: processing SOAP hash...\n";
  open IN, "<$toRemoveFile";
  while (<IN>)
  {
    chomp;

    #s/\/.*$//; # no trailing /1 or /2 from now on
    $hash{$_} = 1;
  }
  close IN;
} else
{
 print STDERR "MOCATScreen_filter.pl: processing BWA hash...\n";
  open my $IN, "$bindir/samtools view $bwa | " or die "ERROR & EXIT: Cannor execute: $bindir/samtools view $bwa | ";
  while ( my $line = <$IN> )
  {
    unless ( $line =~ m/^(\S+)\/[12]\t/ )
    {
      die "INTERNAL ERROR: Please update MOCAT2 to support read names on the form '$line'";
    }
    $hash{$1} = 1;
  }
  close $IN;
}
print STDERR "MOCATScreen_filter.pl: loaded hash ok, processing...\n";

my $EsumLastFile = 0;
my $EsumPair2    = 0;
my $sumLastFile  = 0;
my $sumPair2     = 0;

for my $j ( 0 .. scalar @inputLane - 1 )
{
  print STDERR "MOCATScreen_filter.pl: Processing lane $j:$inputLane[$j]...\n";
  my $inputLane     = $inputLane[$j];
  my $outputLane    = $outputLane[$j];
  my $extractedLane = $extractedLane[$j];

  foreach my $i ( "pair.1.fq.gz", "single.fq.gz", "pair.2.fq.gz" )
  {

    $inputFile = "$inputLane.$i";
    if ($outputLane)    { $outputFile    = "$outputLane.$i"; }
    if ($extractedLane) { $extractedFile = "$extractedLane.$i"; }

    if ( $i eq "pair.1.fq.gz" )
    {
      print STDERR "MOCATScreen_filter.pl: pair.1 : Set sum and Esum to 0\n";
      $EsumLastFile = 0;
      $sumLastFile  = 0;
    }
    if ( $i eq "pair.2.fq.gz" || $i eq "single.fq.gz" )
    {
      print STDERR "MOCATScreen_filter.pl: pair.2 or single : SUM LAST FILE inserts:$sumLastFile e.inserts:$EsumLastFile\n";
      print STDERR "MOCATScreen_filter.pl: pair.2 or single : BEFORE inserts:$inserts e.inserts:$Einserts\n";
      $inserts  += $sumLastFile;
      $Einserts += $EsumLastFile;
      print STDERR "MOCATScreen_filter.pl: pair.2 or single : AFTER  inserts:$inserts e.inserts:$Einserts\n";
      $EsumLastFile = 0;
      $sumLastFile  = 0;
    }

    if ( $outputFile && $extractedFile ) { print STDERR "MOCATScreen_filter.pl: processing $inputFile -> $outputFile (and $extractedFile), deleting entries in $toRemoveFile\n"; }

    my $counter = 1;
    open IN, "gzip -c -d $inputFile |";
    if ($outputFile)
    {
      print STDERR "Writing to $outputFile\n";
      open OUT, "| $ZIP -$ziplevel -c > $outputFile" or die "ERROR & EXIT: Cannot write $outputFile";
      print STDERR "MOCATScreen_filter.pl: opened $outputFile for output\n";
    }
    if ($extractedFile)
    {
      print STDERR "Writing to $extractedFile\n";
      open EX, "| $ZIP -$ziplevel -c > $extractedFile" or die "ERROR & EXIT: Cannot write $extractedFile";
      print STDERR "MOCATScreen_filter.pl: opened $extractedFile for output\n";
    }
    while (<IN>)
    {
      chomp;
      my $line = $_;
      if ( $counter == 1 )
      {
        $line =~ m/^@(.*)\/[12]/;
        if ( $hash{$1} )
        {
          $print = 0;
        } else
        {
          $print = 1;
        }
      }
      if ($print)
      {
        if ($outputFile) { print OUT "$line\n"; }
      } else
      {
        if ($extractedFile) { print EX "$line\n"; }
      }
      if ( $counter == 2 )
      {
        if ($print)
        {
          if ( $max < length($line) )
          {
            $max = length($line);
          }
          $length += length($line);
          $sum++;
          $sumLastFile++;
        } else
        {
          if ( $Emax < length($line) )
          {
            $Emax = length($line);
          }
          $Elength += length($line);
          $Esum++;
          $EsumLastFile++;
        }
      }
      $counter++;
      if ( $counter == 5 )
      {
        $counter = 1;
      }
    }
    close IN;
    if ($outputFile)
    {
      close OUT or die("MOCATScreen_filter: error closing output pipe: $!");
    }
    if ($extractedFile)
    {
      close EX or die("MOCATScreen_filter: error closing extraction pipe: $!");
    }
  }
}

my $avg = 0;
unless ( $sum == 0 )
{
  $avg = int( ( $length / $sum ) + 0.5 );
}
my $kmer;
if ( $avg % 2 == 0 )
{
  $kmer = $avg / 2;
  if ( $kmer % 2 == 0 )
  {
    $kmer = $kmer + 1;
  } else
  {
    $kmer = $kmer + 2;
  }
} else
{
  $kmer = ( $avg + 1 ) / 2;
  if ( $kmer % 2 == 1 )
  {
    $kmer = $kmer + 0;
  } else
  {
    $kmer = $kmer + 1;
  }
}

open STAT, ">$statsFile";
print STAT "Reads\tBases\tMax\tAvg\tKmer\tInserts\tMin % identity\tMin length\tSOAP max mismatches\n$sum\t$length\t$max\t$avg\t$kmer\t$inserts\t$identity\t$minlength\t$soapmaxmm\n";
close STAT;

$avg = 0;
unless ( $Esum == 0 )
{
  $avg = int( ( $Elength / $Esum ) + 0.5 );
}
if ( $avg % 2 == 0 )
{
  $kmer = $avg / 2;
  if ( $kmer % 2 == 0 )
  {
    $kmer = $kmer + 1;
  } else
  {
    $kmer = $kmer + 2;
  }
} else
{
  $kmer = ( $avg + 1 ) / 2;
  if ( $kmer % 2 == 1 )
  {
    $kmer = $kmer + 0;
  } else
  {
    $kmer = $kmer + 1;
  }
}

open STAT, ">$EstatsFile";
print STAT "Reads\tBases\tMax\tAvg\tKmer\tInserts\tMin % identity\tMin length\tSOAP max mismatches\n$Esum\t$Elength\t$Emax\t$avg\t$kmer\t$Einserts\t$identity\t$minlength\t$soapmaxmm\n";
close STAT;

print STDERR "MOCATScreen_filter.pl: Printed stats and done.\n";
print STDERR "MOCATScreen_filter.pl: FINISHED.\n";

exit 0;
