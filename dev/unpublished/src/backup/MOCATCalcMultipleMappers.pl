#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

# %targets, %targetsllByAllHash, @insertTarget, @readTarget keep track of generating the all by all multiple mappers array

print "MultipleMappers.pl :: initializing\n";

my $prev_insert               = 0;
my $prev_read                 = 0;
my $notInCoordFile            = 0;
my $inCoordFile               = 0;
my $counter                   = 0;
my $counterread               = 0;
my $total                     = 0;
my $totalread                 = 0;
my $total_interactions_insert = 0;
my $total_interactions_read   = 0;
my ( %map, $map, $level, %targets, %allByAllHash, @insertTarget, @readTarget, %hash, %hashread, %hashstore, %hashstoreread, @inserts, $total_inserts, $input_file, $insert_stats_file, $input, $sample, $cwd, $data_type, $match, $format, $positionFile, $data_dir, $prog_type, $read );

GetOptions(
	'in=s'     => \$input,
	'format=s' => \$format,
	'sample=s' => \$sample,
	'sum=s'    => \$level,
	'map=s'    => \$map
);

# Load map
if ( $level ne 'none' ) {
	open MAP, "<$map";
	while (<MAP>) {
		chomp;
		my @line = split "\t";
		if ( $level eq 'taxid' ) {
			$map{ $line[0] } = $line[0];
		}
		if ( $level eq 'kingdom' ) {
			$map{ $line[0] } = $line[1];
		}
		if ( $level eq 'phylum' ) {
			$map{ $line[0] } = $line[2];
		}
		if ( $level eq 'class' ) {
			$map{ $line[0] } = $line[3];
		}
		if ( $level eq 'order' ) {
			$map{ $line[0] } = $line[4];
		}
		if ( $level eq 'family' ) {
			$map{ $line[0] } = $line[5];
		}
		if ( $level eq 'genus' ) {
			$map{ $line[0] } = $line[6];
		}
		if ( $level eq 'species' ) {
			$map{ $line[0] } = $line[7];
		}
		if ( $level eq 'curated.species' ) {
			$map{ $line[0] } = $line[8];
		}
	}
	close MAP;
}

#Parse input
unless ( -e "$input" ) {
	die "ERROR & EXIT: Missing input file $input_file\n";
}
if ( $format eq 'BAM' ) {
	open IN, "samtools view $input | ";
}
elsif ( $format eq 'BAM2' ) {
	open IN, "<$input";
}
elsif ( $format eq 'SOAP' ) {
	if ( $input =~ m/gz$/ ) {
		open IN, "zcat $input | ";
	}
	else {
		open IN, "<$input";
	}
}
my $insert;
my $ref_id;

print "MultipleMappers.pl :: processing file\n";
while (<IN>) {
	chomp;
	my @line = split "\t", $_;
	my $first_base;
	my $length;
	my $last_base;
	if ( $format =~ m/BAM/ ) {
		$read = $line[0];
		$line[0] =~ m/(.+)\/[12]$/;    #get insert id
		$insert     = $1;
		$ref_id     = $line[2];
		$first_base = $line[3];
		$length     = length $line[9];
		$last_base  = $first_base + $length - 1;    # read will be revcomp if direction is "-", so $first_base is always left-most base
	}
	elsif ( $format eq 'SOAP' ) {
		$read = $line[0];
		$line[0] =~ m/(.+)\/[12]$/;                 #get insert id
		$insert     = $1;
		$ref_id     = $line[7];
		$first_base = $line[8];
		$length     = $line[5];
		$last_base  = $first_base + $length - 1;    # read will be revcomp if direction is "-", so $first_base is always left-most base
	}

	# Do this when next insert/read
	if ( $insert ne $prev_insert && $prev_insert ne '0' ) {

		# If it's a multiple mapper
		for my $i ( 0 .. scalar @insertTarget - 1 ) {
			for my $j ( $i .. scalar @insertTarget - 1 ) {

				if ( $allByAllHash{ $insertTarget[$i] }{ $insertTarget[$j] }{'i'} ) {
					$allByAllHash{ $insertTarget[$i] }{ $insertTarget[$j] }{'i'}++;
				}
				else {
					$allByAllHash{ $insertTarget[$i] }{ $insertTarget[$j] }{'i'} = 1;
				}
				if ( $insertTarget[$i] ne $insertTarget[$j] ) {
					$total_interactions_insert++;
				}
			}
		}

		@insertTarget = ();

		# Continue with basic stats
		$hashstore{$counter}++;
		$counter = 0;
		%hash    = ();
		$total++;
	}
	if ( $read ne $prev_read && $prev_read ne '0' ) {

		# If it's a multiple mapper
		for my $i ( 0 .. scalar @readTarget - 1 ) {
			for my $j ( $i .. scalar @readTarget - 1 ) {
				if ( $allByAllHash{ $readTarget[$i] }{ $readTarget[$j] }{'r'} ) {
					$allByAllHash{ $readTarget[$i] }{ $readTarget[$j] }{'r'}++;
				}
				else {
					$allByAllHash{ $readTarget[$i] }{ $readTarget[$j] }{'r'} = 1;
				}
				if ( $readTarget[$i] ne $readTarget[$j] ) {
					$total_interactions_read++;
				}
			}
		}

		@readTarget = ();

		# Continue with basic stats
		$hashstoreread{$counterread}++;
		$counterread = 0;
		%hashread    = ();
		$totalread++;
	}

	# If read hasn't been counted, count it
	unless ( $hash{$ref_id} ) {
		$hash{$ref_id} = 1;
		$counter++;
	}
	unless ( $hashread{$ref_id} ) {
		$hashread{$ref_id} = 1;
		$counterread++;
	}

	# Save prev read/insert
	$prev_read   = $read;
	$prev_insert = $insert;

	# Add target to target arrays
	my $ref = $ref_id;
	if ( $level ne 'none' ) {
		$ref_id =~ m/(\d+)\.*/;
		$ref = $map{$1};
	}
	push @insertTarget, $ref;
	push @readTarget,   $ref;

}
close IN;

#Process final block
$hashstore{$counter}++;
$total++;
$hashstoreread{$counterread}++;
$totalread++;

for my $i ( 0 .. scalar @insertTarget - 1 ) {
	for my $j ( $i .. scalar @insertTarget - 1 ) {
		if ( $allByAllHash{ $insertTarget[$j] }{ $insertTarget[$j] }{'i'} ) {
			$allByAllHash{ $insertTarget[$i] }{ $insertTarget[$i] }{'i'}++;
		}
		else {
			$allByAllHash{ $insertTarget[$i] }{ $insertTarget[$i] }{'i'} = 1;
		}
		if ( $insertTarget[$i] ne $insertTarget[$j] ) {
			$total_interactions_insert++;
		}
	}
}

for my $i ( 0 .. scalar @readTarget - 1 ) {
	for my $j ( $i .. scalar @readTarget - 1 ) {
		if ( $allByAllHash{ $readTarget[$j] }{ $readTarget[$j] }{'r'} ) {
			$allByAllHash{ $readTarget[$i] }{ $readTarget[$i] }{'r'}++;
		}
		else {
			$allByAllHash{ $readTarget[$i] }{ $readTarget[$i] }{'r'} = 1;
		}
		if ( $readTarget[$i] ne $readTarget[$j] ) {
			$total_interactions_read++;
		}
	}
}

print "MultipleMappers.pl :: printing results part 1\n";
open OUT, ">$input.multiplemappers.insert.stats";
my $frac  = $hashstore{'1'} / $total;
my $minus = 1 - $frac;
my $mm    = $total - $hashstore{'1'};
print OUT "$sample\tfraction_mapped_to_1\t$frac\n";
print OUT "$sample\tfraction_mapped_to_>1\t$minus\n";
print OUT "$sample\ttotal_multiple_mappers_interactions\t$total_interactions_insert\n";
print OUT "$sample\tmultiple_mappers\t$mm\n";
print OUT "$sample\ttotal\t$total\n";

foreach my $key ( sort { $a <=> $b } keys(%hashstore) ) {
	print OUT "$sample\t$key\t$hashstore{$key}\n";
}
close OUT;

open OUT, ">$input.multiplemappers.read.stats";
$frac  = $hashstoreread{'1'} / $totalread;
$minus = 1 - $frac;
$mm    = $totalread - $hashstoreread{'1'};
print OUT "$sample\tfraction_mapped_to_1\t$frac\n";
print OUT "$sample\tfraction_mapped_to_>1\t$minus\n";
print OUT "$sample\ttotal_multiple_mappers_interactions\t$total_interactions_read\n";
print OUT "$sample\tmultiple_mappers\t$mm\n";
print OUT "$sample\ttotal\t$totalread\n";

foreach my $key ( sort { $a <=> $b } keys(%hashstoreread) ) {
	print OUT "$sample\t$key\t$hashstoreread{$key}\n";
}
close OUT;

print "MultipleMappers.pl :: printing results part 2 ($level)\n";
open OUTr, ">$input.multiplemappers.read.allbyall.$level";
open OUTi, ">$input.multiplemappers.insert.allbyall.$level";
my $sum;

foreach my $key1 ( keys %allByAllHash ) {
	foreach my $key2 ( keys %{ $allByAllHash{$key1} } ) {
		foreach my $letter ( keys %{ $allByAllHash{$key1}{$key2} } ) {
			if ( $allByAllHash{$key2}{$key1}{$letter} && $allByAllHash{$key1}{$key2}{$letter} && ( $key1 ne $key2 ) ) {
				$sum = $allByAllHash{$key1}{$key2}{$letter} + $allByAllHash{$key2}{$key1}{$letter};
			}
			elsif ( $key1 ne $key2 ) {
				$sum = $allByAllHash{$key1}{$key2}{$letter};
				if ( $letter eq 'r' ) {
					print OUTr $key2 . "\t" . $key1 . "\t" . $sum . "\n";
				}
				if ( $letter eq 'i' ) {
					print OUTi $key2 . "\t" . $key1 . "\t" . $sum . "\n";
				}
			}
			else {
				$sum = $allByAllHash{$key1}{$key2}{$letter};
			}
			if ( $letter eq 'r' ) {
				print OUTr $key1 . "\t" . $key2 . "\t" . $sum . "\n";
			}
			if ( $letter eq 'i' ) {
				print OUTi $key1 . "\t" . $key2 . "\t" . $sum . "\n";
			}
		}
	}
}
close OUTr;
close OUTi;
###

exit 0;
