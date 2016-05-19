#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

my $prev_insert            = 0;
my $prev_read              = 0;
my $notInCoordFile         = 0;
my $inCoordFile            = 0;
my $total_bases_covered    = 0;
my $total_length_matched   = 0;
my $total_bases            = 0;
my $seen_this_insert       = 0;
my $total_seen_inserts     = 0;
my $count_sum              = 0;
my $db_average_gene_length = 0;
my $db_gene_length         = 0;
my ( $base1, $base2 ) = 0;
my $missingFile = "";
my ( $print_rownamesTax, %len, $len, $output_tax, %taxBase, %taxInsert, @readTarget1, @readTarget2, @insertTarget, %map, $map, $ZCAT, $bin_dir, %forwardReads, %reverseReads, %hash, %position, %base1, %base2, %inserts, @inserts, $total_inserts, $input_file, $insert_stats_file, $input, $sample, $cwd, $reads, $data_type, $output, $match, $format, $positionFile, $data_dir, $counted_inserts, $file_list, $file, $falen, $coverage_file, $print_rownames, $printOnlyIfSample0 );
my @levels = ( 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'curated.species' );
my $usage = "                                                                                                                                                                                                                     
Can only be used internally in MOCAT. Sorry.
";

GetOptions(
	'i=s'           => \$input_file,
	's=s'           => \$insert_stats_file,
	'sample=s'      => \$sample,
	'cwd=s'         => \$cwd,
	'rcr=s'         => \$reads,
	'dt=s'          => \$data_type,
	'taxout=s'      => \$output_tax,
	'out=s'         => \$output,
	'match=s'       => \$match,
	'pos=s'         => \$positionFile,
	'datadir=s'     => \$data_dir,
	'file_list=s'   => \$file_list,
	'file=s'        => \$file,
	'falen=s'       => \$falen,
	'covfile=s'     => \$coverage_file,
	'rownames=s'    => \$print_rownames,
	'taxrownames=s' => \$print_rownamesTax,
	'counter=i'     => \$printOnlyIfSample0,
	'bin=s'         => \$bin_dir,
	'zcat=s'        => \$ZCAT,
	'map=s'         => \$map,
	'len=s'         => \$len
);

# BAM or SOAP
if ( $match eq "allbest" ) {
	$format = "BAM";
}
else {
	$format = "SOAP";
}

# Check input file
unless ( -e "$input_file" ) {
	die "ERROR & EXIT: Missing input file $input_file\n";
}

# Load map
open MAP, "<$map";
while (<MAP>) {
	chomp;
	my @line = split "\t";
	$map{'kingdom'}{ $line[0] }         = $line[1];
	$map{'phylum'}{ $line[0] }          = $line[2];
	$map{'class'}{ $line[0] }           = $line[3];
	$map{'order'}{ $line[0] }           = $line[4];
	$map{'family'}{ $line[0] }          = $line[5];
	$map{'genus'}{ $line[0] }           = $line[6];
	$map{'species'}{ $line[0] }         = $line[7];
	$map{'curated.species'}{ $line[0] } = $line[8];
}
close MAP;

#load taxonomy length
open LEN, "<$len";
while (<LEN>) {
	chomp;
	my @line = split "\t";
	$len{ $line[0] } = $line[1];
}
close LEN;

# Load Position hash
my $counter = -1;
my @data = split( "_AND_", $positionFile );
if ( $file_list eq 'yes' ) {
	print STDERR "CalculateTaxonomy :: Calculating coverage using a file list, need to make custom .coord file...\n";
	if ( -e "$file.coord" ) {
		print STDERR "CalculateTaxonomy :: $file has a .coord file\n";
	}
	else {
		print STDERR "CalculateTaxonomy :: Creates .coord file for $file:\n";
		print STDERR "CalculateTaxonomy :: EXE $falen -infile $file -outfile $file.tmp\n";
		system "$falen -infile $file -outfile $file.tmp";
		print STDERR "CalculateTaxonomy :: Add middle column...";
		open IN2,  "<$file.tmp"   or die "ERROR & EXIT: Cannot open input file $!\n";
		open OUT2, ">$file.coord" or die "ERROR & EXIT: Cannot open output file $!\n";
		while (<IN2>) {
			chomp;
			my @line = split;
			print OUT2 "$line[0]\t1\t$line[1]\n";
		}
		close IN2;
		close OUT2;
		print " Done.\n";
		print "CalculateTaxonomy :: EXE rm $file.tmp\n";
		system "rm $file.tmp";
	}
	open POS, "<$file.coord";
	print STDERR "CalculateTaxonomy :: Loading positions of reference $file into hash...\n";
	while (<POS>) {
		$counter++;
		chomp;
		my @line = split /\s+/;
		$position{ $line[0] }{$counter} = [ $line[1], $line[2], 0, 0 ];
		$db_gene_length += $line[2] - $line[1] + 1;
	}
}
if ( $file_list ne 'yes' ) {
	foreach my $d (@data) {
		unless ( -e "$data_dir/$d.coord" ) {
			die "ERROR & EXIT: Missing coordinates file $data_dir/$d.coord for database $data_dir/$d\n";
		}
		open POS, '<', "$data_dir/$d.coord";
		print STDERR "CalculateTaxonomy :: Loading positions of reference $data_dir/$d.coord into hash...\n";
		while (<POS>) {
			$counter++;
			chomp;
			unless (m/^\S+\s+\d+\s+\d+$/) {
				die "ERROR & EXIT: $data_dir/$d.coord has line:\n$_\nExpected the line to be of format: NONSPACECHARS\\tDIGIT\\tDIGIT\n";
			}
			my @line = split /\s+/;
			$position{ $line[0] }{$counter} = [ $line[1], $line[2], 0, 0 ];
			$db_gene_length += $line[2] - $line[1] + 1;
		}
	}
}

#Avg DB gene length
$db_average_gene_length = $db_gene_length / ( $counter + 1 );
print STDERR "CalculateTaxonomy :: total db length: $db_gene_length\n";
print STDERR "CalculateTaxonomy :: avg db entry length: $db_average_gene_length\n";

#Get number of total inserts
if ( !( -e "$insert_stats_file.stats" ) ) {
	print STDERR "CalculateTaxonomy :: missing insert file $insert_stats_file.stats, -1 fraction will be incorrect\n";
	$total_bases   = 0;
	$total_inserts = 0;
	$missingFile   = "missing_insert_file";
}
else {
	print STDERR "CalculateTaxonomy :: using stats file $insert_stats_file.stats\n";
	open IN, "$insert_stats_file.stats" or die "ERROR & EXIT: Missing file $insert_stats_file.stats\n";
	my $a;
	$a = <IN>;
	$a = <IN>;
	close IN;
	chomp $a;
	my @line = split /\t/, $a;
	$total_bases = $line[1];
	print STDERR "CalculateTaxonomy :: Total bases $total_bases\n";

	if ( $line[5] ) {
		print STDERR "CalculateTaxonomy :: Loaded inserts from .stats file.\n";
		$total_inserts = $line[5];
		print STDERR "CalculateTaxonomy :: Total inserts $total_inserts\n";
	}
	elsif ( -e "$insert_stats_file.inserts.stats" ) {
		print STDERR "CalculateTaxonomy :: Inserts file exists.\n";
		$total_inserts = `tail -1 $insert_stats_file.inserts.stats`;
		chomp $total_inserts;
		print "CalculateTaxonomy :: Loaded inserts from .inserts.stats file and saves to stats file.\n";
		open IN, "<$insert_stats_file.stats" or die "ERROR & EXIT: Missing file $insert_stats_file.stats\n";
		$a = <IN>;
		chomp $a;
		my $line1 = "$a\tInserts";
		$a = <IN>;
		chomp $a;
		my $line2 = "$a\t$total_inserts";
		close IN;
		open OUT, '>', "$insert_stats_file.stats" or die "ERROR & EXIT: Cannot open $insert_stats_file.stats for ouput.\n";
		print OUT "$line1\n$line2\n";
		close OUT;
	}
	else {
		print STDERR "CalculateTaxonomy :: Inserts file doesn't exist. Creating it in the stats file...\n";
		my $path;
		if ( $reads eq "reads.processed" ) {
			$path = "reads.processed.$data_type";
		}
		else {
			$path = "reads.screened.$reads.$data_type";
		}
		print STDERR "CalculateTaxonomy :: EXE $ZCAT $cwd/$sample/$path/*pair.1*.gz $cwd/$sample/$path/*single*.gz\n | grep -c .";
		chomp( $total_inserts = `$ZCAT $cwd/$sample/$path/*pair.1*.gz $cwd/$sample/$path/*single*.gz | grep -c .` );
		$total_inserts = $total_inserts / 4;
		open IN, "<$insert_stats_file.stats" or die "ERROR & EXIT: Missing file $insert_stats_file.stats\n";
		$a = <IN>;
		chomp $a;
		my $line1 = "$a\tInserts";
		$a = <IN>;
		chomp $a;
		my $line2 = "$a\t$total_inserts";
		close IN;
		open OUT, '>', "$insert_stats_file.stats" or die "ERROR & EXIT: Cannot open $insert_stats_file.stats for ouput.\n";
		print OUT "$line1\n$line2\n";
		close OUT;
	}
}

#Parse input
if ( $format eq 'BAM' ) {
	print STDERR "CalculateTaxonomy :: Format is BAM, IN=samtools view $input_file |\n";
	open IN, "$bin_dir/samtools view $input_file | ";
}
elsif ( $format eq 'SOAP' ) {
	print STDERR "CalculateTaxonomy :: Format is SAOP, IN=$ZCAT $input_file | \n";
	open IN, "$ZCAT $input_file | ";
}
print STDERR "CalculateTaxonomy :: Parsing file: $input_file...\n";
my $insert;
while (<IN>) {
	chomp;
	my @line = split "\t", $_;
	my $ref_id;
	my $first_base;
	my $length;
	my $last_base;
	my $read;
	my $direction;

	if ( $format eq 'BAM' ) {
		$line[0] =~ m/(.+)\/([12])$/;    #get insert id
		$read       = $line[0];
		$insert     = $1;
		$direction  = $2;
		$ref_id     = $line[2];
		$first_base = $line[3];
		$length     = length $line[9];
		$last_base  = $first_base + $length - 1;    # read will be revcomp if direction is "-", so $first_base is always left-most base
	}
	elsif ( $format eq 'SOAP' ) {
		$line[0] =~ m/(.+)\/([12])$/;               #get insert id
		$read       = $line[0];
		$insert     = $1;
		$direction  = $2;
		$ref_id     = $line[7];
		$first_base = $line[8];
		$length     = $line[5];
		$last_base  = $first_base + $length - 1;    # read will be revcomp if direction is "-", so $first_base is always left-most base
	}

	# When next insert
	if ( $insert ne $prev_insert && $prev_insert ne '0' ) {

		if ( $seen_this_insert == 1 ) {
			$total_seen_inserts++;
		}
		$seen_this_insert = 0;

		my $count = scalar keys %hash;              #how many different ref_ids have I seen?
		foreach my $ref_id ( keys %hash ) {
			for my $index ( keys %{ $hash{$ref_id} } ) {
				$count_sum += ( 1 / $count );
				@{ $position{$ref_id}{$index} }[2] += ( 1 / $count );

				#tax
				foreach my $i (@levels) {
					$ref_id =~ m/(\d+)\.*/;
					my $ref    = $map{$i}{$1};
					my $length = $len{$1};
					unless ( $taxInsert{$i}{$ref} ) {
						@{ $taxInsert{$i}{$ref} } = ( 0, 0, 0, 0 );
					}
					@{ $taxInsert{$i}{$ref} }[0] += 1 / $count;
					@{ $taxInsert{$i}{$ref} }[1] += ( 1 / $count ) / $length;
					if ( $count == 1 ) {
						@{ $taxInsert{$i}{$ref} }[2] += 1 / $count;
						@{ $taxInsert{$i}{$ref} }[3] += 1 / $count / $length;
					}
				}
			}
		}
		%hash = ();

		$count = $base1;    #how many different ref_ids have I seen?
		foreach my $ref_id ( keys %base1 ) {
			for my $index ( keys %{ $base1{$ref_id} } ) {
				$total_bases_covered += ( 1 / $count ) * $base1{$ref_id}{$index};
				@{ $position{$ref_id}{$index} }[3] += ( 1 / $count ) * $base1{$ref_id}{$index};

				#tax
				foreach my $i (@levels) {
					$ref_id =~ m/(\d+)\.*/;
					my $ref    = $map{$i}{$1};
					my $length = $len{$1};
					unless ( $taxBase{$i}{$ref} ) {
						@{ $taxBase{$i}{$ref} } = ( 0, 0, 0, 0 );
					}
					@{ $taxBase{$i}{$ref} }[0] += 1 / $count * $base1{$ref_id}{$index};
					@{ $taxBase{$i}{$ref} }[1] += 1 / $count * $base1{$ref_id}{$index} / $length;
					if ( $count == 1 ) {
						@{ $taxBase{$i}{$ref} }[2] += 1 / $count * $base1{$ref_id}{$index};
						@{ $taxBase{$i}{$ref} }[3] += 1 / $count * $base1{$ref_id}{$index} / $length;
					}

				}
			}
		}
		%base1 = ();
		$base1 = 0;

		$count = $base2;    #how many different ref_ids have I seen?
		foreach my $ref_id ( keys %base2 ) {
			for my $index ( keys %{ $base2{$ref_id} } ) {
				$total_bases_covered += ( 1 / $count ) * $base2{$ref_id}{$index};
				@{ $position{$ref_id}{$index} }[3] += ( 1 / $count ) * $base2{$ref_id}{$index};

				#tax
				foreach my $i (@levels) {
					$ref_id =~ m/(\d+)\.*/;
					my $ref    = $map{$i}{$1};
					my $length = $len{$1};
					unless ( $taxBase{$i}{$ref} ) {
						@{ $taxBase{$i}{$ref} } = ( 0, 0, 0, 0 );
					}
					@{ $taxBase{$i}{$ref} }[0] += 1 / $count * $base2{$ref_id}{$index};
					@{ $taxBase{$i}{$ref} }[1] += 1 / $count * $base2{$ref_id}{$index} / $length;
					if ( $count == 1 ) {
						@{ $taxBase{$i}{$ref} }[2] += 1 / $count * $base2{$ref_id}{$index};
						@{ $taxBase{$i}{$ref} }[3] += 1 / $count * $base2{$ref_id}{$index} / $length;
					}

				}
			}
		}
		%base2 = ();
		$base2 = 0;
	}

	# MAIN LOOP for each line
	if ( $position{$ref_id} ) {
		$inCoordFile++;
		for my $index ( keys %{ $position{$ref_id} } ) {
			my $bounds = $position{$ref_id}{$index};
			if ( !( $last_base < $bounds->[0] || $first_base > $bounds->[1] ) ) {
				my $subtract = 0;
				$hash{$ref_id}{$index} = 1;    #i have seen this reference id
				if ( $first_base < $bounds->[0] ) {
					$subtract = $bounds->[0] - $first_base;
				}
				if ( $last_base > $bounds->[1] ) {
					$subtract = $subtract + $last_base - $bounds->[1];
				}
				if ( $direction == 1 ) {
					if ( $base1{$ref_id}{$index} ) {
						$base1{$ref_id}{$index} += $length - $subtract;
						$base1++;
					}
					else {
						$base1{$ref_id}{$index} = $length - $subtract;
						$base1++;
					}
				}
				elsif ( $direction == 2 ) {
					if ( $base2{$ref_id}{$index} ) {
						$base2{$ref_id}{$index} += $length - $subtract;
						$base2++;
					}
					else {
						$base2{$ref_id}{$index} = $length - $subtract;
						$base2++;
					}
				}
				else {
					die "Direction should be 1 or 2. Internal unknown error.\n";
				}
				$seen_this_insert = 1;
				$total_length_matched += $length - $subtract;
			}
		}
	}
	else {
		$notInCoordFile++;
	}
	$prev_insert = $insert;
	$prev_read   = $read;

}
close IN;

#Process final block
print STDERR "CalculateTaxonomy :: Processing final block...\n";

#### last
if ( $seen_this_insert == 1 ) {
	$total_seen_inserts++;
}
my $count = scalar keys %hash;    #how many different ref_ids have I seen?
foreach my $ref_id ( keys %hash ) {
	for my $index ( keys %{ $hash{$ref_id} } ) {
		$count_sum = $count_sum + ( 1 / $count );
		@{ $position{$ref_id}{$index} }[2] += ( 1 / $count );

		#tax
		foreach my $i (@levels) {
			$ref_id =~ m/(\d+)\.*/;
			my $ref    = $map{$i}{$1};
			my $length = $len{$1};
			unless ( $taxInsert{$i}{$ref} ) {
				@{ $taxInsert{$i}{$ref} } = ( 0, 0, 0, 0 );
			}
			@{ $taxInsert{$i}{$ref} }[0] += 1 / $count;
			@{ $taxInsert{$i}{$ref} }[1] += 1 / $count / $length;
			if ( $count == 1 ) {
				@{ $taxInsert{$i}{$ref} }[2] += 1 / $count;
				@{ $taxInsert{$i}{$ref} }[3] += 1 / $count / $length;
			}

		}
	}
}
%hash = ();

$count = $base1;
foreach my $ref_id ( keys %base1 ) {
	for my $index ( keys %{ $base1{$ref_id} } ) {
		$total_bases_covered += ( 1 / $count ) * $base1{$ref_id}{$index};
		@{ $position{$ref_id}{$index} }[3] += ( 1 / $count ) * $base1{$ref_id}{$index};

		#tax
		foreach my $i (@levels) {
			$ref_id =~ m/(\d+)\.*/;
			my $ref    = $map{$i}{$1};
			my $length = $len{$1};
			unless ( $taxBase{$i}{$ref} ) {
				@{ $taxBase{$i}{$ref} } = ( 0, 0, 0, 0 );
			}
			@{ $taxBase{$i}{$ref} }[0] += 1 / $count * $base1{$ref_id}{$index};
			@{ $taxBase{$i}{$ref} }[1] += 1 / $count * $base1{$ref_id}{$index} / $length;
			if ( $count == 1 ) {
				@{ $taxBase{$i}{$ref} }[2] += 1 / $count * $base1{$ref_id}{$index};
				@{ $taxBase{$i}{$ref} }[3] += 1 / $count * $base1{$ref_id}{$index} / $length;
			}

		}
	}
}
%base1 = ();
$base1 = 0;

$count = $base2;
foreach my $ref_id ( keys %base2 ) {
	for my $index ( keys %{ $base2{$ref_id} } ) {
		$total_bases_covered += ( 1 / $count ) * $base2{$ref_id}{$index};
		@{ $position{$ref_id}{$index} }[3] += ( 1 / $count ) * $base2{$ref_id}{$index};

		#tax
		foreach my $i (@levels) {
			$ref_id =~ m/(\d+)\.*/;
			my $ref    = $map{$i}{$1};
			my $length = $len{$1};
			unless ( $taxBase{$i}{$ref} ) {
				@{ $taxBase{$i}{$ref} } = ( 0, 0, 0, 0 );
			}
			@{ $taxBase{$i}{$ref} }[0] += 1 / $count * $base2{$ref_id}{$index};
			@{ $taxBase{$i}{$ref} }[1] += 1 / $count * $base2{$ref_id}{$index} / $length;
			if ( $count == 1 ) {
				@{ $taxBase{$i}{$ref} }[2] += 1 / $count * $base2{$ref_id}{$index};
				@{ $taxBase{$i}{$ref} }[3] += 1 / $count * $base2{$ref_id}{$index} / $length;
			}

		}
	}
}
%base2 = ();
$base2 = 0;
#### last

print STDERR "CalculateTaxonomy :: Number of reads not counted, because their matching ref_id was not in the .coord file: $notInCoordFile\n";
print STDERR "CalculateTaxonomy :: Number of reads counted, with matching ref_id in the .coord file: $inCoordFile\n";
print STDERR "CalculateTaxonomy :: Done parsing. Writing stats...\n";
my $mapped_inserts = $total_seen_inserts;
print STDERR "CalculateTaxonomy :: total bases covered $total_bases_covered\n";
print STDERR "CalculateTaxonomy :: total insert count sum $count_sum\n";
print STDERR "CalculateTaxonomy :: mapped inserts $mapped_inserts\n";
print STDERR "CalculateTaxonomy :: total inserts $total_inserts\n";
print STDERR "CalculateTaxonomy :: stats file is $coverage_file\n";

my $not_mapped_bases;
my $not_mapped_inserts;
my $fraction_mapped_bases = 0;
if ( $total_bases > 0 ) {
	$fraction_mapped_bases = $total_bases_covered / $total_bases;
}

# bfr not correctly calculated. Also not obvious what it'd be useful for.
my $bfr = 0;
if ( $total_length_matched > 0 ) {
	$bfr = $total_bases_covered / $total_length_matched;
}
$not_mapped_bases   = $total_bases - $total_bases_covered;
$not_mapped_inserts = $total_inserts - $mapped_inserts;
print STDERR "CalculateTaxonomy :: not mapped inserts $not_mapped_inserts\n";

my $fr = 0;
if ( $total_inserts > 0 ) {
	$fr = $mapped_inserts / $total_inserts;
}
open STATS, '>', "$coverage_file" or die "ERROR & EXIT: Could not open $coverage_file for output.\n";
print STATS "total_inserts\tmapped_inserts\tfraction_mapped_inserts\ttotal_bases\tmapped_bases\tfraction_mapped_bases\tdb_average_entry_length\n";
print STATS "$total_inserts\t$mapped_inserts\t$fr\t$total_bases\t$total_bases_covered\t$fraction_mapped_bases\t$db_average_gene_length\n";
close STATS;

print STDERR "CalculateTaxonomy :: Print output...\n";
open OUT,      '>', "$cwd/$sample/insert.$output.insert.coverage.count" or die "ERROR & EXIT: Cannot print to output file $cwd/$sample/insert.$output.insert.coverage.count\n";
open BASE,     '>', "$cwd/$sample/base.$output.base.coverage.count"     or die "ERROR & EXIT: Cannot print to output file $cwd/$sample/insert.$output.base.coverage.count\n";
open NORM,     '>', "$cwd/$sample/insert.$output.insert.coverage.norm"  or die "ERROR & EXIT: Cannot print to output file $cwd/$sample/insert.$output.insert.coverage.norm\n";
open BASENORM, '>', "$cwd/$sample/base.$output.base.coverage.norm"      or die "ERROR & EXIT: Cannot print to output file $cwd/$sample/insert.$output.base.coverage.count\n";

my $printrownames = 0;
if ( $print_rownames eq 'yes' ) {
	$printrownames = 1;
	open HEADER,     '>', "$cwd/$sample/rownames.$output.rownames"      or die "ERROR & EXIT: Cannot print to output file $cwd/$sample/rownames.$output.rownames\n";
	open HEADERUNIQ, '>', "$cwd/$sample/rownames.$output.rownames.uniq" or die "ERROR & EXIT: Cannot print to output file $cwd/$sample/rownames.$output.rownames.uniq\n";
}
else {
	if ( $printOnlyIfSample0 == 0 ) {
		$printrownames = 1;
		open HEADER,     '>', "$print_rownames"      or die "ERROR & EXIT: Cannot print to output file $print_rownames\n";
		open HEADERUNIQ, '>', "$print_rownames.uniq" or die "ERROR & EXIT: Cannot print to output file $print_rownames.uniq\n";

	}
}

print OUT "$sample\n";
print BASE "$sample\n";
print BASENORM "$sample\n";
print NORM "$sample\n";

if ($printrownames) {
	print HEADER "ref_id\n";
	print HEADERUNIQ "ref_id\n";
	print HEADER "-1\n";
	print HEADERUNIQ "-1\n";
}
if ( !( -e "$insert_stats_file.stats" ) ) {
	print BASE "missing_insert_file_not_calculated\n";
}
else {
	print BASE "$not_mapped_bases\n";
}

my $fraction = $not_mapped_bases / $db_average_gene_length;

#print BASENORM "$fraction\n";
print BASENORM "NA\n";
if ( !( -e "$insert_stats_file.stats" ) ) {
	print OUT "missing_insert_file_not_calculated\n";
}
else {
	print OUT "$not_mapped_inserts\n";
}
$fraction = $not_mapped_inserts / $db_average_gene_length;

#print NORM "$not_mapped_inserts\n";
print NORM "NA\n";
foreach my $k ( sort keys %position ) {
	for my $bounds ( sort keys %{ $position{$k} } ) {
		my @bounds    = @{ $position{$k}{$bounds} };
		my $count     = $bounds[2];
		my $basecount = $bounds[3];
		my $length    = $bounds[1] - $bounds[0] + 1;
		my $norm      = $count / $length;
		my $basenorm  = $basecount / $length;
		print BASE "$basecount\n";
		print BASENORM "$basenorm\n";
		print OUT "$count\n";
		print NORM "$norm\n";

		if ($printrownames) {
			print HEADER "$k\n";
			print HEADERUNIQ "$k\_$bounds[0]\_$bounds[1]\n";
		}
	}
}
close OUT;
close NORM;
close HEADER;
close HEADERUNIQ;
close BASENORM;
close BASE;

system "ln -fs $print_rownames $cwd/$sample/base.$output.rownames";
system "ln -fs $print_rownames.uniq $cwd/$sample/base.$output.rownames.uniq";
system "ln -fs $print_rownames $cwd/$sample/insert.$output.rownames";
system "ln -fs $print_rownames.uniq $cwd/$sample/insert.$output.rownames.uniq";

# Print tax
print STDERR "CalculateTaxonomy :: Printing taxonomic abundances...\n";
foreach my $i (@levels) {
	open B,   ">$output_tax.base.raw.$i";
	open Bn,  ">$output_tax.base.norm.$i";
	open Bu,  ">$output_tax.base.only.unique.raw.$i";
	open Bun, ">$output_tax.base.only.unique.norm.$i";
	open I,   ">$output_tax.insert.raw.$i";
	open In,  ">$output_tax.insert.norm.$i";
	open Iu,  ">$output_tax.insert.only.unique.raw.$i";
	open Iun, ">$output_tax.insert.only.unique.norm.$i";

	print B "$sample\n";
	print Bn "$sample\n";
	print Bu "$sample\n";
	print Bun "$sample\n";
	print I "$sample\n";	
	print In "$sample\n";
	print Iu "$sample\n";
	print Iun "$sample\n";

	my %taxa;
	foreach my $ref ( keys %{ $map{$i} } ) {
		$taxa{ $map{$i}{$ref} } = 1;
	}
	if ( $printOnlyIfSample0 == 0 ) {
		open TAXHEADER, ">$print_rownamesTax.$i.rownames" or die "ERROR & EXIT: Cannot print to output file $print_rownamesTax.$i.rownames";
		print TAXHEADER "taxa\n-1\n";
	}

	if ( !( -e "$insert_stats_file.stats" ) ) {
		print B "missing_insert_file_not_calculated\n";
	}
	else {
		print B "$not_mapped_bases\n";
	}
	if ( !( -e "$insert_stats_file.stats" ) ) {
		print I "missing_insert_file_not_calculated\n";
	}
	else {
		print I "$not_mapped_inserts\n";
	}
	
	print Bn "NA\n";
	print Bu "NA\n";
	print Bun "NA\n";
	print In "NA\n";
	print Iu "NA\n";
	print Iun "NA\n";
		
	foreach my $taxa ( sort keys %taxa ) {

		unless ( $taxBase{$i}{$taxa} ) {
			@{ $taxBase{$i}{$taxa} } = ( 0, 0, 0, 0 );
		}
		unless ( $taxInsert{$i}{$taxa} ) {
			@{ $taxInsert{$i}{$taxa} } = ( 0, 0, 0, 0 );
		}

		if ( $printOnlyIfSample0 == 0 ) {
			print TAXHEADER "$taxa\n";
		}

		print B "$taxBase{$i}{$taxa}[0]\n";
		print Bn "$taxBase{$i}{$taxa}[1]\n";
		print Bu "$taxBase{$i}{$taxa}[2]\n";
		print Bun "$taxBase{$i}{$taxa}[3]\n";
		print I "$taxInsert{$i}{$taxa}[0]\n";
		print In "$taxInsert{$i}{$taxa}[1]\n";
		print Iu "$taxInsert{$i}{$taxa}[2]\n";
		print Iun "$taxInsert{$i}{$taxa}[3]\n";
	}
	close TAXHEADER;
	close B;
	close Bn;
	close Bu;
	close Bun;
	close I;
	close In;
	close Iu;
	close Iun;
	system "ln -sf $print_rownamesTax.$i.rownames $output_tax.$i.rownames";

}

print STDERR "CalculateTaxonomy :: Done.\n";

if ( $missingFile eq "missing_insert_file" ) {
	die "ERROR & EXIT: INSERT FILE WAS MISSING. -1 FRACTION NOT CALCULATED, BUT OTHER VALUES CORRECT. CALCULATIONS HAVE FINISHED.\n";
}
elsif ( $not_mapped_inserts < 0 ) {
	die "ERROR & EXIT: NUMBER OF MAPPED INSERTS IS $not_mapped_inserts. Something is wrong...\n";
}

exit 0;
