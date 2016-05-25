#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.
my ( @allFiles, $use3files, %saved_qual_cutoff, %saved_first_cycle, $length_cutoff, $qual_cutoff, $paired_end_data, $use_5prime_file, $trim_5prime_end, $src_dir, $bin_dir, $solexaqa_or_fastx, $cwd, $temp_dir, $file_formats_array, $ZCAT, $sample, $trim_5_prime_end );
my $TOTmax            = 0;
my $TOTlength         = 0;
my $TOTsum            = 0;
my $TOTinserts        = 0;
my $fastq_trim_filter = "fastq_trim_filter_v5_EMBL";
my $format_counter    = 0;
GetOptions(
	'length_cutoff=s'      => \$length_cutoff,
	'qual_cutoff=s'        => \$qual_cutoff,
	'paired_end_data=s'    => \$paired_end_data,
	'use_5prime_file=s'    => \$use_5prime_file,
	'trim_5prime_end=s'    => \$trim_5_prime_end,
	'src_dir=s'            => \$src_dir,
	'bin_dir=s'            => \$bin_dir,
	'solexaqa_or_fastx=s'  => \$solexaqa_or_fastx,
	'cwd=s'                => \$cwd,
	'temp_dir=s'           => \$temp_dir,
	'file_formats_array=s' => \$file_formats_array,
	'zcat=s'               => \$ZCAT,
	'sample=s'             => \$sample,
	'use3files=s'          => \$use3files
);

# If running OSX, change binary file
chomp( my $systemType = `uname -s` );
if ( $systemType =~ m/Darwin/ ) {
	$fastq_trim_filter = $fastq_trim_filter . "_OSX";
}

my $ending;
# Get fq or fq.gz files
if ( $use3files == 1 ) {
	@allFiles = `ls -1 $cwd/$sample/*fq $cwd/$sample/*fq.gz 2>/dev/null | grep -v 'trimmed.filtered'  | grep -P '\.single\.|\.pair\.'`;
# this will only be used in case of use3files=1
$ending = $allFiles[0];
$ending =~ s/.*.pair.1.//;
}
else {
	@allFiles = `ls -1 $cwd/$sample/*fq $cwd/$sample/*fq.gz 2>/dev/null | grep -v 'trimmed.filtered' | grep -v '.single.' | grep -v '.pair.'`;
}

for my $i ( 0 .. scalar @allFiles - 1 ) {
	chomp( $allFiles[$i] );
	my @J = split /\//, $allFiles[$i];
	$allFiles[$i] = $J[-1];
}


# Delete old files before starting
print "ReadTrimFilter_aux : SAMPLE $sample\n";
print "ReadTrimFilter_aux : SCRIPT VERSION 4\n";
print "ReadTrimFilter_aux : use3files=$use3files\n";
print "ReadTrimFilter_aux : CLEANUP : rm -f $temp_dir/$sample/temp/*.trimmed.filtered.fq*; mkdir -p $cwd/$sample/stats; rm -fr $cwd/$sample/reads.processed.$solexaqa_or_fastx; mkdir -p $cwd/$sample/reads.processed.$solexaqa_or_fastx\n";

system "rm -f $temp_dir/$sample/temp/*.trimmed.filtered.fq*";
system "mkdir -p $cwd/$sample/stats";
system "rm -fr $cwd/$sample/reads.processed.$solexaqa_or_fastx";
system "mkdir -p $cwd/$sample/reads.processed.$solexaqa_or_fastx";

# get quality statistics
foreach my $j (@allFiles) {
	( my $format, my $sanger ) = get_format();
	my $zip = "cat";
	if ( $j =~ /\.gz$/ )  { $zip = "gunzip -c" }
	if ( $j =~ /\.bz2$/ ) { $zip = "bunzip2 -c" }
	print "ReadTrimFilter_aux : GET QUALITY STATS : mode=$solexaqa_or_fastx sanger=$sanger format=$format sample=$j : EXECUTE $zip $cwd/$sample/$j | $bin_dir/fastx_quality_stats -o $temp_dir/$sample/temp/$j.qual_stats.temp $sanger\n";
	system "$zip $cwd/$sample/$j | $bin_dir/fastx_quality_stats -o $temp_dir/$sample/temp/$j.qual_stats.temp $sanger";
	print "ReadTrimFilter_aux : Get 3 prime end\n";
	$saved_qual_cutoff{"$sample$j"} = get_3prime_end( $sample, $j );
}

# Process all files
$format_counter = -1;
foreach my $j (@allFiles) {
	( my $format, my $sanger ) = get_format();
	my $first_cycle;
	my $qual_cutoff = $saved_qual_cutoff{"$sample$j"};
	if ( $use_5prime_file eq 'yes' ) {
		$first_cycle = get_5prime_from_file( $sample, $j );
		$saved_first_cycle{"$sample$j"} = $first_cycle;
	}
	if ( $use_5prime_file eq 'no' ) {
		$first_cycle = get_quality_statistics( $sample, $j );
		$saved_first_cycle{"$sample$j"} = $first_cycle;
		open X, '>>', "$cwd/MOCAT.cutoff5prime_calculated" or die "ERROR & EXIT: Cannot write to $cwd/MOCAT.cutoff5prime_calculated";
		print X "$sample\t$j\t$first_cycle\n";
		close X;
	}
	if ( $trim_5_prime_end eq "no" ) {
		$first_cycle = 1;
		$saved_first_cycle{"$sample$j"} = $first_cycle;
	}

	# Process non pair-end files here
	if ( $paired_end_data eq "no" && $use3files == 0 ) {
		my $basename = $j;
		$basename =~ s/.fq.bz2$//;
		$basename =~ s/.fq.gz$//;
		$basename =~ s/.fq$//;

		my $cmd = "$bin_dir/$fastq_trim_filter -m $solexaqa_or_fastx -a $cwd/$sample/$j -f $first_cycle -q $qual_cutoff -l $length_cutoff -p 50 $sanger -o $cwd/$sample/reads.processed.$solexaqa_or_fastx/$basename";

		print "ReadTrimFilter_aux : PROCESS SINGLE-END $j : EXECUTE $cmd\n";
		my @out = `$cmd`;
		unless ( scalar @out == 5 ) {
			die "SAMPLE ERROR & EXIT: Sample seem to not have any reads passing read trim filter. Something wrong? Perhaps the cut off value is too high?";
		}
		print "ReadTrimFilter_aux : $fastq_trim_filter OUTPUT: ";
		print $out[0];
		shift @out;
		chomp( my $sum     = $out[0] );
		chomp( my $len     = $out[1] );
		chomp( my $max     = $out[2] );
		chomp( my $inserts = $out[3] );
		$TOTinserts = $TOTinserts + $inserts;
		$TOTsum     = $TOTsum + $sum;
		$TOTlength  = $TOTlength + $len;

		if ( $TOTmax < $max ) {
			$TOTmax = $max;
		}
	}
}

# Post process pair-end files
if ( $paired_end_data eq "yes" && $use3files == 0 ) {
	for ( my $count = 0 ; $count <= scalar @allFiles - 1 ; $count = $count + 2 ) {
		$allFiles[$count] =~ m/^(.+)\.1.fq/;
		my $base_name = $1;
		$format_counter = $count;
		( my $format, my $sanger ) = get_format();

		my $cmd = "$bin_dir/$fastq_trim_filter -m $solexaqa_or_fastx -a $cwd/$sample/" . $allFiles[$count] . " -b $cwd/$sample/" . $allFiles[ $count + 1 ] . " -f " . $saved_first_cycle{ $sample . $allFiles[$count] } . " -2 " . $saved_first_cycle{ $sample . $allFiles[ $count + 1 ] } . " -q $qual_cutoff -l $length_cutoff -p 50 $sanger -o $cwd/$sample/reads.processed.$solexaqa_or_fastx/$base_name";

		print "ReadTrimFilter_aux : PROCESS PAIR-END $base_name : EXECUTE $cmd\n";
		my @out = `$cmd`;
		unless ( scalar @out == 5 ) {
			die "SAMPLE ERROR & EXIT: Sample seem to not have any reads passing read trim filter. Something wrong? Perhaps the cut off value is too high?";
		}
		shift @out;
		chomp( my $sum     = $out[0] );
		chomp( my $len     = $out[1] );
		chomp( my $max     = $out[2] );
		chomp( my $inserts = $out[3] );
		$TOTinserts = $TOTinserts + $inserts;
		$TOTsum     = $TOTsum + $sum;
		$TOTlength  = $TOTlength + $len;

		if ( $TOTmax < $max ) {
			$TOTmax = $max;
		}
	}
}

# new, support for paired end and 3 files
if ( $paired_end_data eq "yes" && $use3files == 1 ) {

	# Both of these checks ahve been done in the .pm file, but why not redo them here...
	unless ( scalar @allFiles % 3 == 0 ) {
		die "ERROR & EXIT: $sample seems to have .pair. and .single. files, but the number of files is not exactly 3 (number is " . ( scalar @allFiles ) . "), please check that no files are mising";
	}
	chomp( my @beginning = `ls -1 $cwd/$sample/*.pair.1.*fq.gz $cwd/$sample/*.pair.1.*.fq $cwd/$sample/*.pair.1.*fq.bz2 2>/dev/null | sed 's/.pair.1.*//'` );
	chomp( my $end       = `ls -1 $cwd/$sample/*.pair.1.*fq.gz $cwd/$sample/*.pair.1.*.fq $cwd/$sample/*.pair.1.*fq.bz2 2>/dev/null | sed 's/.*.pair.1.//' | sort -u` );
	foreach my $beginning (@beginning) {
		unless ( -e "$beginning.single.$end" ) {
			die "ERROR & EXIT: Expected $beginning.single.$end to exist. Do the files have different ending for the sample $sample?";
		}
		unless ( -e "$beginning.pair.1.$end" ) {
			die "ERROR & EXIT: Expected $beginning.single.$end to exist. Do the files have different ending for the sample $sample?";
		}
		unless ( -e "$beginning.pair.2.$end" ) {
			die "ERROR & EXIT: Expected $beginning.single.$end to exist. Do the files have different ending for the sample $sample?";
		}
	}

	my $format_counter = -3;
	foreach my $base_name (@beginning) {
		$base_name =~ s/$cwd\/$sample\///;

		# they should all have the same format within a triple;
		$format_counter = $format_counter + 3;
		( my $format, my $sanger ) = get_format();

		# the order of the files should always be .pair.1 .pair.2 .single
		my $cmd = "$bin_dir/$fastq_trim_filter -m $solexaqa_or_fastx -a $cwd/$sample/$base_name.pair.1.$ending -b $cwd/$sample/$base_name.pair.2.$ending -c $cwd/$sample/$base_name.single.$ending -f " . $saved_first_cycle{"$sample$allFiles[0]"} . " -2 " . $saved_first_cycle{"$sample$allFiles[0]"} . " -3 " . $saved_first_cycle{"$sample$allFiles[0]"} . " -q $qual_cutoff -l $length_cutoff -p 50 $sanger -o $cwd/$sample/reads.processed.$solexaqa_or_fastx/$base_name";

		print "ReadTrimFilter_aux : PROCESS PAIR-END $base_name : EXECUTE $cmd\n";
		my @out = `$cmd`;
		unless ( scalar @out == 5 ) {
			die "SAMPLE ERROR & EXIT: Sample seem to not have any reads passing read trim filter. Something wrong? Perhaps the cut off value is too high?";
		}
		shift @out;
		chomp( my $sum     = $out[0] );
		chomp( my $len     = $out[1] );
		chomp( my $max     = $out[2] );
		chomp( my $inserts = $out[3] );
		$TOTinserts = $TOTinserts + $inserts;
		$TOTsum     = $TOTsum + $sum;
		$TOTlength  = $TOTlength + $len;

		if ( $TOTmax < $max ) {
			$TOTmax = $max;
		}
	}

}

# Nice cleanup of files
print "ReadTrimFilter_aux : SUMMARIZE STATS and write $cwd/$sample/stats/$sample.readtrimfilter.$solexaqa_or_fastx.stats\n";
system "rm -f $temp_dir/$sample/temp/*qual_stats.temp; mv $temp_dir/$sample/temp/*qual_stats $cwd/$sample/reads.processed.$solexaqa_or_fastx;";

# Stats section
my $max     = $TOTmax;
my $length  = $TOTlength;
my $inserts = $TOTinserts;
open STAT, '>', "$cwd/$sample/stats/$sample.readtrimfilter.$solexaqa_or_fastx.stats" or die "ERROR & EXIT: Cannot write to $cwd/$sample/stats/$sample.readtrimfilter.$solexaqa_or_fastx.stats";
if ( $TOTsum == 0 ) {
	die "SAMPLE ERROR & EXIT: Sample $sample seem to not have any reads passing read trim filter. Something wrong?";
}
my $avg = int( ( $TOTlength / $TOTsum ) + 0.5 );
my $kmer;
if ( $avg % 2 == 0 ) {

	#avg_readlen is even
	$kmer = $avg / 2;
	if ( $kmer % 2 == 0 ) {
		$kmer = $kmer + 1;
	}
	else {
		$kmer = $kmer + 2;
	}
}
else {

	#avg_readlen is odd
	$kmer = ( $avg + 1 ) / 2;
	if ( $kmer % 2 == 1 ) {
		$kmer = $kmer + 0;
	}
	else {
		$kmer = $kmer + 1;
	}
}
print STAT "Reads\tBases\tMax\tAvg\tKmer\tInserts\n$TOTsum\t$TOTlength\t$TOTmax\t$avg\t$kmer\t$TOTinserts\n";
close STAT;
print "ReadTrimFilter_aux : COMPLETED READ TRIM FILTER\n";
exit 0;

#-----------
#SUBROUTINES
#-----------
sub get_5prime_from_file {
	my $sample      = $_[0];
	my $j           = $_[1];
	my $first_cycle = "";
	open FPI, '<', "$cwd/MOCAT.cutoff5prime" or die "ERROR &EXIT: Cannot read $cwd/MOCAT.cutoff5prime";
	while (<FPI>) {
		chomp;
		my @line = split;
		if (
			$line[0] eq $sample
			&& (         $line[1] eq $j
				|| "$line[1].fq"     eq $j
				|| "$line[1].fq.gz"  eq $j
				|| "$line[1].fq.bz2" eq $j )
		  )
		{
			$first_cycle = $line[2];
		}
	}
	close FPI;
	if ( $first_cycle eq "" ) {
		die "\nERROR and EXIT: MOCAT.cutoff5prime file incorrect. Last checked sample: $sample. Last checked lane: $j. File should be in format <sample>\\t<lane as complete file name>\\t<base>";
	}
	return $first_cycle;
}

sub get_3prime_end {
	my $sample = $_[0];
	my $j      = $_[1];
	open IN_TEMP,  '<', "$temp_dir/$sample/temp/$j.qual_stats.temp" or die "ERROR & EXIT: Cannot read $temp_dir/$sample/temp/$j.qual_stats.temp";
	open OUT_TEMP, '>', "$temp_dir/$sample/temp/$j.qual_stats"      or die "ERRROR & EXIT: Cannot write to $temp_dir/$sample/temp/$j.qual_stats";
	my ( $raw_read_num, $raw_read_base ) = ( 0, 0 );
	<IN_TEMP>;
	while (<IN_TEMP>) {

		#        unless ($. == 1){
		my @line   = split "\t", $_;
		my $cycle  = $line[0];
		my $mean   = $line[5];
		my $median = $line[7];
		if ( $mean >= $qual_cutoff && $median >= $qual_cutoff ) {
			print OUT_TEMP $_;
		}
		$raw_read_num = $line[1] if ( $raw_read_num < $line[1] );
		$raw_read_base += $line[1];

		#        }
	}
	close IN_TEMP;
	close OUT_TEMP;
	open RAW, "> $cwd/$sample/stats/$j.raw.reads.stats" || die "ERROR and EXIT: Cannot open $cwd/$sample/stats/$j.raw.reads.stats\n";
	print RAW "Reads\tBases\n", $raw_read_num, "\t", $raw_read_base, "\n";
	close RAW;
	return $qual_cutoff;
}

sub get_quality_statistics {
	my $sample = $_[0];
	my $j      = $_[1];
	my @line   = my @cycle_array = my @cycle_number = ();
	my ( %A, %C, %G, %T );
	open IN, '<', "$temp_dir/$sample/temp/$j.qual_stats"
	  || die "ERROR and EXIT: Cannot open $temp_dir/$sample/temp/$j.qual_stats\n";
	while (<IN>) {
		@line = split "\t", $_;
		$A{ $line[0] } = $line[12];
		$C{ $line[0] } = $line[13];
		$G{ $line[0] } = $line[14];
		$T{ $line[0] } = $line[15];
		push( @cycle_number, $line[0] );
	}
	close IN;
	my $read_length = scalar @cycle_number;
	my ( $max_A, $min_A );
	my ( $max_C, $min_C );
	my ( $max_G, $min_G );
	my ( $max_T, $min_T );
	my $warning;
	( $max_A, $min_A, $warning ) = &get_base_range( \%A, $read_length, "$sample/$j" );
	( $max_C, $min_C, $warning ) = &get_base_range( \%C, $read_length, "$sample/$j" );
	( $max_G, $min_G, $warning ) = &get_base_range( \%G, $read_length, "$sample/$j" );
	( $max_T, $min_T, $warning ) = &get_base_range( \%T, $read_length, "$sample/$j" );
	push( @cycle_array, '0trim' );

	if ($warning) {
		print STDERR $warning;
	}

	foreach my $m (@cycle_number) {
		if (         $A{$m} > $max_A
			|| $A{$m} < $min_A
			|| $C{$m} > $max_C
			|| $C{$m} < $min_C
			|| $G{$m} > $max_G
			|| $G{$m} < $min_G
			|| $T{$m} > $max_T
			|| $T{$m} < $min_T )
		{
			my $var = $m . "trim";
			push @cycle_array, $var;
		}
		else {
			my $var = $m . "keep";
			push @cycle_array, $var;
		}
	}
	my $first = 1;
	my $tempstring = join( " ", @cycle_array );
	if ( $tempstring =~ m/^(\d+trim )+(\d+)keep/ ) {
		$first = $2;
	}
	return $first;
}

sub get_base_range {
	my ( $base, $read_length, $file ) = @_;
	my ( $sum_base, $avg_base, $sum_of_squares_base, $var_base, $stdev_base, $warning );
	$sum_base = 0;
	foreach my $i ( keys %{$base} ) {
		$sum_base += $$base{$i};
	}
	if ( $sum_base == 0 ) {
		$warning             = "WARNING: $file had 0 reads passing 5' quality statistics. This means for this lane the 5' cut off defaults to base 1 (no 5' filtering).\n";
		$avg_base            = 0;
		$sum_of_squares_base = 0;
	}
	else {
		$avg_base = $sum_base / $read_length;
		foreach my $i ( keys %{$base} ) {
			$sum_of_squares_base += ( $$base{$i} - $avg_base )**2;
		}
	}
	if ( $read_length == 1 ) {
		$stdev_base = 0;
	}
	else {
		$var_base = $sum_of_squares_base / ( $read_length - 1 );
		$stdev_base = sqrt($var_base);
	}
	return ( $avg_base + ( $stdev_base * 2 ), $avg_base - ( $stdev_base * 2 ), $warning );
}

sub get_format {
	my $format = substr $file_formats_array, $format_counter, 1;
	my $sanger;
	if ( $format eq 's' ) { $format = 'sanger'; $sanger = "-Q 33" }
	if ( $format eq 'x' ) { $format = 'solexa'; $sanger = "" }
	if ( $format eq 'i' ) { $format = 'ia';     $sanger = "" }
	$format_counter++;
	return ( $format, $sanger );
}
