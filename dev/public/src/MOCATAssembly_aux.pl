#!/usr/bin/env perl
use strict;
use warnings;
use MOCATCore;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

my @prefix;
my @sample;
my $folder              = $ARGV[0];
my $db                  = $ARGV[1];
my $max                 = $ARGV[2];
my $bin_dir             = $ARGV[3];
my $kmer                = $ARGV[4];
my $cwd                 = $ARGV[5];
my $sample_file         = $ARGV[6];
my $temp_dir            = $ARGV[7];
my $mode                = $ARGV[8];
my $process_folder      = $ARGV[9];
my $reads               = $ARGV[10];
my $calc_insert_size_by = $ARGV[11];
my $scr_folder          = $ARGV[12];
my $bin_folder          = $ARGV[13];
my $denovo              = $ARGV[14];
my $sample              = $ARGV[15];
my $processors          = $ARGV[16];
my $lengthcutoff        = $ARGV[17];
my $ZCAT                = $ARGV[18];
my $pair_end            = $ARGV[19];
my $LOG                 = $ARGV[20];
my $min_size            = 1;
my $max_size            = 1000;
my $sample_size         = 10000000;    #number of sequences used for estimation of insert size
my $file_name;
my $prefix;
my %info;
my %max_readlen;
my @readlen;
my %pairPath;
my %singlePath;

my @file;
my @F2 = split ',', $folder;
foreach my $folder (@F2) {
	push @file, <$folder/*pair*fq.gz $folder/*single*fq.gz>;
	if ( scalar @file == 0 ) {
		die "ERROR & EXIT: Found no *pair* files in $folder\assembly_db_n";
	}
}

print "<<< RUNNING prepareAssembly >>>\n";
print "rm -f $cwd/$sample/stats/$sample\.assembly.$reads.$mode.K$kmer.inserts.stats\n";
system "rm -f $cwd/$sample/stats/$sample\.assembly.$reads.$mode.K$kmer.inserts.stats";

if ( $pair_end eq 'no' ) {
	my %prefix;
	print "    <<< single read mode >>>\n";
	foreach my $f (@file) {
		print "        <<< Prepare data for $f >>>\n";
		my @line      = split "\/", $f;
		my $file_name = $line[-1];
		my $lane      = $file_name;
		$lane =~ s/\.pair.*//;
		$lane =~ s/\.single.*//;
		my $pre = $f;
		if ( $pre =~ m/\.pair\./ ) {
			$pre =~ s/\.pair.*/\.pair/;
			$prefix{$pre} = 1;
		}
		if ( $f =~ /1\.fq/ ) {
			$pairPath{$sample}{$lane} = $f;
		}
		elsif ( $f =~ /single/ ) {
			$singlePath{$sample}{$lane} = $f;
		}
	}
	foreach my $prefix ( keys %prefix ) {
		print "       <<< Processing $prefix >>>\n";
		$file_name = $prefix;
		my @line = split /\//, $file_name;
		$line[-1] =~ m/(.+)\.pair/;
		my $lane       = $1;
		my $lines      = $sample_size * 4;
		#my $file_name2 = $file_name;
		#$file_name2 =~ s/$cwd/$temp_dir/;
		#$file_name2 =~ s/$sample/$sample\/temp/;
		#$file_name2 =~ s/$process_folder.$mode//;
		$info{$lane} = 1;
	}
}

if ( $pair_end eq 'yes' ) {

	if ( $calc_insert_size_by eq "mapping" ) {
		my %prefix;
		print "    <<< mapping mode >>>\n";
		foreach my $f (@file) {
			print "        <<< Prepare data for $f >>>\n";
			my @line      = split "\/", $f;
			my $file_name = $line[-1];
			my $lane      = $file_name;
			$lane =~ s/\.pair.*//;
			$lane =~ s/\.single.*//;
			my $pre = $f;
			if ( $pre =~ m/\.pair\./ ) {
				$pre =~ s/\.pair.*/\.pair/;
				$prefix{$pre} = 1;
			}
			if ( $f =~ /1\.fq/ ) {
				$pairPath{$sample}{$lane} = $f;
			}
			elsif ( $f =~ /single/ ) {
				$singlePath{$sample}{$lane} = $f;
			}
		}

		print "*** LANES ***\n";
		foreach my $prefix ( keys %prefix ) {
			print "$prefix\n";
		}
		print "*** LANES ***\n";

		foreach my $prefix ( keys %prefix ) {
			print "       <<< Processing $prefix >>>\n";
			$file_name = $prefix;
			my @line = split /\//, $file_name;
			$line[-1] =~ m/(.+)\.pair/;
			my $lane       = $1;
			my $lines      = $sample_size * 4;
			my $file_name2 = "$temp_dir/$sample/temp/$line[-1]";
			print "$ZCAT $file_name.1.fq.gz | head -n $lines > $file_name2.1.temp.fq\n";
			system "$ZCAT $file_name.1.fq.gz | head -n $lines > $file_name2.1.temp.fq";
			print "$ZCAT $file_name.2.fq.gz | head -n $lines > $file_name2.2.temp.fq\n";
			system "$ZCAT $file_name.2.fq.gz | head -n $lines > $file_name2.2.temp.fq";
			print "$bin_dir/soap2.21 -a $file_name2.1.temp.fq -b $file_name2.2.temp.fq -D $db -m $min_size -x $max_size -p $processors -o $file_name2.temp.out -2 $file_name.temp2.unpaired $LOG\n";
			system "$bin_dir/soap2.21 -a $file_name2.1.temp.fq -b $file_name2.2.temp.fq -D $db -m $min_size -x $max_size -p $processors -o $file_name2.temp.out -2 $file_name2.temp.unpaired $LOG";

			unless ( -e "$file_name2.temp.out" ) {
				die "ERROR & EXIT: Missing file $file_name2.temp.out. It seems like this command wasn't executed properly: $bin_dir/soap2.21 -a $file_name2.1.temp.fq -b $file_name2.2.temp.fq -D $db -m $min_size -x $max_size -p $processors -o $file_name2.temp.out -2 $file_name.temp2.unpaired";
			}
			my $filesize = -s "$file_name2.temp.out";
			if ( $filesize eq "0" ) {
				system "rm $file_name2*temp*";
				die "ERROR & EXIT: Insert sizes could not be calculated. Most likely the number of processed reads is too low.\nThere were not enough reads mapping to the mapping db to estimate the distance between two paired reads.\nTo solve this, calculate insert sizes differently.";
			}
			open( FH, "<$file_name2.temp.out" );
			my ( $a, $b, $nextline, $size, $sum, $total );
			while (<FH>) {
				my @line = split "\t", $_;
				if ( $line[6] eq "+" ) {
					$a = $line[8];
				}
				else {
					$a = $line[8] + $line[5];
				}
				$nextline = <FH>;
				my @nextline = split "\t", $nextline;
				if ( $nextline[6] eq "+" ) {
					$b = $nextline[8];
				}
				else {
					$b = $nextline[8] + $nextline[5];
				}
				$size = int( sqrt( ( $a - $b ) * ( $a - $b ) ) );
				$sum++;
				$total += $size;
			}
			close FH;
			$info{$lane} = int( $total / $sum );
			open OUTx, '>>', "$cwd/$sample/stats/$sample\.assembly.$reads.$mode.K$kmer.inserts.stats" or die "ERROR & EXIT: Cannot write to $cwd/$sample/stats/$sample\.assembly.$reads.$mode.K$kmer.inserts.stat";
			print OUTx $lane . "\t" . int( $total / $sum ) . "\n";
			print "Wrote " . $lane . "\t" . int( $total / $sum ) . " to $cwd/$sample/stats/$sample\.assembly.$reads.$mode.K$kmer.inserts.stats\n";
			close OUTx;
			system "rm $file_name2*temp*";
		}
	}

	if ( $calc_insert_size_by eq "assembly" ) {
		print "    <<< assembly mode >>>\n";
		my $outDir  = "$temp_dir/$sample/temp";
		my $name    = "$sample.$mode.K$kmer";
		my $lcutoff = $lengthcutoff;
		`mkdir -p $outDir/index` unless ( -e "$outDir/index" );
		`mkdir -p $outDir/bwa`   unless ( -e "$outDir/bwa" );
		my %lanePath;
		foreach my $f (@file) {
			print "        <<< Prepare data for $f >>>\n";
			$f =~ m/(.+\.(pair|single))/;
			my $prefix = $1;
			push( @prefix, $prefix );
			my @line      = split "\/", $f;
			my $file_name = $line[-1];
			my $lane      = $file_name;
			$lane =~ s/\.pair.*//;
			$lane =~ s/\.single.*//;

			if ( $f =~ /1\.fq/ ) {
				$lanePath{$lane} = $f;
				$pairPath{$sample}{$lane} = $f;
			}
			elsif ( $f =~ /single/ ) {
				$singlePath{$sample}{$lane} = $f;
			}
		}

		my %max;
		my %average;

		my $statsFile;
		if ( $reads eq "reads.processed" ) {
			$statsFile = "$cwd/$sample/stats/$sample.readtrimfilter.$mode.stats";
		}
		else {
			$statsFile = "$cwd/$sample/stats/$sample.screened.$reads.$mode.stats";
		}
		if ( -s "$statsFile" ) {
			open IAL, '<', "$statsFile";
		}
		else {
			die "\nERROR & EXIT: Missing insert file $statsFile.\nHave you specified the correct databse using --assembly 'DATABASE' or 'reads.processed' or 'fastafile'?";
		}

		my $line = <IAL>;
		$line = <IAL>;
		chomp($line);
		my @temp = split( "\t", $line );
		$max{$sample}     = $temp[2];
		$average{$sample} = $temp[3];
		close IAL;

		my %insert;
		unless ( -s "$cwd/MOCAT.preliminary_insert_sizes" ) {
			print "ERROR & EXIT: Missing manually entered preliminary insert sizes in file '$cwd/MOCAT.preliminary_insert_sizes'.\nEither specify assembly_calculate_insert_sizes to be 'mapping' in the config file, or, for each sample specify <SAMPLE><TAB><LANE><TAB><insert size> in the file '$cwd/MOCAT.preliminary_insert_sizes'.";
			die "ERROR & EXIT: Missing manually entered preliminary insert sizes in file '$cwd/MOCAT.preliminary_insert_sizes'.\nEither specify assembly_calculate_insert_sizes to be 'mapping' in the config file, or, for each sample specify <SAMPLE><TAB><LANE><TAB><insert size> in the file '$cwd/MOCAT.preliminary_insert_sizes'.";
		}
		open( II, '<', "$cwd/MOCAT.preliminary_insert_sizes" );
		while (<II>) {
			chomp;
			@temp = split;
			$insert{ $temp[0] }{ $temp[1] } = $temp[2];
		}
		close II;

		`mkdir -p $outDir/$sample` unless ( -e "$outDir/$sample" );
		print "        <<< Create config >>>\n";
		die "average reads length is not exist: $sample\n" unless ( exists $average{$sample} );
		die "max reads length is not exist: $sample\n"     unless ( exists $max{$sample} );
		my $kmer;
		if ( $average{$sample} % 2 == 0 ) {
			$kmer = $average{$sample} / 2;
			if ( $kmer % 2 == 0 ) {
				$kmer += 1;
			}
			else {
				$kmer += 2;
				print OUT "\n\n";

			}
		}
		else {
			$kmer = ( $average{$sample} + 1 ) / 2;
			$kmer += 1 if ( $kmer % 2 == 0 );
		}

		open OC, ">$outDir/$sample/$sample.config" or die "can\'t open config file: $outDir/$sample/$sample.config\n";
		print OC "max_rd_len=$max{$sample}\n";
		foreach my $date ( sort keys %{ $pairPath{$sample} } ) {
			my $date_noscreen = $date;
			$date_noscreen =~ s/\.screened.*//;
			unless ( exists $insert{$sample}{$date_noscreen} ) {
				print "ERROR & EXIT: Missing insert size for sample: $sample and lane: $date_noscreen\nPlease add this information to the file MOCAT.preliminary_insert_sizes.";
				die "ERROR & EXIT: Missing insert size for sample: $sample and lane: $date_noscreen\nPlease add this information to the file MOCAT.preliminary_insert_sizes.";
			}
			my $size = $insert{$sample}{$date_noscreen};
			print OC "[LIB]\n";
			print OC "avg_ins=$size\n";
			print OC "asm_flags=3\n";
			if ( $size >= 800 ) {
				print OC "rank=2\n";
			}
			else {
				print OC "rank=1\n";
			}
			my $reads2 = $pairPath{$sample}{$date};
			$reads2 =~ s/1\.fq/2\.fq/;
			print OC "q1=$pairPath{$sample}{$date}\n";
			print OC "q2=$reads2\n";
			if ( exists $singlePath{$sample}{$date} ) {
				chomp(my $NOL = `$ZCAT $singlePath{$sample}{$date} | wc -l`);
				if ($NOL > 0 ){
					print OC "[LIB]\n";
					print OC "asm_flags=1\n";
					print OC "q=$singlePath{$sample}{$date}\n";
				}
			}
		}
		close OC;
		`mkdir -p $outDir/$sample/K$kmer` unless ( -e "$outDir/$sample/K$kmer" );

		# bit tricky with the -D option. In the end Jens decided to include it in both steps, assuming this will work properly...
		# new in version 1.4.x, I removed the -D 1 option
		print "$denovo pregraph -s $outDir/$sample/$sample.config -K $kmer -o $outDir/$sample/K$kmer/$sample\_K$kmer\_L$max{$sample} > $outDir/$sample/K$kmer/$sample\_K$kmer\_L$max{$sample}.pregraph.log && $denovo contig -D 1 -M 3 -g $outDir/$sample/K$kmer/$sample\_K$kmer\_L$max{$sample} >$outDir/$sample/K$kmer/$sample\_K$kmer\_L$max{$sample}.contig.log\n";
		system "$denovo pregraph -s $outDir/$sample/$sample.config -K $kmer -o $outDir/$sample/K$kmer/$sample\_K$kmer\_L$max{$sample} > $outDir/$sample/K$kmer/$sample\_K$kmer\_L$max{$sample}.pregraph.log && $denovo contig -D 1 -M 3 -g $outDir/$sample/K$kmer/$sample\_K$kmer\_L$max{$sample} >$outDir/$sample/K$kmer/$sample\_K$kmer\_L$max{$sample}.contig.log";
		print "perl $scr_folder/MOCATAssembly_cut_more.pl $outDir/$sample/K$kmer/$sample\_K$kmer\_L$max{$sample}.contig $lcutoff $outDir/index/$name\n";
		`perl $scr_folder/MOCATAssembly_cut_more.pl $outDir/$sample/K$kmer/$sample\_K$kmer\_L$max{$sample}.contig $lcutoff $outDir/index/$name`;
		my $contig = "$outDir/index/$name.more$lcutoff";
		print "$bin_folder/bwa index $contig\n";
		`$bin_folder/bwa index $contig`;
		open OUT2, '>>', "$cwd/$sample/stats/$sample\.assembly.$reads.$mode.K$kmer.inserts.stats";

		foreach my $date ( keys %lanePath ) {
			my $reads1 = $lanePath{$date};
			my $sai1   = ( split /\//, $reads1 )[-1];
			my $reads2 = $lanePath{$date};
			$reads2 =~ s/1\.fq/2\.fq/;
			die "reads2 is not exist: $reads2\n" unless ( -e $reads2 );
			my $sai2 = ( split /\//, $reads2 )[-1];
			print "$bin_folder/bwa aln $contig $reads1 > $outDir/bwa/$sai1.sai\n";
			`$bin_folder/bwa aln $contig $reads1 >$outDir/bwa/$sai1.sai`;
			print "$bin_folder/bwa aln $contig $reads2 > $outDir/bwa/$sai2.sai\n";
			`$bin_folder/bwa aln $contig $reads2 > $outDir/bwa/$sai2.sai`;
			print "$bin_folder/bwa sampe $contig $outDir/bwa/$sai1.sai $outDir/bwa/$sai2.sai $reads1 $reads2 >$outDir/bwa/$sample\_$date.sam 2>$outDir/bwa/$sample\_$date.sam.log\n";
			`$bin_folder/bwa sampe $contig $outDir/bwa/$sai1.sai $outDir/bwa/$sai2.sai $reads1 $reads2 >$outDir/bwa/$sample\_$date.sam 2>$outDir/bwa/$sample\_$date.sam.log`;
			chomp( my $mess = `head -4 $outDir/bwa/$sample\_$date.sam.log | tail -1` );
			my $data = ( split /\:/, $mess )[-1];
			my ( $insert, $sd ) = ( split /\s+/, $data )[ 1, 3 ];
			$insert = int( $insert + 0.5 );
			$sd     = int( $sd + 0.5 );
			print OUT2 "$date\t$insert\t$sd\n";
			$info{$date} = $insert;
		}
		close OUT2;
	}
}

# Write config
print "    <<< Write new config file >>>\n";
print "mkdir -p $cwd/$sample/assembly.$reads.$mode.K$kmer\n";
system "mkdir -p $cwd/$sample/assembly.$reads.$mode.K$kmer";
open OUT, ">", "$cwd/$sample/assembly.$reads.$mode.K$kmer/$sample.config";
print OUT "max_rd_len=$max\n";
foreach my $key ( keys %info ) {
	my $insert_size = $info{$key};
	if ( $pair_end eq 'yes' ) {
		print OUT "[LIB]" . "\n";
		print OUT "avg_ins=$info{$key}\n";
		print OUT "asm_flags=3" . "\n";
		if ( $insert_size >= 500 ) {
			print OUT "rank=2" . "\n";
		}
		else {
			print OUT "rank=1" . "\n";
		}
		my $reads2 = $pairPath{$sample}{$key};
		$reads2 =~ s/1\.fq/2\.fq/;
		print OUT "q1=$pairPath{$sample}{$key}\n";
		print OUT "q2=$reads2\n";
	}
	chomp(my $NOL = `$ZCAT $singlePath{$sample}{$key} | head | grep -c .`);
	if ($NOL > 0 ){
		print OUT "[LIB]" . "\n";
		print OUT "asm_flags=1" . "\n";
		print OUT "q=$singlePath{$sample}{$key}\n";
	}
}
close OUT;
print "<<< Completed prepareAssembly >>>\n";
exit 0;

#SUB
sub uniq {
	return keys %{ { map { $_ => 1 } @_ } };
}
