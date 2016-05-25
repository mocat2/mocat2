package MOCATSampleStatus;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;
use Spreadsheet::WriteExcel;
use Spreadsheet::WriteExcel::Utility qw( xl_range_formula );

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

sub system_ {
	my $cmd = shift;
	system($cmd) == 0 or die("\nERROR & EXIT: system($cmd) failed: $!\n");
}

sub run {
	my $original_sample_file = $sample_file;
	if ($original_sample_file =~ m/\//) {
		my @a = split '/', $original_sample_file;
		$original_sample_file = $a[-1];
	}
	my $sample_file = "$cwd/SUMMARIES/$original_sample_file/$original_sample_file";
	my $sample_file_folder = "$cwd/SUMMARIES/$original_sample_file/";
	system_("mkdir -p $cwd/SUMMARIES/$original_sample_file");
	system_("cp $original_sample_file $sample_file_folder");
	print localtime() . ": PROCESSING SAMPLES\n";
	open OUT, '>', "$sample_file.status" or die "ERROR & EXIT: Cannot write $sample_file.status. Do you have permission to write in $sample_file.status?";

	my ( $p, @c, @s );

	my %FivePrime;
	my @FivePrime = `cat $cwd/MOCAT.cutoff5prime 2>/dev/null| cut -f 1 | sort -u 2>/dev/null`;
	foreach my $fp (@FivePrime) {
		chomp($fp);
		$FivePrime{$fp} = 1;
	}

	my $myMax = 0;
	foreach my $sample (@samples) {
		if ( $myMax < length($sample) ) {
			$myMax = length($sample);
		}
	}
	if ( $myMax < 6 ) {
		$myMax = 6;
	}

	printf OUT "SAMPLE";
	print OUT " " x ( $myMax - 6 );
	print OUT "|FQC|5|R TF|SCREEN|ASMB|R.ASMB|GENEPRED|R.GENEPRED|\n";

	my $counter = 0;

	foreach my $sample (@samples) {

		$counter++;
		print localtime() . ": [$counter/" . scalar @samples . "] Processing $sample...";

		print OUT "$sample";
		print OUT " " x ( $myMax - length($sample) );
		print OUT "|";

		# FastQC files
		$p = 0;
		@c = `ls $cwd/$sample/*_fastqc 2>/dev/null`;
		if ( scalar @c > 0 ) { print OUT "FQC"; $p = 1; }
		unless ($p) {
			print OUT "---";
		}
		print OUT "|";
		$p = 0;

		# Set 5 prime in cutoff5prime file
		$p = 0;
		if ( $FivePrime{$sample} ) {
			print OUT "M";
			$p = 1;
		}
		unless ($p) {
			print OUT "A";
		}
		print OUT "|";
		$p = 0;

		# RTF
		@c = `ls $cwd/$sample/reads.processed.fastx/*pair.1*.gz $cwd/$sample/reads.processed.fastx/*single*.gz $cwd/$sample/reads.processed.fastx/*pair.2*.gz 2>/dev/null`;
		if ( scalar @c % 3 == 0 && scalar @c >= 3 && ( -s "$cwd/$sample/stats/$sample.readtrimfilter.fastx.stats" ) ) { print OUT "XP"; $p++; $p++; }
		@c = `ls $cwd/$sample/reads.processed.fastx/*trimmed*.gz 2>/dev/null`;
		if ( scalar @c > 0 && ( -s "$cwd/$sample/stats/$sample.readtrimfilter.fastx.stats" ) ) { print OUT "XS"; $p++; $p++; }
		for my $i ( $p .. 1 ) {
			print OUT "-";
		}
		$p = 0;
		@c = `ls $cwd/$sample/reads.processed.solexaqa/*pair.1*.gz $cwd/$sample/reads.processed.solexaqa/*single*.gz $cwd/$sample/reads.processed.solexaqa/*pair.2*.gz 2>/dev/null`;
		if ( scalar @c % 3 == 0 && scalar @c >= 3 && ( -s "$cwd/$sample/stats/$sample.readtrimfilter.solexaqa.stats" ) ) { print OUT "QP"; $p++; $p++; }
		@c = `ls $cwd/$sample/reads.processed.solexaqa/*trimmed*.gz 2>/dev/null`;
		if ( scalar @c > 0 && ( -s "$cwd/$sample/stats/$sample.readtrimfilter.solexaqa.stats" ) ) { print OUT "QS"; $p++; $p++; }
		for my $i ( $p .. 1 ) {
			print OUT "-";
		}
		print OUT "|";
		$p = 0;

		# DB/FF Screen
		@s = `ls $cwd/$sample/stats/$sample.screen*.*.fastx.stats 2>/dev/null`;
		@c = `ls $cwd/$sample/reads.screened.*.fastx/*.screened.*.pair.*.fq.gz $cwd/$sample/reads.screened.*.fastx/*.screened.*.single.fq.gz 2>/dev/null | grep -v 'reads.screened.fastafile'`;
		if ( scalar @c % 3 == 0 && scalar @c > 0 && scalar @s > 0 ) { print OUT "SCX"; $p++; $p++; $p++; }
		for my $i ( $p .. 2 ) {
			print OUT "-";
		}
		$p = 0;
		@s = `ls $cwd/$sample/stats/$sample.screen*.*.solexaqa.stats 2>/dev/null`;
		@c = `ls $cwd/$sample/reads.screened.*.solexaqa/*.screened.*.pair.*.fq.gz $cwd/$sample/reads.screened.*.solexaqa/*.screened.*.single.fq.gz 2>/dev/null | grep -v 'reads.screened.fastafile'`;
		if ( scalar @c % 3 == 0 && scalar @c > 0 && scalar @s > 0 ) { print OUT "SCQ"; $p++; $p++; $p++; }
		for my $i ( $p .. 2 ) {
			print OUT "-";
		}
		$p = 0;
		print OUT "|";
		$p = 0;

		#    # Fasta File Screen
		#    @s=`ls $cwd/$sample/stats/$sample.screen*.fastafile.fastx.stats 2>/dev/null`;
		#    @c=`ls $cwd/$sample/reads.screened.fastafile.fastx/*.screened.*.pair.*.fq.gz $cwd/$sample/reads.screened.fastafile.fastx/*.screened.*.single.fq.gz 2>/dev/null`;if(scalar @c % 3 == 0 && scalar @c > 0 && scalar @s > 0){print OUT "FFX";$p++;$p++;$p++;}
		#    for my $i ($p..2) {
		#      print OUT "-";
		#    }
		#    $p=0;
		#    @s=`ls $cwd/$sample/stats/$sample.screen*.fastafile.solexaqa.stats 2>/dev/null`;
		#    @c=`ls $cwd/$sample/reads.screened.fastafile.solexaqa/*.screened.*.pair.*.fq.gz $cwd/$sample/reads.screened.fastafile.solexaqa/*.screened.*.single.fq.gz 2>/dev/null`;if(scalar @c % 3 == 0 && scalar @c > 0 && scalar @s > 0){print OUT "FFQ";$p++;$p++;$p++;}
		#    for my $i ($p..2) {
		#      print OUT "-";
		#    }
		#    ; print OUT "|"; $p=0;

		# Assembly
		@s = `ls $cwd/$sample/stats/$sample.assembly.*.fastx.K*.assembly.stats 2>/dev/null | grep -v 'assembly.revised.'`;
		@c = `ls $cwd/$sample/assembly.*.fastx.K*/$sample.*.scaftig.gz $cwd/$sample/assembly.*.fastx.K*/$sample.*.scaftig 2>/dev/null | grep -v 'assembly.revised.'`;
		if ( scalar @c > 0 && scalar @s > 0 ) { print OUT "AX"; $p++; $p++; }
		for my $i ( $p .. 1 ) {
			print OUT "-";
		}
		$p = 0;
		@s = `ls $cwd/$sample/stats/$sample.assembly.*.solexaqa.K*.scaftig.stats 2>/dev/null | grep -v 'assembly.revised.'`;
		@c = `ls $cwd/$sample/assembly.*.solexaqa.K*/$sample.*.scaftig.gz $cwd/$sample/assembly.*.solexaqa.K*/$sample.*.scaftig 2>/dev/null | grep -v 'assembly.revised.'`;
		if ( scalar @c > 0 && scalar @s > 0 ) { print OUT "AQ"; $p++; $p++; }
		for my $i ( $p .. 1 ) {
			print OUT "-";
		}
		print OUT "|";
		$p = 0;

		# Assembly Revised
		@s = `ls $cwd/$sample/stats/$sample.assembly.revised.*.fastx.K*.baseErrorAndIndelError.stats 2>/dev/null`;
		@c = `ls $cwd/$sample/assembly.revised.*.fastx.K*/$sample.*.scaftig.gz $cwd/$sample/assembly.revised.*.fastx.K*/$sample.*.scaftig 2>/dev/null`;
		if ( scalar @c > 0 && scalar @s > 0 ) { print OUT "ARX"; $p++; $p++; $p++; }
		for my $i ( $p .. 2 ) {
			print OUT "-";
		}
		$p = 0;
		@s = `ls $cwd/$sample/stats/$sample.assembly.revised.*.solexaqa.K*.baseErrorAndIndelError.stats 2>/dev/null`;
		@c = `ls $cwd/$sample/assembly.revised.*.solexaqa.K*/$sample.*.scaftig.gz $cwd/$sample/assembly.revised.*.solexaqa.K*/$sample.*.scaftig 2>/dev/null`;
		if ( scalar @c > 0 && scalar @s > 0 ) { print OUT "ARQ"; $p++; $p++; $p++; }
		for my $i ( $p .. 2 ) {
			print OUT "-";
		}
		print OUT "|";
		$p = 0;

		# Gene Prediction Assembly
		@c = `ls $cwd/$sample/gene.prediction.assembly.*.fastx.K*/$sample.*.tab* 2>/dev/null | grep -v 'assembly.revised'`;
		if ( scalar @c > 0 ) { print OUT "GPAX"; $p++; $p++; $p++; $p++; }
		for my $i ( $p .. 3 ) {
			print OUT "-";
		}
		$p = 0;
		@c = `ls $cwd/$sample/gene.prediction.assembly.*.solexaqa.K*/$sample.*.tab* 2>/dev/null | grep -v 'assembly.revised'`;
		if ( scalar @c > 0 ) { print OUT "GPAQ"; $p++; $p++; $p++; $p++; }
		for my $i ( $p .. 3 ) {
			print OUT "-";
		}
		print OUT "|";
		$p = 0;

		# Gene Prediction Assembly Revised
		@c = `ls $cwd/$sample/gene.prediction.assembly.revised.*.fastx.K*/$sample.*.tab* 2>/dev/null`;
		if ( scalar @c > 0 ) { print OUT "GPARX"; $p++; $p++; $p++; $p++; $p++; }
		for my $i ( $p .. 4 ) {
			print OUT "-";
		}
		$p = 0;
		@c = `ls $cwd/$sample/gene.prediction.assembly.revised.*.solexaqa.K*/$sample.*.tab* 2>/dev/null`;
		if ( scalar @c > 0 ) { print OUT "GPARQ"; $p++; $p++; $p++; $p++; $p++; }
		for my $i ( $p .. 4 ) {
			print OUT "-";
		}
		print OUT "|";
		$p = 0;
		print OUT "\n";
		print " OK!\n";
	}
	close OUT;

	print localtime() . ": Processing statistics files";
	foreach my $sample (@samples) {
		print ".";
		my @files = `ls -1 $cwd/$sample/stats/*.stats 2>/dev/null | grep -v '.multiplemapper.' 2>/dev/null`;
		foreach my $file (@files) {
			chomp($file);
			if ( $file =~ /$sample\.readtrimfilter\.(fastx|solexaqa)\.stats/ ) {
				open OUT, ">", "$sample_file.readtrimfilter.$1.summary";
				print OUT "sample\treads\tbases\tmax\tavg\tkmer\tinserts\n";
				close OUT;
			}
			if ( $file =~ /$sample\.readtrimfilter\.inserts.(fastx|solexaqa)\.stats/ ) {
				open OUT, ">", "$sample_file.readtrimfilter.inserts.$1.summary";
				print OUT "sample\tinserts\n";
				close OUT;
			}
			if ( $file =~ /$sample\.screen(ed)\.(.*).(solexaqa|fastx)\.stats/ ) {
				open OUT, ">", "$sample_file.screened.$2.$3.summary";
				print OUT "sample\treads\tbases\tmax\tavg\tkmer\tinserts\tmin_allowed_identity\tmin_allowed_length\tmax_allowed_SOAP_mismatches\n";
				close OUT;
			}
			if ( $file =~ /$sample\.extracted\.(.*).(solexaqa|fastx)\.stats/ ) {
				open OUT, ">", "$sample_file.extracted.$1.$2.summary";
				print OUT "sample\treads\tbases\tmax\tavg\tkmer\tinserts\tmin_allowed_identity\tmin_allowed_length\tmax_allowed_SOAP_mismatches\n";
				close OUT;
			}
			if ( $file =~ /$sample\.screen(ed)\.(.*)\.inserts\.(solexaqa|fastx)\.stats/ ) {
				open OUT, ">", "$sample_file.screened.$2.inserts.$3.summary";
				print OUT "sample\tinserts\n";
				close OUT;
			}
			
			#Not sure when this is used... But probably is?
			if ( $file =~ /(.*)\.([12])\.screen(ed)\.(.*)\.stats/ ) {
				open OUT, ">", "$sample_file.screened.$4.lanes.summary";
				print OUT "sample\tlane\ttotal_reads\ttotal_hits\tunique_hits\n";
				close OUT;
			}
			
			#PE
			if ( $file =~ /$sample\/stats\/(.*)\.[12]\..*raw\.reads\.stats/ ) {
				open OUT, ">", "$sample_file.raw.reads.lanes.summary";
				print OUT "sample\tlane\treads\tbases\n";
				close OUT;
			}
			
			#SE
			elsif ( $file =~ /$sample\/stats\/(.*)\..*raw\.reads\.stats/ ) {
				open OUT, ">", "$sample_file.raw.reads.lanes.summary";
				print OUT "sample\tlane\treads\tbases\n";
				close OUT;
			}
			
			if ( $file =~ /$sample\.assembly\.(.*)\.inserts.stats/ ) {
				my $n1 = $1;
				$n1 =~ s/\.K\d+//;
				open OUT, ">", "$sample_file.assembly.$n1.inserts.summary";
				print OUT "sample\tlane\tinsert_size\n";
				close OUT;
			}
			if ( $file =~ /$sample\.coverage\.(.*)\.stats/ ) {
				my $n1 = $1;
				$n1 =~ s/\.K\d+//;
				open OUT, ">", "$sample_file.coverage.summary";
				print OUT "sample\tfile\ttotal_inserts\tmapped_inserts\tfraction_mapped_inserts\ttotal_bases\tmapped_bases\tfraction_mapped_bases\tdb_average_gene_length\n";
				close OUT;
			}
			if ( $file =~ /$sample\.assembly\.revised\.(.*).baseErrorAndIndelError\.stats/ ) {
				my $n1 = $1;
				$n1 =~ s/\.K\d+//;
				open OUT, ">", "$sample_file.assembly.revised.$n1.summary";
				print OUT "sample\tno_of_single_base_error\tno_of_short_inserts\tno_of_total_len_inserts\tno_of_short_del\ttot_len_short_del\tno_of_ridiculous_regions\ttot_len_ridiculous_region\n";
				close OUT;
			}
			if ( $file =~ /$sample\.assembly\.(.*)\.(fastx|solexaqa)\.K(.*)\.assembly\.stats/ ) {
				open OUT, ">", "$sample_file.assembly.$1.$2.contig.stats.summary";
				print OUT "sample\tcontig_num\tcontig_length_all\tccontig_n50\tcontig_n90\tcontig_max\tcontig_min\n";
				close OUT;

				open OUT, ">", "$sample_file.assembly.$1.$2.scaffold.stats.summary";
				print OUT "sample\tscaffold_num\tscaffold_length_all\tscaffold_n50\tscaffold_n90\tscaffold_max\tscaffold_min\n";
				close OUT;

				open OUT, ">", "$sample_file.assembly.$1.$2.scaftig.stats.summary";
				print OUT "sample\tscaftig_num\tscaftig_length_all\tscaftig_n50\tscaftig_n90\tscaftig_max\tscaftig_min\n";
				close OUT;
			}
			if ( $file =~ /$sample\.assembly\.revised\.(.*).scaftig\.stats/ ) {
				my $n1 = $1;
				$n1 =~ s/\.K\d+//;
				open OUT, ">", "$sample_file.assembly.revised.$n1.scaftig.stats.summary";
				print OUT "sample\tscaftig_num\tscaftig_length_all\tscaftig_n50\tscaftig_n90\tscaftig_max\tscaftig_min\n";
				close OUT;
			}
		}
	}
	foreach my $sample (@samples) {
		print ".";
		my @files = `ls -1 $cwd/$sample/stats/*.stats 2>/dev/null | grep -v '.multiplemapper.' 2>/dev/null`;
		foreach my $file (@files) {
			chomp($file);
			if ( $file =~ /$sample\.readtrimfilter\.(fastx|solexaqa)\.stats/ ) {
				open IN,  '<',  "$file";
				open OUT, ">>", "$sample_file.readtrimfilter.$1.summary";
				<IN>;
				print OUT "$sample\t" . <IN>;
				close OUT;
				close IN;
			}
			if ( $file =~ /$sample\.readtrimfilter\.inserts.(fastx|solexaqa)\.stats/ ) {
				open IN,  '<',  "$file";
				open OUT, ">>", "$sample_file.readtrimfilter.inserts.$1.summary";
				<IN>;
				print OUT "$sample\t" . <IN>;
				close OUT;
				close IN;
			}
			if ( $file =~ /$sample\.screen(ed)\.(.*).(solexaqa|fastx)\.stats/ ) {
				open IN,  '<',  "$file";
				open OUT, ">>", "$sample_file.screened.$2.$3.summary";
				<IN>;
				print OUT "$sample\t" . <IN>;
				close OUT;
				close IN;
			}
			if ( $file =~ /$sample\.extracted\.(.*).(solexaqa|fastx)\.stats/ ) {
				open IN,  '<',  "$file";
				open OUT, ">>", "$sample_file.extracted.$1.$2.summary";
				<IN>;
				print OUT "$sample\t" . <IN>;
				close OUT;
				close IN;
			}
			if ( $file =~ /$sample\.screen(ed)\.(.*)\.inserts\.(solexaqa|fastx)\.stats/ ) {
				open IN,  '<',  "$file";
				open OUT, ">>", "$sample_file.screened.$2.inserts.$3.summary";
				<IN>;
				print OUT "$sample\t" . <IN>;
				close OUT;
				close IN;
			}
			if ( $file =~ /$cwd\/$sample\/stats\/(.*)\.([12])\.screen(ed)\.(.*)\.stats/ ) {
				open IN,  '<',  "$file";
				open OUT, ">>", "$sample_file.screened.$4.lanes.summary";
				<IN>;
				print OUT "$sample\t$1\t" . <IN>;
				close OUT;
				close IN;
			}
			
			#PE
			if ( $file =~ /$cwd\/$sample\/stats\/(.*\.[12])\..*raw\.reads\.stats/ ) {
				open IN,  '<',  "$file";
				open OUT, ">>", "$sample_file.raw.reads.lanes.summary";
				<IN>;
				print OUT "$sample\t$1\t" . <IN>;
				close OUT;
				close IN;
			}
			
			#SE
			elsif ( $file =~ /$cwd\/$sample\/stats\/(.*)\..*raw\.reads\.stats/ ) {
				open IN,  '<',  "$file";
				open OUT, ">>", "$sample_file.raw.reads.lanes.summary";
				<IN>;
				print OUT "$sample\t$1\t" . <IN>;
				close OUT;
				close IN;
			}
						
			if ( $file =~ /$sample\.assembly\.(.*)\.inserts\.stats/ ) {
				open IN, '<', "$file";
				my $n1 = $1;
				$n1 =~ s/\.K\d+//;
				open OUT, ">>", "$sample_file.assembly.$n1.inserts.summary";
				my %saved;
				while (<IN>) {
					chomp;
					my $save = $_;
					unless ($saved{$_}) {
						print OUT "$sample\t$save\n";
					}
					$saved{$save} = 1;
				}
				close OUT;
				close IN;
			}
			if ( $file =~ /$sample\.coverage\.(.*)\.stats/ ) {
				open IN,  '<',  "$file";
				open OUT, ">>", "$sample_file.coverage.summary";
				<IN>;
				print OUT "$sample\t$1\t" . <IN>;
				close OUT;
				close IN;
			}
			if ( $file =~ /$sample\.assembly\.revised\.(.*).baseErrorAndIndelError\.stats/ ) {
				open IN, '<', "$file";
				my $n1 = $1;
				$n1 =~ s/\.K\d+//;
				open OUT, ">>", "$sample_file.assembly.revised.$n1.summary";
				my $save;
				while (<IN>) {
					chomp;
					$save = $_;
				}
				$save =~ s/^\S+/$sample/;
				print OUT "$save\n";
				close OUT;
				close IN;
			}
			if ( $file =~ /$sample\.assembly\.revised\.(.*).scaftig\.stats/ ) {
				open IN, "<$file";
				my $n1 = $1;
				$n1 =~ s/\.K\d+//;
				open OUT, ">>", "$sample_file.assembly.revised.$n1.scaftig.stats.summary";
				my $save;
				my $header = <IN>;
				my %saved;
				while (<IN>) {
					chomp;
					$save = $_;
					unless ($saved{$_}) {
						print OUT "$sample\t$save\n";
					}
					$saved{$save} = 1;
				}
				close OUT;
				close IN;
			}
			if ( $file =~ /$sample\.assembly\.(.*)\.(fastx|solexaqa)\.K(.*)\.assembly\.stats/ || $file =~ /$sample\.assembly\.revised\.(.*)\.(fastx|solexaqa)\.K(.*)\.assembly\.stats/ ) {
				my $revised = "";
				if ( $file =~ /$sample\.assembly\.revised\.(.*)\.(fastx|solexaqa)\.K(.*)\.assembly\.stats/ ) {
					$revised = "revised.";
				}
				my $n1    = $1;
				my $n2    = $2;
				my $saved = "";
				open IN, "<$file";
				while (<IN>) {
					chomp;
					if ( $_ =~ /contig/ ) {
						close OUT;
						open OUT, ">>", "$sample_file.assembly.$revised$n1.$n2.contig.stats.summary";
						$saved = "";
					}
					if ( $_ =~ /scafold/ ) {
						close OUT;
						open OUT, ">>", "$sample_file.assembly.$revised$n1.$n2.scaffold.stats.summary";
						$saved = "";
					}
					if ( $_ =~ /scaftig/ ) {
						close OUT;
						open OUT, ">>", "$sample_file.assembly.$revised$n1.$n2.scaftig.stats.summary";
						$saved = "";
					}
					if ( $_ =~ /^\d+/ ) {
						my @line = split /\s+/, $_;
						unless ( $line[5] eq $saved ) {
							$saved = $line[5];
							print OUT "$sample\t" . join( "\t", @line ) . "\n";
						}
					}
				}
				close OUT;

			}

		}
	}
	print " DONE!\n";

	print localtime() . ": PROCESSING GENE PREDICTIONS:\n";
	print localtime() . ": Preprocessing";
	foreach my $sample (@samples) {
		print ".";
		my @genepreds = `find $cwd/$sample/ -follow -name '*gene.prediction.*.tab*'`;
		foreach my $genepred (@genepreds) {
			chomp($genepred);
			my @gps = split( "/", $genepred );
			$genepred = $gps[-2];
			$genepred =~ s/\.K\d+//;
			open OUT, ">", "$sample_file.$genepred.summary";
			print OUT "sample\ttotal_genes\tcomplete_genes\tincomplete_genes\tgenes_with_start_codon\tgenes_with_stop_codon\ttotal_bp_sum\ttot_sum_bp_complete_genes\ttot_bp_sum_incomplete_genes\n";
			close OUT;
		}
	}
	print " DONE!\n";
	$counter = 0;
	foreach my $sample (@samples) {
		my @genepreds = `find $cwd/$sample/ -follow -name '*gene.prediction.*.tab*'`;
		$counter++;
		print localtime() . ": [$counter/" . scalar @samples . "] Processing $sample...";
		foreach my $genepred (@genepreds) {
			chomp($genepred);
			my $PRE = '';
			if ($genepred =~ m/.gz$/) {
				$PRE = 'z';
			}
			my ($incomplete, $total, $complete, $has_start, $has_stop, @tot_len);
			if ($PRE eq '') {
				chomp( $incomplete = `grep -c 'incomplete' $genepred` );
				chomp( $total      = `grep -c . $genepred` );
				$complete = $total - $incomplete;
				chomp( $has_start = `cut -f 6 $genepred | grep -c 'yes'` );
				chomp( $has_stop  = `cut -f 7 $genepred | grep -c 'yes'` );
				@tot_len = `awk '{if (\$8==\"complete\"){sum += \$5}; if (\$8==\"incomplete\"){insum += \$5 } } END {print sum \"\\n\" insum}' $genepred`;
			}
			if ($PRE eq 'z') {
				chomp( $incomplete = `zgrep -c 'incomplete' $genepred` );
				chomp( $total      = `zgrep -c . $genepred` );
				$complete = $total - $incomplete;
				chomp( $has_start = `zcat $genepred | cut -f 6 | grep -c 'yes'` );
				chomp( $has_stop  = `zcat $genepred | cut -f 7 | grep -c 'yes'` );
				@tot_len = `zcat $genepred | awk '{if (\$8==\"complete\"){sum += \$5}; if (\$8==\"incomplete\"){insum += \$5 } } END {print sum \"\\n\" insum}'`;
			}
			chomp( my $compl_sum   = $tot_len[0] );
			chomp( my $incompl_sum = $tot_len[1] );

			if ( $compl_sum eq "" ) {
				$compl_sum = "0";
			}
			if ( $incompl_sum eq "" ) {
				$incompl_sum = "0";
			}
			my $tot_sum = $compl_sum + $incompl_sum;
			my @gps = split( "/", $genepred );
			$genepred = $gps[-2];
			$genepred =~ s/\.K\d+//;
			open OUT, ">>", "$sample_file.$genepred.summary";
			print OUT "$sample\t$total\t$complete\t$incomplete\t$has_start\t$has_stop\t$tot_sum\t$compl_sum\t$incompl_sum\n";
			close OUT;
		}
		print " OK!\n";
	}

	print localtime() . ": Processing Excel file...";
	my @files      = `find $sample_file_folder -maxdepth 1 -xtype f -name '$original_sample_file.*' | grep -v '.xls' | grep -v '.sqlite'`;
	my $workbook   = Spreadsheet::WriteExcel->new("$sample_file.xls");

	foreach my $file (@files) {
		chomp(my $name = $file);
		$name =~ s/solexaqa\.lanes/lanes.solexaqa/;
		$name =~ s/fastx\.lanes/lanes.fastx/;
		$name =~ s/solexaqa/S/;
		$name =~ s/fastx/F/;
		$name =~ s/.summary//;
		$name =~ s/$sample_file//;
		$name =~ s/gene.prediction/GP/;
		$name =~ s/assembly.revised/ASS.REV./;
		$name =~ s/assembly/ASS/;
		$name =~ s/^\.\///;
		$name =~ s/^\.//;
		$name =~ s/screened/S/g;
		$name =~ s/extracted/E/g;
		$name =~ s/fastafile/FF/;
		$name =~ s/raw.reads.lane/RAW/;
		$name =~ s/readtrimfilter/RTF/;
		$name =~ s/\.$//;
		my $random_number = rand(999999999);
		$random_number =~ m/(.*)\..*/;
		$name = substr( $name, 0, 20 ) . "_" . $1;
		open IN, "<$file";
		my $worksheet = $workbook->add_worksheet($name);
		my $i         = 0;
		$worksheet->write( $i, 0, $file );
		while (<IN>) {
			$i++;
			chomp;
			my $char;
			if ( $file =~ m/.status/ ) {
				$char = "\\|";
			}
			else {
				$char = "\t";
			}
			my @line = split $char;
			for my $j ( 0 .. scalar @line - 1 ) {
				$worksheet->write( $i, $j, $line[$j] );
			}
		}
		close IN;
	}
	print " OK!\n";

	print localtime() . ": Processing SQLite3 file...";
	#print "\n- - -\nTO CREATE DB, EXECUTE: rm -f $sample_file.sqlite* && $scr_dir/MOCATSampleStatus_DB.pl $sample_file.sqlite $sample_file $bin_dir\n- - -\n";
	#system "rm -f $sample_file.sqlite* && $scr_dir/MOCATSampleStatus_DB.pl $sample_file.sqlite $sample_file $bin_dir";
	print " OK!\n";

	print localtime() . ": COMPLETED CHECKS.\n";
	print localtime() . ": Output saved in $sample_file.xls\n";
	print localtime() . ": Output saved in $sample_file.sqlite3\n";
	print localtime() . ": Output saved in $sample_file.summary.xlsx\n";
	print localtime() . ": COMPLETED CALCULATING.\n";
	print "$sep\n";
	system "cat $sample_file.status";
}


1;
