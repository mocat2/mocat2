package MOCATFunctionalProfiling;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

sub create_job {

	### DEFINE VARIABLES AND OPEN input FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	system "mkdir -p $cwd/logs/$job/jobs/";
	open LOG, '>', "$cwd/logs/$job/jobs/MOCATJob_$job\_$date.log" or die "ERROR & EXIT: Cannot write $cwd/logs/$job/jobs/MOCATJob_$job\_$date.log. Do you have permission to write in $cwd?";
	
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}

	chomp(my $host = `hostname`);
	my $groupByColumnExecutable = "$bin_dir/groupByColumn.$host";
	unless (-e "$groupByColumnExecutable") {
		die "ERROR & EXIT: Missing host specific build of groupByColumn executable '$groupByColumnExecutable'. Please run on a different server or make this version.";
	}


	# Reset arrays used to save file names
	my @functional_profiling_base   = ();
	my @functional_profiling_insert = ();
	my $databases                   = join( "_AND_", @do_functional_profiling );
	$conf{functional_profiling_eggnog_map} = "$data_dir/$databases.functional.map.eggnog";
	$conf{functional_profiling_kegg_map} = "$data_dir/$databases.functional.map.kegg";
	
	chomp( my $user = `whoami` );
	my $memory_folder = "/dev/shm/$user/functional_profiling_$date";
	if ($functional_profiling_temp) {
		$memory_folder = $temp_dir;
	}

	my $output_folder_main  = "$cwd/functional.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
	my $output_folder       = "$cwd/functional.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}/gene.level";
	my $eggnog_base         = "$output_folder/$sample_file_basename.functional.profile.$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.base.raw.cog";
	my $eggnog_insert       = "$output_folder/$sample_file_basename.functional.profile.$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.insert.raw.cog";
	my $kegg_ko_base        = "$output_folder/$sample_file_basename.functional.profile.$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.base.raw.ko";
	my $kegg_module_base    = "$output_folder/$sample_file_basename.functional.profile.$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.base.raw.module";
	my $kegg_pathway_base   = "$output_folder/$sample_file_basename.functional.profile.$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.base.raw.pathway";
	my $kegg_ko_insert      = "$output_folder/$sample_file_basename.functional.profile.$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.insert.raw.ko";
	my $kegg_module_insert  = "$output_folder/$sample_file_basename.functional.profile.$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.insert.raw.module";
	my $kegg_pathway_insert = "$output_folder/$sample_file_basename.functional.profile.$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.insert.raw.pathway";
	print localtime() . ": SAVING TEMPORARY DB FILES AND OUTPUT FILES IN $memory_folder\n";

	### CHECKS ###
	print localtime() . ": Checking map files...";
	unless ( -e "$conf{functional_profiling_eggnog_map}" ) {
		die "\nERROR & EXIT: Missing EggNOG map file $conf{functional_profiling_eggnog_map}";
	}
	unless ( -e "$conf{functional_profiling_kegg_map}" ) {
		die "\nERROR & EXIT: Missing KEGG map file $conf{functional_profiling_kegg_map}";
	}
	print " OK!\n";
	print localtime() . ": Checking database rowname files...";
	unless ( -e "$data_dir/$databases.rownames" ) {
		die "\nERROR & EXIT: Missing database rownames file $data_dir/$databases.rownames";
	}
	print " OK!\n";
	print localtime() . ": Checking database length files...";
	foreach my $db (@do_functional_profiling) {
		unless ( -e "$data_dir/$db.len" ) {
			die "\nERROR & EXIT: Missing $data_dir/$db.len";
		}
	}
	print " OK!\n";

	# Store sample names
	print localtime() . ": Loading and checking sample files...";
	foreach my $sample (@samples) {
		my $insert = "$cwd/$sample/insert.coverage.$databases.$conf{MOCAT_data_type}/$sample.filtered.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.insert.coverage.count";
		my $base   = "$cwd/$sample/base.coverage.$databases.$conf{MOCAT_data_type}/$sample.filtered.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.base.coverage.count";
		unless ( -e $insert ) {
			die "\nERROR & EXIT: Missing $insert";
		}
		unless ( -e $base ) {
			die "\nERROR & EXIT: Missing $base";
		}
		push @taxo_profiling_insert, $insert;
		push @taxo_profiling_base,   $base;
	}
	print " OK!\n";

	# Create structure PLEASE NOTE OUTPUT FOLDER 1 IS .../gene, FURTHER DOWN THIS IS REPLACED WITH A m// WHEN STORING THE FINAL SUMMARIZED PROFILES IN .../
	print localtime() . ": Creating folders...";
	print LOG "mkdir -p $output_folder $memory_folder";
	system "mkdir -p $output_folder $memory_folder";
	print " OK!\n";

	# Make length file
	print localtime() . ": Indexing length files...";
	my $DB;
	foreach my $db (@do_functional_profiling) {
		$DB = $DB . " $data_dir/$db.coord ";
	}
	if ( -e "$memory_folder/$databases.len.ff.cidx" ) {
		print "\nFILE $memory_folder/$databases.len.ff.cidx ALREADY EXISTS. SKIPPIN THIS STEP.\n";
	}
	else {
		print LOG "cat $DB | perl -ane '\$l=\$F[2]-\$F[1]+1;print \">\$F[0]\\n\$l\\n\"' > $memory_folder/$databases.len.ff && $bin_dir/cdbfasta $memory_folder/$databases.len.ff 2>/dev/null >/dev/null\n\n";
		system "cat $DB | perl -ane '\$l=\$F[2]-\$F[1]+1;print \">\$F[0]\\n\$l\\n\"' > $memory_folder/$databases.len.ff && $bin_dir/cdbfasta $memory_folder/$databases.len.ff 2>/dev/null >/dev/null";
		
		print " OK!\n";
	}

	# Make pasted file
	print localtime() . ": Pasting, filtering, converting and indexing sample insert files...";
	if ( -e "$memory_folder/$sample_file_basename.insert.gene.profiles.gt0.cidx" ) {
		print "\nFILE $memory_folder/$sample_file_basename.insert.gene.profiles.gt0.cidx ALREADY EXISTS. SKIPPIN THIS STEP.\n";
	}
	else {
		my %hash;
		my $counter1 = 0;
		my $counter2 = 0;
		for my $k (@taxo_profiling_insert) {
			$counter1++;
			unless ( $hash{$counter2} ) {
				$hash{$counter2} = "";
			}
			$hash{$counter2} = $hash{$counter2} . " $k ";
			if ( $counter1 == 100 ) {
				$counter1 = 0;
				$counter2++;
			}
		}
		print LOG "cp $data_dir/$databases.rownames $memory_folder/$sample_file_basename.insert.all\n\n";
		system "cp $data_dir/$databases.rownames $memory_folder/$sample_file_basename.insert.all";
		foreach my $f ( sort {$a<=>$b} keys %hash ) {
			print LOG "paste $memory_folder/$sample_file_basename.insert.all $hash{$f} > $memory_folder/$sample_file_basename.insert.all.tmp && mv $memory_folder/$sample_file_basename.insert.all.tmp $memory_folder/$sample_file_basename.insert.all\n\n";
			system "paste $memory_folder/$sample_file_basename.insert.all $hash{$f} > $memory_folder/$sample_file_basename.insert.all.tmp && mv $memory_folder/$sample_file_basename.insert.all.tmp $memory_folder/$sample_file_basename.insert.all";
		}
		print LOG "rm -f $memory_folder/$sample_file_basename.insert.all.tmp\n\n";
		system "rm -f $memory_folder/$sample_file_basename.insert.all.tmp";
		
		print LOG " perl -lane 'if( \$. <= 2 ){print}else{\$sum = 0; foreach \$i ( \@F[1..\$#F] ){\$sum += \$i;}; if (\$sum != 0){print}}' $memory_folder/$sample_file_basename.insert.all | sed -e 's/^/>/' -e 's/\\t/\\n/' > $memory_folder/$sample_file_basename.insert.gene.profiles.gt0 && cdbfasta $memory_folder/$sample_file_basename.insert.gene.profiles.gt0 >/dev/null 2>/dev/null\n\n";
		system " perl -lane 'if( \$. <= 2 ){print}else{\$sum = 0; foreach \$i ( \@F[1..\$#F] ){\$sum += \$i;}; if (\$sum != 0){print}}' $memory_folder/$sample_file_basename.insert.all | sed -e 's/^/>/' -e 's/\\t/\\n/' > $memory_folder/$sample_file_basename.insert.gene.profiles.gt0 && cdbfasta $memory_folder/$sample_file_basename.insert.gene.profiles.gt0 >/dev/null 2>/dev/null";
		print " OK!\n";
	}
	print localtime() . ": Pasting, filtering, converting and indexing sample base files...";
	if ( -e "$memory_folder/$sample_file_basename.base.gene.profiles.gt0.cidx" ) {
		print "\nFILE $memory_folder/$sample_file_basename.base.gene.profiles.gt0.cidx ALREADY EXISTS. SKIPPIN THIS STEP.\n";
	}
	else {
		my %hash;
		my $counter1 = 0;
		my $counter2 = 0;
		for my $k (@taxo_profiling_base) {
			$counter1++;
			unless ( $hash{$counter2} ) {
				$hash{$counter2} = "";
			}
			$hash{$counter2} = $hash{$counter2} . " $k ";
			if ( $counter1 == 100 ) {
				$counter1 = 0;
				$counter2++;
			}
		}
		print LOG "cp $data_dir/$databases.rownames $memory_folder/$sample_file_basename.base.all\n\n";
		system "cp $data_dir/$databases.rownames $memory_folder/$sample_file_basename.base.all";
		foreach my $f ( sort {$a<=>$b} keys %hash ) {
			print LOG "paste $memory_folder/$sample_file_basename.base.all $hash{$f} > $memory_folder/$sample_file_basename.base.all.tmp && cp $memory_folder/$sample_file_basename.base.all.tmp $memory_folder/$sample_file_basename.base.all\n\n";
			system "paste $memory_folder/$sample_file_basename.base.all $hash{$f} > $memory_folder/$sample_file_basename.base.all.tmp && cp $memory_folder/$sample_file_basename.base.all.tmp $memory_folder/$sample_file_basename.base.all";
		}
		print LOG "rm -f $memory_folder/$sample_file_basename.base.all.tmp\n\n";
		system "rm -f $memory_folder/$sample_file_basename.base.all.tmp";
		print LOG "perl -lane 'if( \$. <= 2 ){print}else{\$sum = 0; foreach \$i ( \@F[1..\$#F] ){\$sum += \$i;}; if (\$sum != 0){print}}' $memory_folder/$sample_file_basename.base.all | sed -e 's/^/>/' -e 's/\\t/\\n/' > $memory_folder/$sample_file_basename.base.gene.profiles.gt0 && cdbfasta $memory_folder/$sample_file_basename.base.gene.profiles.gt0 >/dev/null 2>/dev/null\n\n";
		system "perl -lane 'if( \$. <= 2 ){print}else{\$sum = 0; foreach \$i ( \@F[1..\$#F] ){\$sum += \$i;}; if (\$sum != 0){print}}' $memory_folder/$sample_file_basename.base.all | sed -e 's/^/>/' -e 's/\\t/\\n/' > $memory_folder/$sample_file_basename.base.gene.profiles.gt0 && cdbfasta $memory_folder/$sample_file_basename.base.gene.profiles.gt0 >/dev/null 2>/dev/null";
		print " OK!\n";
	}

	# EggNOG parts
	my $header = "gene_id\tlength\talignment_fraction\tCOG\t" . join( "\t", @samples ) . "\n";
	open OUT, ">$eggnog_base.gene" or die "ERROR & EXIT: Cannot write to $eggnog_base.gene";
	print OUT $header;
	close OUT;
	open OUT, ">$eggnog_insert.gene" or die "ERROR & EXIT: Cannot write to $eggnog_insert.gene";
	print OUT $header;
	close OUT;
	print localtime() . ": Copying and indexing EggNOG DB file...";
	chomp( my $BN = `basename $conf{functional_profiling_eggnog_map}` );
	print LOG "grep -Pv '^#' $conf{functional_profiling_eggnog_map} |  sed -e 's/^/>/' -e 's/\\t/\\n/' > $memory_folder/$BN.ff && $bin_dir/cdbfasta $memory_folder/$BN.ff >/dev/null 2>/dev/null\n\n";
	system "grep -Pv '^#' $conf{functional_profiling_eggnog_map} |  sed -e 's/^/>/' -e 's/\\t/\\n/' > $memory_folder/$BN.ff && $bin_dir/cdbfasta $memory_folder/$BN.ff >/dev/null 2>/dev/null";
	print " OK!\n";
	print localtime() . ": Extracting EggNOG forward insert entries...";
	print LOG "grep -P '^>' $memory_folder/$BN.ff | sed 's/^>//' | $bin_dir/cdbyank -x $memory_folder/$sample_file_basename.insert.gene.profiles.gt0.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort > $memory_folder/eggnog.gene.insert.out\n\n";
	system "grep -P '^>' $memory_folder/$BN.ff | sed 's/^>//' | $bin_dir/cdbyank -x $memory_folder/$sample_file_basename.insert.gene.profiles.gt0.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort > $memory_folder/eggnog.gene.insert.out";
	print " OK!\n";
	print localtime() . ": Extracting EggNOG forward base entries...";
	print LOG "grep -P '^>' $memory_folder/$BN.ff | sed 's/^>//' | $bin_dir/cdbyank -x $memory_folder/$sample_file_basename.base.gene.profiles.gt0.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort > $memory_folder/eggnog.gene.base.out\n\n";
	system "grep -P '^>' $memory_folder/$BN.ff | sed 's/^>//' | $bin_dir/cdbyank -x $memory_folder/$sample_file_basename.base.gene.profiles.gt0.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort > $memory_folder/eggnog.gene.base.out";
	print " OK!\n";
	print localtime() . ": Extracting EggNOG reverse entries...";
	print LOG "cut -f 1 $memory_folder/eggnog.gene.insert.out | $bin_dir/cdbyank -x $memory_folder/$BN.ff.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort > $memory_folder/eggnog.db.out\n\n";
	system "cut -f 1 $memory_folder/eggnog.gene.insert.out | $bin_dir/cdbyank -x $memory_folder/$BN.ff.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort > $memory_folder/eggnog.db.out";
	print " OK!\n";
	print localtime() . ": Extracting EggNOG gene length entries...";
	print LOG "cut -f 1 $memory_folder/eggnog.gene.insert.out | $bin_dir/cdbyank -x $memory_folder/$databases.len.ff.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort > $memory_folder/eggnog.length.out\n\n";
	system "cut -f 1 $memory_folder/eggnog.gene.insert.out | $bin_dir/cdbyank -x $memory_folder/$databases.len.ff.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort > $memory_folder/eggnog.length.out";
	print " OK!\n";
	print localtime() . ": Joining and processing EggNOG base file...";
	print LOG "join -t \$'\\t' $memory_folder/eggnog.length.out $memory_folder/eggnog.db.out | join -t \$'\\t' - $memory_folder/eggnog.gene.base.out | perl -F\"\\t\" -lane '\$alifrac = (\$F[5] - \$F[4] + 1) / \$F[7]; print \"\$F[0]\\t\$F[1]\\t\$alifrac\\t\$F[3]\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]) unless scalar \@F == 9' >> $eggnog_base.gene\n\n";
	system "join -t \$'\\t' $memory_folder/eggnog.length.out $memory_folder/eggnog.db.out | join -t \$'\\t' - $memory_folder/eggnog.gene.base.out | perl -F\"\\t\" -lane '\$alifrac = (\$F[5] - \$F[4] + 1) / \$F[7]; print \"\$F[0]\\t\$F[1]\\t\$alifrac\\t\$F[3]\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]) unless scalar \@F == 9' >> $eggnog_base.gene";
	print " OK!\n";
	print localtime() . ": Joining and processing EggNOG insert file...";
	print LOG "join -t \$'\\t' $memory_folder/eggnog.length.out $memory_folder/eggnog.db.out | join -t \$'\\t' - $memory_folder/eggnog.gene.insert.out | perl -F\"\\t\" -lane '\$alifrac = (\$F[5] - \$F[4] + 1) / \$F[7]; print \"\$F[0]\\t\$F[1]\\t\$alifrac\\t\$F[3]\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]) unless scalar \@F == 9' >> $eggnog_insert.gene\n\n";
	system "join -t \$'\\t' $memory_folder/eggnog.length.out $memory_folder/eggnog.db.out | join -t \$'\\t' - $memory_folder/eggnog.gene.insert.out | perl -F\"\\t\" -lane '\$alifrac = (\$F[5] - \$F[4] + 1) / \$F[7]; print \"\$F[0]\\t\$F[1]\\t\$alifrac\\t\$F[3]\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]) unless scalar \@F == 9' >> $eggnog_insert.gene";
	print " OK!\n";

	# KEGG parts
	$header = "gene_id\tlength\talignment_fraction\tKO\t" . join( "\t", @samples ) . "\n";
	open OUT, ">$kegg_ko_base.gene" or die "ERROR & EXIT: Cannot write to $kegg_ko_base.gene";
	print OUT $header;
	close OUT;
	open OUT, ">$kegg_ko_insert.gene" or die "ERROR & EXIT: Cannot write to $kegg_ko_insert.gene";
	print OUT $header;
	close OUT;
	$header = "gene_id\tlength\talignment_fraction\tmodule\t" . join( "\t", @samples ) . "\n";
	open OUT, ">$kegg_module_base.gene" or die "ERROR & EXIT: Cannot write to $kegg_module_base.gene";
	print OUT $header;
	close OUT;
	open OUT, ">$kegg_module_insert.gene" or die "ERROR & EXIT: Cannot write to $kegg_module_insert.gene";
	print OUT $header;
	close OUT;
	$header = "gene_id\tlength\talignment_fraction\tpathway\t" . join( "\t", @samples ) . "\n";
	open OUT, ">$kegg_pathway_base.gene" or die "ERROR & EXIT: Cannot write to $kegg_pathway_base.gene";
	print OUT $header;
	close OUT;
	open OUT, ">$kegg_pathway_insert.gene" or die "ERROR & EXIT: Cannot write to $kegg_pathway_insert.gene";
	print OUT $header;
	close OUT;
	print localtime() . ": Copying and indexing KEGG DB file...";
	chomp( $BN = `basename $conf{functional_profiling_kegg_map}` );
	print LOG "grep -Pv '^#' $conf{functional_profiling_kegg_map} |  sed -e 's/^/>/' -e 's/\\t/\\n/' > $memory_folder/$BN.ff && $bin_dir/cdbfasta $memory_folder/$BN.ff >/dev/null 2>/dev/null\n\n";
	system "grep -Pv '^#' $conf{functional_profiling_kegg_map} |  sed -e 's/^/>/' -e 's/\\t/\\n/' > $memory_folder/$BN.ff && $bin_dir/cdbfasta $memory_folder/$BN.ff >/dev/null 2>/dev/null";
	print " OK!\n";
	print localtime() . ": Extracting KEGG forward insert entries...";
	print LOG "grep -P '^>' $memory_folder/$BN.ff | sed 's/^>//' | $bin_dir/cdbyank -x $memory_folder/$sample_file_basename.insert.gene.profiles.gt0.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort > $memory_folder/kegg.gene.insert.out\n\n";
	system "grep -P '^>' $memory_folder/$BN.ff | sed 's/^>//' | $bin_dir/cdbyank -x $memory_folder/$sample_file_basename.insert.gene.profiles.gt0.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort > $memory_folder/kegg.gene.insert.out";
	print " OK!\n";
	print localtime() . ": Extracting KEGG forward base entries...";
	print LOG "grep -P '^>' $memory_folder/$BN.ff | sed 's/^>//' | $bin_dir/cdbyank -x $memory_folder/$sample_file_basename.base.gene.profiles.gt0.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort > $memory_folder/kegg.gene.base.out\n\n";
	system "grep -P '^>' $memory_folder/$BN.ff | sed 's/^>//' | $bin_dir/cdbyank -x $memory_folder/$sample_file_basename.base.gene.profiles.gt0.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort > $memory_folder/kegg.gene.base.out";
	print " OK!\n";
	print localtime() . ": Extracting KEGG reverse entries...";
	print LOG "cut -f 1 $memory_folder/kegg.gene.insert.out $memory_folder/kegg.gene.base.out | $bin_dir/cdbyank -x $memory_folder/$BN.ff.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort -u > $memory_folder/kegg.db.out\n\n";
	system "cut -f 1 $memory_folder/kegg.gene.insert.out $memory_folder/kegg.gene.base.out | $bin_dir/cdbyank -x $memory_folder/$BN.ff.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort -u > $memory_folder/kegg.db.out";
	print " OK!\n";
	print localtime() . ": Extracting KEGG gene length entries...";
	print LOG "cut -f 1 $memory_folder/kegg.gene.insert.out | $bin_dir/cdbyank -x $memory_folder/$databases.len.ff.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort > $memory_folder/kegg.length.out\n\n";
	system "cut -f 1 $memory_folder/kegg.gene.insert.out | $bin_dir/cdbyank -x $memory_folder/$databases.len.ff.cidx | perl -ane 'if (m/>(.*)/){chomp(\$m=\$1)} else {print \"\$m\\t\$_\"}' | sort > $memory_folder/kegg.length.out";
	print " OK!\n";
	print localtime() . ": Joining and processing KEGG base files...";
	print LOG "join -t \$'\\t' $memory_folder/kegg.length.out $memory_folder/kegg.db.out | join -t \$'\\t' - $memory_folder/kegg.gene.base.out | perl -F\"\\t\" -lane '\@ko = split \",\",\$F[6]; foreach \$i (\@ko){ unless (scalar \@F == 9) {print \"\$F[0]\\t\$F[1]\\t\$F[5]\\t\$i\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]); } }' >> $kegg_ko_base.gene\n\n";
	print LOG "join -t \$'\\t' $memory_folder/kegg.length.out $memory_folder/kegg.db.out | join -t \$'\\t' - $memory_folder/kegg.gene.base.out | perl -F\"\\t\" -lane '\@ko = split \",\",\$F[7]; foreach \$i (\@ko){ unless (scalar \@F == 9) {print \"\$F[0]\\t\$F[1]\\t\$F[5]\\t\$i\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]); } }' >> $kegg_module_base.gene\n\n";
	print LOG "join -t \$'\\t' $memory_folder/kegg.length.out $memory_folder/kegg.db.out | join -t \$'\\t' - $memory_folder/kegg.gene.base.out | perl -F\"\\t\" -lane '\@ko = split \",\",\$F[8]; foreach \$i (\@ko){ unless (scalar \@F == 9) {print \"\$F[0]\\t\$F[1]\\t\$F[5]\\t\$i\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]); } }' >> $kegg_pathway_base.gene\n\n";	
	system "join -t \$'\\t' $memory_folder/kegg.length.out $memory_folder/kegg.db.out | join -t \$'\\t' - $memory_folder/kegg.gene.base.out | perl -F\"\\t\" -lane '\@ko = split \",\",\$F[6]; foreach \$i (\@ko){ unless (scalar \@F == 9) {print \"\$F[0]\\t\$F[1]\\t\$F[5]\\t\$i\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]); } }' >> $kegg_ko_base.gene";
	system "join -t \$'\\t' $memory_folder/kegg.length.out $memory_folder/kegg.db.out | join -t \$'\\t' - $memory_folder/kegg.gene.base.out | perl -F\"\\t\" -lane '\@ko = split \",\",\$F[7]; foreach \$i (\@ko){ unless (scalar \@F == 9) {print \"\$F[0]\\t\$F[1]\\t\$F[5]\\t\$i\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]); } }' >> $kegg_module_base.gene";
	system "join -t \$'\\t' $memory_folder/kegg.length.out $memory_folder/kegg.db.out | join -t \$'\\t' - $memory_folder/kegg.gene.base.out | perl -F\"\\t\" -lane '\@ko = split \",\",\$F[8]; foreach \$i (\@ko){ unless (scalar \@F == 9) {print \"\$F[0]\\t\$F[1]\\t\$F[5]\\t\$i\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]); } }' >> $kegg_pathway_base.gene";
	print " OK!\n";
	print localtime() . ": Joining and processing KEGG insert files...";
	print LOG "join -t \$'\\t' $memory_folder/kegg.length.out $memory_folder/kegg.db.out | join -t \$'\\t' - $memory_folder/kegg.gene.insert.out | perl -F\"\\t\" -lane '\@ko = split \",\",\$F[6]; foreach \$i (\@ko){ unless (scalar \@F == 9) {print \"\$F[0]\\t\$F[1]\\t\$F[5]\\t\$i\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]); } }' >> $kegg_ko_insert.gene\n\n";
	print LOG "join -t \$'\\t' $memory_folder/kegg.length.out $memory_folder/kegg.db.out | join -t \$'\\t' - $memory_folder/kegg.gene.insert.out | perl -F\"\\t\" -lane '\@ko = split \",\",\$F[7]; foreach \$i (\@ko){ unless (scalar \@F == 9) {print \"\$F[0]\\t\$F[1]\\t\$F[5]\\t\$i\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]); } }' >> $kegg_module_insert.gene\n\n";
	print LOG "join -t \$'\\t' $memory_folder/kegg.length.out $memory_folder/kegg.db.out | join -t \$'\\t' - $memory_folder/kegg.gene.insert.out | perl -F\"\\t\" -lane '\@ko = split \",\",\$F[8]; foreach \$i (\@ko){ unless (scalar \@F == 9) {print \"\$F[0]\\t\$F[1]\\t\$F[5]\\t\$i\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]); } }' >> $kegg_pathway_insert.gene\n\n";
	system "join -t \$'\\t' $memory_folder/kegg.length.out $memory_folder/kegg.db.out | join -t \$'\\t' - $memory_folder/kegg.gene.insert.out | perl -F\"\\t\" -lane '\@ko = split \",\",\$F[6]; foreach \$i (\@ko){ unless (scalar \@F == 9) {print \"\$F[0]\\t\$F[1]\\t\$F[5]\\t\$i\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]); } }' >> $kegg_ko_insert.gene";
	system "join -t \$'\\t' $memory_folder/kegg.length.out $memory_folder/kegg.db.out | join -t \$'\\t' - $memory_folder/kegg.gene.insert.out | perl -F\"\\t\" -lane '\@ko = split \",\",\$F[7]; foreach \$i (\@ko){ unless (scalar \@F == 9) {print \"\$F[0]\\t\$F[1]\\t\$F[5]\\t\$i\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]); } }' >> $kegg_module_insert.gene";
	system "join -t \$'\\t' $memory_folder/kegg.length.out $memory_folder/kegg.db.out | join -t \$'\\t' - $memory_folder/kegg.gene.insert.out | perl -F\"\\t\" -lane '\@ko = split \",\",\$F[8]; foreach \$i (\@ko){ unless (scalar \@F == 9) {print \"\$F[0]\\t\$F[1]\\t\$F[5]\\t\$i\\t\".join(\"\\t\", \@F[9..(scalar \@F-1)]); } }' >> $kegg_pathway_insert.gene";
	print " OK!\n";

	# Summarize
	print localtime() . ": Summarizing files...";
	for my $file ( $eggnog_base, $eggnog_insert, $kegg_ko_base, $kegg_ko_insert, $kegg_module_base, $kegg_module_insert, $kegg_pathway_base, $kegg_pathway_insert ) {
		my $norm = $file;
		$norm =~ s/\/gene.level\//\//;
		$norm =~ s/\.raw\./.raw.sum./;
		
		my $sum = $norm;
		my $median = $norm;
		my $mean = $norm;
		$median =~ s/\.raw\.sum\./.raw.median./;
		$mean =~ s/\.raw\.sum\./.raw.mean./;
		
		my $mean2 = $mean;
		$mean =~ s/\.raw\.mean\./.raw.non.zero.values.mean./;
		
		# New version
		print LOG "$groupByColumnExecutable -g 4 -f 5 -h 1 -m $median -a $mean -b $mean2 $file.gene > $sum 2>/dev/null\n\n";
		system "$groupByColumnExecutable -g 4 -f 5 -h 1 -m $median -a $mean -b $mean2 $file.gene > $sum 2>/dev/null";
		
		# Old version
		#system "$conf{MOCAT_DEV_DIR}/MOCATTaxoProfiling_groupByColumn.pl -i $file.gene -g 4 -f 5 -h 1 > $norm";
		#print ".";
		#system "$conf{MOCAT_DEV_DIR}/MOCATTaxoProfiling_groupByColumn_medians.pl -i $file.gene -g 4 -f 5 -h 1 > $norm";
		print ".";
	}
	print " OK!\n";


	#norm (gene length normalization)
	print localtime() . ": Gene length normalizing files...";
	for my $file ( $eggnog_base, $eggnog_insert, $kegg_ko_base, $kegg_ko_insert, $kegg_module_base, $kegg_module_insert, $kegg_pathway_base, $kegg_pathway_insert ) {
		my $norm = $file;
		my $alg  = $file;
		$norm =~ s/(base|insert).raw/$1.length.norm/;
		$alg  =~ s/(base|insert).raw/$1.length.algfrac.norm/;
		print LOG "perl -F\"\\t\" -lane 'if (\$. >= 2){\$len = \$F[1]; \$alg = \$F[2]; \@b = map {\$_ / \$len * \$alg} (\@F[4..\$#F]); \@a = map {\$_ / \$len} (\@F[4..\$#F]); print STDOUT \"\$F[0]\\t\$F[3]\\t\" . join \"\\t\", \@a; print STDERR \"\$F[0]\\t\$F[3]\\t\" . join \"\\t\", \@b} else{print STDOUT \"\$F[0]\\t\$F[3]\\t\" . join \"\\t\", \@F[4..\$#F]; print STDERR \"\$F[0]\\t\$F[3]\\t\" . join \"\\t\", \@F[4..\$#F]}' $file.gene > $norm.gene 2> $alg.gene\n\n";
		system "perl -F\"\\t\" -lane 'if (\$. >= 2){\$len = \$F[1]; \$alg = \$F[2]; \@b = map {\$_ / \$len * \$alg} (\@F[4..\$#F]); \@a = map {\$_ / \$len} (\@F[4..\$#F]); print STDOUT \"\$F[0]\\t\$F[3]\\t\" . join \"\\t\", \@a; print STDERR \"\$F[0]\\t\$F[3]\\t\" . join \"\\t\", \@b} else{print STDOUT \"\$F[0]\\t\$F[3]\\t\" . join \"\\t\", \@F[4..\$#F]; print STDERR \"\$F[0]\\t\$F[3]\\t\" . join \"\\t\", \@F[4..\$#F]}' $file.gene > $norm.gene 2> $alg.gene";
		print ".";
	}
	print " OK!\n";

	# Summarize
	print localtime() . ": Summarizing normlized files...";
	for my $file ( $eggnog_base, $eggnog_insert, $kegg_ko_base, $kegg_ko_insert, $kegg_module_base, $kegg_module_insert, $kegg_pathway_base, $kegg_pathway_insert ) {
		my $norm = $file;
		my $alg  = $file;
		$norm =~ s/(base|insert).raw/$1.length.norm/;
		$alg  =~ s/(base|insert).raw/$1.length.algfrac.norm/;
		my $out  = $norm;
		my $out2 = $alg;
		$out  =~ s/\/gene.level\//\//;
		$out2 =~ s/\/gene.level\//\//;
		$out =~ s/\.length\.norm\./.length.norm.sum./;
		$out2 =~ s/\.length\.algfrac\.norm\./.length.algfrac.norm.sum./;
		
		my $sum = $out;
		my $sum_alg = $out2;
		my $median = $out;
		my $median_alg = $out2;
		my $mean = $out;
		my $mean_alg = $out2;
		
		$median =~ s/\.length\.norm\.sum\./.length.norm.median./;
		$median_alg =~ s/\.length\.algfrac\.norm\.sum\./.length.algfrac.norm.median./;
		$mean =~ s/\.length\.norm\.sum\./.length.norm.mean./;
		$mean_alg =~ s/\.length\.algfrac\.norm\.sum\./.length.algfrac.norm.mean./;
		
		my $mean2 = $mean;
		my $mean2_alg = $mean_alg;
		$mean2 =~ s/\.length\.norm\.mean\./.length.norm.non.zero.values.mean./;
		$mean2_alg =~ s/\.length\.algfrac\.norm\.mean\./.length.algfrac.norm.non.zero.values.mean./;
		
		
		# new version
		print LOG "$groupByColumnExecutable -g 2 -f 3 -h 1 -m $median -a $mean -b $mean2 $norm.gene > $sum 2>/dev/null\n\n";
		system "$groupByColumnExecutable -g 2 -f 3 -h 1 -m $median -a $mean -b $mean2 $norm.gene > $sum 2>/dev/null";
		print ".";
		print LOG "$groupByColumnExecutable -g 2 -f 3 -h 1 -m $median_alg -a $mean_alg -b $mean2_alg $alg.gene > $sum_alg 2>/dev/null\n\n";
		system "$groupByColumnExecutable -g 2 -f 3 -h 1 -m $median_alg -a $mean_alg -b $mean2_alg $alg.gene > $sum_alg 2>/dev/null";
		print ".";
		
		# old version, note that to make this work (again), you have to rename out and out2 before the third system command
		#system "$conf{MOCAT_DEV_DIR}/MOCATTaxoProfiling_groupByColumn.pl -i $norm.gene -g 2 -f 3 -h 1 > $out";
		#print ".";
		#system "$conf{MOCAT_DEV_DIR}/MOCATTaxoProfiling_groupByColumn.pl -i $alg.gene -g 2 -f 3 -h 1 > $out2";
		#print ".";
		#system "$conf{MOCAT_DEV_DIR}/MOCATTaxoProfiling_groupByColumn_medians.pl -i $norm.gene -g 2 -f 3 -h 1 > $out";
		#print ".";
		#system "$conf{MOCAT_DEV_DIR}/MOCATTaxoProfiling_groupByColumn_medians.pl -i $alg.gene -g 2 -f 3 -h 1 > $out2";
		#print ".";
	}
	print " OK!\n";

	# Delete files
	print localtime() . ": Removing temporary files...";

	print LOG "rm -rf $memory_folder\n\n";
	system "rm -rf $memory_folder";
	print " OK!\n";

	# Zip gene files
	print localtime() . ": Compressing gene files...";
	print LOG "$ZIP -$ziplevel -f $output_folder/$sample_file_basename*\n\n";
	system "$ZIP -$ziplevel -f $output_folder/$sample_file_basename*";
	print " OK!\n";
	
	# Zip main files
	print localtime() . ": Compressing main files...";
	print LOG "$ZIP -$ziplevel -f $output_folder_main/$sample_file_basename*\n\n";
	system "$ZIP -$ziplevel -f $output_folder_main/$sample_file_basename*";
	print " OK!\n";
	print localtime() . ": RESULTS SAVED IN $output_folder_main\n";
	

}

1;
