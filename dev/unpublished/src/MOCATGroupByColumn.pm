package MOCATGroupByColumn;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

sub run {
	### DEFINE VARIABLES AND OPEN OUTPUT FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	chomp( my $user = `whoami` );
	my $temp_dir = "/dev/shm/$user/group_by_column_$date";
	if ($functional_profiling_temp) {
		$temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
	}
	if ($temp_dir eq "/dev/shm/$user/group_by_column_$date") {
		system "mkdir -p $temp_dir";
	}
	
	chomp(my $host = `hostname`);
	my $groupByColumnExecutable = "$bin_dir/groupByColumn.$host";
	unless (-e "$groupByColumnExecutable") {
		die "ERROR & EXIT: Missing host specific build of groupByColumn executable '$groupByColumnExecutable'. Please run on a different server or make this version.";
	}
	
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}

	my $uniq = "";
	if ($pcf_uniq) {
		$uniq = ".uniq";
	}
	
	my $map_file = $gbc_map;
	unless ($map_file) {
		die "ERROR & EXIT: Missing -map option";
	}
	unless (-e $map_file) {
		die "ERROR & EXIT: Missing file $map_file";
	}

	my $databases = join( "_AND_", @do_group_by_column );
	my $assembly_type = 'assembly';
	my $end;
	my $assembly = 0;
	chomp( my $bn = `basename $sample_file` );
	chomp( my $bn_map_file = `basename $map_file` );
	if (         $do_group_by_column[0] eq 's'
		|| $do_group_by_column[0] eq 'c'
		|| $do_group_by_column[0] eq 'f'
		|| $do_group_by_column[0] eq 'r' )
	{
		if ( $do_group_by_column[0] eq 's' ) { $end = 'scaftig' }
		if ( $do_group_by_column[0] eq 'c' ) { $end = 'contig' }
		if ( $do_group_by_column[0] eq 'f' ) { $end = 'scafSeq' }
		if ( $do_group_by_column[0] eq 'r' ) {
			$assembly_type = 'assembly.revised';
			$end           = 'scaftig';
		}
		$assembly  = 1;
		$databases = "$assembly_type.$end";
	}

	system "mkdir -p $cwd/abundance.tables/$databases/group.by.column/$bn_map_file";
	print localtime() . ": Loading mapping file...";
	open IN, "<$map_file";
	chomp( my $line = <IN>);
	my @header = split "\t", $line;
	my %hash;
	my $columns = -1;
	while (<IN>){
		chomp;
		my @line = split "\t";
		unless ($columns == scalar @line || $columns == -1) {
			die "\nERROR & EXIT: Mapping file has different number of columns on different rows";
		}
		$columns = scalar @line;
		$hash{$line[0]} = join("\t", @line[1 .. scalar @line-1]);
	}
	close IN;
	print " OK!\n";
	
	foreach my $boi ( ( 'base', 'insert' ) ) {
		foreach my $norm ( ( 'raw', 'norm' ) ) {
			print localtime() . ": Pasting $boi:$norm file to $temp_dir...";
			my $input_file = "$cwd/abundance.tables/$databases/$bn.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$boi.coverage.$norm";
			my $output_file = "$temp_dir/$bn.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$boi.coverage.$norm.tmp";
			open IN, "<$input_file" or die "\nERROR & EXIT: Missing $input_file";
			open OUT, ">$output_file" or die "\nERROR & EXIT: Cannot write to $output_file";
			
			# process 2 fist lines
			chomp(my $line1 = <IN>);
			my @line1 = split "\t", $line1;
			chomp(my $line2 = <IN>);
			my @line2 = split "\t", $line2;
			print OUT join("\t", @header) . "\t" . join("\t", @line1[1 .. scalar @line1-1]) . "\n";
			print OUT "$line2[0]\t". "-1\t"x($columns-1) . join("\t", @line2[1 .. scalar @line2-1]) . "\n";
			
			while (<IN>) {
				chomp;
				my @line = split "\t";
				if ($hash{$line[0]}) {
					print OUT "$line[0]\t$hash{$line[0]}\t" . join("\t", @line[1 .. scalar @line-1]) . "\n";
				} else {
					print OUT "$line[0]\t". "unassigned\t"x($columns-1) . join("\t", @line[1 .. scalar @line-1]) . "\n";
				}
			}
			close IN;
			close OUT;
			print " OK!\n";
			my $index = 1;
			my $columns2 = scalar @header + 1;
			for my $name (@header[1 .. scalar @header-1]) {
				$index++;
				print localtime() . ": Generating $boi:$norm:$name files...";
				my $final_output_sum_file = "$cwd/abundance.tables/$databases/group.by.column/$bn_map_file/$bn.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$boi.$norm.sum.$name";
				my $final_output_median_file = "$cwd/abundance.tables/$databases/group.by.column/$bn_map_file/$bn.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$boi.$norm.median.$name";
				my $final_output_mean_file = "$cwd/abundance.tables/$databases/group.by.column/$bn_map_file/$bn.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$boi.$norm.mean.$name";
				my $final_output_only_non_zero_mean_file = "$cwd/abundance.tables/$databases/group.by.column/$bn_map_file/$bn.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$boi.$norm.only.non.zero.values.mean.$name";
				system "$groupByColumnExecutable -g $index -f $columns2 -h 2 -a $final_output_mean_file -b $final_output_only_non_zero_mean_file -m $final_output_median_file $output_file > $final_output_sum_file 2> /dev/null";
				system "$ZIP -$ziplevel -f $final_output_mean_file $final_output_only_non_zero_mean_file $final_output_median_file $final_output_sum_file";
				print " OK!\n";
			}
			system "rm -fr $output_file";
		}
	}
	if ($temp_dir eq "/dev/shm/$user/group_by_column_$date") {
		system "rm -r $temp_dir";
	}
	print localtime() . ": Completed.\n";
	print localtime() . ": RESULTS SAVED IN $cwd/abundance.tables/$databases/group.by.column/$bn_map_file/*\n";
}

1;
