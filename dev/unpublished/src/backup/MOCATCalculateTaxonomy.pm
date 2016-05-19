package MOCATCalculateTaxonomy;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.
sub create_job {
	### DEFINE VARIABLES AND OPEN OUTPUT FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?\n";
	my $warnings = 0;
	print localtime() . ": Creating $job jobs...";
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
	my $print_rownames_file;
	my $read_type = 'screened';

	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}

	$taxo_profiling_map = $conf{taxo_profiling_map};
	if ($taxo_profiling_map_in) {
		$taxo_profiling_map = $taxo_profiling_map_in;
	}
	if ( $multiple_mappers_summarize_level ne 'none' ) {
		unless ( -e "$taxo_profiling_map" ) {
			die "\nERROR & EXIT: Missing map file $taxo_profiling_map";
		}
	}
	else {
		$taxo_profiling_map = "not_used";
	}

	$taxo_profiling_cog      = $conf{taxo_profiling_len};
	unless ( -e "$taxo_profiling_cog" ) {
		die "\nERROR & EXIT: Missing len file $taxo_profiling_cog\n";
	}
	my @tmp = `awk -F\"\\t\" '{print NF}' $taxo_profiling_cog | sort -u`;
	foreach my $tmp (@tmp) {
		chomp($tmp);
		unless ( $tmp eq '2' ) {
			die "\nERROR & EXIT: Length file doesn't have exactly 2 columns on all lines.\n";
		}
	}
	
	my @mtfl;
	my $map_reads_to_filelist = 'no';
	my $assembly_type         = 'assembly';
	my $end;
	if (         $databases[0] eq 's'
		|| $databases[0] eq 'c'
		|| $databases[0] eq 'f'
		|| $databases[0] eq 'r' )
	{
		$print_rownames_file   = 'yes';
		$map_reads_to_filelist = 'yes';
		if ( $databases[0] eq 's' ) { $end = 'scaftig' }
		if ( $databases[0] eq 'c' ) { $end = 'contig' }
		if ( $databases[0] eq 'f' ) { $end = 'scafSeq' }
		if ( $databases[0] eq 'r' ) {
			$assembly_type = 'assembly.revised';
			$end           = 'scaftig';
		}
		foreach my $sample (@samples) {
			my $file;
			( my $max, my $avg, my $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-r" );
			if ( $reads eq 'reads.processed' ) {
				$file = "$cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end";
			}
			else {
				$file = "$cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end";
			}
			unless ( -e "$file.gz" ) {
				die "\nERROR & EXIT: Missing $end file: $file.gz\n";
			}
			print JOB "$ZIP -dc $file.gz > $file && ";
			push @mtfl, $file;
		}
	}
	### JOB SPECIFIC ###
	# Check .coord file, if it doesn't exist, create it from the .len file. If the .len file doesn't exists, create that and the .coord file.
	if ( $map_reads_to_filelist ne 'yes' ) {
		foreach my $databases (@databases) {
			unless ( -e "$data_dir/$databases.coord" ) {
				unless ( -e "$data_dir/$databases.len" ) {
					unless ( -e "$data_dir/$databases" ) {
						die "\nERROR & EXIT: $data_dir/$databases does not exist. Cannot create .len and .coord files\n";
					}
					print "\n" . localtime() . ": Creating length file $data_dir/$databases.len...";
					system "$scr_dir/MOCATFilter_falen.pl -infile $data_dir/$databases -outfile $data_dir/$databases.len";
					print " OK!";
				}
				print "\n" . localtime() . ": Creating $data_dir/$databases.coord, and then continue...";

				if ( $systemType =~ m/Darwin/ ) {
					system "sed \'s/[[:space:]]/	1	/\' $data_dir/$databases.len > $data_dir/$databases.coord";
				}
				else {
					system "sed \'s/\\t/\\t1\\t/\' $data_dir/$databases.len > $data_dir/$databases.coord";
				}

				unless ( -e "$data_dir/$databases.coord" ) {
					print localtime() . ": Error creating .coord file\n";
					die "ERROR & EXIT: Could not create $data_dir/$databases.coord from $data_dir/$databases.len\nDoes the system have write permission?\n";
				}
			}
		}
	}
	my $counter = -1;
	foreach my $sample (@samples) {
		my ( $max, $avg, $kmer );
		if ( $map_reads_to_filelist eq 'yes' ) {
			( $max, $avg, $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-mr" );
		}
		$counter++;
		my $folder_filtered;
		my $LOG = " 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.log >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.log ";
		my ( $output, $outputTax, $outputTaxRownames, $databases );
		if ( $map_reads_to_filelist ne 'yes' ) {
			$databases           = join( "_AND_", @databases );
			$print_rownames_file = $databases;
			$output              = "$sample.filtered.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
			$outputTax           = "$sample.taxonomic.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
			$outputTaxRownames   = "$sample_file_basename.taxonomic.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
			$folder_filtered     = "$cwd/$sample/reads.filtered.$databases.$conf{MOCAT_data_type}";
		}
		my $fullname = "NOFILE";
		if ( $map_reads_to_filelist eq 'yes' ) {
			$fullname = $mtfl[$counter];
			my @parts = split '/', $fullname;
			$databases       = $parts[-1];
			$output          = "$sample.filtered.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
			$outputTax       = "$sample.taxonomic.profile.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
			$outputTaxRownames = "$sample_file_basename.taxonomic.profile.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.rownames";
			$folder_filtered = "$cwd/$sample/reads.filtered.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}";
		}

		#    if ($coverage_version_1) {
		#      $folder_filtered = "$cwd/$sample/reads.filtered.$reads.$conf{MOCAT_data_type}";
		#    } else {

		#    }
		my $folder_coverage = "coverage.$databases.$conf{MOCAT_data_type}";
		my $folder_tax      = "taxonomic.profiles.$databases.$conf{MOCAT_data_type}";
		if ( $map_reads_to_filelist eq 'yes' ) {
			$folder_coverage = "coverage.$databases";
			$folder_tax      = "taxonomic.profiles.$databases";
		}
		my $input;
		if ( $conf{MOCAT_mapping_mode} eq "allbest" ) {
			$input = "$folder_filtered/$output.bam";
		}
		elsif (      $conf{MOCAT_mapping_mode} eq "unique"
			|| $conf{MOCAT_mapping_mode} eq "random" )
		{
			$input = "$folder_filtered/$output.soap.gz";
		}
		my $inserts;
		if ( $reads eq 'reads.processed' ) {
			$inserts = "$cwd/$sample/stats/$sample.readtrimfilter.$conf{MOCAT_data_type}";
		}
		else {
			$inserts = "$cwd/$sample/stats/$sample.$read_type.$reads.$conf{MOCAT_data_type}";
		}
		unless ( -e $input ) {

			# Support for older versions
			$folder_filtered = "$cwd/$sample/reads.filtered.$reads.$conf{MOCAT_data_type}";
			if ( $conf{MOCAT_mapping_mode} eq "allbest" ) {
				$input = "$folder_filtered/$output.bam";
			}
			elsif (      $conf{MOCAT_mapping_mode} eq "unique"
				|| $conf{MOCAT_mapping_mode} eq "random" )
			{
				$input = "$folder_filtered/$output.soap.gz";
			}
			unless ( -e $input ) {
				die "\nERROR & EXIT: Missing filtered mapping results file: $input\n";
			}
			else {
				$warnings = $input;
			}
		}
		my $covfile = "$cwd/$sample/stats/$sample.coverage.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.stats";
		if ( $map_reads_to_filelist eq 'yes' ) {
			$covfile = "$cwd/$sample/stats/$sample.coverage.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.stats";
		}
		print JOB
		  "mkdir -p $cwd/$sample/$folder_tax && mkdir -p $cwd/$sample/base.$folder_coverage && mkdir -p $cwd/$sample/insert.$folder_coverage && $conf{MOCAT_DEV_DIR}/MOCATCalculateTaxonomy.pl -zcat $ZCAT -bin $bin_dir -i $input -s $inserts -sample $sample -cwd $cwd -rcr $reads -dt $conf{MOCAT_data_type} -taxrownames $data_dir/$outputTaxRownames -rownames $data_dir/$databases.rownames -out $folder_coverage/$output -taxout $cwd/$sample/$folder_tax/$outputTax -match $conf{MOCAT_mapping_mode} -datadir $data_dir -pos $databases -file_list $map_reads_to_filelist -file $fullname -falen $scr_dir/MOCATFilter_falen.pl -covfile $covfile -counter $counter -map $taxo_profiling_map -len $taxo_profiling_cog $LOG";
		if ( $map_reads_to_filelist eq 'yes' ) {
			print JOB "rm -f $fullname\n";
		}
		else {
			print JOB "\n";
		}

	}
	close JOB;
	print " OK!\n";
	if ($warnings) {
		print "\n\nWARNING! WARNING! WARNING!\nDid not find new version of filtered results, but found old version (one example):\n$warnings\nSO, CURRENTLY USING THIS OLD VERSION FOR CALCULATING COVERAGES!\n\n";
		print "Please read warnings and press enter to continue.\n";
		chomp( my $key = <STDIN> );
	}
}

sub post_check_files {
	print localtime() . ": Checking files... ";
	print " OK!\n";
}
1;
