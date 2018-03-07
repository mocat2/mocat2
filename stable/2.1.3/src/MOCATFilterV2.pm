package MOCATFilterV2;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

sub system_ {
	my $cmd = shift;
	system($cmd) == 0 or die("\nERROR & EXIT: system($cmd) failed: $!\n");
}

sub create_job {
	### DEFINE VARIABLES AND OPEN OUTPUT FILE ###
	my $job        = shift;
	my $processors = shift;
	$ZIP =~ s/pigz.*/pigz -p 2/;
	my $jobfile = "$jobdir/MOCATJob_$job\_$date";
	open JOB, '>', $jobfile or die "ERROR & EXIT: Cannot open $jobfile (for writing). Do you have permission to write in $jobdir?";

	my $OSX = "";
	if ( $systemType =~ m/Darwin/ ) {
		$OSX = "_OSX";
	}

	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}
	print localtime() . ": Creating $job jobs...";

	# Loop over samples
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
	if ( $SHM ) {
		$temp_dir = "/dev/shm/$username/MOCAT_temp/";
	}

	foreach my $sample (@samples) {
	    MOCATCore::mkdir_or_die("$temp_dir/$sample/temp");

		my $output_folder;
		my $output_file;
		my $input_file = "";
		my $input_folder;
		my $len_file;
		my $full_name   = "";

		my $stats;
		if ( $reads eq 'reads.processed' ) {
			$stats = "$cwd/$sample/stats/$sample.readtrimfilter.$conf{MOCAT_data_type}.stats";
		}
		else {
			$stats = "$cwd/$sample/stats/$sample.$read_type.$reads.$conf{MOCAT_data_type}.stats";
		}

		# Define
		my $LOG = " 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log ";

		print JOB "exec 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log && ";
		print JOB "mkdir -p $temp_dir/$sample/temp && ";
		
		# If filtering a scaftig, contig, scafSeq or revised scaftig
		if (         $databases[0] eq 's'
			|| $databases[0] eq 'c'
			|| $databases[0] eq 'f'
			|| $databases[0] eq 'r' )
		{
			my $assembly_type = 'assembly';
			my $end;
			if ( $databases[0] eq 's' ) {
				$end = 'scaftig';
			}
			if ( $databases[0] eq 'c' ) {
				$end = 'contig';
			}
			if ( $databases[0] eq 'f' ) {
				$end = 'scafSeq';
			}
			if ( $databases[0] eq 'r' ) {
				$assembly_type = 'assembly.revised';
				$end           = 'scaftig';
			}
			( my $max, my $avg, my $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-r" );
			$input_folder  = "$cwd/$sample/reads.mapped.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}";
			$input_file    = "$input_folder/$sample.mapped.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.soap.gz";
			$output_folder = "$cwd/$sample/reads.filtered.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}";
			$output_file   = "$output_folder/$sample.filtered.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
			$full_name     = "$end.$assembly_type.K$kmer";

			#if ( $conf{filter_make_unique_sorted} eq 'yes' ) {

				# Get input file, and check if it exists, also create .len file
				my $assembly_file = "$cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end";
				$len_file = "$assembly_file.len";
				( -e "$assembly_file.gz" ) or die "\nERROR & EXIT: Missing $end file: $assembly_file.gz";
				( -e $input_file ) or die "\nERROR & EXIT: Missing mapping file: $input_file";
				unless ( -e $len_file ) {
					print JOB "$scr_dir/MOCATFilter_falen.pl -infile $assembly_file.gz -outfile $len_file -zip && ";
				}
			#}
		}    # End scaftig

		# Filtering against a DB
		else {
			$len_file = "$temp_dir/$sample/temp/lengths.$date";
			if ( defined $database_name ) {
				$full_name = $database_name;
			}
			else {
				#$full_name = join( "_AND_", @databases );
				$full_name = MOCATCore::checkAndReturnDB ( \@databases );
			}
			$output_folder = "$cwd/$sample/reads.filtered.$full_name.$conf{MOCAT_data_type}";
			$output_file   = "$output_folder/$sample.filtered.$read_type.$reads.on.$full_name.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
			if ( $conf{filter_make_unique_sorted} eq 'yes' ) {
				foreach my $databases (@databases) {
					if ( -e "$data_dir/$databases" ) {
						unless ( -e "$data_dir/$databases.len" ) {
							print "\n" . localtime() . ": Creating length file $data_dir/$databases.len...";
							system "$scr_dir/MOCATFilter_falen.pl -infile $data_dir/$databases -outfile $data_dir/$databases.len";
							print " OK!\n";
							print localtime() . ": Continuing creating $job jobs...";
						}
					}
					else {
						die "\nERROR & EXIT: Missing database $databases in $data_dir";
					}
				}
			}
			foreach my $databases (@databases) {
				unless ( -e "$data_dir/$databases.len" ) {
					print "\n" . localtime() . ": $data_dir/$databases.len is missing. Creating it...";
					system_("$scr_dir/MOCATFilter_falen.pl -infile $data_dir/$databases -outfile $data_dir/$databases.len");
				}
				if ( $reads eq 'reads.processed' ) {
					$input_folder = "$cwd/$sample/reads.mapped.$databases.$conf{MOCAT_data_type}";
				}
				else {
					$input_folder = "$cwd/$sample/reads.mapped.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}";
				}
				my $file;
				if ( $reads eq 'reads.processed' ) {
					$file = "$sample.mapped.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.soap.gz";
				}
				else {
					$file = "$sample.mapped.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.soap.gz";
				}
				unless ( -e "$input_folder/$file" ) {
					die "\nERROR & EXIT: Missing $input_folder/$file";
				}
				$input_file = "$input_folder/$file " . $input_file;
			}
			print JOB "cat $data_dir/" . join( ".len $data_dir/", @databases ) . ".len | sort -u > $len_file && ";

		}

		MOCATCore::mkdir_or_die($output_folder);

		# Print job
		chomp( my $pipefail = `set -o | grep 'pipefail' | grep . -c` );
			if ( $pipefail == 1 ) {
				$pipefail = 'set -o pipefail;';
			}
			else {
				$pipefail = '';
			}
			my $screen_source;
			if ( $reads eq 'reads.processed' ) {
				$screen_source = "reads.processed";
			}
			else {
				$screen_source = "reads.$read_type.$reads";
			}
			chomp( my @pairs   = `ls $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/*pair*gz 2>/dev/null` );
			chomp( my @singles = `ls $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/*single*gz 2>/dev/null` );
			if (scalar @pairs == 0 || scalar @singles == 0) {
				die "\nERROR & EXIT: It seems like one of these two commands failed to report any files:\nls $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/*pair*gz\nls $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/*single*gz";
			}
			
			print JOB "set -e; mkdir -p $temp_dir/$sample/temp && $pipefail $scr_dir/MOCATFilter.p1.pl -zip '$ZIP -2' -input_files $input_file -references ". join (" ", @pairs) . " ". join (" ", @singles) . " -temp_output_file $temp_dir/$sample/temp/$date.tmp.filtering $LOG ";
			print JOB " | $scr_dir/MOCATFilter.p2.pl -temp_input_file $temp_dir/$sample/temp/$date.tmp.filtering $LOG ";
			my $coords = "";
			foreach my $database (@databases) {
				if ( -e "$data_dir/$database.coord" ) {
					$coords = "$coords $data_dir/$database";
				}
			}
			if ( $coords ne "" ) {
				print JOB " | perl $scr_dir/MOCATFilter_remove_in_padded.pl -db $coords $LOG";
			}
		
		if ( $conf{MOCAT_mapping_mode} eq "random" || $conf{MOCAT_mapping_mode} eq "unique" ) {
			print JOB " | perl -F\"\\t\" -lane '\$len=\$F[5];\$mm = \$F[-1] =~ tr/[a-zA-Z]/[a-zA-Z]/;\$as = 100-(\$mm/\$len)*100; if (\$as >= $conf{filter_percent_cutoff} && \$len >= $conf{filter_length_cutoff}){\$read{\$F[0]}++;print \$_}' "
			. " | $ZIP -$ziplevel -c > $output_file.soap.gz.tmp "
			. " | $scr_dir/MOCATFilter_stats.pl -db $full_name -format SOAP -length $conf{filter_length_cutoff} -identity $conf{filter_percent_cutoff} -stats $stats $LOG"
			. " | $ZIP -$ziplevel -c > $output_file.soap.gz.tmp"
#			  . " && export LC_ALL=C"
#			  . " && $bin_dir/psort --parallel=$processors -S $conf{filter_psort_buffer} -T $temp_dir/$sample/temp/psort $tosort_file $LOG"
#			  . " | $scr_dir/MOCATFilter_stats.pl -db $full_name -format SOAP -length $conf{filter_length_cutoff} -identity $conf{filter_percent_cutoff} -stats $stats $LOG"
#			  . " | $ZIP -$ziplevel -c > $output_file.soap.gz.tmp"
#			  . " && rm -f $tosort_file"
			  . " && sync"    # Important for NFS
			  . " && test -e $output_file.soap.gz.tmp" . " && mv $output_file.soap.gz.tmp $output_file.soap.gz && test -e $output_file.soap.gz";
		}
		if ( $conf{MOCAT_mapping_mode} eq "allbest" ) {
			print JOB " | $scr_dir/MOCATFilter_soap2sam.awk -v MIN_LEN=$conf{filter_length_cutoff} -v MIN_AS=$conf{filter_percent_cutoff} $LOG ";
			if ( $conf{filter_paired_end_filtering} eq 'yes' ) {
				print JOB "| $scr_dir/MOCATFilter_filterPE.pl $LOG ";
			}
			print JOB " | $scr_dir/MOCATFilter_besthit.pl $LOG ";
			print JOB " | $scr_dir/MOCATFilter_stats.pl -db $full_name -format SAM -length $conf{filter_length_cutoff} -identity $conf{filter_percent_cutoff} -stats $stats $LOG ";
			print JOB " | $bin_dir/msamtools$OSX -Sb -m merge -t $len_file - $LOG > $output_file.bam.tmp";
			print JOB " && sync  && test -e $output_file.bam.tmp && mv $output_file.bam.tmp $output_file.bam && test -e $output_file.bam ";
#			print JOB " | $bin_dir/msamtools$OSX " . "-S "
#			  . "-m filter "
#			  . "-l $conf{filter_length_cutoff} "
#			  . "-p $conf{filter_percent_cutoff} " . "-z 0 "
#			  . "--besthit "
#			  . "-t $len_file "
#			  . "- $LOG"
#			  . " | $scr_dir/MOCATFilter_stats.pl "
#			  . " -db $full_name "
#			  . "-format SAM "
#			  . "-length $conf{filter_length_cutoff} "
#			  . "-identity $conf{filter_percent_cutoff} "
#			  . "-stats $stats " . "$LOG "
#			  . " | $bin_dir/msamtools$OSX " . "-Sb "
#			  . "-m merge "
#			  . "-t $len_file "
#			  . "- $LOG > $output_file.bam.tmp"
#			  . " && rm -f $tosort_file"
#			  . " && sync "    # Wait for the system to equilibrate [sometimes important in network filesystems]
#			  . " && test -e $output_file.bam.tmp" . " && mv $output_file.bam.tmp $output_file.bam" . " && test -e $output_file.bam ";
			if ( $conf{filter_make_unique_sorted} eq 'yes' ) {
				print JOB " && $bin_dir/samtools view $output_file.bam" . "| $bin_dir/msamtools$OSX" . " -Sb" . " -m filter " . " -l $conf{filter_length_cutoff} " . " -p $conf{filter_percent_cutoff} " . " -z 0 " . " --uniqhit " . " -t $len_file - $LOG " . " | $bin_dir/samtools sort -m $conf{filter_samtools_memory} - > $output_file.unique.sorted" . " && rm -f $len_file";
			}
			else {
				print JOB " && rm -f $len_file";
			}
		}
		if ( $conf{filter_remove_mapping_files} eq 'yes' ) {
			print JOB " && rm $input_file\n";
		}
		else {
			print JOB "\n";
		}
	}    # End loop samples
	close JOB;
	print " OK!\n";
}
1;
