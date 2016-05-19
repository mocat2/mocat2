package MOCATSummarizeCarmaCoverages;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

sub run {

	print "MOCATSummarizeCarmaCoverages : starting\n";

	### DEFINE VARIABLES AND OPEN input FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
		
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}
	my %taxa;

	### JOB SPECIFIC ###
	my $db = join( "_AND_", @do_summarize_carma_coverages );
	foreach my $sample (@samples) {
		print "MOCATSummarizeCarmaCoverages : processing $sample\n";
		foreach my $start ( "base.coverage", "insert.coverage" ) {
			my $input = "$sample.filtered.$read_type.$reads.on.$db.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$start";
			foreach my $end ( "norm", "count" ) {
				#print "MOCATSummarizeCarmaCoverages : processing $sample - $start - $end\n";
				my $inputfolder = "$cwd/$sample/taxonomic.profiles/$input/$start.$end";
				unless ( -e "$inputfolder/$sample.$end.taxonomic.composition.tab" ) {
					die "ERROR & EXIT: Missing $inputfolder/$sample.$end.taxonomic.composition.tab";
				}
				open IN, "<$inputfolder/$sample.$start.$end.taxonomic.composition.tab";
				my $counter = 0;
				while (<IN>) {
					$counter++;
					chomp;
					my @line = split /\t/, $_;
					$line[1] =~ s/\"//g;
					$taxa{$start}{$end}{ $line[0] }{ $line[1] }{$sample}{"count"}    = $line[2];
					$taxa{$start}{$end}{ $line[0] }{ $line[1] }{$sample}{"fraction"} = $line[3];
				}
				#print "MOCATSummarizeCarmaCoverages : processing $sample - $start - $end - loaded 2 x $counter entries\n";
			}
		}
	}
	print "MOCATSummarizeCarmaCoverages : all files processed\n";
	print "MOCATSummarizeCarmaCoverages : creating folders\n";
	my $folder = "$cwd/taxonomic.profiles/carma.annotations/$carma_file/$read_type.$reads.on.$db.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
	system "mkdir -p $folder/base.coverage.norm";
	system "mkdir -p $folder/base.coverage.count";
	system "mkdir -p $folder/insert.coverage.norm";
	system "mkdir -p $folder/insert.coverage.count";

	foreach my $start ( "base.coverage", "insert.coverage" ) {		
		foreach my $end ( "norm", "count" ) {
			print "MOCATSummarizeCarmaCoverages : processing $start.$end\n";
			foreach my $genera ( keys %{$taxa{$start}{$end}} ) {
				#print "MOCATSummarizeCarmaCoverages : printing $start - $end - $genera\n";
				open COUNT, ">$folder/$start.$end/$sample_file.$read_type.$reads.on.$db.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.taxo_profile.$start.$end.$genera";
				print COUNT "$genera\t" . join( "\t", @samples ) . "\n";
				open FRAC, ">$folder/$start.$end/$sample_file.$read_type.$reads.on.$db.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.taxo_profile.$start.$end.$genera.fraction";
				print FRAC "$genera\t" . join( "\t", @samples ) . "\n";
				foreach my $entry ( sort keys %{$taxa{$start}{$end}{$genera}} ) {
					my $line_to_print_count = "";
					my $line_to_print_fraction = "";
					foreach my $sample (@samples) {
						foreach my $type ( "count", "fraction" ) {
							unless ( $taxa{$start}{$end}{$genera}{$entry}{$sample}{$type} ) {
								$taxa{$start}{$end}{$genera}{$entry}{$sample}{$type} = 0;
							}
							if ( $type eq "count" ) {
								$line_to_print_count = $line_to_print_count . "$taxa{$start}{$end}{$genera}{$entry}{$sample}{$type}\t";
							}
							if ( $type eq "fraction" ) {
								$line_to_print_fraction = $line_to_print_fraction . "$taxa{$start}{$end}{$genera}{$entry}{$sample}{$type}\t";
							}
						}
					}
					$line_to_print_count =~ s/\t$/\n/;
					$line_to_print_fraction =~ s/\t$/\n/;
					print COUNT "$entry\t$line_to_print_count";
					print FRAC  "$entry\t$line_to_print_fraction";
				}
			}
		}
	}
	print "MOCATSummarizeCarmaCoverages : completed\n";
}

1;
