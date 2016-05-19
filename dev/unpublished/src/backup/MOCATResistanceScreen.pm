package MOCATResistanceScreen;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

sub create_job {
	### DEFINE VARIABLES AND OPEN input FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?\n";

	print localtime() . ": Creating $job jobs...\n";
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}
	my $read_type2 = 'screened';
	if ($use_extracted_reads2) {
		$read_type2 = 'extracted';
	}

	# Define variables
	my $databases = join( "_AND_", @do_resistance_screen );
	my $cur_sp_tax = "$cwd/taxonomic.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}/$sample_file_basename.taxonomic.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.base.raw.curated.species";;
	my $folder = "$data_dir";
	my $output_folder = "$cwd/resistance.potential/$read_type2.$resistance_screen_reads2.on.263RefGeneCatalog.padded/$read_type.$reads.on.$databases";
	$output_file = "$output_folder/$sample_file_basename.resistance.potential.tab";
	system "mkdir -p $output_folder";
	
	# NOTE THIS IS HARD CODED FOR NOW AS IT REQUIRES 263REFGENECATALOG
	my $base_cov   = "$cwd/abundance.tables/263RefGeneCatalog.padded/$sample_file_basename.$read_type2.$resistance_screen_reads2.on.263RefGeneCatalog.padded.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.base.coverage.raw";
	
	# Check files
	print localtime() . ": Checking taxonomic curated.species file...";
	unless (-e $cur_sp_tax) {
		die "\nERROR & EXIT: Missing $cur_sp_tax\nHave you run -taxo_profiling using this specific sample file ($sample_file_basename)?\n";
	}
	print " OK!\n";
	if ($resistance_screen_reads2 eq "") {
		die "ERROR & EXIT: -r2 [READS ORIGIN] has to be the reads you used when screening against 263RefGeneCatalog.padded\n";
	}
	print localtime() . ": Checking base coverage file...";
	unless (-e $base_cov) {
		die "\nERROR & EXIT: Missing $base_cov\nHave you mapped reads against 263RefGeneCatalog.padded and then run -pcf?\n";
	}
	print " OK!\n";
	foreach my $file ("ResistanceScreenData.drug2species","ResistanceScreenData.genes","ResistanceScreenData.genomesize","ResistanceScreenData.marker.lengths","ResistanceScreenData.output","ResistanceScreenData.readlengths","ResistanceScreenData.referencepopulation","ResistanceScreenData.refgenelengths") {
		unless (-e "$folder/$file") {
			die "ERROR & EXIT: Missing $folder/$file\n";
		}
	}
	
	print localtime() . ": Creating read length file...";
	system "mkdir -p $temp_dir/resistance_screen_temp/";
	open OUT, ">$temp_dir/resistance_screen_temp/avg.readlengths" or die "\nERROR & EXIT: Cannot write to $temp_dir/resistance_screen_temp/avg.readlengths\n";
	foreach my $sample (@samples) {
		( my $max, my $avg, my $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-r" );
		print OUT "$sample\t$avg\n";
	}
	close OUT;
	print " OK!\n";
	print localtime() . ": Printing job file...";
	
	# Static and 263RefGeneCatalog dependent
	print JOB "$conf{MOCAT_DEV_DIR}/ResistanceScreen/screenSamplesForResistance.pl -refgenelengthfile=$folder/ResistanceScreenData.refgenelengths -refgene2resgenefile=$folder/ResistanceScreenData.genes";
	print JOB " -resgene2speciesfile=$folder/ResistanceScreenData.drug2species -genomesizefile=$folder/ResistanceScreenData.genomesize -markergenelengthfile=$folder/ResistanceScreenData.marker.lengths -referencepopulationfile=$folder/ResistanceScreenData.referencepopulation";
	print JOB " -genomesizefile=$folder/ResistanceScreenData.genomesize -markergenelengthfile=$folder/ResistanceScreenData.marker.lengths -referencepopulationfile=$folder/ResistanceScreenData.referencepopulation -rownamefile=$folder/ResistanceScreenData.rownames";
	
	# Variable
	print JOB " -speciesbasecovfile=$cur_sp_tax";
	print JOB " -genebasecovfile=$base_cov";
	print JOB " -readlengthfile=$temp_dir/resistance_screen_temp/avg.readlengths";
	print JOB " -outputfile=$output_file\n";
	print " OK!\n";
	print localtime() . ": Temp directory is $temp_dir/SAMPLE/temp\n";	
}


1;
