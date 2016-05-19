package MOCATResistanceScreen;
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
	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?";
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

	if ($identity2 eq "") {
		$identity2 = $conf{filter_percent_cutoff};
	}
	
	my $folder = "$data_dir";
	if ($resistance_screen_suffix ne '') {
		$resistance_screen_suffix = ".$resistance_screen_suffix";
		unless (-e "$folder/ResistanceScreenData.genes$resistance_screen_suffix") {
			die "ERROR & EXIT: Missing file $folder/ResistanceScreenData.genes$resistance_screen_suffix for suffix $resistance_screen_suffix"
		}
		unless (-e "$folder/ResistanceScreenData.drug2clusters$resistance_screen_suffix") {
			die "ERROR & EXIT: Missing file $folder/ResistanceScreenData.drug2clusters$resistance_screen_suffix for suffix $resistance_screen_suffix"
		}
		unless (-e "$folder/ResistanceScreenData.referencepopulation$resistance_screen_suffix") {
			die "ERROR & EXIT: Missing file $folder/ResistanceScreenData.referencepopulation$resistance_screen_suffix for suffix $resistance_screen_suffix"
		}
	}
	


	# Define variables
	my $databases = join( "_AND_", @do_resistance_screen );
	my $cur_sp_tax = "$cwd/taxonomic.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}/$sample_file_basename.taxonomic.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.base.norm.specI_clusters";
	my $output_folder = "$cwd/resistance.potential/$read_type2.$resistance_screen_reads2.on.$resistance_screen_refgenecatalog/$read_type.$reads.on.$databases";
	$output_file = "$output_folder/$sample_file_basename.resistance.potential$resistance_screen_suffix";
	my $output_dataframe = "$output_folder/$sample_file_basename.resistance.potential$resistance_screen_suffix.dataframe";
	system "mkdir -p $output_folder";
	
	# NOTE THIS IS HARD CODED FOR NOW AS IT REQUIRES 263REFGENECATALOG
	my $base_cov   = "$cwd/abundance.tables/$resistance_screen_refgenecatalog/$read_type2.$resistance_screen_reads2.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$identity2/$sample_file_basename.$read_type2.$resistance_screen_reads2.on.$resistance_screen_refgenecatalog.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$identity2.base.coverage.raw";
	
	# Check files
	print localtime() . ": Checking rownames file...";
	unless (-e "$data_dir/$resistance_screen_refgenecatalog.rownames") {
		die "\nERROR & EXIT: Missing $data_dir/$resistance_screen_refgenecatalog.rownames";
	}
	print " OK!\n";
	print localtime() . ": Checking taxonomic specI_clusters file...";
	unless (-e "$cur_sp_tax.gz") {
		die "\nERROR & EXIT: Missing $cur_sp_tax.gz\nHave you run -taxo_profiling using this specific sample file ($sample_file_basename)?";
	}
	print " OK!\n";
	if ($resistance_screen_reads2 eq "") {
		die "ERROR & EXIT: -r2 [READS ORIGIN] has to be the reads you used when screening against $resistance_screen_refgenecatalog";
	}
	print localtime() . ": Checking base coverage file...";
	unless (-e "$base_cov.gz") {
		die "\nERROR & EXIT: Missing $base_cov.gz\nHave you mapped reads against $resistance_screen_refgenecatalog and then run -pcf?";
	}
	print " OK!\n";
	foreach my $file ("ResistanceScreenData.drug2clusters","ResistanceScreenData.genes","ResistanceScreenData.genomesize","ResistanceScreenData.marker.lengths","ResistanceScreenData.readlengths","ResistanceScreenData.referencepopulation","ResistanceScreenData.refgenelengths") {
		unless (-e "$folder/$file") {
			die "ERROR & EXIT: Missing $folder/$file";
		}
	}
	
	print localtime() . ": Creating read length file...";
	system "mkdir -p $cwd/resistance_screen_temp/";
	open OUT, ">$cwd/resistance_screen_temp/avg.readlengths" or die "\nERROR & EXIT: Cannot write to $cwd/resistance_screen_temp/avg.readlengths";
	foreach my $sample (@samples) {
		( my $max, my $avg, my $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-r" );
		print OUT "$sample\t$avg\n";
	}
	my $LOG       = "2>> $cwd/logs/$job/samples/MOCATJob_$job.MAIN_OUPUT.log >> $cwd/logs/$job/samples/MOCATJob_$job.MAIN_OUPUT.log";
	close OUT;
	print " OK!\n";
	print localtime() . ": Printing job file...";
	
	# Static and 263RefGeneCatalog dependent
	print JOB "$conf{MOCAT_DEV_DIR}/ResistanceScreen/screenSamplesForResistance.pl -refgenelengthfile=$folder/ResistanceScreenData.refgenelengths -refgene2resgenefile=$folder/ResistanceScreenData.genes$resistance_screen_suffix";
	print JOB " -resgene2speciesfile=$folder/ResistanceScreenData.drug2clusters$resistance_screen_suffix -genomesizefile=$folder/ResistanceScreenData.genomesize -markergenelengthfile=$folder/ResistanceScreenData.marker.lengths ";
	print JOB " -genomesizefile=$folder/ResistanceScreenData.genomesize -referencepopulationfile=$folder/ResistanceScreenData.referencepopulation$resistance_screen_suffix -rownamefile=$data_dir/$resistance_screen_refgenecatalog.rownames";
	
	# Variable
	# omitted in new version : -markergenelengthfile=$folder/ResistanceScreenData.marker.lengths
	print JOB " -speciesbasecovfile=$cur_sp_tax.gz";
	print JOB " -genebasecovfile=$base_cov.gz";
	print JOB " -readlengthfile=$cwd/resistance_screen_temp/avg.readlengths";
	print JOB " -outputfile=$output_file";
	print JOB " -outputdataframe=$output_dataframe $LOG\n";
	print " OK!\n";
	print localtime() . ": Temp directory is $temp_dir/resistance_screen_temp\n";	
}


1;
