package MOCATPasteTaxonomyCoverageFiles;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

sub run {
	### DEFINE VARIABLES AND OPEN OUTPUT FILE ###
	my $job        = $_[0];
	my $processors = $_[1];

	my $read_type = 'screened';
	if ($use_extracted_reads) {
		$read_type = 'extracted';
	}

	my $uniq = "";
	if ($pcf_uniq) {
		$uniq = ".uniq";
	}

	my $databases = join( "_AND_", @databases );
	my $assembly_type = 'assembly';
	my $end;
	my $assembly = 0;
	chomp( my $bn = `basename $sample_file` );
	if (         $databases[0] eq 's'
		|| $databases[0] eq 'c'
		|| $databases[0] eq 'f'
		|| $databases[0] eq 'r' )
	{
		if ( $databases[0] eq 's' ) { $end = 'scaftig' }
		if ( $databases[0] eq 'c' ) { $end = 'contig' }
		if ( $databases[0] eq 'f' ) { $end = 'scafSeq' }
		if ( $databases[0] eq 'r' ) {
			$assembly_type = 'assembly.revised';
			$end           = 'scaftig';
		}
		$assembly  = 1;
		$databases = "$assembly_type.$end";
	}

	system "mkdir -p $cwd/taxonomic.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
	my @levels = ( 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'curated.species' );
	foreach my $boi ( ( 'base', 'insert' ) ) {
		foreach my $norm ( ( 'raw', 'norm', 'only.unique.raw', 'only.unique.norm' ) ) {
			foreach my $i (@levels) {

				print localtime() . ": Preparing $boi:$norm:$i";
				my $to_paste = '';
				my $rownames;
				my $folder;
				my $file;
				my %hash;
				my $counter1 = 0;
				my $counter2 = 0;
				foreach my $sample (@samples) {
					$counter1++;
					if ( !$assembly ) {
						$folder = "$cwd/$sample/taxonomic.profiles.$databases.$conf{MOCAT_data_type}";

						#$rownames = "$folder/$sample.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.rownames$uniq";
						#$file     = "$folder/$sample.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$boi.coverage.$norm";
						$file     = "$folder/$sample.taxonomic.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$boi.$norm.$i";
						$rownames = "$folder/$sample.taxonomic.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$i.rownames";

					}
					elsif ($assembly) {
						( my $max, my $avg, my $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-r" );
						$folder = "$cwd/$sample/$boi.coverage.$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end";

						#$rownames = "$folder/$sample.filtered.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.rownames$uniq";
						#$file     = "$folder/$sample.filtered.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$boi.coverage.$norm";
						$file     = "$folder/$sample.taxonomic.profile.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$boi.$norm.$i";
						$rownames = "$folder/$sample.taxonomic.profile.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$i.rownames";

					}
					unless ( -e $file ) {
						die "\nERROR & EXIT: Missing $file\n";
					}
					unless ( $hash{$counter2} ) {
						$hash{$counter2} = "";
					}
					$hash{$counter2} = $hash{$counter2} . " $file ";
					if ( $counter1 == 100 ) {
						$counter1 = 0;
						$counter2++;
					}
					print ".";
				}
				my $norm2 = $norm;
				if ( $norm eq 'count' ) {
					$norm2 = 'raw';
				}
				unless ( -e "$rownames" ) {
					die "\nERROR & EXIT: Missing $rownames\n";
				}
				print " OK!\n";
				print localtime() . ": Processing $boi:$norm:$i...";

				my $name = "$cwd/taxonomic.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}/$sample_file_basename.taxonomic.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$boi.$norm.$i";
				system "cp $rownames $name";
				foreach my $f ( sort { $a <=> $b } keys %hash ) {
					system "paste $name $hash{$f} > $name.tmp && mv $name.tmp $name";
					print ".";
				}
				system "rm -f $name.tmp";
				system "$ZIP -f -$ziplevel $name";
				print " OK!\n";
			}
		}
	}
	print localtime() . ": Completed.\n";
	print localtime() . ": RESULTS SAVED IN $cwd/taxonomic.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}\n";
}

1;
