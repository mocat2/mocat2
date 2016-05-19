package MOCATImportOld;
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
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	my $db = $_[3];
	
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );

	my $mode = $conf{MOCAT_data_type};

	### JOB SPECIFIC ###

	my $counter;

	print "Processing samples:\n";
	foreach my $sample (@samples) {
		unless ( $sample eq "" ) {

			$counter++;
			print localtime() . ": [$counter/" . scalar @samples . "] Processing $sample...";

			# Get Kmer
			my $statsFile;
			$statsFile = "$sample/$sample.stats_read_trim_filter_$mode";
			if ( -s "$statsFile" ) {
				open IL, '<', "$statsFile";
			}
			else {
				die "\nERROR & EXIT: Missing insert file $statsFile";
			}
			my $temp;
			while (<IL>) {
				chomp;
				$temp = $_;
			}
			close IL;
			my @temp = split( /\t/, $temp );
			my $kmer = $temp[4];

			chomp( my $file = `ls -1 $sample/assembly.$mode.K$kmer/$sample.scaftig` );
			$file =~ s/.*$sample\/assembly.$mode/\.\.\/assembly.$mode/;
			chomp( my $file2 = `ls -1 $sample/assembly.$mode.K$kmer/$sample.scafSeq` );
			system "mkdir -p $sample/assembly.FROMOLDi8.$mode.K$kmer/";
			system "mkdir -p $sample/stats";
			system "mkdir -p $sample/reads.screened.FROMOLDi8.$mode";

			if ( $db ne "NO" ) {
				system "ln -fs $file $sample/assembly.FROMOLDi8.$mode.K$kmer/$sample.assembly.FROMOLDi8.$mode.K$kmer.scaftig.old";
				system "perl $scr_dir/MOCATAssembly_getScaf.pl $db $file2 $sample > $sample/assembly.FROMOLDi8.$mode.K$kmer/$sample.assembly.FROMOLDi8.$mode.K$kmer.scaftig";
			}
			else {
				system "ln -fs $file $sample/assembly.FROMOLDi8.$mode.K$kmer/$sample.assembly.FROMOLDi8.$mode.K$kmer.scaftig";
			}

			open IN,  '<', "$sample/$sample.stats_assembly_insert_sizes";
			open OUT, '>', "$sample/stats/$sample.assembly.FROMOLDi8.$mode.K$kmer.inserts.stats";
			while (<IN>) {
				chomp $_;
				my @line = split( "\t", $_ );
				print OUT "$line[0].screened.FROMOLDi8\t$line[1]\n";
			}
			close IN;
			close OUT;

			system "ln -fs ../$sample.stats_read_trim_filter_$mode $sample/stats/$sample.screen.FROMOLDi8.$mode.stats";
			system " for f in $sample/reads.processed.$mode/*gz; do h=`echo \$f | sed 's/.*reads.processed.$mode\\///'`; g=`echo \$f | sed 's/.*reads.processed.$mode\\///' | sed 's/all.fq.gz/screened.FROMOLDi8.all.fq.gz/' | sed 's/.single/.screened.FROMOLDi8.single/' | sed 's/.pair/.screened.FROMOLDi8.pair/'`; ln -fs ../reads.processed.$mode/\$h $sample/reads.screened.FROMOLDi8.$mode/\$g; done";
			print " OK!\n";
		}
	}
	print "All samples processed.\n";
}

1;
