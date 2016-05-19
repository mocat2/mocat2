package MOCATExtractGenes;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.




##### 2012 12 05
##### Note this version doesn not function, because stored in the final DB file, are the scaftig IDs and not the GENE IDs.





sub run {

	print localtime() . ": EXECUTING:\n";
	my $genesFile        = $do_extract_genes;
	my $genesContigsFile = $extract_genes_scaftig;

	my @genesFile_name = split "/", $genesFile;
	unless ( -e $genesContigsFile ) {
		die "ERROR & EXIT: Missing $genesContigsFile";
	}
	unless ( -e $genesFile ) {
		die "ERROR & EXIT: Missing $genesFile";
	}
	unless ( -e "$genesContigsFile.cidx" ) {
		print localtime() . ": Indexing contig files...";
		system "$bin_dir/cdbfasta $genesContigsFile 2>/dev/null";
		print " OK!\n";
	}

	print localtime() . ": Grepping genes...";
	system "mkdir -p $cwd/genes.extracted";
	open IN, '<', "$genesFile";
	system "rm -fr $cwd/genes.extracted/$genesFile_name[-1]";
	open TABLE, ">$cwd/genes.extracted/$genesFile_name[-1].coord";
	open PM, ">$cwd/genes.extracted/$genesFile_name[-1].gene2contig";
	print PM "#gene_id\tstart|stop\tcontig_id\tstart-100bp\tstop+100bp\n";
	my $counter = 0;
	my @to_run  = ();
	while (<IN>) {
		$counter++;
		chomp;
		my $line = $_;
		my $name;
		my $start;
		my $stop;
		my $strand;
		my $gene;

		if (m/^>MC\d+\.MG.*/) {
			my @F = split /\s+/, $line;
			$F[0] =~ s/^>//;
			$gene   = $F[0];
			$line   = m/>(\S+)\s+.*start:(\d+)\s+end:(\d+)\s+strand (.)/;
			$name   = $1;
			$start  = $2;
			$stop   = $3;
			$strand = $4;
		}
		elsif (m/^>.*complete$/) {
			$line =~ m/>(\S+)_gene(\S+) strand:([+-]) start:(\d+) stop:(\d+) length:(\d+) start_codon:([a-z]+) stop_codon:([a-z]+) gene_type:([a-z]*)/;
			$name   = $1;
			$gene   = "$1_gene$2";
			$start  = $4;
			$stop   = $5;
			$strand = $3;
		}
		else {
			next;
		}

		my $newstart = $start - 100;
		if ( $newstart < 1 ) { $newstart = 1; }
		my $newstop   = $stop + 100;
		my $pos_start = $start - $newstart + 1;
		my $pos_stop  = $stop - $start + $pos_start;
		$name =~ s/GP\d+\.//;
		$name =~ s/\.G\d+//;

		push @to_run, "$name $gene\t$newstart\t$newstop";
		print TABLE "$gene\t$pos_start\t$pos_stop\n";
		print PM "$gene\t$pos_start\t$pos_stop\t$name\t$newstart\t$newstop\n";

		if ( $counter == 100000 ) {
			$counter = 0;
			system "rm -fr $cwd/genes.extracted/tmp";
			open OUT, ">$cwd/genes.extracted/tmp";
			foreach my $line (@to_run) {
				print OUT "$line\n";
			}
			close OUT;
			system "cat $cwd/genes.extracted/tmp | cdbyank -Q -R $genesContigsFile.cidx | sed -e 's/\\t>.*//' -e 's/^%/>/' -e s'/%\$//' >> $cwd/genes.extracted/$genesFile_name[-1]";
			@to_run = ();
		}
	}
	system "rm -fr $cwd/genes.extracted/tmp";
	open OUT, ">$cwd/genes.extracted/tmp";
	foreach my $line (@to_run) {
		print OUT "$line\n";
	}
	close OUT;
	system "cat $cwd/genes.extracted/tmp | cdbyank -Q -R $genesContigsFile.cidx >> $cwd/genes.extracted/$genesFile_name[-1]";
	system "rm -fr $cwd/genes.extracted/tmp";
	print " DONE!\n";
	close IN;
	close PM;
	close TABLE;

}

1;
