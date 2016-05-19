package MOCATExtractGenes_old;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

sub run {

	print localtime() . ": EXECUTING:\n";
	my $genesFile        = $do_extract_genes;
	my $genesContigsFile = $extract_genes_scaftig;

	my @genesFile_name = split "/", $genesFile;
	unless ( -e $genesContigsFile ) {
		#die "ERROR & EXIT: Missing $genesContigsFile";
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
	open OUT,   ">$cwd/genes.extracted/$genesFile_name[-1].padded";
	open TABLE, ">$cwd/genes.extracted/$genesFile_name[-1]..padded.coord";
	open PM, ">$cwd/genes.extracted/$genesFile_name[-1].gene2contig";
	print PM "#gene_id\tstart|stop\tcontig_id\tstart-100bp\tstop+100bp\n";
	my $counter = 0;
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
		elsif (m/^>.*complete/) {
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

		chomp( my $seq = `$bin_dir/cdbyank $genesContigsFile.cidx -a \"$name\t$newstart\t$newstop\" -R` );

		#my $seq="$name\t$newstart\t$newstop";
		$seq =~ s/^>(.*)/>$gene\t$pos_start\t$pos_stop/;
		### Additional info :: gene:$gene cds_start:$start cds_stop:$stop extracted_start:$newstart extracted_stop:$newstop here_start:$pos_start here_stop:$pos_stop strand: $strand
		if ( $counter % 100 == 0 ) {
			print ".";
		}
		print OUT "$seq\n";
		print TABLE "$gene $pos_start $pos_stop\n";
		print PM "$gene\t$pos_start\t$pos_stop\t$name\t$newstart\t$newstop\n";
	}

	print " DONE!\n";
	close IN;
	close OUT;
	close TABLE;
	close PM;
}

1;
