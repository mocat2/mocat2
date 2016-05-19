#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

################################
#USAGE
my $usage = "\nUSAGE: GenerateSummary.pl -m|mode (RefMG | mOTU) OPTIONS\n
mode = RefMG
OPTIONS:
\t-c <count.data> (required)                               tab-delimieted file: <GeneID>\\t<Count>
\t-t <taxo.annotation.map> (required)                      tab-delimited file: <TaxID>\\t<Kingdom>\\t<Phylum>\\t<Class>\\t<Order>\\t<Family>\\t<Genus>\\t<Species>\\t<CuratedSpecies>
\t-l <lengths.concatenated.cogs> (required)                tab-delimited file: <TaxID>\\t<Length_of_concatenated_COGs(in_bp)>
EXAMPLE: GenerateSummary.pl -m RefMG -c sample.count -t RefMGv9.taxid.annot.curated -l RefMGv9.taxid.len\n
mode = mOTU
OPTIONS:
\t-c <count.data> (required)                                tab-delimieted file: <GeneID>\\t<Count>
\t-g <gene.annotation.map> (required)                       tab-delimited file:<GeneID>\\t<GeneLength>\\t<COGID>\\t<OTU>
EXAMPLE: GenerateSummary.pl -m mOTU -c sample.count -g RefMGv9.map\n
\n";

################################
#Defaults

################################
#Global variables
my ( $mode, $count, $taxomap, $concatlen, $genemap, $ogmap );
my $all_tax_avg_len = 0;
my ( @cog,    @abund,        @sample );
my ( %GeneID, %GeneIDLength, %COGID, %OTU, %OGID, %OG );
my ( %TaxID,  %Kingdom,      %Phylum, %Class, %Order, %Family, %Genus, %Species, %CuratedSpecies );
my %tax_avg_len;
my %cog_avg_len;

################################
#Get options
GetOptions(
	'm=s' => \$mode,
	'c=s' => \$count,
	't=s' => \$taxomap,
	'l=s' => \$concatlen,
	'g=s' => \$genemap
);

################################
#Parse input
unless ( $mode eq "RefMG" || $mode eq "mOTU" ) { system "clear"; die "set mode to either RefMG or mOTU" . $usage }
if ( $mode eq "RefMG" ) {
	unless ( -e $count && -e $taxomap && -e $concatlen ) { system "clear"; die $usage }
}
elsif ( $mode eq "mOTU" ) {
	unless ( -e $count && -e $genemap ) { system "clear"; die $usage }
}
elsif ( $mode eq "cog" ) {
	unless ( -e $count && -e $ogmap ) { system "clear"; die $usage }
}
else {
	system "clear";
	die $usage;
}

################################
#Start
#system "clear";

#Mode=RefMG
if ( $mode eq "RefMG" ) {

	#Load taxonomic annotations
	print "Loading taxonomic annotations...";
	open( TAXO, "<$taxomap" ) or die "MISSING $taxomap\n";
	while (<TAXO>) {
		my @line = split "\t", $_;
		chomp @line;
		$TaxID{ $line[0] }          = 0;
		$Kingdom{ $line[0] }        = $line[1];
		$Phylum{ $line[0] }         = $line[2];
		$Class{ $line[0] }          = $line[3];
		$Order{ $line[0] }          = $line[4];
		$Family{ $line[0] }         = $line[5];
		$Genus{ $line[0] }          = $line[6];
		$Species{ $line[0] }        = $line[7];
		$CuratedSpecies{ $line[0] } = $line[8];
	}
	close TAXO;
	print "done\n";

	#Generate annotation file if it does not exist yet
	print "Generating rownames file if it does not exist...";
	if ( -e "$taxomap.rownames" ) {
		print "$taxomap.rownames already exists.\n";
	}
	else {
		open( ROW, ">$taxomap.rownames" );
		my $header = "TaxID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tCuratedSpecies\n";
		$header .= "-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\n";
		print ROW $header;
		foreach my $k ( sort keys %TaxID ) {
			print ROW "$k\t$Kingdom{$k}\t$Phylum{$k}\t$Class{$k}\t$Order{$k}\t$Family{$k}\t$Genus{$k}\t$Species{$k}\t$CuratedSpecies{$k}\n";
		}
		print "done.\n";
	}
	close ROW;

	#Load concatenated gene lengths for tax ids
	print "Loading concatenated gene lengths for taxonomy identifiers...";
	open( AVGLEN, "<$concatlen" ) or die "MISSING $concatlen\n";
	while (<AVGLEN>) {
		my @line = split "\t", $_;
		chomp $_;
		$tax_avg_len{ $line[0] } = $line[1];
		$all_tax_avg_len += $line[1];
	}
	close AVGLEN;
	$all_tax_avg_len = $all_tax_avg_len / scalar( keys %tax_avg_len );
	print "done\n";

	#Generate output files
	print "Parsing count data...";
	open( COUNT, "<$count" );
	while (<COUNT>) {
		my @line = split "\t", $_;
		chomp @line;
		if ( $_ =~ /^-1/ ) {
			$TaxID{ $line[0] } = $line[1];    #Add -1
		}
		elsif ( $_ !~ /^ref_id/ ) {
			my @tax = split /\./, $line[0];
			$TaxID{ $tax[0] } += $line[1];
		}
	}
	close COUNT;
	print "done.\n";

	print "Generating output files...";
	my @out = split /\./, $count;
	pop @out;
	my $output = join( "\.", @out );
	my $sample;
	if ( $output =~ /\// ) {
		my @a    = split /\//, $output;
		my $file = $a[-1];
		my @b    = split /\.filtered/, $file;
		$sample = $b[0];
	}
	else {
		$sample = $out[0];
	}

	open( NORM, ">$output.by.taxonomy.norm" )   or die "Can't open output file\n";
	open( RAW,  ">$output.by.taxonomy.raw" )    or die "Can't open output file\n";
	open( SCA,  ">$output.by.taxonomy.scaled" ) or die "Can't open output file\n";
	print RAW "$sample\n";
	print NORM "$sample\n";
	print SCA "$sample\n";
	my $weighted_sum = 0;
	my $total_sum    = 0;
	my $weighted_avg = 0;
	foreach my $k ( sort keys %TaxID ) {
		my $coverage_raw = $TaxID{$k};
		if ( $k =~ /-1/ ) {
			$tax_avg_len{$k} = $all_tax_avg_len;
		}
		unless ( defined( $tax_avg_len{$k} ) ) {
			die "ERROR & EXIT: DATABASE ENTRY $k HAS MAPPED TAXA ID $TaxID{$k}.\nThis indicates that perhaps not all entries in the database (rownames file) are on the format TAXAID.xxx.\nIf This isn't the case, mapping to TAXA ID isn't possible and this whole pipeline won't work.\n";
		}
		my $coverage_norm = $TaxID{$k} / $tax_avg_len{$k};
		$weighted_sum += $TaxID{$k};
		$total_sum += $coverage_norm;
		print RAW "$coverage_raw\n";
		print NORM "$coverage_norm\n";
		$coverage_raw  = 0;
		$coverage_norm = 0;
	}
	$weighted_avg = $weighted_sum / $total_sum;
	foreach my $k ( sort keys %TaxID ) {
		my $coverage_norm = $TaxID{$k} / $tax_avg_len{$k};
		my $coverage_sca  = $coverage_norm * $weighted_avg;
		print SCA "$coverage_sca\n";
	}
	print "done.\n";
	close NORM;
	close RAW;
	close SCA;
}

#Mode=mOTU
if ( $mode eq "mOTU" ) {

	#Load gene annotations
	print "Loading gene annotations...";
	open( MAP, "<$genemap" );
	while (<MAP>) {
		my @line = split "\t", $_;
		chomp @line;
		$GeneID{ $line[0] }       = 0;
		$GeneIDLength{ $line[0] } = $line[1];
		$COGID{ $line[0] }        = $line[2];
		$OTU{ $line[0] }          = "$line[2].$line[3]";
	}
	close MAP;
	print "done\n";

	#Generate annotation file if it does not exist yet
	print "Generating rownames file if it does not exist...";
	if ( -e "$genemap.rownames" ) {
		print "$genemap.rownames already exists.\n";
	}
	else {
		open( ROW, ">$genemap.rownames" ) or die "MISSING $genemap.rownames\n";
		my $header = "GeneID\tCOGID\tOTU\n";
		print ROW $header;
		foreach my $k ( sort keys %GeneID ) {
			print ROW "$k\t$COGID{$k}\t$OTU{$k}\n";
		}
		print "done.\n";
	}
	close ROW;

	print "Parsing count data...";
	open( COUNT, "<$count" ) or die "MISSING $count\n";
	while (<COUNT>) {
		my @line = split "\t", $_;
		chomp @line;
		if ( $_ =~ /^-1/ ) {
			next;
		}
		elsif ( $_ !~ /^ref_id/ ) {
			$GeneID{ $line[0] } += $line[1];
		}
	}
	close COUNT;
	print "done.\n";

	print "Generating output files...";
	my @out = split /\./, $count;
	pop @out;
	my $output = join( "\.", @out );
	my $sample;
	if ( $output =~ /\// ) {
		my @a    = split /\//, $output;
		my $file = $a[-1];
		my @b    = split /\.filtered/, $file;
		$sample = $b[0];
	}
	else {
		$sample = $out[0];
	}

	#Generate output files
	open( NORM, ">$output.by.mOTU.norm" )   or die "Can't open output file\n";
	open( RAW,  ">$output.by.mOTU.raw" )    or die "Can't open output file\n";
	open( SCA,  ">$output.by.mOTU.scaled" ) or die "Can't open output file\n";
	print RAW "$sample\n";
	print NORM "$sample\n";
	print SCA "$sample\n";
	my $weighted_sum = 0;
	my $total_sum    = 0;
	my $weighted_avg = 0;
	foreach my $k ( sort keys %GeneID ) {
		my $coverage_raw  = $GeneID{$k};
		my $coverage_norm = $GeneID{$k} / $GeneIDLength{$k};
		$weighted_sum += $GeneID{$k};
		$total_sum += $coverage_norm;
		print RAW "$coverage_raw\n";
		print NORM "$coverage_norm\n";
		$coverage_raw  = 0;
		$coverage_norm = 0;
	}
	$weighted_avg = $weighted_sum / $total_sum;
	foreach my $k ( sort keys %GeneID ) {
		my $coverage_norm = $GeneID{$k} / $GeneIDLength{$k};
		my $coverage_sca  = $coverage_norm * $weighted_avg;
		print SCA "$coverage_sca\n";
	}
	print "done.\n";
	close NORM;
	close RAW;
	close SCA;
}

#Mode=OG
if ( $mode eq "OG" ) {

	#Load gene annotations
	print "Loading gene annotations...";
	open( MAP, "<$genemap" ) or die "MISSING $genemap";
	while (<MAP>) {
		my @line = split "\t", $_;
		chomp @line;
		$GeneID{ $line[0] }       = 0;
		$GeneIDLength{ $line[0] } = $line[1];
		$OGID{ $line[0] }         = $line[2];
	}
	close MAP;
	print "done\n";

	#Generate annotation file if it does not exist yet
	print "Generating rownames file if it does not exist...";
	if ( -e "$genemap.rownames" ) {
		print "$genemap.rownames already exists.\n";
	}
	else {
		open( ROW, ">$genemap.rownames" ) or die "MISSING $genemap.rownames\n";
		my $header = "GeneID\tOGID\n";
		print ROW $header;
		foreach my $k ( sort keys %GeneID ) {
			print ROW "$k\t$OGID{$k}\n";
		}
		print "done.\n";
	}
	close ROW;

	print "Parsing count data...";
	open( COUNT, "<$count" ) or die "MISSING $count\n";
	while (<COUNT>) {
		my @line = split "\t", $_;
		chomp @line;
		if ( $_ =~ /^-1/ ) {
			next;
		}
		elsif ( $_ !~ /^ref_id/ ) {
			$GeneID{ $line[0] } += $line[1];
		}
	}
	close COUNT;
	print "done.\n";

	print "Generating output files...";
	my @out = split /\./, $count;
	pop @out;
	my $output = join( "\.", @out );
	my $sample;
	if ( $output =~ /\// ) {
		my @a    = split /\//, $output;
		my $file = $a[-1];
		my @b    = split /\.filtered/, $file;
		$sample = $b[0];
	}
	else {
		$sample = $out[0];
	}

	#Generate output files
	open( NORM, ">$output.by.OG.norm" ) or die "Can't open output file\n";
	open( RAW,  ">$output.by.OG.raw" )  or die "Can't open output file\n";

	#    open (SCA, ">$output.by.OG.scaled") or die "Can't open output file\n";
	print RAW "$sample\n";
	print NORM "$sample\n";

	#    print SCA "$sample\n";
	#    my $weighted_sum = 0;
	#    my $total_sum = 0;
	#    my $weighted_avg = 0;
	foreach my $k ( sort keys %GeneID ) {
		my $coverage_raw  = $GeneID{$k};
		my $coverage_norm = $GeneID{$k} / $GeneIDLength{$k};

		#		$weighted_sum += $GeneID{$k};
		#		$total_sum += $coverage_norm;
		print RAW "$coverage_raw\n";
		print NORM "$coverage_norm\n";
		$coverage_raw  = 0;
		$coverage_norm = 0;
	}

	#    $weighted_avg = $weighted_sum / $total_sum;
	#    foreach my $k (sort keys %GeneID){
	#		my $coverage_norm = $GeneID{$k} / $GeneIDLength{$k};
	#		my $coverage_sca = $coverage_norm * $weighted_avg;
	#		print SCA "$coverage_sca\n";
	#   }
	print "done.\nPART 1/2 complete, data still has to be summarized by column.\n";
	close NORM;
	close RAW;

	#    close SCA;
}

exit 0;
