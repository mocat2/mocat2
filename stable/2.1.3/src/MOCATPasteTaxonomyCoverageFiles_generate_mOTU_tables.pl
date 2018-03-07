use warnings;
use strict;
use POSIX;
use Getopt::Long;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

my $wd;
my $annot;
my $motu;
my $prefix;
my $full_output = 0;

my $PREVALENCE_MINIMUM = 2;
GetOptions(
	'wd:s'         => \$wd,
	'map:s'        => \$annot,
	'table:s'      => \$motu,
	'prefix:s'     => \$prefix,
	'prevalence:s' => \$PREVALENCE_MINIMUM,
	'full-output'  => \$full_output,
);

if ( scalar @ARGV > 4 ) {
	$PREVALENCE_MINIMUM = $ARGV[4];
	( $PREVALENCE_MINIMUM > 0 )
	  or die "Prevalence minimum should be greater than zero (got $PREVALENCE_MINIMUM)";
}

sub sum_row {
	my $res = 0;
	foreach my $v (@_) { $res += $v; }
	return $res;
}

sub mean_row_if_valid {
	my $len = scalar @_;
	if ( $len < $PREVALENCE_MINIMUM ) { return 0.; }
	my $m = ( sum_row @_ ) / $len;
	if ( $m < 1. ) { return 0.; }
	return $m;
}

chdir $wd;

open( ANNOT, '<', $annot )
  or die("Cannot open annotation file '$annot': $!");

my %annotation;
while (<ANNOT>) {
	unless (m/^#/) {
		chomp;
		my @tokens  = split /\t/;
		my $motu_id = shift @tokens;
		$annotation{$motu_id} = \@tokens;
	}
}
close(ANNOT);
my @colnames_mOTU_taxo_annot = ( "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Cluster", "annotated.mOTU.clusters", "mOTU.clusters" );

#open( MOTU, "gunzip -c $motu|" )
#  or die("Cannot $motu file: $!");
open( MOTU, "<", $motu ) or die("Cannot $motu file: $!");


my @samples;
my %motu_data;
while (<MOTU>) {
	if (/^#/) { next }
	chomp;
	my @tokens  = split /\t/;
	my $motu_id = shift @tokens;
	unless (@samples) {
		@samples = @tokens;
		next;
	}
	@tokens = map { floor($_); } @tokens;
	if ( sum_row(@tokens) ) {
		$motu_data{$motu_id} = \@tokens;
	}
}

close(MOTU);

foreach my $i ( 8, 9 ) {
	my $level_name = $colnames_mOTU_taxo_annot[$i];
	my %level_count;
	my %activemotus;
	foreach my $motu_id ( keys %annotation ) {
		my $group = $annotation{$motu_id}[$i];
		if ( $group eq "NA" ) { next }
		if ($full_output) {
			$activemotus{$group} = 1;
		}
		if ( not exists $motu_data{$motu_id} ) { next; }
		my @row = @{ $motu_data{$motu_id} };
		foreach my $i ( 0 .. ( scalar @row - 1 ) ) {
			my $sample = $samples[$i];
			$level_count{"$group.$sample"} = [] unless exists $level_count{"$group.$sample"};
			if ( $row[$i] ) {
				push( @{ $level_count{"$group.$sample"} }, $row[$i] );
				$activemotus{$group} = 1;
			}
		}
	}
	my %total;
	my %valid;
	my %printed;
	foreach my $sample (@samples) {
		$total{$sample} = 0;
	}
	open( MOTU_COUNT_MEAN, '>', "$prefix.$level_name.tab" )
	  or die("Cannot open output file: $!");
	print MOTU_COUNT_MEAN "\t";
	print MOTU_COUNT_MEAN join( "\t", @samples );
	print MOTU_COUNT_MEAN "\n";
	foreach my $group ( sort keys %activemotus ) {
		my $maxn = 0;
		foreach my $sample (@samples) {
			my $values_ref = $level_count{"$group.$sample"};
			if ( defined $values_ref ) {
				my @values = @{$values_ref};
				$maxn = scalar @values if ( scalar @values > $maxn );
			}
		}
		if ( $maxn < $PREVALENCE_MINIMUM && !$full_output ) { next; }
		print MOTU_COUNT_MEAN "$group";
		foreach my $sample (@samples) {
			my $values_ref = $level_count{"$group.$sample"};
			if ( defined $values_ref ) {
				my @values = @{$values_ref};
				my $mean   = mean_row_if_valid(@values);
				print MOTU_COUNT_MEAN "\t$mean";
				$total{$sample} += $mean;
				$valid{"$group.$sample"} = $mean;
			}
			else {
				print MOTU_COUNT_MEAN "\t0";
			}
			$printed{$group} = 1;
		}
		print MOTU_COUNT_MEAN "\n";
	}
	close(MOTU_COUNT_MEAN);

	open( MOTU_FRACTION_MEAN, '>', "$prefix.$level_name.fraction.tab" )
	  or die("Cannot open output file: $!");
	print MOTU_FRACTION_MEAN "\t";
	print MOTU_FRACTION_MEAN join( "\t", @samples );
	print MOTU_FRACTION_MEAN "\n";
	foreach my $group ( sort keys %printed ) {
		print MOTU_FRACTION_MEAN "$group";
		foreach my $sample (@samples) {
			my $t = $total{$sample};
			my $v = $valid{"$group.$sample"};
			if ( defined $t && defined $v && $t > 0 ) {
				my $normed = $v / $t;
				print MOTU_FRACTION_MEAN "\t$normed";
			}
			else {
				print MOTU_FRACTION_MEAN "\t0";
			}
		}
		print MOTU_FRACTION_MEAN "\n";
	}
	close(MOTU_FRACTION_MEAN);

}
