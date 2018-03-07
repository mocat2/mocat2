#! /usr/bin/env perl

package Smash::Analyses::Annotator;
use strict;
use warnings;
use Smash::Core;
use Smash::Utils::BlastParser;

our @ISA = qw(Smash::Core);

sub init {
	my $this = shift;
	$this->SUPER::init();
	my $name = $this->name;
	my $installed_version = $this->get_conf_value("Current Version", lc($name));
	if ($this->{VERSION} ne $installed_version) {
		printf STDERR "ERROR: Incompatible $name version.\n";
		printf STDERR "\tRequested version: %s\n", $this->{VERSION};
		printf STDERR "\tInstalled version: %s\n", $installed_version;
		die "Will quit now!\n";
	}
}

sub blast_flavor {
	shift->{FLAVOR};
}

1;

package Smash::Analyses::Annotator::Functional;
use strict;
use warnings;
use Smash::Core;

our @ISA = qw(Smash::Analyses::Annotator);

sub name {shift->{NAME}}
sub type {shift->{TYPE}}

sub new {
	my $class = shift;
	my %param = @_;
	my @available = qw(KEGG eggNOG);
	my ($name) = grep {lc($param{REF_DB}) eq lc($_)} @available;
	if (!$name) {
		warn "Unknown functional database: ".$param{REF_DB}."\n";
		die  "Known: eggnog, kegg\n";
	}
	$class .= "::$name";
	my $this  = bless {%param}, $class;
	select(STDERR); $| = 1; select(STDOUT);
	$this->{PROGRESS} = \*STDERR;
	$this->{TYPE} = "functional";
	$this->{NAME} = $name;
	return $this;
}

1;

package Smash::Analyses::Annotator::Functional::KEGG;
use strict;
use warnings;
use Smash::Core;
#use Smash::Utils::MatrixIO qw(write_R_matrix zero_fill_matrix);

our @ISA = qw(Smash::Analyses::Annotator::Functional Exporter);
#our @EXPORT_OK = qw($KeggKO $KeggModule $KeggPathway);
#our %EXPORT_TAGS = ('all' => [@EXPORT_OK]);

sub annotate {
	my $this         = shift;
	my ($kegg_blast, $kegg_map) = @_;

	my $Length     = {};
	my $Protein2KO = {};
	my $KO2Module  = {};
	my $KO2Pathway = {};

	print STDERR "# Getting information from database\n";

	my $dbh = $this->get_refproteindb_handle;

	# Get protein lengths

	{ 
		my $kegg_length_sth = $dbh->prepare("SELECT protein, length FROM kegg_details");
		$kegg_length_sth->execute();
		my ($kegg, $length);
		$kegg_length_sth->bind_columns(\$kegg, \$length);
		while ($kegg_length_sth->fetch()) {
			$Length->{$kegg} = $length;
		}
	}

	# Get protein2ko

	{
		my ($kegg, $ko);
		my $sth = $dbh->prepare("SELECT protein, ko FROM protein2ko");
		$sth->execute();
		$sth->bind_columns(\$kegg, \$ko);
		while ($sth->fetch()) {
			$Protein2KO->{$kegg}->{$ko} = 1;
		}
	}

	# Get ko2module

	{
		my ($ko, $module);
		my $sth = $dbh->prepare("SELECT ko, module FROM ko2module");
		$sth->execute();
		$sth->bind_columns(\$ko, \$module);
		while ($sth->fetch()) {
			$KO2Module->{$ko}->{$module} = 1;
		}
	}

	# Get ko2pathway

	{
		my ($ko, $pathway);
		my $sth = $dbh->prepare("SELECT ko, pathway FROM ko2pathway");
		$sth->execute();
		$sth->bind_columns(\$ko, \$pathway);
		while ($sth->fetch()) {
			$pathway =~ s/^(map|ko)//;
			$KO2Pathway->{$ko}->{$pathway} = 1;
		}
	}

	$this->close_refproteindb_handle();

	print STDERR "# Performing mapping\n";

	# Read the kegg_maps

	open(KEGG, "<$kegg_blast") || die "Cannot open $kegg_blast: $!";
	open(MAP, ">$kegg_map") || die "Cannot open $kegg_map: $!";

	my $parser = new Smash::Utils::BlastParser(FILE => $kegg_blast, FLAVOR => $this->blast_flavor);
	GENE:while (my $report = $parser->nextReport) {
		my $gene = $report->queryName;
		HIT:while (my $sbjct = $report->nextSbjct) {

			my $kegg = $sbjct->sbjctName;
			next HIT unless $Protein2KO->{$kegg};   # If no KO annotation, then go to the next one

			# KO

			my @kos = keys (%{$Protein2KO->{$kegg}});

			# Modules and Pathways

			my $Modules  = {};
			my $Pathways = {};
			foreach my $ko (@kos) {
				map {$Modules->{$_}  = 1} keys %{$KO2Module->{$ko}};
				map {$Pathways->{$_} = 1} keys %{$KO2Pathway->{$ko}};
			}

			my @modules  = keys %$Modules;
			my @pathways = keys %$Pathways;

			my @hsps;
			while (my $group = $sbjct->nextGroup) {
				while (my $hsp = $group->nextHSP) {
					push(@hsps, $hsp);
				}
			}

			# Estimate covered region of KEGG protein

			# This is protein blast, so dont worry about strand. 
			# Assume qb < qe, sb < se.

			my $Covered = [];
			my $covered = 0;
			my $bits;
			my $sim;
			if (@hsps == 1) {
				my ($hsp) = @hsps;
				$covered = $hsp->se-$hsp->sb+1;
				$bits = $hsp->bits;
				$sim  = $hsp->percent;
			} else {
				my @bits;
				my @sim;
				foreach my $hsp (@hsps) {
					map {$Covered->[$_] = 1} ($hsp->sb..$hsp->se);
					push(@bits, $hsp->bits);
					push(@sim, $hsp->percent);
				}
				map {$covered++ if ($Covered->[$_])} (0..$#$Covered);
				$bits = join(",", @bits);
				$sim  = join(",", @sim);
			}
			$covered /= $Length->{$kegg};

			printf MAP "%s\t%s\t%s\t%s\t%.6f\t%s\t%s\t%s\n", $gene, $kegg, $sim, $bits, $covered, join(",", @kos), join(",", @modules), join(",", map {"ko$_"} @pathways);
			last HIT;
		}
	}
=begin OLD
	LINE:while (my ($gene, $kegg, $gstart, $gend, $kstart, $kend) = $this->parse_next_best_hit_line(\*KEGG)) {
		last LINE unless $gene;

		# KO

		my @kos = sort keys (%{$Protein2KO->{$kegg}});

		# Modules and Pathways

		my $Modules  = {};
		my $Pathways = {};
		foreach my $ko (@kos) {
			map {$Modules->{$_}  = 1} keys %{$KO2Module->{$ko}};
			map {$Pathways->{$_} = 1} keys %{$KO2Pathway->{$ko}};
		}

		my @modules  = sort keys %$Modules;
		my @pathways = sort keys %$Pathways;

		print MAP join("\t", $gene, $kegg, $gstart, $gend, $kstart, $kend, join(",", @kos), join(",", @modules), join(",", map {"ko$_"} @pathways));
		print MAP "\n";
	}
=cut

	close(MAP);
	close(KEGG);
}

sub parse_next_best_hit_line {
	my $this = shift;
	my $FH   = shift;
	my $line = $this->get_next_uncommented_line($FH);
	return 0 unless $line;
	my @w = split(/\s+/, $line);
	my ($gene, $kegg, $gstart, $gend, $kstart, $kend) = @w;
	return ($gene, $kegg, $gstart, $gend, $kstart, $kend);
}

1;
