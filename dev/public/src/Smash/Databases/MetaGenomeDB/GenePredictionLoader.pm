package Smash::Databases::MetaGenomeDB::GenePredictionLoader;

use strict;
use warnings;

our @ISA = qw(Smash::Databases::MetaGenomeDB::Loader);

use File::Path;
use Smash::Utils::GFF;
use Smash::Utils::GTF;

my $PROGRESS;

sub label {shift->{LABEL}}

sub init {
	my $this       = shift;
	$this->SUPER::init();

	if (!$this->{LABEL}) {  # fasta file is given, so use it!
		$this->{LABEL} = $this->{GENEPRED};
	}
	$PROGRESS = \*STDERR;
	select($PROGRESS); $| = 1; select(STDOUT);
}

sub unload_db {
	my $this = shift;

	# Locations etc

	my $genepred    = $this->genepred;
	my $genepred_id = $this->get_id_by_name("gene_prediction", $genepred);

	# Remove from MC* database
	print "Removing entries from gene and SmashDB.gene_prediction\n";

	my $dbh = $this->get_db_handle;
	{
		my $sth = $dbh->prepare('DELETE FROM gene WHERE gene_prediction_id=?');
		$sth->execute($genepred_id);
		$sth->finish();
	}
	$dbh->commit();
	$this->close_db_handle();
}

sub wipe_out {
	my $this         = shift;
	my $genepred     = $this->genepred;
	my $genepred_dir = $this->genepred_dir($genepred);
	my $genepred_id  = $this->get_id_by_name("gene_prediction", $genepred);
	my $dbh;

	print "Removing entry from SmashDB.gene_prediction\n";
	$dbh = $this->get_smashdb_handle;
	{
		my $sth = $dbh->prepare('DELETE FROM gene_prediction WHERE gene_prediction_id=?');
		$sth->execute($genepred_id);
		$sth->finish();
	}
	$this->close_smashdb_handle();

	print "Removing files from $genepred_dir\n";
	rmtree($genepred_dir);

}

############################################
# Actual execution of the object
############################################

sub run {
	my $this         = shift;
	my $assembly     = $this->assembly;
	my $genepred     = $this->genepred;
	my $genepred_dir = $this->genepred_dir($genepred);
	my $label        = $this->label;


	my $assembly_id  = $this->get_id_by_name("assembly", $assembly)        || "Assembly $assembly does not exist in Smash!";
	my $genepred_id  = $this->get_id_by_name("gene_prediction", $genepred) || "Gene prediction $genepred does not exist in Smash!";

	my %ContigExternal2Internal; # for contig external to internal mapping

	# Read the GFF file

	my $gff_file;
	if ($label eq $genepred) { # full run
		$gff_file = "$genepred_dir/$label.contig2gene.gff";
	} else {                  # parallelized run, split
		$gff_file = "$genepred_dir/$label/$label.contig2gene.gff";
	}

	####
	# Get contig external2internal maps from DB
	####

	print $PROGRESS "Getting contig information ...";
	my $dbh = $this->get_db_handle;
	{
		my $contig_sth   = $dbh->prepare('SELECT contig_id, external_id FROM contig WHERE assembly_id=? ORDER BY contig_id');
		$contig_sth->execute($assembly_id);
		my ($bind_id, $bind_name);
		$contig_sth->bind_columns( \( $bind_id, $bind_name ) );
		while ($contig_sth->fetchrow_arrayref()) {
			$ContigExternal2Internal{$bind_name} = $bind_id;
		}
	}
	$this->close_db_handle();
	print $PROGRESS " done\n";

	####
	# process genes
	####

	print $PROGRESS "Parsing GFF file ...";
	my $genes = Smash::Utils::GTF::parse_gtf($gff_file);
	print $PROGRESS " done\n";

	print $PROGRESS "Loading genes into database ...";
	{
		$dbh = $this->get_db_handle;
		my $gene_sth     = $dbh->prepare('INSERT INTO gene(external_id, gene_prediction_id, contig_id, length, start, end, strand, start_codon, stop_codon, gc) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?);');
		foreach my $gene (@$genes) {
			my $contig_id   = $ContigExternal2Internal{$gene->seqname};
			my $external_id = $gene->name;
			my $strand      = $gene->strand;
			my $length      = 0;
			my $start_codon = 0;
			my $stop_codon  = 0;
			foreach my $feature (@{$gene->features}) {
				my $type = $feature->type;

				# do we have start or stop?

				# Have to handle new GTF as well as old GFF with start_codon and stop_codon as attributes!

				$start_codon = 1 if ($type eq "start_codon") or Smash::Utils::GFF::text2flag($feature->get_attribute("start_codon"));
				$stop_codon  = 1 if ($type eq "stop_codon")  or Smash::Utils::GFF::text2flag($feature->get_attribute("stop_codon"));

				# start and stop dont count; unknown could be intron etc...

				$length += ($feature->end-$feature->start+1) if $type eq "CDS" or $type eq "rRNA" or $type eq "miscRNA" or $type eq "tRNA";
			}
			$gene_sth->execute($external_id, $genepred_id, $contig_id, $length, $gene->start, $gene->end, $strand, $start_codon, $stop_codon, 0);
		}
	}
	$dbh->commit();
	$this->close_db_handle();
	print $PROGRESS " done\n";
}

1;
