package Smash::Databases::MetaGenomeDB::AssemblyLoader;

use strict;
use warnings;

our @ISA = qw(Smash::Databases::MetaGenomeDB::Loader);

use File::Path;
use Smash::Utils::GFFlite;

############################################
# Constructor
############################################

my $PROGRESS;

sub assembler {shift->{ASSEMBLER}}

sub init {
	my $this = shift;
	$this->parse_config();
	my $assembly = $this->assembly;
	my ($collection, $metagenome) = $this->parse_assembly_id($assembly);
	$this->{COLLECTION} = $collection;
	$this->{METAGENOME} = $metagenome;
	$this->SUPER::init();

	# Database stuff

	my $dbh = $this->get_smashdb_handle;
	{
		my $sth = $dbh->prepare("SELECT name FROM assembly INNER JOIN program ON assembler=program_id WHERE external_id=?");
		$sth->execute($assembly);
		my ($program_name) = $sth->fetchrow_array();
		$this->{ASSEMBLER} = $program_name;
	}
	$this->close_smashdb_handle();

	die "Assembly $assembly not found in Smash!\n" unless $this->assembler;

	$PROGRESS = \*STDERR;
	select($PROGRESS); $| = 1; select(STDOUT);
}

sub finish {
	my $this = shift;
	$this->SUPER::finish();
}

sub unload_db {
	my $this = shift;

	# Locations etc

	my $assembly    = $this->assembly;
	my $assembly_id = $this->get_id_by_name("assembly", $assembly);

	# Remove from MC* database
	print "Removing entries from scaffold, scaffold2contig, contig and contig2read\n";

	my $dbh = $this->get_db_handle();
	{
		####
		# The first two deletes are skipped, since I added an "ON DELETE CASCADE" constraing on scaffold2contig and contig2read.
		# UPDATE on 05.11.2010. SQLite does not respect ON DELETE CASCADE through DBD::SQLite. It only works on command line sqlite3.
		# So I re-enabled them and implemented the deletion of contig2read and scaffold2contig rows the old classic way without a JOIN.
		####

		my $sths = $dbh->prepare('DELETE FROM scaffold2contig WHERE scaffold_id IN (SELECT scaffold_id FROM scaffold WHERE assembly_id=?)');
		$sths->execute($assembly_id);
		my $sthc = $dbh->prepare('DELETE FROM contig2read WHERE contig_id IN (SELECT contig_id FROM contig WHERE assembly_id=?)');
		$sthc->execute($assembly_id);

		my $sth1 = $dbh->prepare('DELETE FROM scaffold WHERE assembly_id=?');
		$sth1->execute($assembly_id);
		#$sth1->fetchrow_array(); # perl DBI shush
		$sth1 = undef;
		my $sth2 = $dbh->prepare('DELETE FROM contig WHERE assembly_id=?');
		$sth2->execute($assembly_id);
		#$sth2->fetchrow_array(); # perl DBI shush
		$sth2 = undef;
	}
	$dbh->commit();
	$this->close_db_handle();
}

sub wipe_out {
	my $this         = shift;
	my $assembly     = $this->assembly;
	my $assembly_dir = $this->assembly_dir($assembly);
	my $assembly_id  = $this->get_id_by_name("assembly", $assembly);

	print "Removing entry from SmashDB.assembly\n";
	my $dbh = $this->get_smashdb_handle;
	{
		my $sth = $dbh->prepare('DELETE FROM assembly WHERE assembly_id=?');
		$sth->execute($assembly_id);
		#$sth->fetchrow_array(); # perl DBI shush
	}
	$dbh->commit();
	$this->close_smashdb_handle();

	print "Removing files from $assembly_dir\n";
	rmtree($assembly_dir);

	my $workspace_dir = sprintf("%s/%s/%s/%s/%s", $this->data_dir, $this->get_smash_conf_value("workspace_dir"), "Assembler", $this->assembler, $assembly);
	print "Removing files from $workspace_dir\n";
	rmtree($workspace_dir);
}

############################################
# Load reads and xml information into db
############################################

sub run {
	use Fcntl qw(:seek);
	my $this = shift;
	
	# Locations etc

	my $assembly     = $this->assembly;
	my $assembly_dir = $this->assembly_dir($assembly);
	my $contigfasta  = "$assembly_dir/$assembly.contigs.fa";
	my $scaffasta    = "$assembly_dir/$assembly.scaffolds.fa";
	my $contig2read  = "$assembly_dir/$assembly.contig2read.gff";
	my $scaf2contig  = "$assembly_dir/$assembly.scaf2contig.gff";

	my $assembly_id  = $this->get_id_by_name("assembly", $assembly);
	if (!$assembly_id) {
		die "Assembly $assembly does not exist in Smash. Please check your assembly run!";
	}

	my %SeqLength; # Track contig/scaf length
	my %SeqCount; # Track contig/scaf featurecount
	my %ContigExternal2Internal; # for contig external to internal mapping
	my %ScaffoldExternal2Internal; # for scaffold external to internal mapping

	# Parse contig fasta

	print $PROGRESS "Contigs:\n";

	print $PROGRESS "\tReading FASTA file ...";

	open(CONTIGS, "<$contigfasta") || die "Cannot open $contigfasta: $!";
	my $fasta = new FAlite(\*CONTIGS);
	while (my $entry = $fasta->nextEntry) {
		my $def = $entry->def;
		$def =~ s/^>//;
		$def =~ s/\s+.*//;
		$SeqLength{$def} = length($entry->seq);
	}
	close(CONTIGS);
	print $PROGRESS " done\n";

	# Parse contig2read GFF

	print $PROGRESS "Contig2Read:\n";

	print $PROGRESS "\tReading GFF file ...";

	open(CONTIG2READ, "<$contig2read") || die "Cannot open $contig2read: $!";
	my $gff   = new Smash::Utils::GFFlite(\*CONTIG2READ);
	while (my $f = $gff->nextFeature) {
		my $name = $f->seqname; 
		$SeqCount{$name}++;
	}
	print $PROGRESS " done\n";

	# Process contig2read maps

	print $PROGRESS "\tLoading into database ...";
	my $dbh = $this->get_db_handle();
	{
		my $contig_sth = $dbh->prepare('INSERT INTO contig(external_id, read_count, length, assembly_id) VALUES(?, ?, ?, ?);');
		foreach my $external_id (keys %SeqLength) {
			$contig_sth->execute($external_id, $SeqCount{$external_id}, $SeqLength{$external_id}, $assembly_id); # We use score field for length
			$ContigExternal2Internal{$external_id} = $this->last_db_insert_id;
		}
	}
	print $PROGRESS ".";

	# read the contig2read gff file again. Don't reopen, since it is a waste. Just rewind
	# but rebless a new GFFlite object!

	{
		my $c2r_sth = $dbh->prepare('INSERT INTO contig2read(contig_id, read_id, start, end, strand) VALUES(?, ?, ?, ?, ?);');
		seek(CONTIG2READ, 0, SEEK_SET);
		$gff   = new Smash::Utils::GFFlite(\*CONTIG2READ);
		while (my $feature = $gff->nextFeature) {
			my $read_id = $feature->get_attribute("read");
			$c2r_sth->execute($ContigExternal2Internal{$feature->seqname}, $read_id, $feature->location->[0]->start, $feature->location->[0]->end, $feature->strand);
		}
		close(CONTIG2READ);
	}
	print $PROGRESS " done\n";
	$dbh->commit();
	$this->close_db_handle();

	%SeqLength = ();
	%SeqCount = ();

	##############
	# Arachne does not generate scaffold fasta files. Therefore, to get the lengths of the
	# scaffold sequences, I use the <score> field of GFF. It is used only when length of
	# a scaffold is not set. Loading external assemblies will usually make a fake scaffold
	# and write it to $scaffasta, so that would be ok.
	##############

	if (-f $scaffasta) {
		# Parse scaffold fasta

		print $PROGRESS "Scaffolds:\n";

		print $PROGRESS "\tReading FASTA file ...";

		open(SCAFFOLDS, "<$scaffasta") || die "Cannot open $scaffasta: $!";
		my $fasta = new FAlite(\*SCAFFOLDS);
		while (my $entry = $fasta->nextEntry) {
			my $def = $entry->def;
			$def =~ s/^>//;
			$def =~ s/\s+.*//;
			$SeqLength{$def} = length($entry->seq);
		}
		close(SCAFFOLDS);
		print $PROGRESS " done\n";
	}

	# Parse scaf2contig GFF

	print $PROGRESS "Scaf2Contig:\n";

	print $PROGRESS "\tReading GFF file ...";

	open(SCAF2CONTIG, "<$scaf2contig") || die "Cannot open $scaf2contig: $!";
	$gff   = new Smash::Utils::GFFlite(\*SCAF2CONTIG);
	while (my $f = $gff->nextFeature) {
		my $name = $f->seqname; 
		$SeqCount{$name}++;
		$SeqLength{$name} = $f->score unless $SeqLength{$name};
	}
	print $PROGRESS " done\n";

	# Process scaf2contig maps

	print $PROGRESS "\tLoading into database ...";
	$dbh = $this->get_db_handle();
	{
		my $scaffold_sth = $dbh->prepare('INSERT INTO scaffold(external_id, contig_count, length, assembly_id) VALUES(?, ?, ?, ?);');
		foreach my $external_id (keys %SeqLength) {
			$scaffold_sth->execute($external_id, $SeqCount{$external_id}, $SeqLength{$external_id}, $assembly_id); # We use score field for length
			$ScaffoldExternal2Internal{$external_id} = $this->last_db_insert_id;
		}
	}
	print $PROGRESS ".";

	{
		my $s2c_sth = $dbh->prepare('INSERT INTO scaffold2contig(scaffold_id, contig_id, start, end, strand) VALUES(?, ?, ?, ?, ?);');
		seek(SCAF2CONTIG, 0, SEEK_SET);
		$gff   = new Smash::Utils::GFFlite(\*SCAF2CONTIG);
		while (my $feature = $gff->nextFeature) {
			my $contig_id = $ContigExternal2Internal{$feature->get_attribute("contig")};
			$s2c_sth->execute($ScaffoldExternal2Internal{$feature->seqname}, $contig_id, $feature->start, $feature->end, $feature->strand);
		}
	}
	$dbh->commit();
	$this->close_db_handle();
	print $PROGRESS " done\n";
}

1;
