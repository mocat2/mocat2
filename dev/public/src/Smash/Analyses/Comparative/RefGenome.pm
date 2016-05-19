#! /usr/bin/env perl

package Smash::Analyses::Comparative::RefGenome;
use strict;
use warnings;
use Smash::Core qw(:all);
use Smash::Global qw(:all);
use Smash::Utils::Taxonomy qw(:all);
use Smash::Utils::MatrixIO qw(read_R_matrix write_R_matrix zero_fill_matrix get_column_labels);

our @ISA = qw(Smash::Analyses::Comparative);

my $MIN_ALIGN_LEN = 100;

my $READ_PHYLO_MAP;
my $AVERAGE_GENOME_SIZE;
my %Species2Tax = ();
my %Sequence2Tax = ();
my %GenomeSize   = ();

sub init {
	my $this = shift;
	$this->SUPER::init();
	$this->{TYPE} = "refgenome";
	Smash::Utils::Taxonomy::init("BERGEY");
	Smash::Utils::Taxonomy::init("NCBI.DIV");

	# if localdb is used, then dont connect to database

# There is a problem in estimating genome size. Some times some genomes
# are sequenced multiple times, and they even use the same definition field.
# But this procedure used to just add all the sizes together for each genome,
# which is plain wrong. However, since some have multiple chromosomes, you
# should not use average of all either! The compromise is to only use one
# value per definition, so if the same sequence was sequenced twice then only
# one is counted. To make it reproducible, the first query is sorted by 
# sequence_id.
# CAVEAT: If there are two WGS genomes with the same tax_id and using different
#         assemblies, and hence different contig/scaffold names, this cannot be
#         solved. This must be handled manually!

	my $smash = new Smash::Core();
	$smash->init();
	my $dbh   = $smash->get_refgenomedb_handle();

	# init seq2tax
	# init genome size

	my %Size;
	my $seq_sth = $dbh->prepare_cached("SELECT taxonomy_id, sequence_id, definition, length FROM sequence ORDER BY sequence_id");
	$seq_sth->execute();
	while (my ($tax_id, $seq_id, $definition, $length) = $seq_sth->fetchrow_array()) {
		my $genome_id = $tax_id;
		$Sequence2Tax{$seq_id} = $genome_id;
		$Size{$genome_id}{$seq_id} = $length;
	}
	foreach my $genome_id (keys %Size) {
		$GenomeSize{$genome_id} = 0;
		foreach my $seq_id (keys %{$Size{$genome_id}}) {
			$GenomeSize{$genome_id} += $Size{$genome_id}{$seq_id};
		}
	}

	$smash->close_refgenomedb_handle();
	$smash->finish();

	# estimate average genome size
	# used to be 3539003

	my $total = 0;
	my @sizes = values %GenomeSize;
	$total += $_ foreach @sizes;
	$AVERAGE_GENOME_SIZE = int($total/@sizes+0.5);
	print STDERR "Avg. genome size = $AVERAGE_GENOME_SIZE\n";
}

sub attach_sample {
	my $this       = shift;
	my $metagenome = shift || die "attach_sample needs metagenome";
	my $label      = shift || die "attach_sample needs label";
	my $smash      = $this->smash;

	my %options;

	# if localdb is used, then dont connect to database

	if (!$this->is_local) {
		%options = (METAGENOME => $metagenome);
	}
	$smash->finish() if $smash;
	$smash = new Smash::Core(%options);
	$smash->init();
	$this->{SMASH} = $smash;

	# for local dbfile, we skipped setting METAGENOME
	# if we didnt, then SMASH will try to connect to the database.
	# now we fix this by explicitly setting METAGENOME so that the next steps work

	$smash->{METAGENOME} = $metagenome unless $smash->metagenome; 

	$this->{LABEL} = $label;
	push(@{$this->{LABELS}}, $label);
}

sub get_data_count {
	my $this  = shift;
	my $metagenome = $this->smash->metagenome;

	my $count;
	my $dbh = $this->smash->get_db_handle();
	{
		my $sth   = $dbh->prepare_cached("SELECT COUNT(DISTINCT template_id) FROM readinfo r INNER JOIN library l USING (library_id) INNER JOIN sample s USING (sample_id) WHERE metagenome_id=? AND r.length >= ?");
		$sth->execute($metagenome, $MIN_ALIGN_LEN);
		($count) = $sth->fetchrow_array();
		$sth->fetchrow_array();
	}
	$this->smash->close_db_handle();
	return $count;
}

# convert to abundance by normalizing by genome size

sub convert_to_abundance {

	my $this = shift;
	my $feature_count = shift;

	my $samples = get_column_labels($feature_count);

	zero_fill_matrix($feature_count);

	# record raw count in case of replicate

	if (defined($this->replicate)) {
		my $file_name = sprintf("%s.feature_vector.rawcount.%s.%s.%s", $this->name, $this->type, $this->refdb, $this->replicate);
		write_R_matrix($file_name, $feature_count);
	}

	# convert mapped counts to coverage, which is a proxy for number of individuals

	foreach my $genome_id (grep {$_ >= 0} keys %$feature_count) {
		my $genome_size = $GenomeSize{$genome_id} || die "Genome size unknown for $genome_id!";
		foreach my $label (@$samples) {
			my $c = $feature_count->{$genome_id}{$label} || 0;
			$feature_count->{$genome_id}{$label} = $c/$genome_size;
		}
	}

	# unmapped counts to coverage using avg genome length in refgenome set

	map {$feature_count->{-1}{$_} /= $AVERAGE_GENOME_SIZE if $feature_count->{-1}{$_};} @$samples;

	return $feature_count;
}

sub analyze_one_sample {
	use DBD::SQLite;
	use File::Temp;
	my $this              = shift;
	my $label             = $this->label;
	my $identity          = $this->identity;
	my $feature_count     = $this->feature_count;
	my $refdb             = $this->refdb;
	my $smash             = $this->smash;
	my $metagenome        = $smash->metagenome;
	my $analyses_dir      = $smash->analyses_dir($metagenome);
	my $blast_file        = "$analyses_dir/$metagenome.$refdb.blastn";

	my $multilevel_threshold;
	if ($this->multilevel) {
	}

	# local alignment files, not stored in SMASH repository

	if ($this->is_local) {
		$blast_file   = "$metagenome.$refdb.blastn";
	}

	my $BLAST = safe_open_file($blast_file);

	my %Template;

	# if this is not local, check for paired end data from the database

	if (!$this->is_local) {
		my ($id, $template);
		my $dbh = $smash->get_db_handle();
		{
			my $sth = $dbh->prepare("SELECT read_id, template_id FROM readinfo WHERE read_id LIKE ?");
			$sth->execute("$metagenome.%");
			$sth->bind_columns(\$id, \$template);
			while ($sth->fetch()) {
				$Template{$id} = $template;
			}
		}
		$smash->close_db_handle();
	}

	my $db = new File::Temp();
	my $dbname = $db->filename;
	
	my $dbh = DBI->connect("DBI:SQLite:dbname=$dbname.sqlite", "", "", {AutoCommit => 0, RaiseError => 1}) || die "Couldn't connect to database: " .DBI->errstr;

	$dbh->do("CREATE TABLE templates(template_id VARCHAR(255) NOT NULL, taxonomy_id INTEGER NOT NULL)");
	$dbh->do("CREATE TABLE read_ref(read_id VARCHAR(255) NOT NULL, taxonomy_id INTEGER NOT NULL)");

	# Read the refgenome blasts

	{
		my $sth = $dbh->prepare("INSERT INTO templates(template_id, taxonomy_id) VALUES(?, ?)");
		my $sth2 = $dbh->prepare("INSERT INTO read_ref(read_id, taxonomy_id) VALUES(?, ?)");

		# WU-BLAST or NCBI-BLAST?

		my $line = <$BLAST>;
		chomp($line);
		my @dummy  = split(/\s+/, $line);
		my $flavor = "WU";
		   $flavor = "NCBI" if (@dummy == 12);
		seek $BLAST, 0, 0;

		my ($pi_field, $qb_field, $qe_field, $bits_field);
		if ($flavor eq "WU") {
			$pi_field = 10;
			$qb_field = 17;
			$qe_field = 18;
			$bits_field = 4;
		} else {
			$pi_field = 2;
			$qb_field = 6;
			$qe_field = 7;
			$bits_field = 11;
		}


		# Track best scores for each read

		my %BestScore;
		HSP:while (<$BLAST>) {
			chomp;
			my @words = split(/\s+/);
			my ($read, $bits) = @words[0, $bits_field];
			$BestScore{$read} = $bits if not $BestScore{$read};
			$BestScore{$read} = $bits if $bits > $BestScore{$read};
		}

		# Re-parse the BLAST output and keep only the best hits

		seek $BLAST, 0, 0;
		HSP:while (<$BLAST>) {
			chomp;
			my @words = split(/\s+/);
			my ($read, $ref_sequence, $pi, $qb, $qe, $bits) = @words[0,1,$pi_field, $qb_field, $qe_field, $bits_field];

			next HSP unless $pi > $identity;
			next HSP if abs($qe - $qb) < $MIN_ALIGN_LEN;
			next HSP if $bits < $BestScore{$read};

			my ($tax_id) = $Sequence2Tax{$ref_sequence};
			defined($tax_id) || die "No match for $ref_sequence\n";

			$sth->execute($Template{$read} || $read, $tax_id); # for non paired-end, template is the same as read
			$sth2->execute($read, $tax_id);
		}
		$BLAST->close();
	}

	# create index to make it faster

	$dbh->do("CREATE INDEX idx_t1 ON templates(template_id, taxonomy_id)");

	# we have a redundant list of (template, tax) now.

	# this is skippped now.
	# uncommenting this will treat single end and paired end mappings to a taxon the same. (fwd+rev, fwd only, rev only all get the same weight).
	#$dbh->do("CREATE TABLE templates_filtered(template_id VARCHAR(255), taxonomy_id INTEGER, CONSTRAINT pk_tf1 PRIMARY KEY(template_id, taxonomy_id))");
	#$dbh->do("INSERT INTO templates_filtered(template_id, taxonomy_id) SELECT DISTINCT template_id, taxonomy_id FROM templates");

	# Let's get the number of mappings per template, since we need that for sharing across taxa later

	$dbh->do("CREATE TABLE template_mapping_count(template_id VARCHAR(255) PRIMARY KEY, map_count INTEGER NOT NULL)");
	$dbh->do("INSERT INTO template_mapping_count(template_id, map_count) SELECT template_id, COUNT(*) AS map_count FROM templates GROUP BY template_id");

	# Let us now get a (template, tax, fraction) from them all
	# each such k-tuple give the fraction of template assigned to that tax. sum(fraction) over tax for each template must be 1.

	$dbh->do("CREATE TABLE template_tax_map(template_id VARCHAR(255) NOT NULL, taxonomy_id INT NOT NULL, count FLOAT(6,4), CONSTRAINT pk_tm1 PRIMARY KEY(template_id, taxonomy_id))");
	$dbh->do("INSERT INTO template_tax_map(template_id, taxonomy_id, count) SELECT t.template_id, taxonomy_id, SUM(1.0/map_count) FROM templates t INNER JOIN template_mapping_count tmc ON t.template_id=tmc.template_id GROUP BY t.template_id, taxonomy_id");

	$dbh->do("CREATE INDEX idx_r1 ON read_ref(read_id, taxonomy_id)");
	$dbh->do("CREATE INDEX idx_r2 ON read_ref(taxonomy_id)");
	$dbh->do("CREATE TABLE read_ref_map(read_id VARCHAR(255), taxonomy_id INTEGER, confidence FLOAT, CONSTRAINT pk_tf1 PRIMARY KEY(read_id, taxonomy_id))");
	{
		my %TaxMaps;
		my $sth = $dbh->prepare("SELECT DISTINCT taxonomy_id FROM read_ref");
		$sth->execute();
		while (my ($tax_id) = $sth->fetchrow_array()) {
			my $new_id = $tax_id;
			$TaxMaps{$tax_id} = $new_id;
		}

		my %ReadRef;
		my $sth2 = $dbh->prepare("SELECT read_id, taxonomy_id FROM read_ref");
		$sth2->execute();
		while (my ($read, $tax_id) = $sth2->fetchrow_array()) {
			$ReadRef{$read}{$TaxMaps{$tax_id}} = 1;
		}

		my $read_refgenome_map = sprintf("%s.%s.read_refgenome_map.txt", $this->name, $metagenome);
		open(READ_REFGENOME_MAP, ">$read_refgenome_map") || die "Cannot open $read_refgenome_map: $!";
		foreach my $read (keys %ReadRef) {
			my @taxa = keys %{$ReadRef{$read}};
			foreach my $tax_id (@taxa) {
				printf READ_REFGENOME_MAP "%s\t%d\t%f\n", $read, $tax_id, 1/@taxa;
			}
		}
		close(READ_REFGENOME_MAP);
	}

	$dbh->commit();

	my $assigned = 0;
	{
		my ($tax, $count);
		my $sth = $dbh->prepare("SELECT taxonomy_id, SUM(count) FROM template_tax_map GROUP BY taxonomy_id");
		$sth->execute();
		$sth->bind_columns(\$tax, \$count);
		while ($sth->fetch()) {
			$feature_count->{$tax}->{$label} = $count;
			$assigned += $count;
		}
	}
	#$sth->finish();
	#$dbh->do("DROP TABLE taxonomy_map");

	$dbh->do("DROP TABLE template_tax_map");
	$dbh->do("DROP TABLE template_mapping_count");
	$dbh->do("DROP TABLE templates");
	$dbh->commit();
	$dbh->disconnect();
	undef $dbh;
	unlink "$dbname.sqlite";
	#print STDERR "multihits in $label: $multihits\n";

	return $assigned;
}

1;
