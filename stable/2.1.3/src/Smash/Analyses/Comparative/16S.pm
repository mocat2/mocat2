#! /usr/bin/env perl

package Smash::Analyses::Comparative::16S;
use strict;
use warnings;
use Smash::Core;
use Smash::Global qw(:all);
use Smash::Utils::Taxonomy qw(:all);
use Smash::Utils::MatrixIO qw(read_R_matrix write_R_matrix get_column_labels);

our @ISA = qw(Smash::Analyses::Comparative);

my $READ_PHYLO_MAP;

sub init {
	my $this = shift;
	$this->SUPER::init();
	$this->{TYPE} = "16S";
	Smash::Utils::Taxonomy::init("bergey");
	$this->{READ_LENGTH_THRESHOLD} = 250 unless $this->{READ_LENGTH_THRESHOLD};
	$this->{CONFIDENCE_THRESHOLD}  = 0.6 unless $this->{CONFIDENCE_THRESHOLD};
	if ($this->experiment eq "16S") {
		$this->{READ_LENGTH_THRESHOLD} = 250;
		$this->{CONFIDENCE_THRESHOLD}  = 0.5;
	}
	return $this;
}

sub attach_sample {
	my $this       = shift;
	my $metagenome = shift || die "attach_sample needs metagenome";
	my $label      = shift || die "attach_sample needs label";
	my $smash      = $this->smash;

	$smash->finish() if $smash;
	$smash = new Smash::Core(METAGENOME => $metagenome);
	$smash->init();
	$this->{SMASH} = $smash;

	$this->{LABEL} = $label;
	push(@{$this->{LABELS}}, $label);
}

sub get_data_count {
	my $this  = shift;
	my $metagenome = shift;
	my $dbh   = $this->smash->get_db_handle();
	my $count;
	{
		my $sth   = $dbh->prepare_cached("SELECT COUNT(DISTINCT template_id) FROM readinfo r INNER JOIN library l USING (library_id) INNER JOIN sample s USING (sample_id) WHERE metagenome_id=? AND r.length >= ?");
		$sth->execute($metagenome, 100);
		($count) = $sth->fetchrow_array();
		$sth->fetchrow_array();
	}
	$this->smash->close_db_handle();
	return $count;
}

########################
# Analyze one sample
########################

sub analyze_one_sample {

	my $this = shift;
	my $label = $this->label;
	my $smash = $this->smash;
	my $feature_count = $this->feature_count;
	my $metagenome = $smash->metagenome;
	my $class_file = $smash->analyses_dir($metagenome)."/$metagenome.16S.@{[$this->refdb]}.class";

	# local alignment files, not stored in SMASH repository

	if ($this->is_local) {
		$class_file   = "$metagenome.16S.@{[$this->refdb]}.class";
	}

	my $required_level = $Rank2Num{$this->level};

	my %Template;
	my %ReadLength;

	# if this is not local, check for paired end data from the database

	if (!$this->is_local) {
		my $dbh = $smash->get_db_handle;
		{
			my $sth = $dbh->prepare("SELECT read_id, template_id, length FROM readinfo WHERE read_id LIKE ?");
			$sth->execute("$metagenome.%");
			my ($id, $template, $length);
			$sth->bind_columns(\$id, \$template, \$length);
			while ($sth->fetch()) {
				$Template{$id} = $template;
				$ReadLength{$id} = $length;
			}
		}
		$smash->close_db_handle();
	}

	my $db = new File::Temp();
	my $dbname = $db->filename;
	
	my $dbh = DBI->connect("DBI:SQLite:dbname=$dbname.sqlite", "", "", {AutoCommit => 0, RaiseError => 1}) || die "Couldn't connect to database: " .DBI->errstr;

	$dbh->do("CREATE TABLE templates(template_id VARCHAR(255) NOT NULL, taxonomy_id INTEGER NOT NULL)");
	$dbh->do("CREATE TABLE templates_filtered(template_id VARCHAR(255) NOT NULL, taxonomy_id INTEGER NOT NULL)");
	$dbh->do("CREATE TABLE template_occurrences(template_id VARCHAR(255) PRIMARY KEY, occurrence INTEGER NOT NULL)");
	$dbh->do("CREATE TABLE template_map(template_id VARCHAR(255) NOT NULL, taxonomy_id INT NOT NULL, count FLOAT(6,4))");
	#$dbh->do("CREATE TABLE taxonomy_map(taxonomy_id INT NOT NULL, count FLOAT(12,4))");

	{ # to keep $sth within limited scope
		my $sth = $dbh->prepare("INSERT INTO templates(template_id, taxonomy_id) VALUES(?, ?)");
		open(CLASS, "<$class_file") || die "Cannot open $class_file: $!";
		LINE:while (<CLASS>) {
			chomp();

			my @pairs  = split(/\t/);
			my $read   = shift(@pairs);
			shift(@pairs); # dont need strand now
			next LINE if ( $ReadLength{$read} && $ReadLength{$read} < $this->{READ_LENGTH_THRESHOLD});

			my $phylum;
			my $tax_id;
			my $prev_name;
			my $saved_node = undef;
			RANK: while (scalar(@pairs) > 1) {
				my $name = shift(@pairs); 
				my $rank = shift(@pairs); 
				my $prob = shift(@pairs);

				# If the probability is too low, we quit!

				if ($prob < $this->{CONFIDENCE_THRESHOLD}) { last RANK; }

				# Fix names: complete Incertae Sedis, and remove quotes

				$name   =~ s/^['"]//;
				$name   =~ s/['"]$//;
				if ($name =~ /^Incertae/i) {
					$name = "$prev_name $name";
				}
				$tax_id = $BergeyTree->get_id_for_ranked_name($rank, $name);

				my $level = $Rank2Num{$rank};

				# Save the tax_id if we reach the level we want
				# If we passed it without seeing what we want, try
				# going back up until we pass it the other way

				if ($level == $required_level) { # right there!
					$saved_node = $BergeyTree->nodes->{$tax_id};
					last RANK;
				} elsif ($level > $required_level && !defined($saved_node)) { # passed it without saving!
					my $node = $BergeyTree->nodes->{$tax_id};
					PARENT: while (1) {
						my $parent       = $node->parentlink;
						my $parent_level = $Rank2Num{$parent->rank};
						if ($parent_level < $required_level) {
							$saved_node = $node;
							last RANK;
						}
						$node = $parent;
					}
				}
				$prev_name = $name;
			}

			# At this point, we have $tax_id = the final prediction, $saved_tax_id = prediction closest to the required level

			if (!defined($saved_node)) { # didnt reach the required level
				$saved_node = $BergeyTree->nodes->{$tax_id};
				if (!$saved_node) {
					print "$tax_id\n";
				}
			}

			my $rank = $saved_node->rank;
			my $name = $saved_node->name;

			# $saved_tax_id now has the Bergey tax id for the right candidate at approximately $required_level

			# Want at least $opt_level, others dont matter (buffer 4)
			# For genus, this is suborder, which already includes family!

			if ($Rank2Num{$rank} >= $required_level-4) { 
				$sth->execute($Template{$read}, $tax_id);
			}

			# Note the final mapping down

			if (1 == 0) { # disabled for now
				my $node = $BergeyTree->nodes->{$tax_id};
				$rank = $node->rank;
				if ($Rank2Num{$rank} >= $Rank2Num{"domain"}) {
					print $READ_PHYLO_MAP "$read\t$tax_id\n";
				}
			}
		}
		close(CLASS);
		#$sth->finish();
	}

	$dbh->do("INSERT INTO templates_filtered(template_id, taxonomy_id) SELECT DISTINCT template_id, taxonomy_id FROM templates");
	$dbh->do("INSERT INTO template_occurrences(template_id, occurrence) SELECT template_id, COUNT(*) AS occurrence FROM templates GROUP BY template_id");
	$dbh->do("INSERT INTO template_map(template_id, taxonomy_id, count) SELECT template_id, taxonomy_id, 1.0/occurrence FROM templates INNER JOIN template_occurrences USING (template_id)");
	#$dbh->do("CREATE TABLE taxonomy_map AS SELECT taxonomy_id, SUM(count) AS count FROM template_map GROUP BY taxonomy_id");
	$dbh->commit();

	my $assigned = 0;
	{
		my ($tax, $count);
		my $sth = $dbh->prepare("SELECT taxonomy_id, SUM(count) FROM template_map GROUP BY taxonomy_id");
		$sth->execute();
		$sth->bind_columns(\$tax, \$count);
		while ($sth->fetch()) {
			$feature_count->{$tax}->{$label} = $count;
			$assigned += $count;
		}
	}
	#$dbh->do("DROP TABLE taxonomy_map");

	$dbh->do("DROP TABLE template_map");
	$dbh->do("DROP TABLE template_occurrences");
	$dbh->do("DROP TABLE templates_filtered");
	$dbh->do("DROP TABLE templates");
	$dbh->commit();
	$dbh->disconnect();
	undef $dbh;
	unlink "$dbname.sqlite";

	return $assigned;
}

sub get_attribute {
	my ($string, $attribute) = @_;
	my (undef, $value) = $string =~ m#${attribute}=(["'])([^\1]*?)\1#i;
	return $value;
}

sub parse_rrndb_info {
	my $this = shift;

	my ($rdp_classifier_dir) = $this->smash->software_dir("rdp_classifier", "current");

	# Estimate the 16S copy number for the whole tree.

	my %CopyNumberSum;
	my %Count;
	my %Genera;
	my %Species;
	my %Strains;
	my %Species2Genus;
	my %Strain2Species;
	my %SpCount;

	# As of 28.09.2009, rrnDB has copy number info for 979 genomes.
	# As of 17.06.2010, rrnDB has 1068 genomes.

	open(COPYNUMBER, "<$rdp_classifier_dir/rrn_copy_numbers.txt") || die "Cannot open rrn copy number file: $!";
	while (<COPYNUMBER>) {
		if (!m/^#/) {
			my ($genus, $species, $strain, $copies16S) = (split(/\t/))[0..3];
			next if ($copies16S =~ /NA/);

			if ($species eq "(no species)" || $species eq "sp.") {
				$species .= (($SpCount{"$genus $species"} || 0)+1);
			}
			$species = "$genus $species";
			$strain  = "$species $strain";
			$CopyNumberSum{$strain} = $copies16S;

			# Mark sp/ge/str

			$Genera{$genus} = 1;
			$Species{$species} = 1;
			$Strains{$strain} = 1;

			# Map level-up

			$Species2Genus{$species} = $genus;
			$Strain2Species{$strain} = $species;

			# Set counters

			$CopyNumberSum{$species} = 0;
			$CopyNumberSum{$genus} = 0;
			$Count{$species} = 0;
			$Count{$genus} = 0;
		}

	}
	close(COPYNUMBER);

	# Agglomorate at species level

	foreach my $strain (keys %Strains) {
		my $species = $Strain2Species{$strain};
		$CopyNumberSum{$species} += $CopyNumberSum{$strain};
		$Count{$species}++;
	}

	# Agglomorate at genus level

	foreach my $species (keys %Species) {
		my $genus = $Species2Genus{$species};
		$CopyNumberSum{$genus} += ($CopyNumberSum{$species}/$Count{$species});
		$Count{$genus}++;
	}

	# Summarize at genus level

	foreach my $genus (keys %Genera) {
		my $tax_id = $BergeyTree->get_id_for_ranked_name("genus", $genus);
		if ($tax_id) {
			$BergeyTree->nodes->{$tax_id}->{DATA} = $CopyNumberSum{$genus}/$Count{$genus};
		}
	}

	# Calculate average of the whole phylo tree

	$BergeyTree->propagate_data_average_of_children_to_root();
	$BergeyTree->propagate_data_to_leaves();

	# For debugging purposes

	#$BergeyTree->print_tree();
}

sub merge_feature_vectors {
	my $this = shift;
	my $feature_count = shift;
	my $merge_level   = shift;

	my $MergedFeatureCount;
	my $labels = get_column_labels($feature_count);

	my @ids = grep {$_ >= 0} keys %$feature_count; # skip unmapped = -1 here

	TAX_ID:foreach my $tax_id (@ids) {
		my $rdp_tax_id = $tax_id;
		my $rank = get_taxonomic_rank($BergeyTree, $rdp_tax_id, $merge_level);
		foreach my $label (@$labels) {
			$MergedFeatureCount->{$rank}->{$label} += ($feature_count->{$tax_id}{$label}||0);
		}
	}

	# handle unmapped
	foreach my $label (@$labels) {
		$MergedFeatureCount->{-1}->{$label} += ($feature_count->{-1}{$label}||0);
	}

	return $MergedFeatureCount;
}

sub convert_to_abundance {

	my $this = shift;
	my $feature_count = shift;

	my @samples = @{$this->labels};
	my $data_count = $this->data_count;

	# record raw count in case of replicate

	if (defined($this->replicate)) {
		my $file_name = sprintf("%s.feature_vector.rawcount.%s.%s.%s", $this->name, $this->type, $this->refdb, $this->replicate);
		write_R_matrix($file_name, $feature_count);
	}

	# Parse rrnDB info file and store copy numbers at different nodes in tree

	$this->parse_rrndb_info();

	# Normalize 16S counts by copy number, now we have ~number of individuals

	my @features = keys %$feature_count;
	foreach my $label (@samples) {
		foreach my $tx (@features) {
			my $c = ($feature_count->{$tx}{$label} || 0);
			if ($c) {

				# Get the copy number at the lowest common ancestor
				# Unassigned or unknown 16S gets copy number of root!

				my $copy_number = $BergeyTree->rootlink->data;
				if ($tx != -1) {
					my @lineage = $BergeyTree->get_full_lineage($tx);
					LINEAGE: foreach my $ancestor (@lineage) {
						$copy_number = $BergeyTree->nodes->{$ancestor}->data;
						if ($copy_number) {
#if ($tx == 572511 || $tx == 3017) {
#	printf "Blautia: $copy_number from %s\n", $Tax2Name{$ancestor};
#}
							last LINEAGE;
						}
					}
				}
				$c = $c/$copy_number;
				$feature_count->{$tx}{$label} = $c;
			}
		}
	}

	return $feature_count;
}

1;
