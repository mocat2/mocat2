#! /usr/bin/env perl

package Smash::Analyses::Comparative::Functional;
use strict;
use warnings;
use Smash::Core;
use Smash::Utils::MatrixIO qw(write_R_matrix zero_fill_matrix);

our @ISA = qw(Smash::Analyses::Comparative Exporter);
our @EXPORT_OK = qw($KeggKO $KeggModule $KeggPathway);
our %EXPORT_TAGS = ('all' => [@EXPORT_OK]);

sub normalize {shift->{NORMALIZE}}

# Should be implemented by subclasses

sub get_popup_info {}
sub parse_annotations {}

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
	return $this;
}

sub init {
	my $this = shift;
	if ($this->normalize eq "data") {
		warn "functional annotation normally uses hit-normalization and not data-normalization.";
	}
	$this->SUPER::init();
}

sub get_functional_hierarchy {
	my $this  = shift;
	my $level = $this->level;
	return $this->{TREE}->{uc($level)};
}

sub attach_sample {
	my $this     = shift;
	my $genepred = shift || die "attach_sample needs genepred";
	my $label    = shift || die "attach_sample needs label";
	my $smash    = $this->smash;

	$smash->finish() if $smash;
	$smash = new Smash::Core(GENEPRED => $genepred);
	$smash->init();
	$this->{SMASH} = $smash;

	$this->{LABEL} = $label;
	push(@{$this->{LABELS}}, $label);
}

sub get_data_count {
	my $this  = shift;
	my $genepred = shift;
	my $dbh   = $this->smash->get_db_handle();
	my $count;
	{
		my $sth   = $dbh->prepare_cached("SELECT gene_unit_count FROM summary_gene_stats WHERE gene_prediction_id=?");
		$sth->execute($genepred);
		($count) = $sth->fetchrow_array();
		$sth->fetchrow_array();
	}
	$this->smash->close_db_handle();
	return $count;
}

sub analyze_one_sample {
	my $this         = shift;
	my $label        = $this->label;
	my $smash        = $this->smash;
	my $metagenome   = $smash->metagenome;
	my $genepred     = $smash->genepred;
	my $genepred_id  = $smash->get_id_by_name("gene_prediction", $genepred);
	my $genepred_dir = $smash->genepred_dir($genepred);
	my $map_file     = "$genepred_dir/$genepred.@{[$this->refdb]}mapping.txt";

	# hack
	# warn "This includes a hack. Please fix the code soon!\n";
	# $map_file     = "$genepred.@{[$this->refdb]}mapping.txt";
	# end hack

	my $Coverage;
	my $Length;
	my $dbh = $smash->get_db_handle();
	{
		my ($gene, $length, $cov);
		my $sth = $dbh->prepare("SELECT gene_name, coverage FROM gene_details WHERE gene_prediction_id=?");
		$sth->execute($genepred_id);
		$sth->bind_columns(\$gene, \$cov);
		while ($sth->fetch()) {
			$Coverage->{$gene} = $cov;
		}
		$sth = $dbh->prepare("SELECT external_id, length FROM gene WHERE gene_prediction_id=?");
		$sth->execute($genepred_id);
		$sth->bind_columns(\$gene, \$length);
		while ($sth->fetch()) {
			$Length->{$gene} = $length;
		}
	}
	$smash->close_db_handle();

	$this->analyze_file($map_file, $Coverage, $Length);
}

sub analyze_one_genome {
	my $this         = shift;
	my $label        = $this->label;
	my $smash        = $this->smash;
	my $kegg_map     = "$label.@{[$this->refdb]}mapping.txt";

	my $Coverage = {};
	my $Length = {};
	{
		my ($gene, $length, $cov);
		my $dbh = $smash->get_refgenomedb_handle;
		my $sth = $dbh->prepare("SELECT gene.external_id, gene.length FROM gene WHERE external_id like ?");
		$sth->execute("$label.%");
		$sth->bind_columns(\$gene, \$length);
		while ($sth->fetch()) {
			$Length->{$gene} = $length;
			$Coverage->{$gene} = 1;
		}
	}
	$this->analyze_file($kegg_map, $Coverage, $Length);
}

sub convert_to_abundance {
	my $this          = shift;
	my $feature_count = shift;

	zero_fill_matrix($feature_count);

	# record raw count in case of replicate

	if (defined($this->replicate)) {
		my $file_name = sprintf("%s.feature_vector.rawcount.%s.%s.%s", $this->name, $this->type, $this->refdb, $this->replicate);
		write_R_matrix($file_name, $feature_count);
	}

	return $feature_count;
}

1;

package Smash::Analyses::Comparative::Functional::KEGG;
use strict;
use warnings;
use Smash::Utils::MatrixIO qw(read_R_matrix write_R_matrix transpose_matrix write_two_column_hash);
use Smash::Utils::Tree;

our @ISA = qw(Smash::Analyses::Comparative::Functional);

sub ko2gene                {shift->{KO2GENE}}
sub size_normalize         {shift->{SIZE_NORMALIZE}}

sub init {
	my $this = shift;
	$this->{SIZE_NORMALIZE} = 0;
	$this->SUPER::init();
}

sub analyze_file {
	my $this         = shift;
	my $label        = $this->label;
	my $level        = $this->level;
	my $FeatureCount = $this->feature_count;
	my ($kegg_map, $Coverage, $Length) = @_;

	my $WRITE_HASH = 0;

	my %Assigned;

	# Read the kegg_maps

	open(KEGG, "<$kegg_map") || die "Cannot open $kegg_map: $!";

	my $assigned=0;

	# Process the whole dataset, coverage needs to be used here, since this is the original data

	# Process the kegg map file by calculating coverage of every kegg protein that predicted genes map to

	my $Gene2Kegg = {};	# array per kegg protein that tracks positions covered by sample
	my $KeggCoverage = {};
	my $KeggGeneInheritance = {};
	my $prev_gene = "dummy007";
	LINE:while (my ($gene, $kegg, $gstart, $gend, $kstart, $kend, $ko, $module, $pathway) = $this->parse_next_kegg_map_line(\*KEGG)) {
		last LINE unless $gene;
		if ($gene ne $prev_gene) { 

			# update running sums for all kegg proteins that this gene maps to
			# input here is a BEST-blast hit file. Multiple KEGG proteins per query is 
			# possible when the scores have a tie.
			# In those cases, we should avoid double-counting of two proteins from the same KO. 
			# So we divide by the number 
			# of KEGG genes each query maps to.

			my $keggs = scalar(keys %$Gene2Kegg);
			foreach my $k (keys %$Gene2Kegg) {
				my $overlap = 0;
				map {$overlap++ if ($Gene2Kegg->{$k}->[$_]);} (0..$#{$Gene2Kegg->{$k}});
				my $transferred = ($overlap*$Coverage->{$prev_gene}/$keggs);
				$KeggGeneInheritance->{$k}->{$prev_gene} += $transferred;
			}
			$Gene2Kegg = {};
		}
		map {$Gene2Kegg->{$kegg}->[$_] = 1} ($kstart..$kend);
		$prev_gene = $gene;
	}

	# process the final group after reaching end-of-file

	my $keggs = scalar(keys %$Gene2Kegg);
	foreach my $k (keys %$Gene2Kegg) {
		my $overlap = 0;
		map {$overlap++ if ($Gene2Kegg->{$k}->[$_]);} (0..$#{$Gene2Kegg->{$k}});
		my $transferred = ($overlap*$Coverage->{$prev_gene}/$keggs);
		$KeggGeneInheritance->{$k}->{$prev_gene} += $transferred;
	}

	# normalize coverage by kegg protein length

	# Get protein lengths

	my $dbh = $this->smash->get_refproteindb_handle;

	# scoping for $kegg_length_sth
	{ 
		my $kegg_length_sth = $dbh->prepare("SELECT protein, length FROM kegg_details");
		$kegg_length_sth->execute();
		my ($kegg, $length);
		$kegg_length_sth->bind_columns(\$kegg, \$length);
		while ($kegg_length_sth->fetch()) {
			if ($KeggGeneInheritance->{$kegg}) {
				foreach my $gene (keys %{$KeggGeneInheritance->{$kegg}}) {
					my $transferred = $KeggGeneInheritance->{$kegg}->{$gene} / $length;
					$KeggGeneInheritance->{$kegg}->{$gene} = $transferred;
					#print "$gene\t$kegg\t$transferred\n";
					$KeggCoverage->{$kegg} += $transferred;
				}
			}
		}
	}

=begin DEBUG

write_two_column_hash("$label.kegg_gene_transferred.txt", transpose_matrix($KeggGeneInheritance));

=cut

	if ($this->level ne "kegg") {
		my ($kegg, $ko);
		my ($count, $length);
		my $KOCoverage;

		# transfer coverage from kegg protein to ko

		my $sth = $dbh->prepare("SELECT protein, ko FROM protein2ko");
		$sth->execute();
		$sth->bind_columns(\$kegg, \$ko);
		while ($sth->fetch()) {
			if ($KeggCoverage->{$kegg}) {
				$KOCoverage->{$ko} += $KeggCoverage->{$kegg};
			}
		}

		# create feature vectors in KO

		$assigned = 0;
		while (my ($ko, $coverage) = each %$KOCoverage) {
			$FeatureCount->{$ko}->{$label} = $coverage;
			$assigned += $coverage;
		}

		# record feature vector if required

		if ($WRITE_HASH == 1) {
			my $genepred = $this->smash->genepred;
			open(KO, ">$genepred.sample_ko.txt");
			while (my ($ko, $coverage) = each %$KOCoverage) {
				print KO "$ko\t$label\t$coverage\n";
			}
			close(KO);
		}

		# track gene2function maps

	} else {

		# create feature vectors in KEGG proteins

		$assigned = 0;
		while (my ($kegg, $coverage) = each %$KeggCoverage) {
			$FeatureCount->{$kegg}->{$label} = $coverage;
			$assigned += $coverage;
		}

	}

	$this->smash->close_refproteindb_handle();
	close(KEGG);
	return $assigned;
}

sub parse_next_kegg_map_line {
	my $this = shift;
	my $FH   = shift;
	my $line = $this->get_next_uncommented_line($FH);
	return 0 unless $line;
	my @w = split(/\s+/, $line);
	my ($gene, $kegg, $gstart, $gend, $kstart, $kend, $ko, $module, $pathway) = @w;
	return ($gene, $kegg, $gstart, $gend, $kstart, $kend, $ko, $module, $pathway);
}

sub convert_to_abundance_wrong {

	# convert from KO counts to KO abundance
	# I used to divide the KO abundance by the number of proteins in each KO, thinking EVERY hit above 60bits
	# was used. But apparently it is a best-blast-hit procedure, so I shouldnt do it. The only thing I should
	# do is to divide by the number of KEGG genes each predicted gene hits with the same score, to avoid 
	# double counting. So I removed the normalization of KO abundances by the KO protein_count. Mani, 18.11.2009

}

sub parse_annotations {
	my $this = shift;

	return if ($this->{KO2NAME});

	my $dbh  = $this->smash->get_refproteindb_handle();
	{
		my ($ko, $name);
		my $sth = $dbh->prepare("SELECT ko, name FROM ko2name");
		$sth->execute();
		$sth->bind_columns(\$ko, \$name);
		while ($sth->fetch()) {
			$this->{KO2NAME}->{$ko} = $name;
		}
	}
	{
		my ($ko, $module);
		my $sth = $dbh->prepare("SELECT ko, module FROM ko2module");
		$sth->execute();
		$sth->bind_columns(\$ko, \$module);
		while ($sth->fetch()) {
			$this->{KO2MODULE}->{$ko}->{$module} = 1;
		}
	}
	{
		my ($ko, $pathway);
		my $sth = $dbh->prepare("SELECT ko, pathway FROM ko2pathway");
		$sth->execute();
		$sth->bind_columns(\$ko, \$pathway);
		while ($sth->fetch()) {
			$this->{KO2PATHWAY}->{$ko}->{$pathway} = 1;
		}
	}
	{
		my ($module, $name);
		my $sth = $dbh->prepare("SELECT module, name FROM module2name");
		$sth->execute();
		$sth->bind_columns(\$module, \$name);
		while ($sth->fetch()) {
			$this->{MODULE2NAME}->{$module} = $name;
		}
	}
	{
		my ($module, $pathway);
		my $sth = $dbh->prepare("SELECT module, pathway FROM module2pathway");
		$sth->execute();
		$sth->bind_columns(\$module, \$pathway);
		while ($sth->fetch()) {
			$this->{MODULE2PATHWAY}->{$module}->{$pathway} = 1;
		}
	}
	{
		my ($pathway, $name);
		my $sth = $dbh->prepare("SELECT pathway, name FROM pathway2name");
		$sth->execute();
		$sth->bind_columns(\$pathway, \$name);
		while ($sth->fetch()) {
			$this->{PATHWAY2NAME}->{$pathway} = $name;
		}
	}
	$this->smash->close_refproteindb_handle();

	my $version = $this->smash->get_conf_value("Current Version", "kegg");
	$version =~ s/kegg//;
	my $root;
	my ($KeggKO, $KeggModule, $KeggPathway);

	$KeggKO  = new Smash::Utils::Tree(NAME => "KEGG v$version", TYPE => "functional");
	$root    = new Smash::Utils::Tree::Node(ID => "0", NAME => "root", RANK => "root");
	$KeggKO->add_node($root);
	$KeggKO->set_root($root);
	my $KO2NAME = $this->{KO2NAME};
	foreach my $ko (keys %$KO2NAME) {
		my $name = $KO2NAME->{$ko} || "Unknown";
		my $node = new Smash::Utils::Tree::Node(ID => $ko, NAME => "$ko: $name", RANK => "KO");
		$KeggKO->add_node($node);
		$root->add_child($node);
	}
	$this->{TREE}->{KO} = $KeggKO;

	$KeggModule  = new Smash::Utils::Tree(NAME => "KEGG v$version", TYPE => "functional");
	$root    = new Smash::Utils::Tree::Node(ID => "0", NAME => "root", RANK => "root");
	$KeggModule->add_node($root);
	$KeggModule->set_root($root);
	my $MODULE2NAME = $this->{MODULE2NAME};
	foreach my $module (keys %$MODULE2NAME) {
		my $name = $MODULE2NAME->{$module} || "Unknown";
		my $node = new Smash::Utils::Tree::Node(ID => $module, NAME => "$module: $name", RANK => "Module");
		$KeggModule->add_node($node);
		$root->add_child($node);
	}
	$this->{TREE}->{MODULE} = $KeggModule;

	$KeggPathway  = new Smash::Utils::Tree(NAME => "KEGG v$version", TYPE => "functional");
	$root    = new Smash::Utils::Tree::Node(ID => "0", NAME => "root", RANK => "root");
	$KeggPathway->add_node($root);
	$KeggPathway->set_root($root);
	my $PATHWAY2NAME = $this->{PATHWAY2NAME};
	foreach my $pathway (keys %$PATHWAY2NAME) {
		my $name = $PATHWAY2NAME->{$pathway} || "Unknown";
		my $node = new Smash::Utils::Tree::Node(ID => $pathway, NAME => "$pathway: $name", RANK => "Pathway");
		$KeggPathway->add_node($node);
		$root->add_child($node);
	}
	$this->{TREE}->{PATHWAY} = $KeggPathway;

}

sub get_popup_info {
	my $this     = shift;
	my $features = shift;
	my $level    = lc($this->level);
	my $popup    = {};
	foreach my $query (keys %$features) {
		my $description = $this->{uc($level)."2NAME"}->{$query} || "Unknown";
		my $title = "$query: $description";
		my $text  = "";
		if ($level eq "ko") {
			my $ko = $query;
			$text .= "<H2><B>KEGG Orthology:</B></H2>";
			$text .= "<TEXTFORMAT TABSTOPS=\"[42]\"><A TARGET=\"_itol_kegg\" HREF=\"http://www.kegg.jp/entry/$ko\">$ko</A>\\t: $description</TEXTFORMAT><BR>";
			$text .= "<BR>";
			my @modules  = keys %{$this->{KO2MODULE}->{$ko}};
			my @pathways = keys %{$this->{KO2PATHWAY}->{$ko}};
			if (@modules) {
				$text .= "<H2><B>KEGG Modules:</B></H2>";
				foreach my $module (@modules) {
					my $name = $this->{MODULE2NAME}->{$module} || "Missing";
					$text .= "<TEXTFORMAT TABSTOPS=\"[42]\"><A TARGET=\"_itol_kegg\" HREF=\"http://www.kegg.jp/entry/$module\">$module</A>\\t: $name</TEXTFORMAT><BR>";
				}
				$text .= "<BR>";
			}
			if (@pathways) {
				$text .= "<H2><B>KEGG Pathways:</B></H2>";
				foreach my $pathway (@pathways) {
					my $name = $this->{PATHWAY2NAME}->{$pathway} || "Missing";
					$text .= "<TEXTFORMAT TABSTOPS=\"[42]\"><A TARGET=\"_itol_kegg\" HREF=\"http://www.genome.jp/kegg-bin/show_pathway?ko$pathway+$query\">map$pathway</A>\\t: $name</TEXTFORMAT><BR>";
				}
				$text .= "<BR>";
			}
		} elsif ($level eq "module") {
			my $module = $query;
			$text .= "<H2><B>KEGG Module:</B></H2>";
			$text .= "<TEXTFORMAT TABSTOPS=\"[42]\"><A TARGET=\"_itol_kegg\" HREF=\"http://www.kegg.jp/entry/$module\">$module</A>\\t: $description</TEXTFORMAT><BR>";
			$text .= "<BR>";
			my @pathways = keys %{$this->{MODULE2PATHWAY}->{$module}};
			if (@pathways) {
				$text .= "<H2><B>KEGG Pathways:</B></H2>";
				foreach my $pathway (@pathways) {
					my $name = $this->{PATHWAY2NAME}->{$pathway} || "Missing";
					$text .= "<TEXTFORMAT TABSTOPS=\"[42]\"><A TARGET=\"_itol_kegg\" HREF=\"http://www.genome.jp/kegg-bin/show_pathway?ko$pathway+$query\">map$pathway</A>\\t: $name</TEXTFORMAT><BR>";
				}
				$text .= "<BR>";
			}
		} elsif ($level eq "pathway") {
			my $pathway = $query;
			$text .= "<H2><B>KEGG Pathway:</B></H2>";
			$text .= "<TEXTFORMAT TABSTOPS=\"[42]\"><A TARGET=\"_itol_kegg\" HREF=\"http://www.kegg.jp/entry/ko$pathway\">ko$pathway</A>\\t: $description</TEXTFORMAT><BR>";
		}

		my $html = "<BOX><BR>$text</BOX>";
		$html =~ s/\n//sg;
		$popup->{$query} = "$title\t$html";
	}
	return $popup;
}

sub merge_feature_vectors {
	my $this          = shift;
	my $feature_count = shift;
	my $level         = shift;

	# merge if necessary

	if ($level eq "module" || $level eq "pathway") {
		return $this->merge_feature_vectors_internal($feature_count, "ko", $level);
	} elsif ($level eq "ko") {
		return $feature_count;
	} else {
		die "$this does not know how to summarize at level=$level!";
	}
}

sub merge_feature_vectors_internal {
	my $this          = shift;
	my $feature_count = shift;
	my $from          = shift;
	my $to            = shift;

	my $dbh           = $this->smash->get_refproteindb_handle();

	my $merged_feature_count   = {};

	my $WRITE_HASH = 0;

	# transfer coverage from ko to module

	if ($from eq "ko" && ($to eq "module" || $to eq "pathway") ||
	   ($from eq "protein" && $to eq "ko")) {
		my $function2gene = {};
		my $query = "SELECT $to FROM ${from}2$to WHERE $from=?";
		my $sth = $dbh->prepare($query);
		foreach my $k (keys %$feature_count) {
			$sth->execute($k);
			my @samples = keys %{$feature_count->{$k}};
			while (my ($function) = $sth->fetchrow_array()) {

				# propagate sample2ko mapping to this level
				foreach my $sample (@samples) {
					$merged_feature_count->{$function}->{$sample} += $feature_count->{$k}->{$sample};
				}
			}
		}

		if ($this->size_normalize == 1) {
			my %Count;
			my ($function, $count);
			$query = "SELECT $to, COUNT($from) FROM ${from}2$to GROUP BY $to";
			$sth = $dbh->prepare($query);
			$sth->execute();
			$sth->bind_columns(\$function, \$count);
			while ($sth->fetch()) {
				$Count{$function} = $count;
			}
			foreach my $function (keys %$merged_feature_count) {
				my $count = $Count{$function};
				foreach my $sample (keys %{$merged_feature_count->{$function}}) {
					$merged_feature_count->{$function}->{$sample} /= $count;
				}
				foreach my $gene (keys %{$function2gene->{$function}}) {
					$function2gene->{$function}->{$gene} /= $count;
				}
			}
		}
	}

	$this->smash->close_refproteindb_handle();
	return $merged_feature_count;
}

1;

package Smash::Analyses::Comparative::Functional::eggNOG;
use strict;
use warnings;
use Smash::Utils::MatrixIO qw(read_R_matrix write_R_matrix transpose_matrix read_two_column_hash write_two_column_hash read_multi_column_matrix);

our @ISA = qw(Smash::Analyses::Comparative::Functional);

sub analyze_file {
	my $this         = shift;
	my $label        = $this->label;
	my $level        = $this->level;
	my $FeatureCount = $this->feature_count;
	my $genepred     = $this->smash->genepred;
	my ($eggnog_map, $Coverage, $Length) = @_;

	my $WRITE_HASH = 0;

	my %Assigned;

	# Read the eggnog maps

	open(STRING, "<$eggnog_map") || die "Cannot open $eggnog_map: $!";

	my $assigned=0;

	# Process the whole dataset, coverage needs to be used here, since this is the original data

	# Process the eggNOG map file by calculating coverage of every eggNOG protein and OG that predicted genes map to

	my $Gene2String = {};
	my %StringLength;
	my $StringCoverage = {};
	my %OverlapScore;
	my $prev_gene = "dummy007";
	LINE:while (my ($gene, $string, $sstart, $send, $cog, $length) = parse_next_eggNOG_map_line(\*STRING)) {

		last unless $gene; # strange cases where last line is comment line will not set EOF at the last valid line
				   # it has to be caught this way!

		$StringLength{$string} = $length; # record string protein's length. might do it multiple times, but it's ok!

		if (eof(STRING) || $gene ne $prev_gene) { # last line will match eof

			# if this was the last line, then include the last line to $Gene2String

			if (eof(STRING)) {
				map {$Gene2String->{$cog}->{$string}->[$_] = 1} ($sstart..$send); 
			}

			##################
			# each metagenome gene is mapped to one or more COG(s) via a string protein.
			# we keep track of the mappings to the COG and the protein by a 3-D array $Gene2String.
			# we then transfer the coverage from this gene to string proteins.
			# since multiple proteins can be used to map to the same COG in principle, 
			# we divide by the number of proteins this gene mapped to under this COG.
			##################

			my $coverage = $Coverage->{$prev_gene} || 1;        # If coverage is missing (sometimes Celera does this), make it 1.
			foreach my $c (keys %$Gene2String) {                             # for every COG this gene maps to
				my $strings = scalar(keys %{$Gene2String->{$c}});        #   count number of string proteins under this COG that the gene maps to
				foreach my $s (keys %{$Gene2String->{$c}}) {             #   for every string protein under this COG
					my $overlap = 0;
					map {$overlap++ if ($Gene2String->{$c}->{$s}->[$_]);} (0..$#{$Gene2String->{$c}->{$s}});  # count the number of overlapping positions.
					$StringCoverage->{$c}->{$s} += ($overlap*$coverage/$strings);                             # if there are multiple string proteins, share coverage.
					                                                                                          # when we add them together for the COG later, this will make sense.
					$OverlapScore{$s}{$prev_gene} = $overlap;
				}
			}
			$Gene2String = {};
			if (eof(STRING)) {
				last LINE;
			}
		}
		map {$Gene2String->{$cog}->{$string}->[$_] = 1} ($sstart..$send);
		$prev_gene = $gene;
	}

	# create feature vectors in COG

	########
	# each COG will probably have multiple proteins.
	# each COG gets coverage of each protein under it, normalized by the protein length.
	########

	$assigned = 0;
	foreach my $cog (keys %$StringCoverage) {
		my $sum = 0;
		foreach my $string (keys %{$StringCoverage->{$cog}}) {
			$sum += $StringCoverage->{$cog}->{$string}/$StringLength{$string};
		}
		$FeatureCount->{$cog}->{$label} = $sum;
		$assigned += $sum;
	}

	# record feature vector if required

	if ($WRITE_HASH == 1) {
		open(COG, ">$genepred.sample_cog.txt");
		close(COG);
	}

	close(STRING);
	return $assigned;
}

sub get_popup_info {
	my $this     = shift;
	my $features = shift;
	my $level    = lc($this->level);
	$level       = "og" if $level eq "cog";
	my $popup    = {};
	return $popup if $level eq "funcat";
	foreach my $query (keys %$features) {
		my $description = $this->{uc($level)."2NAME"}->{$query} || "Unknown";
		my $title = "$query: $description";
		my $fun_cat_text = "None";
		if ($this->{OG2FUNCAT}->{$query}) {
			my @fun_cat = split('', $this->{OG2FUNCAT}->{$query});
			if (@fun_cat) {
				$fun_cat_text = join("<BR>", map {"<B>$_</B> - ".$this->{FUNCAT2NAME}->{$_}} @fun_cat);
			}
		}
		my @kos     = keys %{$this->{OG2KO}->{$query}};
		my $ko_text = "";
		if (@kos) {
			$ko_text  = "<BR><TEXTFORMAT TABSTOPS=\"[42]\"><B>KEGG</B>\\t: ";
			$ko_text .= join(",&nbsp;", map {"<A TARGET=\"_itol_kegg\" HREF=\"http://www.genome.jp/dbget-bin/www_bget?ko:$_\">$_</A>"} @kos);
			$ko_text .= "</TEXTFORMAT>";
		}
		my $html =<<EOF;
<BOX>
<H1><B>$description</B></H1>
<BR>
<H2><B>Functional categories: </B></H2>
$fun_cat_text
<BR>
<BR>
<H2><B>External links: </B></H2>
<TEXTFORMAT TABSTOPS="[42]">
<B>eggNOG</B>\\t: <A TARGET="_itol_eggnog" HREF="http://eggnog.embl.de/cgi_bin/display_multi_clusters.pl?number_of_entries=1&species=auto_detect&1=$query">$query</A>
</TEXTFORMAT>
<BR>
<TEXTFORMAT TABSTOPS="[42]">
<B>STRING</B>\\t: <A TARGET="_itol_string" HREF="http://string.embl.de/version_8_3/newstring_cgi/show_network_section.pl?identifier=$query&all_channels_on=1&interactive=yes&network_flavor=evidence&targetmode=cogs">$query</A>
</TEXTFORMAT>
$ko_text
</BOX>
EOF
		$html =~ s/\n//sg;
		$popup->{$query} = "$title\t$html";
	}
	return $popup;
}

sub parse_annotations {
	my $this = shift;

	return if ($this->{OG2NAME});

	my $version         = $this->smash->get_conf_value("Current Version", "eggnog");
	   $version =~ s/eggnog//;
	my $eggnog          = "eggnog$version";
	my $external_dir    = $this->smash->get_smash_conf_value("data_dir")."/external";
	my $og_annot_file   = "$external_dir/${eggnog}_og_annotations.txt";
	my $fun_cat_file    = "$external_dir/${eggnog}_fun_cat_descriptions.txt";
	my $og_fun_cat_file = "$external_dir/${eggnog}_og_fun_cat.txt";
	my $og2ko_file      = "$external_dir/${eggnog}_og_ko_map.txt";
	$this->{OG2KO}      = read_multi_column_matrix($og2ko_file);
	$this->{OG2NAME}    = read_two_column_hash($og_annot_file);
	$this->{OG2FUNCAT}  = read_two_column_hash($og_fun_cat_file);
	my $parser          = new Smash::Config::ConfigParser($fun_cat_file);
	my $conf            = $parser->parse();

	my $root;
	my ($eggNogOG, $eggNogFunCat);

	$eggNogOG  = new Smash::Utils::Tree(NAME => "eggNOG v$version", TYPE => "functional");
	$root    = new Smash::Utils::Tree::Node(ID => "0", NAME => "root", RANK => "root");
	$eggNogOG->add_node($root);
	$eggNogOG->set_root($root);
	my $OG2NAME = $this->{OG2NAME};
	foreach my $og (keys %$OG2NAME) {
		my $name = $OG2NAME->{$og} || "Unknown";
		my $node = new Smash::Utils::Tree::Node(ID => $og, NAME => "$og: $name", RANK => "Orthologous Group");
		$eggNogOG->add_node($node);
		$root->add_child($node);
	}
	$this->{TREE}->{OG}  = $eggNogOG;
	$this->{TREE}->{COG} = $eggNogOG;


	$eggNogFunCat  = new Smash::Utils::Tree(NAME => "eggNOG v$version", TYPE => "functional");
	$root    = new Smash::Utils::Tree::Node(ID => "0", NAME => "root", RANK => "root");
	$eggNogFunCat->add_node($root);
	$eggNogFunCat->set_root($root);
	my $id = 1;
	foreach my $section (keys %$conf) {
		my $high_node = new Smash::Utils::Tree::Node(ID => $id++, NAME => $section, RANK => "Cellular Function");
		$eggNogFunCat->add_node($high_node);
		$root->add_child($high_node);
		foreach my $funcat (keys %{$conf->{$section}}) {
			my $name = $conf->{$section}->{$funcat}; #." ($section)";
			my $node = new Smash::Utils::Tree::Node(ID => $funcat, NAME => "$funcat: $name", RANK => "COG Functional Category");
			$eggNogFunCat->add_node($node);
			$high_node->add_child($node);
			$this->{FUNCAT2NAME}->{$funcat} = $name;
		}
	}
	$this->{TREE}->{FUNCAT}  = $eggNogFunCat;
}

sub merge_feature_vectors {
	my $this          = shift;
	my $feature_count = shift;
	my $level         = shift;

	# merge if necessary

	if ($level eq "funcat") {
		return $this->merge_feature_vectors_internal($feature_count, "og", $level);
	} elsif ($level eq "og" || $level eq "cog") {
		return $feature_count;
	} else {
		die "$this does not know how to summarize at level=$level!";
	}
}

sub merge_feature_vectors_internal {
	my $this          = shift;
	my $feature_count = shift;
	my $from          = shift;
	my $to            = shift;

	my $merged_feature_count   = {};

	# transfer coverage from og to funcat

	if ($from eq "og" && $to eq "funcat") {
		my $function2gene = {};
		OG:foreach my $og (keys %$feature_count) {
			next OG unless $this->{OG2FUNCAT}->{$og};
			my @fun_cat = split('', $this->{OG2FUNCAT}->{$og});
			my @samples = keys %{$feature_count->{$og}};
			foreach my $funcat (@fun_cat) {
				foreach my $sample (@samples) {
					$merged_feature_count->{$funcat}->{$sample} += $feature_count->{$og}->{$sample};
				}
			}
		}
	}

	return $merged_feature_count;
}

#####################
# eggNOG mapping output file looks as follows:
#
# gene, string_protein, COG, string_protein_start, string_protein_end, bit_score, string_protein_length
#
#####################

sub parse_next_eggNOG_map_line {
	my $fh = shift;
	my $line;

	return 0 if eof $fh;

	# get the next non-comment line

	do {
		$line = <$fh>;
	} while (!eof($fh) && $line =~ /^#/);

	return 0 unless $line;
	return 0 if ($line =~ /^#/);
	chomp($line);

	my ($gene, $string, $cog, $sstart, $send, $score, $length) = split(/\s+/, $line);
	return ($gene, $string, $sstart, $send, $cog, $length);
}

1;
