package Smash::Utils::Taxonomy;

use strict;
use warnings;
use Smash::Global qw(:all);
use Smash::Core;
use Smash::Utils::Tree;
use Smash::Utils::HTML qw(get_html_page post_html_page strip_html get_attribute retrieve_remote_file);

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(%Rank2Num %Num2Rank);

# NCBI related stuff

push(@EXPORT_OK, qw(get_taxonomy_for_id get_ordered_taxonomy_for_id));
push(@EXPORT_OK, qw(get_taxonomic_rank get_ncbi_taxonomic_rank www_get_taxonomy_for_complete_name www_get_taxonomy_for_token_set www_get_taxonomy_for_id));

# RDP / Bergey related stuff

push(@EXPORT_OK, qw($NCBITree NCBIEukTree));
push(@EXPORT_OK, qw($BergeyTree));
push(@EXPORT_OK, qw($SmashBergeyTree $SmashNCBITree));
push(@EXPORT_OK, qw(bergey2ncbi bergey2ncbi_name_rank ncbi2bergey ncbi2bergey_name_rank));
push(@EXPORT_OK, qw(%Smash2RDP));

our %EXPORT_TAGS = ('taxonomy' => [grep {$_ =~ /taxonomy/} @EXPORT_OK], 'all' => [@EXPORT_OK]);

# Ivica's tree parsing routine:

my %names;

# General

our %Rank2Num;
our %Num2Rank;

# Trees

our $NCBITree;
our $NCBIEukTree;
our $BergeyTree;
our $SmashNCBITree;
our $SmashBergeyTree;

# Smash

our %Smash2RDP;

=head1 NAME

Smash::Utils::Taxonomy - NCBI and RDP taxonomy related utility functions

=head1 SYNOPSIS

	use Smash::Utils::Taxonomy qw(:all);

	# Parse the NCBI tree from the NCBI taxonomy dump files

	Smash::Utils::Taxonomy::init("NCBI");
	print $NCBITree->root->newick_name;

	# Parse the RDP tree distributed with RDP classifier
	# Only works when rdp_classifier is installed

	Smash::Utils::Taxonomy::init("RDP");
	print $BergeyTree->root->newick_name;

=head1 DESCRIPTION

C<Smash::Utils::Taxonomy> provides several useful functions that are related to NCBI and RDP taxonomic trees.
When C<init()> is called, it populates two variables of type L<Smash::Utils::Tree|Utils::Tree>. These are 
C<$NCBITree> and C<$BergeyTree> for NCBI and RDP trees, respectively. These tree objects can
be manipulated or queried using all the methods from L<Smash::Utils::Tree|Utils::Tree> as well as special
methods implemented here.

=cut

sub init {

	my $parse_flag = shift;

	my $BERGEY = 0;
	my $NCBI   = 0; # 1 - prok, euk separately; 2 - prok, euk with Chloroplast hack

	# parse smash.conf
	# Set global variables

	Smash::Core->new()->init();

	if ($parse_flag) {
		if ($parse_flag =~ /RDP_REGROUP/i) {
			$BERGEY = 2;
			$NCBI = 2;
		} elsif ($parse_flag =~ /BERGEY|RDP/i) {
			$BERGEY = 1;
		} elsif ($parse_flag =~ /NCBI.DIV/i) {
			$NCBI = 2;
			$BERGEY = 1;
		} elsif ($parse_flag =~ /NCBI/i) {
			$NCBI = 1;
		} 
	}

	if ($BERGEY) {
		if (!$BergeyTree) {
			# This also populates $Smash2RDP->{}
			$BergeyTree = get_rdp_taxonomy_tree();
			# Add regrouped genera
			if ($BERGEY == 2) {
				my @new_nodes =("10001:genus:Lachnospiraceae_g1:family:Lachnospiraceae",
						"10002:genus:Lachnospiraceae_g2:family:Lachnospiraceae",
						"10003:genus:Lachnospiraceae_g3:family:Lachnospiraceae",
						"10004:genus:Lachnospiraceae_g4:family:Lachnospiraceae",
						"10005:genus:Lachnospiraceae_g5:family:Lachnospiraceae",
						"10006:genus:Lachnospiraceae_g6:family:Lachnospiraceae",
						"10007:genus:Lachnospiraceae_g7:family:Lachnospiraceae",
						"10008:genus:Lachnospiraceae_g8:family:Lachnospiraceae",
						"10009:genus:Lachnospiraceae_g9:family:Lachnospiraceae",
						"10010:genus:Lachnospiraceae_g10:family:Lachnospiraceae",
						"10011:genus:Lachnospiraceae_g11:family:Lachnospiraceae",
						"10101:genus:Ruminococcaceae_g1:family:Ruminococcaceae",
						"10102:genus:Ruminococcaceae_g2:family:Ruminococcaceae",
						"10201:genus:Clostridiales_g1:order:Clostridiales",
						);
				foreach my $n (@new_nodes) {
					my ($id, $rank, $name, $parent_rank, $parent_name) = split(":", $n);
					my $pid = $BergeyTree->get_id_for_ranked_name($parent_rank, $parent_name);
					my $parent = $BergeyTree->nodes->{$pid};
					my $node = Smash::Utils::Tree::Node->new(ID => $id, NAME => $name, RANK => $rank);
					$BergeyTree->add_node($node);
					$parent->add_child($node);
				}
			}
		}
	}

	if ($NCBI) {
		if (!$TAXONOMY_LOCAL_LOCATION) {
			die "Cannot find the local taxonomy repository! Did you forget to set 'local_repository' in the [Taxonomy] section of smash.conf?\n";
		}
		if (!(-e "$TAXONOMY_LOCAL_LOCATION/nodes.dmp" && -e "$TAXONOMY_LOCAL_LOCATION/names.dmp" && -e "$TAXONOMY_LOCAL_LOCATION/merged.dmp")) {
			print STDERR "Missing NCBI Taxonomy dump files in $TAXONOMY_LOCAL_LOCATION\n";
			update_files();
		}

		if (!$NCBITree) {
			($NCBITree, $NCBIEukTree) = get_ncbi_taxonomy_tree($parse_flag);
			if ($NCBI == 2) { # Chloroplast hack

				# Add new nodes for "Cyanobacteria <class>" and "Chloroplast"
				# Move genera listed under Chloroplast in RDP to the new node

				my $id                  = $NCBITree->get_id_for_ranked_name("phylum", "Cyanobacteria");
				my $ncbi_cyanobacteria  = $NCBITree->nodes->{$id};

				my $ncbi_cyano_class    = Smash::Utils::Tree::Node->new(NAME => "Cyanobacteria <class>", ID => 1234567890, RANK => "class", PARENT => $ncbi_cyanobacteria->id);
				$NCBITree->add_node($ncbi_cyano_class);

				my $ncbi_chloroplast    = Smash::Utils::Tree::Node->new(NAME => "Chloroplast", ID => 1234567891, RANK => "class", PARENT => $ncbi_cyanobacteria->id);
				$NCBITree->add_node($ncbi_chloroplast);

				my $rdp_chl_id          = $BergeyTree->get_id_for_ranked_name("family", "Chloroplast");
				my $rdp_chloroplast     = $BergeyTree->nodes->{$rdp_chl_id};
				foreach my $rdp_node (values %{$rdp_chloroplast->children}) {
					my $id;
					foreach my $rank ("no rank", qw(genus family class phylum)) {
						$id = $NCBIEukTree->get_id_for_ranked_name($rank, $rdp_node->name);
						last if $id;
					}
					if ($id) {
						my $node = $NCBIEukTree->nodes->{$id};
						$node->parentlink->remove_child($node);
						$node->{RANK} = $rdp_node->rank;
						$node->{PARENT} = $ncbi_chloroplast->id;
						$NCBITree->add_node($node);
					}
				}
			}
		}
		#$NCBITree->print_newick();
	}

	# ranks of the levels

	############################################################################
	#
	# Here's the hierarchy of ranks from:
	#   International Code of Nomenclature of Bacteria
	#   Bacteriological Code, 1990 Revision
	#   Chapter 3, Rules of Nomenclature with Recommendations
	#   http://www.ncbi.nlm.nih.gov/books/NBK8808/#A191
	#
	# my @ranks    = qw(domain superkingdom kingdom superphylum phylum class subclass order suborder family subfamily tribe subtribe supergenus genus subgenus speciesgroup speciessubgroup species subspecies norank);
	#
	# We practically care only about phylum, class, order, family and genus. 
	# The subranks subdomain, subphylum, ..., etc are parsed, but these will 
	# be moved up one level if they are the terminal ranks.
	# Also, we dont care about Incertae Cedis, since it is no information we can use.
	#
	############################################################################

	my @ranks    = qw(domain superkingdom kingdom superphylum phylum class subclass order suborder family subfamily tribe subtribe supergenus genus subgenus speciesgroup speciessubgroup species subspecies norank);
	%Rank2Num = map {$ranks[$_] => $_} 0..$#ranks;
	$Rank2Num{"no rank"} = $Rank2Num{"norank"};
	$Rank2Num{"no_rank"} = $Rank2Num{"norank"};
	$Rank2Num{"species group"} = $Rank2Num{"speciesgroup"};
	$Rank2Num{"species_group"} = $Rank2Num{"speciesgroup"};
	$Rank2Num{"species subgroup"} = $Rank2Num{"speciessubgroup"};
	$Rank2Num{"species_subgroup"} = $Rank2Num{"speciessubgroup"};
	%Num2Rank = map {$_ => $ranks[$_]} 0..$#ranks;
}

sub get_smash_tree {

	# Get pruned trees:
	# For Bergey, get the subtree with all the RDP ranks in genomes
	# For NCBI, copy the bergey tree, and transform each node into an NCBI node, creating a novel node if necessary (e.g. Group II)

	# Hack: NCBI and RDP don't agree on the taxa in the tree. These will be fixed here.
	# I use a hash called %Remap that maps things back and forth
	# 1. NCBI merged Silicibacter to Ruegeria

	my %Remap = (
			#$BergeyTree->get_id_for_ranked_name("genus", "Silicibacter") => $BergeyTree->get_id_for_ranked_name("genus", "Ruegeria"),
		);

	# end hack

	my %bergey_keep;
	foreach my $tax (keys %Smash2RDP) {
		my $rdp_parent = $Smash2RDP{$tax};
		#if ($Remap{$rdp_parent}) {
		#	$rdp_parent = $Remap{$rdp_parent};
		#}
		$bergey_keep{$rdp_parent} = 1;
	}
	$SmashBergeyTree = $BergeyTree->get_pruned_tree([keys %bergey_keep]);

	# Now add the genomes to the pruned Bergey tree.
	# Each genome goes as a child to the RDP classified parent.

	foreach my $tax (keys %Smash2RDP) {
		my $rdp_parent = $Smash2RDP{$tax};
		#if ($Remap{$rdp_parent}) {
		#	$rdp_parent = $Remap{$rdp_parent};
		#}
		my $node = $NCBITree->nodes->{$tax};
		if (!$node) {
			die "ERROR: taxonomy id $tax not found in NCBI taxonomy dump!\n";
		}
		my $n = $node->copy_node();
		$n->{ID} = "NCBI.".$n->id;
		$n->{PARENTLINK} = undef;
		$n->{PARENT} = undef;
		$SmashBergeyTree->add_node($n);
		$SmashBergeyTree->nodes->{$rdp_parent}->add_child($n);
	}

	# now copy bergey tree and convert into NCBI ids

	$SmashNCBITree = $SmashBergeyTree->copy_tree();
	$SmashNCBITree->{NAME} = "NCBI";
	delete $SmashNCBITree->{RANKEDNAME2ID};

	my $NewNodes =  {};
	my $NewRankedName2Id = {};
	my $new_ncbi_tax_id = 2200000000;

	# Remap each node to NCBI id
	# And fix the rankedname2tax map

	foreach my $node (values %{$SmashNCBITree->nodes}) {
		my $id = $node->id;
		my $new_id;
		# Leaves already have NCBI ids, so no need to convert bergey2ncbi
		if ($id =~ /^NCBI/) {
			$new_id = $id;
			$new_id =~ s/^NCBI\.//;
		} else {
			if ($Remap{$id}) {
				$id = $Remap{$id};
			}
			# NCBI merged Silicibacter into Ruegeria, but RDP keeps them separate
			# So make a new node for Silicibater in SmashNCBITree, with Ruegeria's lineage
			if ($node->name eq "Silicibacter") {
				my $x = $SmashBergeyTree->get_id_for_ranked_name("genus", "Ruegeria");
				my $n = $SmashNCBITree->nodes->{$x};
				$new_id = $new_ncbi_tax_id;
				$new_ncbi_tax_id++;
				$node->parentlink->remove_child($node);
				$n->parentlink->add_child($node);
			} else {
				$new_id = bergey2ncbi($id);
				if (!$new_id) {
					$new_id = $new_ncbi_tax_id;
					$new_ncbi_tax_id++;
				}
			}
		}
		$node->{ID} = $new_id;
		$NewNodes->{$new_id} = $node;
		$NewRankedName2Id->{$node->rank}->{$node->name} = $new_id;
	}

	# Fix parents
	foreach my $node (values %{$SmashNCBITree->nodes}) {
		$node->{PARENT} = $node->parentlink->id if $node->parentlink;
	}

	delete $SmashNCBITree->{NODES};
	$SmashNCBITree->{NODES} = $NewNodes;
	$SmashNCBITree->{RANKEDNAME2ID} = $NewRankedName2Id;
	$SmashNCBITree->{ROOT} = $SmashNCBITree->rootlink->id;
}

########################
# Get the Bergey taxonomy information
########################

sub get_rdp_taxonomy_tree {

	my $tree = new Smash::Utils::Tree(NAME => "RDP", TYPE => "phylogenetic");

	# parse tree

	my $smash     = new Smash::Core();
	   $smash->init();
	my ($rdp_dir) = $smash->software_dir("rdp_classifier", "current");
	   $smash->finish();
	my $tree_file = "$rdp_dir/bergeyTrainingTree.xml";
	open(TREE, "<$tree_file") || die "Cannot open $tree_file: $!";
	LINE: while (<TREE>) {
		chomp();
		if (!m/<TreeNode/) {
			next LINE;
		}

		my $name   = get_attribute($_, "name");
		my $rank   = get_attribute($_, "rank");
		my $tax    = get_attribute($_, "taxid");
		my $parent = get_attribute($_, "parentTaxid");

		# remove quotes

		$name =~ s/^['"]//;
		$name =~ s/['"]$//;
		$name =~ s/^&quot;//;
		$name =~ s/&quot;$//;


		# RDP uses -1 for global root, but we use it for unknown, so change -1 to -10000

		$parent = -10000 if $parent == -1;

		# Make this node

		my $node   = new Smash::Utils::Tree::Node(ID => $tax, NAME => $name, RANK => $rank, PARENT => $parent);

		# Add this node to the tree

		$tree->add_node($node);

		# Store root if this is root

		if (lc($name) eq "root") {
			$node->{NAME} = "root";
			$node->{RANK} = "no rank";
			$tree->set_root($node);
		}
	}
	close(TREE);

	# Prefix ambiguous names
	# E.g., If something's called just "Incertae Cedis XIII" or so, we prepend the previous rank to resolve ambiguity

	my @ids = sort {$a <=> $b} keys %{$tree->nodes};
	foreach my $id (@ids) {
		my $node = $tree->nodes->{$id};
		my $name = $node->name;
		if ($name =~ /^\bGp[0-9IVX]+[a-z]?\b/ || $name =~ /\bFamily [0-9IVX]+\b/ || $name =~ /\bIncertae Sedis [0-9IVX]+\b/) {
			my $parent = $node->parentlink->name;
			$name   = "$parent $name";
			$node->{NAME} = $name;
			$tree->{RANKEDNAME2ID}->{$node->rank}->{$name} = $id;
		}
	}

	# parse SMASH2RDP

	my $remap_file = "$rdp_dir/smash2rdp.txt";
	open(SMASH2RDP, "<$remap_file") || die "Cannot find Smash-to-RDP remapping file at $remap_file: $!";
	while (<SMASH2RDP>) {
		chomp();
		next if m/^\s*#/;
		my ($smash_id, $ncbi_id, $rdp_id) = split(/\t/);
		$Smash2RDP{$smash_id} = $rdp_id;
	}
	close(SMASH2RDP);

	#$tree->print_tree();
	return $tree;
}

########################
# Get the NCBI taxonomy information
# Since parent of a node can appear after the child in the tree file,
# we have to parse them separately and add_child() in the second round
########################

sub get_ncbi_taxonomy_tree {

	use Fcntl qw(:seek);

	# parse the ncbi dump files

	my $nodes_file = "$TAXONOMY_LOCAL_LOCATION/nodes.dmp";
	my $names_file = "$TAXONOMY_LOCAL_LOCATION/names.dmp";
	my $merge_file = "$TAXONOMY_LOCAL_LOCATION/merged.dmp";

	my $prok_tree = new Smash::Utils::Tree(NAME => "NCBI", TYPE => "phylogenetic");

	open(NODES, "<$nodes_file") || die "Cannot open $nodes_file: $!";
	while (<NODES>) {
		chomp();
		s/\t\|$//;
		my ($tax_id, $parent, $rank) = split(/\t\|\t/);

		# Make this node

		my $node   = new Smash::Utils::Tree::Node(ID => $tax_id, RANK => $rank, PARENT => $parent);
		$prok_tree->nodes->{$tax_id} = $node;
	}
	close(NODES);

	open(MERGED, "<$merge_file") || die "Cannot open $merge_file: $!";
	while (<MERGED>) {
		chomp();
		s/\t\|$//;
		my ($old_id, $new_id) = split(/\t\|\t/);

		# Also give names to the old taxonomy_id's

		$prok_tree->nodes->{$old_id} = $prok_tree->nodes->{$new_id};
	}
	close(MERGED);

	# Make Eukaryote tree, but just copy the nodes. Don't duplicate

	my $euk_tree = new Smash::Utils::Tree(NAME => "NCBI Eukaryotes");
	$euk_tree->{NODES} = $prok_tree->nodes;

	# parse names.dmp


	my @extensions = qw(.prok .euk);
	my @trees      = ($prok_tree, $euk_tree);
	foreach my $i (0..$#extensions) {
		my $tree = $trees[$i];
		my $file = $names_file.$extensions[$i];
		make_prokaryote_eukaryote_tree() unless -f $file;
		open(NAMES, "<$file") || die "Cannot open $file: $!";
		NAME: while (<NAMES>) {
			chomp();
			s/\t\|$//;
			my ($tax_id, $name, $uname, $class) = split(/\t\|\t/);
			$tree->{VALIDNODE}->{$tax_id} = 1;
			my $node = $tree->nodes->{$tax_id};
			if (!defined $node) {
				die "$tax_id in $file not found in tree $tree!";
			}
			my $rank = $node->rank;
			if ($class eq "synonym") {
				$tree->{RANKEDNAME2ID}->{$rank}->{$name} = $tax_id;
			} elsif ($class eq "scientific name") {
				$tree->{RANKEDNAME2ID}->{$rank}->{$name} = $tax_id;
				if ($uname) {
					warn "$rank: $name\n" if $tree->{RANKEDNAME2ID}->{$rank}->{$uname};
					$tree->{RANKEDNAME2ID}->{$rank}->{$uname} = $tax_id;
				}
				$node->{NAME} = $name;

				# Store root
				# Root's parent is also root, so dont add a child to its parent.

				if (lc($name) eq "root") {
					$tree->set_root($node);
				} else {
					my $parent = $tree->nodes->{$node->parent};
					$parent->add_child($node);
				}
			}

		}
		close(NAMES);
	}

	return ($prok_tree, $euk_tree);
}

sub make_prokaryote_eukaryote_tree {

	use Fcntl qw(:seek);

	my $ROOT;
	my $CELL;
	my $BACTERIA;
	my $ARCHAEA;
	my $EUKARYOTA;
	my %Tax2Parent;

	# parse the ncbi dump files

	my $nodes_file = "$TAXONOMY_LOCAL_LOCATION/nodes.dmp";
	my $names_file = "$TAXONOMY_LOCAL_LOCATION/names.dmp";

	open(NODES, "<$nodes_file") || die "Cannot open $nodes_file: $!";
	while (<NODES>) {
		chomp();
		s/\t\|$//;
		my ($tax_id, $parent, $rank) = split(/\t\|\t/);
		$Tax2Parent{$tax_id} = $parent;
	}
	close(NODES);

	# parse names.dmp

	open(NAMES, "<$names_file") || die "Cannot open $names_file: $!";

	# get $BACTERIA and $ARCHAEA

	PRESCAN: while (<NAMES>) {
		chomp();
		s/\t\|$//;
		my ($tax_id, $name, $uname, $class) = split(/\t\|\t/);
		if ($uname eq "Bacteria <prokaryote>" && $class eq "scientific name") {
			$BACTERIA = $tax_id;
		} elsif ($name eq "Archaea" && $class eq "scientific name") {
			$ARCHAEA = $tax_id;
		} elsif ($name eq "Eukaryota" && $class eq "scientific name") {
			$EUKARYOTA = $tax_id;
		} elsif ($name eq "root" && $class eq "scientific name") {
			$ROOT = $tax_id;
		} elsif ($name eq "cellular organisms" && $class eq "scientific name") {
			$CELL = $tax_id;
		}
		last PRESCAN if ($BACTERIA && $ARCHAEA && $CELL);
	}

	die "Could not find ids for Bacteria and Archaea!\n" if (!$BACTERIA && !$ARCHAEA);

	# reparse

	open(PROK, ">$names_file.prok") || die "Cannot open $names_file.prok";
	open(EUK,  ">$names_file.euk" ) || die "Cannot open $names_file.euk";
	seek(NAMES, 0, SEEK_SET);
	NAME: while (my $line = <NAMES>) {
		chomp($line);
		$line =~ s/\t\|$//;
		my ($tax_id, $name, $uname, $class) = split(/\t\|\t/, $line);
		if ($tax_id == $ROOT || $tax_id == $CELL) {
			print PROK "$line\n";
			print EUK  "$line\n";
			next NAME;
		}
		my $node = $tax_id;
		do {
			if ($node == $BACTERIA || $node == $ARCHAEA) {
				print PROK "$line\n";
				next NAME;
			} elsif ($node == $EUKARYOTA) {
				print EUK  "$line\n";
				next NAME;
			}
		} while (($node = $Tax2Parent{$node}) != 1);
	}
	close(NAMES);
	close(EUK);
}

=head1 Functions accessing the local tree objects

=over 4

=item B<init($type)>

Initializes C<$NCBITree> and/or C<$BergeyTree> objects based on C<$type>. When
C<$type> is C<"RDP">, it parses the RDP/Bergey tree from C<rdp_classifier>.
When C<$type> is C<"NCBI">, it parses the NCBI tree from the NCBI taxonomy dump
files.

=item B<update_files()>

Updates the NCBI taxonomy dump files by retrieving the latest from NCBI website.

=item B<bergey2ncbi($id)>

Returns the NCBI taxonomy id corresponding to the relevant taxon from RDP tree.

=item B<ncbi2bergey($id)>

Returns the RDP taxonomy id corresponding to the relevant taxon from NCBI tree.

=cut

sub update_files {
	use File::Temp qw(tempdir);
	use File::Copy;
	use File::Path;
	select(STDERR); $| = 1; select(STDOUT);
	my $filename = "taxdump.tar.gz";
	my $tmpdir = tempdir(CLEANUP => 1);
	print STDERR "NCBI Taxonomy repository:\n";
	print STDERR "Remote location: $NCBI_TAXONOMY_REMOTE_LOCATION\n";
	print STDERR "Local  location: $TAXONOMY_LOCAL_LOCATION\n";
	print STDERR "Retrieving remote files from remote location ...";
	system("wget --progress=dot:mega -O $tmpdir/$filename $NCBI_TAXONOMY_REMOTE_LOCATION/$filename") == 0 || die "download failed: $!"; 
	print STDERR " done\n";
	print STDERR "Extracting files ...";
	system("cd $tmpdir && tar xfz $filename") == 0 || die "File extraction failed: $!";
	print STDERR " done\n";
	if (! -d "$TAXONOMY_LOCAL_LOCATION") {
		mkpath $TAXONOMY_LOCAL_LOCATION;
	}
	copy "$tmpdir/nodes.dmp", "$TAXONOMY_LOCAL_LOCATION/";
	copy "$tmpdir/names.dmp", "$TAXONOMY_LOCAL_LOCATION/";
	copy "$tmpdir/merged.dmp", "$TAXONOMY_LOCAL_LOCATION/";
	print STDERR "Making prokaryote/eukaryote trees ...";
	make_prokaryote_eukaryote_tree();
	print STDERR " done\n";
	print STDERR "Local repository updated.\n";
}

# map from bergey taxonomy to ncbi taxonomy

sub bergey2ncbi {
	my $id = shift;
	my $node = $BergeyTree->nodes->{$id};
	die "Node for $id not found in RDP tree!\n" unless $node;
	my $new_id = bergey2ncbi_name_rank($node->name, $node->rank);
	return $new_id;
}

sub bergey2ncbi_name_rank {
	my ($name, $rank) = @_;

	# hacks

	if ($name eq "Escherichia/Shigella" && $rank eq "genus") {
		$name = "Escherichia";
	} elsif ($name eq "GpIIa" && $rank eq "genus") {
		$name = "Prochlorococcus";
	} elsif ($name eq "Parasutterella" && $rank eq "genus") {
		$rank = "subtribe";
	} elsif ($name eq "Marvinbryantia" && $rank eq "genus") {
		$name = "Bryantella";
	} elsif ($name eq "Gracilibacteraceae" && $rank eq "family") {
		$name = "Graciibacteraceae";
	} elsif ($name eq "Archaea" && $rank eq "domain") {
		$rank = "superkingdom";
	} elsif ($name eq "Bacteria" && $rank eq "domain") {
		$rank = "superkingdom";
	} elsif ($name eq "Cyanobacteria" && $rank eq "phylum") {
		$name = "Cyanobacteria/Chloroplast";
	} elsif ($name eq "TM7" && $rank eq "phylum") {
		$rank = "no rank";
		$name = "candidate division TM7";
	#} elsif ($name eq "Cyanobacteria" && $rank eq "class") {
	#	$rank = "phylum";
	} elsif ($name eq "Root" && $rank eq "norank") {
		$rank = "no rank";
		$name = "root";
	}

	my $ncbi_tax_id = $NCBITree->get_id_for_ranked_name($rank, $name) 
			|| $NCBITree->get_id_for_ranked_name($rank, "$name ($rank)") 
			|| $NCBITree->get_id_for_ranked_name("no rank", $name);

	return $ncbi_tax_id;
}

sub ncbi2bergey {
	my $id = shift;
	my $node = $NCBITree->nodes->{$id};
	my $new_id = ncbi2bergey_name_rank($node->name, $node->rank);
	warn "WARNING: NCBI_id=$id\n" unless $new_id;
	return $new_id;
}

sub ncbi2bergey_name_rank {
	my ($name, $rank) = @_;

	# hacks

	if (($name eq "Escherichia" || $name eq "Shigella") && $rank eq "genus") {
		$name = "Escherichia/Shigella";
	} elsif ($name eq "Prochlorococcus" && $rank eq "genus") {
		$name = "GpIIa";
	} elsif ($name eq "Parasutterella" && $rank eq "subtribe") {
		$rank = "genus";
	} elsif ($name eq "Bryantella" && $rank eq "genus") {
		$name = "Marvinbryantia";
	} elsif ($name eq "Graciibacteraceae" && $rank eq "family") {
		$name = "Gracilibacteraceae";
	} elsif ($name eq "Archaea" && $rank eq "superkingdom") {
		$rank = "domain";
	} elsif ($name eq "Bacteria" && $rank eq "superkingdom") {
		$rank = "domain";
	} elsif ($name eq "Cyanobacteria/Chloroplast" && $rank eq "phylum") {
		$name = "Cyanobacteria";
	} elsif ($name eq "candidate division TM7" && $rank eq "no rank") {
		$rank = "phylum";
		$name = "TM7";
	} elsif ($name eq "root" && $rank eq "no rank") {
		$rank = "norank";
		$name = "Root";
	}

	my $rdp_tax_id = $BergeyTree->get_id_for_ranked_name($rank, $name);
	warn "WARNING: No RDP/Bergey match for $rank=$name\n" unless $rdp_tax_id;

	return $rdp_tax_id;
}

=item B<get_taxonomy_for_id>

Same as C<www_get_taxonomy_for_id>, only gets the information from local NCBI taxonomy
dump files. Called as:

	$NCBITree->get_taxonomy_for_id(435590);
	$BergeyTree->get_taxonomy_for_id(443);

=item B<get_ordered_taxonomy_for_id>

Same as C<get_taxonomy_for_id>, but prepends the ordinal rank followed by
underscore so that the hash can be sorted using the ordinal rank.

=item B<get_ncbi_taxonomic_rank($tax_id, $rank)>

returns the NCBI taxonomy id of the ancestor of C<$tax_id> at 
C<$rank>. For example, C<get_ncbi_taxonomic_rank(435590, 'genus')> returns 816
(435590 is for B<"Bacteroides vulgatus"> and 816 is for B<"Bacteroides">).
If you want the name, then use C<$NCBITree-E<gt>nodes-E<gt>{816}-E<gt>name> to 
get B<"Bacteroides"> back.

=back

=cut

sub get_phylum {
	my $tree = shift;
	my $tax_id = shift;
	my $candidate = $tree->nodes->{tax_id};
	SEARCH:while (1) {
		my $rank = $candidate->rank;
		if ($rank =~ /phylum$/) {
			return $candidate->id;
		} elsif ($candidate->id == $tree->root) {
			return undef;
		}
		$candidate = $candidate->parentlink;
	}
}

sub get_taxonomic_rank {
	my $tree = shift;
	my ($in_tax_id, $required_rank) = @_;
	my $node = $tree->nodes->{$in_tax_id};
	my $required_level = $Rank2Num{$required_rank};
	do {
		if (ref($node) !~ /Smash/) {
			warn "No node for $in_tax_id, $required_rank\n";
			die  "Perhaps you should update the $tree database?\n";
		}
		if (!defined($Rank2Num{$node->rank})) {
			warn "WARNING: $node has no rank!\n";
			$node->print_node(0, \*STDERR);
		} elsif ($Rank2Num{$node->rank} <= $required_level) {
			return $node->id;
		}
	} while (($node = $node->parentlink) && ($node->id != $tree->root));
	return $in_tax_id;
}

sub get_ordered_taxonomy_for_id {
	my $tree = shift;
	my $tax_id = shift;
	my %Taxonomy;
	my $node = $tree->nodes->{$tax_id};
	my $level = 0;
	SEARCH:while ($node->id > 1) {
		my $rank = $node->rank;
		$rank =~ s/ /_/g;
		$rank = "${level}__$rank";
		$Taxonomy{$rank} = $node->name;
		$node = $node->parentlink;
		$level++;
	}
	return %Taxonomy;
}

sub get_taxonomy_for_id {
	my $tree = shift;
	my $id = shift;
	my %taxonomy = get_ordered_taxonomy_for_id($tree, $id);
	%taxonomy = map {my $k=$_; s/^[0-9]+__//; $_ => $taxonomy{$k}} keys %taxonomy;
	$taxonomy{ncbi_name} = $tree->nodes->{$id}->name;
	$taxonomy{tax_id}    = $id;
	return %taxonomy;
}

sub get_ncbi_taxonomic_rank {
	my ($tax_id, $required_rank) = @_;

	$required_rank = "no rank" if $required_rank eq "norank";

	# ranks of the levels

	my $required_level = $Rank2Num{$required_rank};

	# if you dont understand what they want, just give the input back

	if (!defined($required_level)) {
		warn "Do not understand $required_rank";
		return $tax_id;
	}

	my @lineage = $NCBITree->get_full_lineage($tax_id);
	
	if ($required_rank eq "species") {

		# hack for Clostridium botulinum
		# since there are many clades under it, i have a priority list. if you hit this, you assign it
		my $hit;
		foreach my $problem (591968, 36826, 36827, 36828, 36829, 36830, 36831, 447213) {
			if (scalar(grep {$problem == $_} @lineage) > 0) {
#warn "$tax_id\t$problem\n";
				return $problem;
			}
		}
		# end


		# hack for wrong species

		my %Remap = (
				# E. coli species
				556266 => 562, # Shigella sp. D9
				457400 => 562, # Escherichia sp. 1_1_43
				469598 => 562, # Escherichia sp. 3_2_53FAA
				457401 => 562, # Escherichia sp. 4_1_40B

				1765   => 1773, # Mycobacterium bovis

				264    => 263, # http://www.bacterio.cict.fr/f/francisella.html#novicida
				234    => 29459,# http://www.bacterio.cict.fr/b/brucella.html#melitensis
				);
		foreach my $problem (keys %Remap) {
			if (scalar(grep {$problem == $_} @lineage) > 0) {
#warn "$tax_id\t".$Remap{$problem}."\n";
				return $Remap{$problem};
			}
		}

		# end
	}


	my $new_id;

	# try to get the right level
	TAXONOMY:foreach my $candidate (@lineage) {
		my $rank = $NCBITree->nodes->{$candidate}->rank;
		my $level = $Rank2Num{$rank};
		if (defined($level) && $level == $required_level) {
			$new_id = $candidate;
			last TAXONOMY;
		}
	}

	# failed to get the right level
	# so get the closest higher than required-1
	# e.g.,
	# cellular organisms; Bacteria; Firmicutes; Clostridia; Clostridiales; unclassified Clostridiales; unclassified Clostridiales (miscellaneous); butyrate-producing bacterium SR1/5
	# genus returns "unclassified Clostridiales (miscellaneous)"

	if (!$new_id) {
		my $attempted_id;
		TAXONOMY:foreach my $i (0..($#lineage-1)) {
			my $candidate = $lineage[$i];
			my $rank = $NCBITree->nodes->{$candidate}->rank;
			my $level = $Rank2Num{$rank};
			if (defined($level)) {
				if ($level == $required_level+1) {
					return $lineage[$i+1];
				}
			}
		}
	}

	# so get the closest lower than required+1
	# e.g.,
	# cellular organisms; Bacteria; Firmicutes; Clostridia; Clostridiales; unclassified Clostridiales; unclassified Clostridiales (miscellaneous); butyrate-producing bacterium SR1/5
	# genus returns "unclassified Clostridiales"

	if (!$new_id) {
		my $attempted_id;
		TAXONOMY:foreach my $candidate (@lineage) {
			my $rank = $NCBITree->nodes->{$candidate}->rank;
			my $level = $Rank2Num{$rank};
			if (defined($level)) {
				if ($level < $required_level) {
					$new_id = $attempted_id;
					last TAXONOMY;
				} else {
					$attempted_id = $candidate;
				}
			}
		}
	}

	if (!$new_id) {
		my $name = $NCBITree->nodes->{$tax_id}->rank;
		printf STDERR "%d (%s) does not have %s\n", $tax_id, $name || "name unknown", $required_rank;
		$new_id = $tax_id;
	}
	return $new_id;
}

=head1 Functions that query NCBI Taxonomy website directly

C<Smash::Utils::Taxonomy> provides several useful functions that are related to NCBI taxonomy ids.
These functions return a hash with keys are NCBI ranks and values as the rank values. E.g.,

	%taxonomy = (   
			tax_id  => 435590,
			name    => "Bacteroides vulgatus ATCC 8482", 
			species => "Bacteroides vulgatus",
			genus   => "Bacteroides",
			family  => "Bacteroidaceae",
			order   => "Bacteroidales",
			class   => "Bacteroidia",
			phylum  => "Bacteroidetes",
		   superkingdom => "Bacteria");

These functions are NOT object oriented, so you would call them as:

	my %taxonomy = www_get_taxonomy_for_id(435590);

=over 4

=item B<www_get_taxonomy_for_complete_name>

Performs a search using "complete name" mode on NCBI taxonomy website and
returns the taxonomy if found.

=item B<www_get_taxonomy_for_token_set>

Performs a search using "token set" mode on NCBI taxonomy website and
returns the taxonomy if found.

=item B<www_get_taxonomy_for_id>

Gets the taxonomy for given taxonomy id from NCBI taxonomy website.

=back

=cut
	
sub www_get_taxonomy_for_complete_name {
	www_get_taxonomy_by_searching(shift, 1);
}

sub www_get_taxonomy_for_token_set {
	www_get_taxonomy_by_searching(shift, 3);
}

sub www_get_taxonomy_for_id {
	www_get_taxonomy_by_searching(shift, 5);
}

sub www_get_taxonomy_by_searching {
	my $srchtext = shift;
	my $srchmode = shift;
	my %taxonomy;
	my $name;
	my $tax_id;
	my $page;
	my $rank;
	my ($begin, $end);
	my $buf; 
	my %params = (name => $srchtext, srchmode => $srchmode, keep => 0, lvl => 3);

	###############################################
	# Post the request to NCBI Taxonomy Page
	###############################################

	$page = post_html_page("http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi", \%params);

	###############################################
	# Get the required section
	###############################################

	if ($page !~ m#<em>Taxonomy ID:#) {
		($page) = ($page =~ m^<LI TYPE=circle>(<A.*?</A>)^);
		$page || return;
		my ($url) = get_attribute($page, "href");
		$url   = "http://www.ncbi.nlm.nih.gov$url";
		$page   = get_html_page($url);
	}
	($page) = ($page =~ m#<!--  the contents   -->(.*)</dl>#s);

	###############################################
	# Get NCBI taxonomy ID
	###############################################

	$buf   = $page;
	($tax_id) = ($buf =~ m^<em>Taxonomy ID:.*?</em>(.*?)<br>^s);
	($rank)   = ($buf =~ m^<em>Rank:.*?</em>(.*?)<br>^s);

	###############################################
	# Get NCBI standard name
	###############################################

	$buf   = $page;
	if ($buf =~ m#<h2>(.*)</h2>#s) {
		$name = strip_html($1);
	}

	###############################################
	# Get the full NCBI lineage
	###############################################

	my @pieces = ($buf =~ m#(<a ALT=.*?</a>)#sg);

	###############################################
	# Some pages have a slightly different format
	###############################################

	if (!@pieces) {
		@pieces = ($buf =~ m#(<A HREF="Taxonomy".*?</A>);#sg);
	}

	###############################################
	# Parse the lineage into a hash
	###############################################

	my $level = 1;
	%taxonomy = map {
		my (undef, $key) = m#ALT=(["'])([^\1]*?)\1#;
		my $val = strip_html($_);
		$key =~ s/ /_/g;
		if ($key eq "no_rank") {
			$key .= "_$level";
		}
		$level++;
		$key => $val
	} @pieces;
	$taxonomy{ncbi_name} = $name;
	$taxonomy{tax_id} = $tax_id;
	if ($rank ne "no rank" && !defined($taxonomy{$rank})) {
		$taxonomy{$rank} = $name;
	}
	return %taxonomy;
}

1;
