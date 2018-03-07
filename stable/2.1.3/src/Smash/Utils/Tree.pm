#!/usr/bin/env perl

package Smash::Utils::Tree;

use strict;
use warnings;
use Smash::Utils::PhyloXML qw(make_tree_from_phyloxml);
use overload '""' => '_overload';

=head1 NAME

Smash::Utils::Tree - Object encapsulating a tree structure and associated functions.

=head1 SYNOPSIS

	use Smash::Utils::Tree;

	# make tree

	my $tree = new Smash::Utils::Tree(NAME => "NCBI", TYPE => "phylogenetic");

	# add root

	my $root = new Smash::Utils::Tree::Node(ID => 1, \
						NAME => "root", \
						RANK => "no rank", \
						PARENT => -10000);
	$tree->add_node($root);
	$tree->set_root($root);

	# add a child to the root

	my $bacteria = new Smash::Utils::Tree::Node(ID => 2, \
						NAME => "Bacteria", \
						RANK => "superkingdom");
	$tree->add_node($bacteria);
	$root->add_child($bacteria);

	my $bac_chlor = new Smash::Utils::Tree::Node(ID => 3, \
						NAME => "Bacteroidetes/Chlorobi", \
						RANK => "superphylum");
	$tree->add_node($bac_chlor);
	$bacteria->add_child($bac_chlor);

	my $bacteroidetes = new Smash::Utils::Tree::Node(ID => 4, \
						NAME => "Bacteroidetes", \
						RANK => "phylum");
	$tree->add_node($bacteroidetes);
	$bac_chlor->add_child($bacteroidetes);

	# get information for a node

	my $node = $tree->nodes->{3};
	$node->id;             # 3
	$node->rank;           # superphylum
	$node->name;           # Bacteroidetes/Chlorobi
	$node->newick_name;    # Bacteroidetes_Chlorobi

=head1 DESCRIPTION

C<Smash:Utils::Tree> provides functions to manipulate phylogenetic
trees as objects.

=head1 FUNCTIONS

=head2 Properties

=over 4

=item C<name>

returns the name of the tree.

=item C<type>

returns the type of the tree (commonly used are B<"phylogenetic"> and B<"functional">).

=item C<root>

returns the id of the root of the tree.

=item C<rootlink>

returns the link to the root of the tree.

=item C<nodes>

returns the hash containing all nodes, indexed by the id's.

=item C<get_id_for_ranked_name($rank, $name)>

returns the id of the node for a given name at the given rank. E.g.,

	$tree->get_id_for_ranked_name("genus", "Bacteroides"); # returns 816 in NCBI tree

=back

=head2 Tree manipulation functions

=over 4

=item C<add_node($node)>

adds the node to this tree.

=item C<remove_node($node)>

removes the node from this tree.

=item C<count_leaves>

Returns the number of leaf nodes.

=item C<delete_data>

Deletes data at every node.

=item C<propagate_data_to_root($node_id, $value)>

Propagates $value at DATA all the way upto root from the given node.

=item C<prune_tree>

Prunes the tree so that nodes without values for the DATA field will
be removed. Needs DATA to be set in the nodes to remain, otherwise
the whole tree will disappear and it cannot be undone.
For example, if you have a list of node-ids that
you want to keep, then set DATA to 1 in each of these nodes, 
propagate it all the way up to ROOT, then prune the tree.

=item C<prune_tree_recursive>

Internal recursive routine called by C<prune_tree>.

=item C<remove_straight_links>

Removes internal nodes with just one child. For example,

	    |
	____|____
	|       |
	|       |
	A       |
	|       |
	|       |
	B       |
	|       |
	|       |
	C       D

In this case, C is connected to A directly and becomes

	    |
	____|____
	|       |
	|       |
	|       |
	|       |
	|       |
	|       |
	|       |
	|       |
	C       D

=item C<propagate_data_average>

Propagate values set in the DATA field bottom-up by setting the value
at a node as the average of its immediate children.

=item C<propagate_data_average_recursive>

Internal recursive routine called by C<propagate_data_average>.

=item C<print_newick>

Prints the tree in Newick format. This prints the whole tree. If
you want a subtree to be printed, prune the tree beforehand using
the DATA field. 
For example, if you have a list of node-ids that
you want to keep, then set DATA to 1 in each of these nodes, 
propagate it all the way up to ROOT, then prune the tree.

=item C<print_newick_recursive>

Internal recursive routine called by C<print_newick>.

=item C<print_tree>

Prints the tree in XML format.

=item C<print_tree_recursive>

Internal recursive routine called by C<print_tree>.

=item C<print_data>

Prints the tree in XML format, only for nodes with DATA being set.
Similar to C<print_tree>.

=item C<print_data_recursive>

Internal recursive routine called by C<print_data>.

=item C<dump_tree>

Dumps two files: F<names.dmp> and F<nodes.dmp> in NCBI taxonomy
dump format corresponding to the tree.

=back

=cut

sub name          {shift->{NAME}}
sub type          {shift->{TYPE}}
sub root          {shift->{ROOT}}
sub rootlink      {shift->{ROOTLINK}}
sub nodes         {shift->{NODES}}
sub validnode     {shift->{VALIDNODE}}
sub rankedname2id {shift->{RANKEDNAME2ID}}
sub _overload     {return shift->{NAME}}

sub new {
	my $class = shift;
	my %params = @_;
	my $this;
	if ($params{PHYLOXML}) {
		$this = make_tree_from_phyloxml($params{PHYLOXML}, $params{NAME}, $params{TYPE});
	} else {
		$this = bless {%params}, $class;
		if (!$this->{NODES}) {
			$this->{NODES} = {};
		}
	}
	return $this;
}

sub set_root {
	my $this = shift;
	my $node = shift;

	$this->{ROOT}     = $node->id;
	$this->{ROOTLINK} = $node;
	$node->{PARENT}   = -10000;
	$node->{PARENTLINK} = undef;
}

sub add_unknown {
	my $this = shift;
	my $node = Smash::Utils::Tree::Node->new(ID => -1, NAME => "Unknown", RANK => "no rank", PARENT => $this->root);
	$this->add_node($node);
}

sub get_id_for_ranked_name {
	my $this = shift;
	my ($rank, $name) = @_;
	return $this->{RANKEDNAME2ID}->{$rank}->{$name};
}

sub add_node {
	my $this = shift;
	my $node = shift;

	# Set values in the tree for this node

	$this->{NODES}->{$node->id} = $node;

	# Add this node as a child to its parent

	if (defined($node->parent) && ($node->parent ne "-10000") && $this->{NODES}->{$node->parent}) {
		$this->{NODES}->{$node->parent}->add_child($node);
	}

	# Add name mapping if information exists

	if ($node->rank && $node->name) {
		$this->{RANKEDNAME2ID}->{$node->rank}->{$node->name} = $node->id;
	}
}

sub remove_node {
	my $this = shift;
	my $node = shift;
	$node->parentlink->remove_child($node) if ($node->parentlink);
	delete $this->{NODES}->{$node->id};
}

sub copy_tree {
	my $this = shift;
	my $copy = new Smash::Utils::Tree(NAME => $this->name, TYPE => $this->type);
	$this->copy_nodes_recursive($copy, $this->rootlink);
	my %rankedname2id = %{$this->{RANKEDNAME2ID}};
	$copy->{RANKEDNAME2ID} = \%rankedname2id;
	return $copy;
}

sub copy_nodes_recursive {
	my $this = shift;
	my $copy = shift;
	my $curr = shift;
	my $parent = shift;

	# Copy the current node, and make it the child of new parent node

	my $node = $curr->copy_node();
	$copy->add_node($node);

	# If parent is not sent through args, this is the root

	if ($parent) {
		$parent->add_child($node);
	} else {
		$copy->set_root($node);
	}

	# Traverse through each child, and recursively call with child and copy of curr
	foreach my $child (values %{$curr->children}) {
		$this->copy_nodes_recursive($copy, $child, $node);
	}
}

sub get_leaf_ids {
	my $this = shift;
	my $list = [];
	$this->get_leaf_ids_recursive($this->rootlink, $list);
	return $list;
}

sub get_leaf_ids_recursive {
	my $this = shift;
	my $curr = shift;
	my $list = shift;

	if ($curr->nchildren == 0) {
		push(@$list, $curr->id);
	} else {
		foreach my $child (values %{$curr->children}) {
			$this->get_leaf_ids_recursive($child, $list);
		}
	}
}

sub reset_ids {
	my $this = shift;
	return $this->reset_ids_recursive($this->rootlink);
}

sub reset_ids_recursive {
	my $this = shift;
	my $key  = shift;
	my $curr = shift;

		my $id = $curr->id;
		if ($id ne $key) {
			$this->nodes->{$id}   = $this->nodes->{$key};
			delete $this->nodes->{$key};
		}

	if ($curr->nchildren == 0) {
		return 1;
	} else {
		foreach my $key (keys %{$curr->children}) {
			$this->reset_ids_recursive($key, $curr->children->{$key});
		}
	}
}

sub count_leaves {
	my $this = shift;
	return $this->count_leaves_recursive($this->rootlink);
}

sub count_leaves_recursive {
	my $this = shift;
	my $curr = shift;

	if ($curr->nchildren == 0) {
		return 1;
	} else {
		my $count = 0;
		foreach my $child (values %{$curr->children}) {
			$count += $this->count_leaves_recursive($child);
		}
		return $count;
	}
}

sub delete_data {
	my $this = shift;
	$this->delete_data_recursive($this->rootlink);
}

sub delete_data_recursive {
	my $this = shift;
	my $curr = shift;

	delete $curr->{DATA};
	foreach my $child (values %{$curr->children}) {
		$this->delete_data_recursive($child);
	}
}

sub get_pruned_tree {
	my $this  = shift;
	my $list  = shift;

	my $PrunedTree = new Smash::Utils::Tree(NAME => $this->name, TYPE => $this->type);
	my $root = $this->rootlink->copy_node();
	$PrunedTree->add_node($root);
	$PrunedTree->set_root($root);

	foreach my $id (@$list) {
		my $node = $this->nodes->{$id};
		die "Cannot find node for $id in tree" unless defined($node);
		my $node_copy = $PrunedTree->nodes->{$node->id};
		if (!defined($node_copy)) {
			$node_copy = $node->copy_node();
			$PrunedTree->add_node($node_copy);
			# Re-id'd by NCBI
			if ($node_copy->id ne $id) {
				$PrunedTree->nodes->{$id} = $node_copy;
			}
		}
		while ($node_copy->id ne $PrunedTree->root) {
			my $parent = $node->parentlink;
			my $parent_copy = $PrunedTree->nodes->{$parent->id};
			if (defined($parent_copy)) {
				$parent_copy->add_child($node_copy);
				last;
			} else {
				$parent_copy = $parent->copy_node();
				$PrunedTree->add_node($parent_copy);
				$parent_copy->add_child($node_copy);
				$node = $parent;
				$node_copy = $parent_copy;
			}
		}
	}

	return $PrunedTree;
}

sub prune_tree {
	my $this = shift;
	$this->prune_tree_recursive($this->rootlink);
}

sub prune_tree_recursive {
	my $this = shift;
	my $curr = shift;

	foreach my $child (values %{$curr->children}) {
		$this->prune_tree_recursive($child);
	}
	if (!$curr->data) {
		$this->remove_node($curr);
	}
}

sub remove_straight_links {
	my $this = shift;
	$this->remove_straight_links_recursive($this->rootlink);
}

sub remove_straight_links_recursive {
	my $this = shift;
	my $curr = shift;

	foreach my $child (values %{$curr->children}) {
		$this->remove_straight_links_recursive($child);
	}
	return if ($curr->id eq $this->root);
	my $parent = $curr->parentlink;
	if ($curr->nchildren == 1) {
		my ($child) = values %{$curr->children};
		$curr->remove_child($child);
		$this->remove_node($curr);
		$parent->add_child($child);
	}
}

sub propagate_data_to_root {
	my $this = shift;
	my $id   = shift;
	my $val  = shift;
	my $node = $this->nodes->{$id};
	do {
		$node->{DATA} = $val;
	} while ($node->id != $this->root && ($node = $node->parentlink));

	# The order of the && clause is important, since we want to set DATA
	# at root as well. If you switch the members of the clause, then
	# DATA will not be set for ROOT, and that will be a problem.
}

sub propagate_data_to_leaves {
	my $this = shift;
	die "DATA undefined at root" unless defined($this->rootlink->data);
	$this->propagate_data_to_leaves_recursive($this->rootlink);
}

sub propagate_data_to_leaves_recursive {
	my $this = shift;
	my $curr = shift;

	foreach my $child (values %{$curr->children}) {
		if (!defined($child->data)) {
			$child->{DATA} = $curr->data;
		}
		$this->propagate_data_to_leaves_recursive($child);
	}
}

sub propagate_data_average_of_leaves_to_root {
	my $this = shift;
	$this->propagate_data_average_of_leaves_to_root_recursive($this->rootlink, 0);
}

sub propagate_data_average_of_leaves_to_root_recursive {
	my $this = shift;
	my $curr = shift;
	my $sum  = shift;

	if (!defined($curr->data)) {
		my $sum = 0;
		my $count = 0;
		foreach my $child (values %{$curr->children}) {
			my $child_data = $this->propagate_data_average_of_leaves_to_root_recursive($child, $sum);
			if ($child_data) {
				$sum += $child_data;
				$count++;
			}
		}
		$curr->{DATA} = $sum/$count if $count > 0;
	}

	return $curr->data;
}

sub propagate_data_average_of_children_to_root {
	my $this = shift;
	$this->propagate_data_average_of_children_to_root_recursive($this->rootlink);
}

sub propagate_data_average_of_children_to_root_recursive {
	my $this = shift;
	my $curr = shift;

	if (!defined($curr->data)) {
		my $sum = 0;
		my $count = 0;
		foreach my $child (values %{$curr->children}) {
			my $child_data = $this->propagate_data_average_of_children_to_root_recursive($child);
			if ($child_data) {
				$sum += $child_data;
				$count++;
			}
		}
		$curr->{DATA} = $sum/$count if $count > 0;
	}

	return $curr->data;
}

sub print_newick {
	my $this = shift;
	$this->print_newick_name(@_);
}

sub print_newick_id {
	my $this = shift;
	my $FH   = shift || \*STDOUT;
	$this->print_newick_recursive($this->rootlink, 0, $FH);
	print $FH ";\n";
}

sub print_newick_name {
	my $this = shift;
	my $FH   = shift || \*STDOUT;
	$this->print_newick_recursive($this->rootlink, 1, $FH);
	print $FH ";\n";
}

sub print_newick_recursive {
	my $this   = shift;
	my $curr   = shift;
	my $name   = shift;
	my $FH     = shift;

	if ($curr->nchildren > 0) {
		print $FH "(";
		my @children = sort {$a->id cmp $b->id} values %{$curr->children};
		my $last_child = pop(@children);
		foreach my $child (@children) {
			$this->print_newick_recursive($child, $name, $FH);
			print $FH ",";
		}
		$this->print_newick_recursive($last_child, $name, $FH);
		print $FH ")";
		my $newick_name = $curr->newick_name;
		print $FH ($name)?$newick_name:$curr->id;
		if ($curr->{"branch_length"}) {
			print $FH ":".$curr->{branch_length};
		}
		if ($curr->{"bootstrap"}) {
			print $FH "[".$curr->{"bootstrap"}."]";
		}
	} else {
		my $newick_name = $curr->newick_name;
		print $FH ($name)?$newick_name:$curr->id;
		if ($curr->{"branch_length"}) {
			print $FH ":".$curr->{"branch_length"};
		}
		if ($curr->{"bootstrap"}) {
			print $FH "[".$curr->{"bootstrap"}."]";
		}
	}
}

sub print_phyloxml {
	my $this = shift;
	my $FH   = shift || \*STDOUT;

	printf $FH <<EOL;
<?xml version="1.0" encoding="UTF-8"?>
<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd" xmlns="http://www.phyloxml.org">
EOL

	printf $FH '<phylogeny rooted="false" name="%s" source="%s">', $this->name, "SmashCommunity v1.6";
	print $FH "\n";
	$this->print_phyloxml_recursive($this->rootlink, 0, $FH);
	print $FH "</phylogeny>\n";
	print $FH "\n</phyloxml>\n";
}

sub print_phyloxml_recursive {
	my $this   = shift;
	my $curr   = shift;
	my $indent = shift;
	my $FH     = shift;

	$curr->print_node_phyloxml_open ($indent, $FH) unless ($curr->id eq $this->root);
	$curr->print_node_phyloxml($indent, $FH);
	foreach my $child (sort {$a->name cmp $b->name} values %{$curr->children}) {
		$this->print_phyloxml_recursive($child, $indent+1, $FH);
	}
	$curr->print_node_phyloxml_close($indent, $FH) unless ($curr->id eq $this->root);
}

sub print_tree {
	my $this = shift;
	my $FH   = shift || \*STDOUT;
	printf $FH '<Tree name="%s" source="%s">', $this->name, "SmashCommunity v1.6";
	print $FH "\n";
	$this->print_tree_recursive($this->rootlink, 1, $FH);
	print $FH "</Tree>\n";
}

sub print_tree_recursive {
	my $this   = shift;
	my $curr   = shift;
	my $indent = shift;
	my $FH     = shift;

	$curr->print_node($indent, $FH);

	foreach my $child (sort {$a->id cmp $b->id} values %{$curr->children}) {
		$this->print_tree_recursive($child, $indent+1, $FH);
	}
}

sub print_data {
	my $this = shift;
	my $FH   = shift || \*STDOUT;
	printf $FH '<Tree name="%s" source="%s">', $this->name, "SmashCommunity v1.6";
	print $FH "\n";
	$this->print_data_recursive($this->rootlink, 1, $FH);
	print $FH "</Tree>\n";
}

sub print_data_recursive {
	my $this   = shift;
	my $curr   = shift;
	my $indent = shift;
	my $FH     = shift;

	$curr->print_node($indent, $FH) if $curr->data;

	foreach my $child (sort {$a->id cmp $b->id} values %{$curr->children}) {
		$this->print_data_recursive($child, $indent+1, $FH);
	}
}

sub dump_tree {
	my $this = shift;
	my ($NAMES_FH, $NODES_FH);
	open($NAMES_FH, ">names.dmp") || die "Cannot open $NAMES_FH: $!";
	open($NODES_FH, ">nodes.dmp") || die "Cannot open $NODES_FH: $!";
	$this->dump_tree_recursive($this->rootlink, $NAMES_FH, $NODES_FH);
	close($NAMES_FH);
	close($NODES_FH);
}

sub dump_tree_recursive {
	my $this     = shift;
	my $curr     = shift;
	my $NAMES_FH = shift;
	my $NODES_FH = shift;

	printf $NAMES_FH join("\t|\t", $curr->id, $curr->name, "", "scientific name\t|\n");
	if ($curr->id eq $this->root) {
		printf $NODES_FH join("\t|\t", $curr->id, $curr->id,     $curr->rank, "", 8, 0,  1, 0, 0, 0, 0, 0, "\t|\n");
	} else {
		printf $NODES_FH join("\t|\t", $curr->id, $curr->parent, $curr->rank, "", 0, 1, 11, 1, 0, 1, 0, 0, "\t|\n");
	}

	foreach my $child (sort {$a->id cmp $b->id} values %{$curr->children}) {
		$this->dump_tree_recursive($child, $NAMES_FH, $NODES_FH);
	}
}

=head2 Phylogenetic tree-specific functions

=over 4

=item C<get_full_lineage()>

returns list containing the ids of all the ancestral nodes
in order (bottom-to-top).

=item C<get_full_lineage_string_for_id()>

returns the taxonomic lineage as string in the following
format:

	<lineage superkingdom="Bacteria" phylum="Bacteroidetes">

=item C<get_common_ancestor()>

returns the lowest common ancestor (LCA) for a given list
of node ids (calls one of the two versions specified below).

=item C<get_common_ancestor_for_two()>

returns the lowest common ancestor (LCA) for a given list
of exactly two node ids.

=item C<get_common_ancestor_generic()>

returns the lowest common ancestor (LCA) for a given list
of node ids.

=back

=cut

sub get_full_lineage_string_for_id {
	my $this   = shift;
	my $tax_id = shift;
	my @lineage = reverse $this->get_full_lineage($tax_id);
	my $string = "";
	foreach my $n (@lineage) {
		my $node = $this->nodes->{$n};
		my $key = $node->rank;
		my $val = $node->name;
		$key =~ s/ /_/g;
		$string .= "$key=\"$val\" ";
	}
	return "<lineage $string/>";
}

sub get_full_lineage {
	my $this   = shift;
	my $tax_id = shift;

	if (!defined($this->nodes->{$tax_id}->{PARENT})) {
		die "Cannot find match for $tax_id. Please update $this taxonomy information.";
	}

	my $candidate = $this->nodes->{$tax_id};
	my @lineage = ();
	SEARCH:do {
		push(@lineage, $candidate->id);
	} while (($candidate->id ne $this->root) && ($candidate = $candidate->parentlink));
	return @lineage;
}

sub get_common_ancestor_for_two {
	my $this = shift;
	my ($a, $b) = @_;
	my @la = get_full_lineage($a);
	my @lb = get_full_lineage($b);
	foreach my $anode (@la) {
		foreach my $bnode (@lb) {
			if ($anode eq $bnode) {
				return $anode;
			}
		}
	}
	return -1;
}

sub get_common_ancestor_generic {
	my $this = shift;
	my @x    = @_;
	my $rca  = shift @x; #recent common ancestor
	while (my $node = shift(@x)) {
		$rca = $this->get_common_ancestor_for_two($rca, $node);
	}
	return $rca;
}

sub get_common_ancestor {
	my $this = shift;
	my @x    = @_;
	my $count = scalar(@x);
	if ($count == 0) {
		return -1;
	} elsif ($count == 1) {
		return $x[0];
	} elsif ($count == 2) {
		return $this->get_common_ancestor_for_two(@x);
	} else {
		return $this->get_common_ancestor_generic(@x);
	}
}

1;

package Smash::Utils::Tree::Node;

use strict;
use warnings;
use overload '""' => '_overload';

=head1 NAME

Smash::Utils::Tree::Node - Object encapsulating the node of a tree structure.

=head1 DESCRIPTION

=head1 FUNCTIONS

Properties of the object and methods.

=head2 Properties

=over 4

=item C<id>

returns taxonomy id.

=item C<name>

returns the name of the organism at the node.

=item C<newick_name>

returns the name that is newick compatible.

=item C<rank>

returns the rank of the organism at the node.

=item C<parent>

returns the id of the parent node.

=item C<parentlink>

returns the parent node object.

=item C<nchildren>

returns the number of children.

=item C<children>

returns the hash containing all children nodes.

=item C<data>

returns the data stored at the given node.

=back

=cut

sub id        {shift->{ID}}
sub name      {shift->{NAME}}
sub rank      {shift->{RANK}}
sub parent    {shift->{PARENT}}
sub parentlink{shift->{PARENTLINK}}
sub data      {shift->{DATA}}
sub children  {shift->{CHILDREN}}
sub nchildren {shift->{NCHILDREN}}
sub _overload {shift->{NAME}}

sub new {
	my $class = shift;
	my %params = @_;
	my $this = bless {%params}, $class;
	$this->{NCHILDREN} = 0;
	$this->{CHILDREN} = {};
	return $this;
}

sub newick_name {
	my $this = shift;
	my $name = $this->name;
	$name =~ tr/\(\) :,/<>___/;

	###
	# Hack for showing only species names
	# Variation paper Fig S1.
	###

	 $name = join("_", (split("_", $name))[0,1]) if ($name !~ /_sp\._|obeum|torques|prausnitzii|butyrate-producing|Clostridiales_bacterium|Erysipelotrichaceae_bacterium/ && $name =~ /_/);

	###
	# End Hack
	###

	return $name;
}

=head2 Tree manipulation functions

=over 4

=item C<add_child($child)>

add the given node to the children of this node.

=item C<remove_child($child)>

removes the given from the children of this node.

=item C<copy_node()>

returns a new C<Smash::Utils::Tree::Node> object that is
an exact copy of this node.

=back

=cut

sub add_child {
	my $this = shift;
	my $baby = shift;
	my $children = $this->children;
	if (!$children->{$baby->id}) {
		$children->{$baby->id} = $baby;
		$this->{NCHILDREN}++;
		$baby->{PARENTLINK} = $this;
		$baby->{PARENT} = $this->id;
	}
}

sub remove_child {
	my $this = shift;
	my $baby = shift;
	my $children = $this->children;
	if ($children->{$baby->id}) {
		delete $children->{$baby->id};
		$this->{NCHILDREN}--;
		delete $baby->{PARENTLINK};
		delete $baby->{PARENT};
	}
}

sub remove_all_children {
	my $this = shift;
	my $children = $this->children;
	my @children = keys %$children;
	foreach my $child (@children) {
		delete $children->{$child};
	}
	$this->{NCHILDREN} = 0;
}

sub copy_node {
	my $this = shift;
	my $buffer = {};
	foreach my $key (keys %$this) {
		$buffer->{$key} = $this->{$key};
	}
	my $copy = Smash::Utils::Tree::Node->new(%$buffer);
	return $copy;
}

sub print_node {
	my $this   = shift;
	my $indent = shift || 0;
	my $FH     = shift || \*STDOUT;
	printf $FH ("\t" x $indent);
	printf $FH ('<Node id="%s" name="%s" rank="%s" parent="%s" nchildren="%d"', $this->id, $this->name, $this->rank, ($this->parent > 0)?$this->parent:"none", $this->nchildren);
	printf $FH (' data="%.2f"', $this->data) if $this->data;
	print  $FH (" />\n");
}

sub print_node_phyloxml_open {
	my $this   = shift;
	my $indent = shift || 0;
	my $FH     = shift || \*STDOUT;
	printf $FH ("  " x $indent);
	printf $FH "<clade>";
	printf $FH "\n";
}

sub print_node_phyloxml_close {
	my $this   = shift;
	my $indent = shift || 0;
	my $FH     = shift || \*STDOUT;
	printf $FH ("  " x $indent);
	printf $FH "</clade>";
	printf $FH "\n";
}

sub print_node_phyloxml {
	my $this   = shift;
	my $indent = shift || 0;
	my $FH     = shift || \*STDOUT;
	my $leader = ("  " x ($indent+1));
	if ($this->{branch_length}) {
		printf $FH ("%s<branch_length>%s</branch_length>\n", $leader, $this->{branch_length});
	}
	if ($this->{bootstrap}) {
		printf $FH ("%s<confidence type=\"bootstrap\">%d</confidence>\n", $leader, $this->{bootstrap});
	}
	if ($this->{taxonomy}) {
		printf $FH "%s<taxonomy>\n", $leader;
		printf $FH "%s  <code>%d</code>\n", $leader, $this->{ID};
		printf $FH "%s</taxonomy>\n", $leader;
	}
}

1;
