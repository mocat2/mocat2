#!/usr/bin/env perl

package Smash::Utils::iTOL;

use strict;
use warnings;
use POSIX;
use File::Basename;
use Math::Round;

use Smash::Core qw(:all);
use Smash::Global qw(:all);
use Smash::Utils::Taxonomy qw(:all);
use Smash::Utils::MatrixIO qw(:all);

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(SEPARATOR);
our %EXPORT_TAGS = ('all' => [@EXPORT_OK]);

=head1 NAME

Smash::Utils::iTOL - Wrapper module for interfacing with iTOL

=head1 SYNOPSIS

	use Smash::Utils::iTOL;

	# Make the tree

	my $tree = Smash::Utils::Tree->new(NAME => "NCBI");
	...

	# Populate data to plot
	# $node_id in $data should correspond to the node ids in $tree

	my $data = {};
	$data->{$node_id}->{sample1} = 0.1;
	$data->{$node_id}->{sample2} = 0.2;
	...

	# Make iTOL object

	my $itol = Smash::Utils::iTOL->new(
					NAME => "MetaHIT", 
					TREE => $tree, 
					DESCRIPTION => "MetaHIT phylo composition",
					PICTURE_SIZE => 2000,
					LABELS => ["sample1", "sample2"],
					UPLOAD_ID => "7wztysd"
					);
	$itol->init();
	$itol->add_dataset(
			NAME => "Danish",
			DATA => $data,
			GRAPH_TYPE => "barplot",
			SCALE => 100
			);
	$itol->finish();

=head1 DESCRIPTION

Smash:Utils::iTOL provides a wrapper class to interface with the iTOL webtool
in order to upload and download trees along with datasets. 

=head1 FUNCTIONS

=head2 Object-oriented functions

=over 4

=item C<new()>

makes a new iTOL object.

=item C<init()>

initializes the object by opening the relevant files.

=item C<add_dataset(dataset_name, $dataset_hash_reference)>

Adds a dataset to the existing iTOL upload job. Passing the reference to the dataset hash could be destructive.
If you need it intact, send a copy.

=item C<finish()>

winds up by finalizing the files, uploading the tree and downloading it.

=item C<upload()>

uploads the tree and the associated datasets in the iTOL object.

=item C<download()>

downloads the tree from iTOL.

=back

=cut

my $PROGRESS = \*STDERR;
select($PROGRESS); $| = 1; select(STDOUT);

sub SEPARATOR {":::"}

sub name             {shift->{NAME}}
sub uploadID         {shift->{UPLOAD_ID}}
sub tree_id          {shift->{TREE_ID}}
sub tree             {shift->{TREE}}
sub description      {shift->{DESCRIPTION}}
sub labels           {shift->{LABELS}}
sub dataset_name     {shift->{DATASET_NAME}}
sub dataset_count    {shift->{DATASET_COUNT}}
sub dataset_colors   {shift->{DATASET_COLORS}}
sub log_scale        {shift->{LOG_SCALE}}
sub precision        {shift->{PRECISION}}
sub picture_size     {shift->{PICTURE_SIZE}}
sub graph_width      {shift->{GRAPH_WIDTH}}
sub display_threshold{shift->{DISPLAY_THRESHOLD}}
sub ncbi_taxonomy_id {shift->{NCBI_TAXONOMY_ID}}
sub print_scale      {shift->{PRINT_SCALE}}
sub tree_file        {shift->{TREE_FILE}}
sub ncbi_file        {shift->{NCBI_FILE}}
sub merge_level      {shift->{MERGE_LEVEL}}
sub color_file       {shift->{COLOR_FILE}}
sub config_file      {shift->{CONFIG_FILE}}
sub config_fh        {shift->{CONFIG_FH}}
sub post_fix         {shift->{POST_FIX}}
sub display_mode     {shift->{DISPLAY_MODE} || "normal"}
sub font_size        {shift->{FONT_SIZE} || 20}
sub line_width       {shift->{LINE_WIDTH} || 3}

sub new {
	my $class = shift;
	my %params = @_;
	my $this  = bless {%params}, $class;

	my $name  = $this->name;
	$this->{DATASET_COUNT} = 0;
	$this->{CONFIG_FILE} = "$name.cfg";
	$this->{TREE_FILE}   = "$name.tree";
	$this->{NCBI_FILE}   = "$name.taxids";
	$this->{COLOR_FILE}  = "$name.colors";
	$this->{POPUP_INFO_FILE} = "$name.popup";

	# fill in defaults if missing

	my %defaults = (
			LOG_SCALE => 0,
			PRINT_SCALE => 0,
			DISPLAY_THRESHOLD => 0,
			NCBI_TAXONOMY_ID => 0,
			NODES_FOR_TREE => {},
			DATASET_COLORS => {},
			PRECISION => 1,
			GRAPH_WIDTH => 100,
			POPUP_INFO => undef,
			TREE => 0,
			);

	foreach my $property (keys %defaults) {
		if (!defined($this->{$property})) {
			$this->{$property} = $defaults{$property};
		}
	}

	return $this;
}

sub init {
	my $this = shift;
	my $name = $this->name;
	my $desc = $this->description;
	my $tree_file = $this->tree_file;
	my $ncbi_file = $this->ncbi_file;
	my $color_file = $this->color_file;
	my $uploadID = $this->uploadID;


	my $config_file = $this->config_file;
	my $treeName = "$desc";
	$treeName .= sprintf("; Time:%s", scalar localtime);

	my $FH;
	open($FH, ">$config_file") || die "Cannot open $config_file: $!";

	print $FH <<EOF;
treeName    = $treeName
uploadID    = $uploadID
projectName = itol_uploader
colorDefinitionFile = $color_file
EOF

	if ($this->tree) {
		print $FH <<EOF;
treeFile    = $tree_file
treeFormat  = newick
showInternalIDs = 1
EOF
	} elsif ($this->ncbi_taxonomy_id) {
		Smash::Utils::Taxonomy::init("ncbi");
		print $FH <<EOF;
ncbiFile    = $ncbi_file
ncbiFormat  = idsCollapsed
assignTaxonomy = 1
EOF
	} else {
		print $FH <<EOF;
treeFile    = $tree_file
treeFormat  = newick
showInternalIDs = 1
EOF
	}
	$this->{CONFIG_FH} = $FH;

	return $this;
}

sub get_ncbi_link {
	my $id = shift;
	return "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=$id";
}

sub get_rdp_link {
	my $id = shift;
	return "http://rdp.cme.msu.edu/genome/main.spr?currentRoot=$id&displayDepth=0&openNode=null&taxonomy=rdpHome&type=false&reponly=false";
}

sub print_popup_info_recursive {
	my $tree = shift;
	my $lineage_tree = shift;
	my $node = shift;
	my $FH   = shift || \*STDOUT;

	foreach my $child (values %{$node->children}) {
		print_popup_info_recursive($tree, $lineage_tree, $child, $FH);
	}

	my $id   = $node->id;
	my $name = $node->name;
	my $rank = $node->rank;
	my $newick_name = $node->newick_name;
	my $title = "Leaf information:";
	   $title = "Clade information:" if ($node->nchildren > 0);
	my @lineage = reverse $lineage_tree->get_full_lineage($id);
	my $lineage;
	my $html;
	my $node_link;
	my $node_text;
	if ($tree->type eq "functional") {
		shift @lineage;
		if ($tree->name =~ /eggNOG/i) {
			$node_text = $id;
			$lineage = join("; ", map {$lineage_tree->nodes->{$_}->name} @lineage);
		} elsif ($tree->name =~ /KEGG/i) {
		}
		$html =<<EOF;
<BOX>
<H1><B>$name</B></H1>
<BR>
<H2><B>$tree classification:</B></H2>
<TEXTFORMAT TABSTOPS="[28,35]"><B>Entry\\t:</B>\\t$node_text</TEXTFORMAT>
<BR>
<TEXTFORMAT TABSTOPS="[28,35]"><B>Level\\t:</B>\\t$rank</TEXTFORMAT>
<BR>
<B>Hierarchy:</B>
<BR>
&nbsp;&nbsp;$lineage
<BR>
</BOX>
EOF
	} else { # I assume it is phylogenetic if not "functional"
		if ($tree->name =~ /NCBI/i) {
			$node_link = get_ncbi_link($id);
			$node_text = "<A TARGET='_itol_taxonomy' HREF='$node_link'>$id</A>";
			$lineage = join("; ", map {"<A TARGET='_itol_taxonomy' HREF='".get_ncbi_link($_)."'>".$lineage_tree->nodes->{$_}->name."</A>"} @lineage);
		} elsif ($tree->name =~ /RDP/i) {
			$node_link = get_rdp_link($id);
			$node_text = "<A TARGET='_itol_taxonomy' HREF='$node_link'>$id</A>";
			$lineage = join("; ", map {"<A TARGET='_itol_taxonomy' HREF='".get_rdp_link($_)."'>".$lineage_tree->nodes->{$_}->name."</A>"} @lineage);
		} else {
			$node_text = $id;
			$lineage = join("; ", map {$lineage_tree->nodes->{$_}->name} @lineage);
		}
		$html =<<EOF;
<BOX>
<H1><B>$name</B></H1>
<BR>
<TEXTFORMAT TABSTOPS="[42]"><B>$tree Node</B>\\t: $node_text</TEXTFORMAT>
<BR>
<TEXTFORMAT TABSTOPS="[42]"><B>$tree Rank</B>\\t: $rank</TEXTFORMAT>
<BR>
<B>Lineage:</B>
<BR>
&nbsp;&nbsp;$lineage
<BR>
</BOX>
EOF
	}
	$html =~ s/\n//sg;
	print $FH "$newick_name\t$title\t$html\n";
}

sub finish {
	use File::Copy;
	my $this = shift;
	my $name = $this->name;
	my $Tree = $this->tree;

	my $keep_ids = 0;

	my $nodes_hash = {};

	# Track the nodes with data in them

	foreach my $i (1..$this->dataset_count) {
		my $data = $this->{DATALINE}->[$i];
		if ($this->{GRAPH_TYPE}->[$i] eq "connections") {
			apply_to_matrix_keys($data, sub {$nodes_hash->{$_[0]} = 1; $nodes_hash->{$_[1]} = 1;});
			if (1==0) { # Just in case apply_to_matrix doesnt work
				map {$nodes_hash->{$_} = 1} keys %$data;
				foreach my $node (keys %$data) {
					$nodes_hash->{$node} = 1;
					foreach my $node2 (keys %{$data->{$node}}) {
						$nodes_hash->{$node2} = 1;
					}
				}
			}
		} else {
			map {$nodes_hash->{$_} = 1} keys %$data;
		}
	}

	# If there are datasets, show only the nodes (leaves) corresponding to the datasets
	# If not, show the whole damn tree... Someone must really want to see the tree!

	my @nodes = keys %$nodes_hash;
	if ($this->dataset_count == 0) {
		@nodes = @{$Tree->get_leaf_ids()};
		map {$nodes_hash->{$_} = 1} @nodes;
	}

	# Pop-up Information
	# Since we show lineage information in pop-up, pop-up should be written
	# before we remove straight links. Therefore, it's split into two parts.

	my $popup = $this->{POPUP_INFO};
	print {$this->config_fh} "popupInfoFile = $name.popup\n";;
	open(POPUP, ">$name.popup") || die "Cannot open $name.popup: $!";
	print POPUP<<EOF;
CSS	box { font-family: Arial,Helvetica,sans-serif; font-size: 10px;}
CSS	h1 { font-size: 13px; }
CSS	h2 { font-size: 11px; }
CSS	a { color: #0000ff;}
CSS	a:hover { color: #0000ff;  text-decoration: underline;}
EOF


	# Let iTOL make a tree of NCBI tax_ids

	my $PrunedTree = $Tree;

	################
	# First section is when trees are given as objects:
	# 1. reading in phyloxml file
	# 2. created by parsing ncbi taxdump
	# 3. created by parsing RDP taxonomy file
	# For case 1, if ncbi_taxonomy_id() is set, then assign NCBI names to nodes with valid ids.
	################

	if ($Tree) {
		my $tree_file = $this->tree_file;

		if ($this->dataset_count != 0) {
			$PrunedTree = $Tree->get_pruned_tree(\@nodes);
			$PrunedTree->remove_straight_links();
		}

		foreach my $node (@nodes) {
			delete $nodes_hash->{$node} unless $PrunedTree->nodes->{$node};
		}
		@nodes = keys %$nodes_hash;

		if ($this->ncbi_taxonomy_id) {
			foreach my $node (@nodes) {
				if ($node =~ /^\d+$/) { # valid ncbi id?
					# Copy relevant info from NCBI tree
					map {$PrunedTree->nodes->{$node}->{$_} = $NCBITree->nodes->{$node}->{$_}} qw(NAME RANK);
					map {$Tree->nodes->{$node}->{$_} = $NCBITree->nodes->{$node}->{$_}} qw(NAME RANK);
				}
			}
		}

		if (!defined($popup)) {
			print_popup_info_recursive($PrunedTree, $Tree, $PrunedTree->rootlink, \*POPUP); # write the pop-up info now before removing straight links
		} else {
			foreach my $node (@nodes) {
				print POPUP $PrunedTree->nodes->{$node}->newick_name."\t".$popup->{$node}."\n" if $popup->{$node};
			}
		}

		open(LIST, ">$tree_file") || die "Cannot open $tree_file: $!";
		if ($keep_ids) {
			$PrunedTree->print_newick_id(\*LIST);
		} else {
			$PrunedTree->print_newick(\*LIST);
		}
		close(LIST);
	} elsif ($this->ncbi_taxonomy_id) {
		my $ncbi_file = $this->ncbi_file;
		open(LIST, ">$ncbi_file") || die "Cannot open $ncbi_file: $!";
		foreach my $node (@nodes) {
			print LIST "$node\n";
		}
		close(LIST);
	} else {
		my $tree_file = $this->tree_file;
		@nodes = sort {$a cmp $b} @nodes;
		open(LIST, ">$tree_file") || die "Cannot open $tree_file: $!";
		print LIST "(".join(',', @nodes).");\n";
		close(LIST);
		foreach my $node (@nodes) {
			print POPUP "$node\t".$popup->{$node}."\n" if $popup->{$node};
		}
	}
	close(POPUP);

	# Write dataset files
	# In a tree, if you supply multibar data for internal nodes, it is plotted
	# by iTOL as pie chart on the internal node. So I have to remove the entries
	# for the internal nodes only when it is a Tree object and not NCBI tax-ids.
	# We dont have to do this for connections, since it can go to internal nodes
	# as well.

	foreach my $i (1..$this->dataset_count) {
		my $datafile = $this->{DATASETFILE}->[$i];
		my $dataline = $this->{DATALINE}->[$i];
		open(DATASET, ">$datafile") || die "Cannot open file $datafile: $!";
		if ($this->{GRAPH_TYPE}->[$i] =~ /barplot|piechart|heatmap/) {
			print DATASET "LABELS\t".$this->{LABELLINE}->[$i]."\n";
		}
		print DATASET "COLORS\t".$this->{COLORLINE}->[$i]."\n" if $this->{COLORLINE}->[$i];
		if ($this->{GRAPH_TYPE}->[$i] eq "connections") {
			if ($PrunedTree) {
				apply_to_matrix_keys($dataline, sub {print DATASET join("\t", (map {$PrunedTree->nodes->{$_}->newick_name()} @_[0,1]), $_[2]); print DATASET "\n";});
			} else {
				apply_to_matrix_keys($dataline, sub {print DATASET join("\t", @_); print DATASET "\n";});
			}
		} else {
			if ($PrunedTree) {
				my @deleted_nodes = ();
				foreach my $node (@nodes) {
					if ($PrunedTree->nodes->{$node}->nchildren > 0) {
						push(@deleted_nodes, $node);
						delete $nodes_hash->{$node};
					}
				}
				@nodes = keys %$nodes_hash;
				if (@deleted_nodes) {
					print "WARNING: In dataset $i, data corresponding to the following internal nodes will not be displayed:\n";
					map {$PrunedTree->nodes->{$_}->print_node(1)} @deleted_nodes;
				}
			}
			foreach my $node (@nodes) {
				if ($dataline->{$node}) {
					if ($PrunedTree) {
						if ($keep_ids) {
							printf DATASET "%s\t%s\n", $PrunedTree->nodes->{$node}->id(), $dataline->{$node};
						} else {
							printf DATASET "%s\t%s\n", $PrunedTree->nodes->{$node}->newick_name(), $dataline->{$node};
						}
					} else {
						printf DATASET "%s\t%s\n", $node, $dataline->{$node};
					}
				}
			}
		}
		close(DATASET);
	}

	# COLOR_CLADES tells iTOL to color the leaf nodes according to clades.
	# By default, it will color the nodes according to the given tree.
	# However, there are cases where the given tree is built from scratch
	# and hence does not have useful internal node information. In this 
	# case, COLOR_CLADES=2 will transfer NCBI's taxonomy information
	# to the leaf nodes and color them. Note that the leaf node id's must
	# be NCBI taxonomy id's for this to work

	if ($this->{COLOR_CLADES}) {
		my $color_file = $this->color_file;
		my %Color = (
				# http://mudcu.be/sphere/ 39,92,81,251,220,162

				Archaea               => "#fabcac",
				Actinobacteria        => "#fcddae",
				Bacteroidetes         => "#aefdae",
				Proteobacteria        => "#affdce",
				Clostridiales         => "#aaaaf8",
				Erysipelotrichales    => "#daabf9", 

				# colorbrewer, 9, qualitative, pastel1

				Actinobacteria        => "#FBB4AE",
				Clostridia            => "#B3CDE3",
				Archaea               => "#CCEBC5",
				Erysipelotrichia      => "#DECBE4", 
				Selenomonadales       => "#FED9A6",
				Bacteroidetes         => "#FFFFCC",
				Bacilli               => "#E5D8BD",
				Proteobacteria        => "#FDDAEC",
				Bacteria              => "#F2F2F2",

				#Archaea               => "#d9fff2",
				#Actinobacteria        => "#ffc4dd",   # Lachnospiraceae
				#Bacteroidetes         => "#ffdec7",   # Erysipelotrichaceae
				#Proteobacteria        => "#e4ffa8",   # Ruminococcaceae
				#Clostridiales         => "#d9cfff",
				#Erysipelotrichales    => "#b5e0ff", 
				#Lactobacillales       => "#ff6721",
				#Gammaproteobacteria   => "#b5e0ff",
				#Betaproteobacteria   => "#e4ffa8",
				#Alphaproteobacteria  => "#ff6721",
				#Firmicutes            => "#d9cfff",   # Clostridiaceae
			# eukaryotes
				Eukaryota             => "#ff8373",
				Eurotiomycetes        => "#d9fff2",
				Dothideomycetes       => "#ffdec7",
				Sordariamycetes       => "#d9cfff",
				);

		open(COLORS, ">$color_file") || die "Cannot open $color_file: $!";

		my $tree = $Tree;
		if ($this->{COLOR_CLADES} == 2) {
			$tree = $NCBITree;
		} elsif ($this->ncbi_taxonomy_id) {
			$tree = $NCBITree;
		}
		foreach my $n (@nodes) {
			my $node = $n;
			if ($tree->nodes->{"NCBI.$node"}) {
				$node = "NCBI.$node";
			} elsif (!$tree->nodes->{$node}) {
				next;
			}
			if (my @lineage = $tree->get_full_lineage($node)) {
				foreach my $ancestor (map {$tree->nodes->{$_}->name} reverse @lineage) {
					if (my $color = $Color{$ancestor}) {
						if ($this->{COLOR_CLADES} == 2) {
							print COLORS $tree->nodes->{$node}->newick_name();
						} elsif ($this->ncbi_taxonomy_id || $keep_ids) {
							print COLORS $node;
						} else {
							print COLORS $tree->nodes->{$node}->newick_name();
						}
						print COLORS "\trange\t$color\t$ancestor\n";
					}
				} 
			}
		}
		close(COLORS);
	}

	close($this->config_fh);
	$this->upload();
	my $download = $this->download("svg");
	$this->fix_svg_file($download) if $this->post_fix && $this->display_mode eq "normal";
}

sub fix_svg_file {
	my $this     = shift;
	my $download = shift;

	open(SVG, "<$download") || die "Cannot open $download: $!";
	my $svg = "";
	while (<SVG>) {$svg .= $_};
	close(SVG);

	# remove data labels since we will add it later

	$svg =~ s/<g [^>]*id=\"Rdataset\d+\">.*?<\/g>//sgi;
	$svg =~ s/<g [^>]*id=\"Tdataset\d+\">.*?<\/g>//sgi;
	$svg =~ s/<g id=\'legend\'.*?<\/g>//si;

	# move the color ranges legend up

	$svg =~ s/<g id=\'ranges_legend\' transform=\"translate\(([-\d\.]+),[-\d\.]+\)\">/<g id=\'ranges_legend\' transform=\"translate\(@{[1.5*$1]},0) scale(1.5)\">/si; # scale(2), translate(2*$1,0)

	# add the sample labels above the plot

	my @labels = @{$this->labels};
	my $picture_size = $this->picture_size;
	foreach my $i (1..$this->dataset_count) {
		my $type = $this->{GRAPH_TYPE}->[$i];

		my $scale_width   = 10;
		my $font_size     = 30;

		if ($type eq "barplot") {
			# calculate original graph width before iTOL resized it

			my $original_graph_width = 0;
			foreach my $label (@labels) {
				$original_graph_width += $this->{FINAL_GRAPH_WIDTH}->[$i]->{$label};
			}

			# position the new label on top after getting the resized graph widths

			# since we rotate the text, lets use half the space at most, so they are not too close.
			# that's 10*@labels*$font_size in symmetric fonts
			# @labels*2*$font_size will share $picture_size space

			$font_size = sprintf("%.4f", $picture_size/scalar(@labels)/4);
			if ($font_size > 30) {
				$font_size = 30;
			}

			# find out where the graph starts in the svg

			#my ($graph) = $svg =~ /<g id="dataset0*${i}" display="block"><line\s+stroke="black" stroke-width="[\d\.]+" stroke-dasharray="[\d\.]+,[\d\.]+" x1="[\d\.]+" y1="[\d\.]+" x2="[\d\.]+" y2="[\d\.]+" style="stroke-width: 1; opacity: [\d\.]+" \/>(.*?)<\/g>/si;
			my ($graph) = $svg =~ /<g id="dataset0*${i}" display="block">(.*?)<\/g>/si;
			my ($graph_start_x, $graph_end_y) = (INT_MAX, INT_MIN);
			my @rects = $graph =~ /(<rect x=.*? y=.*? height=.*? width=.*? style=.*? \/>)/gsi;
			foreach my $rect (@rects) {
				my ($x, $y, $height) = $rect =~ /<rect x=['"]([\d\.\-]+)['"] y=['"]([\d\.\-]+)['"] height=['"]([\d\.\-]+)['"] width=.*? style=.*? \/>/si;
				if ($x < $graph_start_x) {
					$graph_start_x = $x;
				}
				if ($y+$height > $graph_end_y) {
					$graph_end_y = $y+$height;
				}
			}

			# make the panel

			my $longest_label = (sort {length($b) <=> length($a)} @labels)[0];
			my $panel_width   = 0.75*length($longest_label)*$font_size;

			($graph_start_x, $panel_width) = map {sprintf("%.8f", $_)} ($graph_start_x, $panel_width);

			# we make the panel by stacking the text and rectangles one on top of another and then rotate them
				# calculate the proportion of this graph's width in the original graph's width and normalize by iTOL picture size

			my $panel_bg = "";
			my $prev_color = "dummy";
			my $panel_color_start = 0;
			my $cumulative_graph_width = 0;
			my $buffer = "";
			foreach my $label (@labels) {
				my $color       = $this->get_color_for_label($label);
				my $graph_width = sprintf("%.8f", $picture_size * $this->{FINAL_GRAPH_WIDTH}->[$i]->{$label} / $original_graph_width);

				if ($color ne $prev_color) {
					$buffer .= "<rect x=\"0\" y=\"$panel_color_start\" width=\"$panel_width\" height=\"@{[$cumulative_graph_width-$panel_color_start]}\" style=\"fill:$prev_color;fill-opacity:0.2;\" \/>";
					$panel_color_start = $cumulative_graph_width;
				}
				$cumulative_graph_width += $graph_width;
				$prev_color = $color;
			}
			$buffer .= "<rect x=\"0\" y=\"$panel_color_start\" width=\"$panel_width\" height=\"@{[$cumulative_graph_width-$panel_color_start]}\" style=\"fill:$prev_color;fill-opacity:0.5;\" \/>";
			$panel_bg  = "<g id=\"scale_legend\" transform=\"translate(@{[$graph_start_x-20]}, -$scale_width)\"><text fill=\"#000000\" style=\"font-size:@{[0.75*$font_size]}px; text-anchor:end;\">Scale</text></g>";
			$panel_bg .= "<g id=\"panel_background_$i\" transform=\"translate($graph_start_x, @{[-5*$scale_width]}) rotate(270 0 0)\">".$buffer."</g>\n";

			# add label

			my $panel_text = "";
			$cumulative_graph_width = 0;
			foreach my $label (@labels) {
				my $color       = $this->get_color_for_label($label);
				my $graph_width = sprintf("%.8f", $picture_size * $this->{FINAL_GRAPH_WIDTH}->[$i]->{$label} / $original_graph_width);
				$panel_text .= "<text fill=\"$color\" x=\"-2\" y=\"@{[-$cumulative_graph_width-2]}\" style=\"font-size:${font_size}px; font-weight:bold; text-anchor:end;\">$label<\/text>";

				$cumulative_graph_width += $graph_width;
			}
			$panel_text = "<g id=\"panel_text_$i\" transform=\"translate($graph_start_x, @{[-5*$scale_width]}) rotate(90 0 0)\">".$panel_text."</g>\n";

			# add scale

			my $panel_scale = "";
			$cumulative_graph_width = 0;
			foreach my $label (@labels) {
				my $color       = $this->get_color_for_label($label);
				my $graph_width = sprintf("%.8f", $picture_size * $this->{FINAL_GRAPH_WIDTH}->[$i]->{$label} / $original_graph_width);
				my $scale       = 0.80*$this->{FINAL_GRAPH_WIDTH}->[$i]->{$label};
				SCALE:foreach my $try (100000, 10000, 1000, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001) {
					if ($try*int($scale/$try) > 0) {
						$scale = $try*int($scale/$try);
						last SCALE;
					}
				}
				my $scale_px = sprintf("%.8f", $picture_size * $scale / $original_graph_width);

				# rescale for display, since FINAL_GRAPH_WIDTH where $scale comes from has been multiplied by SCALE in add_dataset()

				$scale /= $this->{SCALE}->[$i];

				# make a scale bar with $scale_width width and value centered on top of it.

				$panel_scale .= "<rect x=\"@{[-2*$scale_width]}\" y=\"@{[-$cumulative_graph_width-$scale_px]}\" width=\"$scale_width\" height=\"$scale_px\" style=\"fill:$color;fill-opacity:0.9;\" \/>";
				my $scale_text = $scale;
				if (defined($this->{UNIT}->[$i])) {
					if ($this->{UNIT}->[$i] eq "%") {
						$scale_text = (100*$scale)."%";
					} else {
						$scale_text .= $this->{UNIT}->[$i];
					}
				}
				$panel_scale .= "<g transform=\"translate(@{[-2*$scale_width-0.2*$scale_width]}, -$cumulative_graph_width) rotate(270 0 0)\"><text fill=\"$color\" x=\"0\" y=\"0\" style=\"font-size:@{[0.75*$font_size]}px; font-weight:bold;\">$scale_text<\/text></g>";

				$cumulative_graph_width += $graph_width;
			}
			$panel_scale = "<g id=\"panel_scale_$i\" transform=\"translate($graph_start_x, 0) rotate(90 0 0)\">".$panel_scale."</g>\n";

			$svg =~ s/(<g id="dataset0*${i}" display="block">)/${panel_bg}${panel_text}${panel_scale}$1/si;

		} elsif ($type eq "boxplot") {

			my $font_size = 30;
			my $min = nearest(0.1, $this->{MIN_VALUE}->[$i]*$this->{SCALE}->[$i]);
			my $max = nearest(0.1, $this->{MAX_VALUE}->[$i]*$this->{SCALE}->[$i]);
			# find out where the graph starts in the svg

			my ($graph) = $svg =~ /<g id="dataset0*${i}" display="block"><line\s+stroke="black" stroke-width="[\d\.]+" stroke-dasharray="[\d\.]+,[\d\.]+" x1="[\d\.]+" y1="[\d\.]+" x2="[\d\.]+" y2="[\d\.]+" style="stroke-width: 1; opacity: [\d\.]+" \/>(.*?)<\/g>/si;
			my ($graph_start_x, $graph_end_y) = (INT_MAX, INT_MIN);
			my @lines = $graph =~ /(<line style=.*? x1=.*? y1=.*? x2=.*? y2=.*? \/>)/gsi;
			foreach my $line (@lines) {
				my ($x1, $x2, $y2) = $line =~ /<line style=.*? x1=['"]([\d\.\-]+)['"] y1=.*? x2=['"]([\d\.\-]+)['"] y2=['"]([\d\.\-]+)['"] \/>/si;
				if ($x1 < $graph_start_x) {
					$graph_start_x = $x1;
				}
				if ($y2 > $graph_end_y) {
					$graph_end_y = $y2;
				}
			}
			my $panel_bg  = "<g id=\"scale_legend\" transform=\"translate(@{[$graph_start_x-20]}, -$scale_width)\"><text fill=\"#000000\" style=\"font-size:@{[0.75*$font_size]}px; text-anchor:end;\">Scale</text></g>\n";

			my $half_scale = int(0.5+$scale_width/2);
			my $quarter_scale = int(0.5+$scale_width/4);
			my $panel_scale  = "<line x1=\"0\" y1=\"0\" x2=\"0\" y2=\"-$picture_size\" style=\"stroke: #000000; stroke-width: $quarter_scale;\" />\n";
			   $panel_scale .= "<line x1=\"-$half_scale\" y1=\"0\" x2=\"$half_scale\" y2=\"0\" style=\"stroke: #000000; stroke-width: 2;\" />\n";
			   $panel_scale .= "<line x1=\"-$half_scale\" y1=\"-$picture_size\" x2=\"$half_scale\" y2=\"-$picture_size\" style=\"stroke: #000000; stroke-width: 2;\" />\n";
			   $panel_scale .= "<text fill=\"#000000\" x=\"@{[-6*$scale_width]}\" y=\"0\" style=\"font-size:@{[0.75*$font_size]}px; font-weight:bold;\">$min<\/text>\n";
			   $panel_scale .= "<text fill=\"#000000\" x=\"@{[-6*$scale_width]}\" y=\"-$picture_size\" style=\"font-size:@{[0.75*$font_size]}px; font-weight:bold;\">$max<\/text>\n";

			$panel_scale = "<g id=\"panel_scale_$i\" transform=\"translate($graph_start_x, 0) rotate(90 0 0)\">\n$panel_scale</g>\n";
			$svg =~ s/(<g id="dataset0*${i}" display="block">)/${panel_bg}${panel_scale}$1/si;
		}
	}

	copy($download, "$download.original");
	open(SVG, ">$download") || die "Cannot open $download: $!";
	print SVG $svg;
	close(SVG);
}

sub upload {
	my $this        = shift;
	my $config_file = $this->config_file;
	my $tree_id;
	my $success  = 0;
	my $smash    = Smash::Core->new()->init();
	my $uploader = $smash->get_script_path("iTOL_uploader.pl");
	my $command  = "$SMASH_PERL $uploader --config $config_file";
	print "Running: $command\n";
	open(RESPONSE, "$command |") || die "Cannot open pipe to iTOL_uploader: $!";
	while (<RESPONSE>) {
		chomp();
		if (m/^Upload failed/) {
			print "ITOL: $_\n";
			while (<RESPONSE>) {
				print "ITOL: $_";
			}
			return undef;
		} elsif (m/^Upload successful\./) {
			$success = 1;
		} elsif (m/^(\d+)$/ && $success) {
			$tree_id = $1;
		}
		print "ITOL: $_\n";
	}
	$this->{TREE_ID} = $tree_id;
	print "Upload successful. You can view this tree on iTOL using this URL:\n";
	print "http://itol.embl.de/external.cgi?tree=".$this->tree_id."\n";
	return $tree_id;
}

sub download {
	my $this        = shift;
	my $format      = shift;
	my $name        = $this->name;
	my $tree_id     = $this->tree_id;
	my $display_mode = $this->display_mode;
	my $font_size   = $this->font_size;
	my $line_width  = $this->line_width;
	my $output_file = "$name.downloaded.$format";
	my $smash       = Smash::Core->new()->init();
	my $downloader  = $smash->get_script_path("iTOL_downloader.pl");

	my $datasets    = join(",", map {"dataset$_"} 1..$this->dataset_count);
	my $command = "$SMASH_PERL $downloader --tree $tree_id --format $format --outputFile $output_file --lineWidth $line_width --fontSize $font_size --displayMode $display_mode --colorBranches 1 --alignLabels 1 --scaleFactor 1 --showBS 1 --BSdisplayValue M0 --BSdisplayType text";
	$command .= " --datasetList $datasets" if $this->dataset_count;
	print "Running: $command\n";
	system($command);
	return $output_file;
}

sub merge_taxonomy_ids {
	my $this    = shift;
	my $dataset = shift;
	my $labels  = $this->labels;
	my $merge_level = $this->merge_level;
	my $tree    = $this->tree;

	my %MergedDataset;

	if ($tree) {
		my @ids = grep {$_ >= 0} keys %$dataset;
		foreach my $old_id (@ids) {
			my $new_id = get_taxonomic_rank($tree, $old_id, $merge_level);
			foreach my $label (@$labels) {
				$MergedDataset{$new_id}{$label} += ($dataset->{$old_id}{$label}||0);
			}
		}
		return \%MergedDataset;
	} elsif ($this->ncbi_taxonomy_id) {
		Smash::Utils::Taxonomy::init("prokaryotes");
		my @ids = grep {$_ >= 0} keys %$dataset;
		foreach my $old_id (@ids) {
			my $new_id = get_ncbi_taxonomic_rank($old_id, $merge_level);
			foreach my $label (@$labels) {
				$MergedDataset{$new_id}{$label} += ($dataset->{$old_id}{$label}||0);
			}
		}
		return \%MergedDataset;
	}
}

####
# print the colors
# if there's a color, that will be used. default will be red.
# colors will be matched by looking for a prefix of the graph label that matches with one of the available color labels
# e.g., AM-AD-1 will match A, AM, AM-, AM-A, AM-AD, AM-AD-, AM-AD-1. If more than one match is present, the longest will be used!
####

sub get_color_for_label {
	my $this   = shift;
	my $label  = shift;
	my $colors = $this->dataset_colors;
	my @colored_prefixes = sort {length($b) <=> length($a)} keys %$colors;
	foreach my $available_color (@colored_prefixes) {
		if ($available_color eq substr($label, 0, length($available_color))) {
			return $colors->{$available_color};
		}
	}
	return "#8E191B";
}

# REMEMBER TO USE THIS LINE LATER
#print DATASET ($tree)?$tree->nodes->{$node}->newick_name():$node;

sub add_dataset {
	# All keys %Label are sorted in reverse so that $opt_tag.dataset and $opt_tag.piechart match in their order

	my $this           = shift;
	my %params         = @_;
	my $dataset_name   = $params{NAME};
	my $in_dataset_ref = $params{DATA};
	my $scale          = $params{SCALE} || 1;
	my $unit           = $params{UNIT};
	my $graph_type     = $params{GRAPH_TYPE};
	my $keep_top_n     = $params{KEEP_TOP_N};
	my $popup_info     = $params{POPUP_INFO};
	my $boxplot_raw    = $params{BOXPLOT_RAW};

	my %dataset_hash   = %$in_dataset_ref;
	my $dataset        = \%dataset_hash;

	my $new_dataset    = {};

	my $name           = $this->name;
	my @labels         = @{$this->labels};
	my $dataset_colors = $this->dataset_colors;
	my $log_scale      = $this->log_scale;
	my $precision      = $this->precision;
	my $picture_size   = $this->picture_size;
	my $graph_width    = $this->graph_width;
	my $spacer         = $graph_width / 10;
	my $display_threshold = $this->display_threshold;
	my $dataset_index  = $this->dataset_count+1;
	my $dataset_file   = "$name.$dataset_name.dataset";

	if ($popup_info) {
		map {$this->{POPUP_INFO}->{$_} = $popup_info->{$_}} keys %$popup_info;
	}

	# print iTOL config files

	print {$this->config_fh} <<EOF;
dataset${dataset_index}Label = $dataset_name
dataset${dataset_index}Separator = tab
EOF


	#delete $dataset->{-1};

	# merge tax_ids if necessary

	if (defined($this->merge_level)) {
		$dataset = $this->merge_taxonomy_ids($dataset);
	}

	# Keep only the top N per sample
	# WARNING: It keeps only topN in the graph... it should keep the topN, but also others that are topN in other samples
	# FIXED: 27.10.2009

	if ($keep_top_n) {
		my %Keep;
		foreach my $label (@labels) {
			my @nodes = sort {$dataset->{$b}{$label} <=> $dataset->{$a}{$label}} grep {defined($dataset->{$_}{$label})} keys %$dataset;
			if (scalar(@nodes) > $keep_top_n) {
				map {$Keep{$_} = 1} @nodes[0..($keep_top_n-1)];
			} else {
				map {$Keep{$_} = 1} @nodes;
			}
		}
		foreach my $node (keys %$dataset) {
			if (!defined($Keep{$node})) {
				delete $dataset->{$node};
			}
		}
	}

	########################
	# Track the max values for each graph to set the width of each graph
	########################

	my %GRAPH_WIDTH = map {$_ => $graph_width} @labels;

	# Throw away very low numbers
	# Find max value of a node in all datasets and delete if it is lower than threshold

	foreach my $node (keys %$dataset) {
		next if ($this->tree && $node eq "-1"); # Dont include Unknown (-1) in counting MAX
		my $max_value_in_sample = 0;
		foreach my $label (@labels) {
			my $value = $dataset->{$node}->{$label} || 0;
			if ($value >= 0) { 
				if ($value > $GRAPH_WIDTH{$label}) {$GRAPH_WIDTH{$label} = $value;};
				if ($value > $max_value_in_sample) {$max_value_in_sample = $value;};
			}
		}
		if (!$keep_top_n) {
			if ($max_value_in_sample < $display_threshold) { delete $dataset->{$node};}
		}
	}

	# because of the zero_fill_matrix() call, no need to check for presence of value in cell from now on!

	if ($graph_type ne "connections") {
		zero_strip_matrix($dataset);
		zero_fill_matrix($dataset);
	}

	if ($this->{RE_SORT}) { # re-sort the abundances and forget the actual features. Just keep the ordinal positions
		my @tnodes = keys %$dataset;
		foreach my $label (@labels) {
			my @tmp_nodes = sort {$dataset->{$b}->{$label} <=> $dataset->{$a}->{$label}} grep {defined($dataset->{$_}->{$label})} @tnodes;
			map {$new_dataset->{1+$_}->{$label} = $dataset->{$tmp_nodes[$_]}->{$label}} (0..$#tmp_nodes);
		}
		$dataset = $new_dataset;
	}

	# reassess nodes to display on tree after tree-pruning and editing

	my @nodes = sort keys %$dataset;

	# set dataset filename

	$this->{DATASETFILE}->[$dataset_index] = $dataset_file;

	if ($graph_type eq "barplot") {

		########################
		# set the labels and colors
		########################

		$this->{LABELLINE}->[$dataset_index] = join("\t", map {"$_\t${_}_blank"} @labels);
		$this->{COLORLINE}->[$dataset_index] = join("\t", map {$this->get_color_for_label($_)."\t#FFFFFF"} @labels);

		########################
		# set the dataset lines
		########################

		# round up the maxvalue for display using precision and add the spacer as well

		foreach my $label (@labels) {
			$GRAPH_WIDTH{$label} = $precision*POSIX::ceil(($GRAPH_WIDTH{$label}+$spacer)*$scale/$precision);
		}

		# print the data. for each value to be a bar chart, it needs to be padded upto the graphwidth.

		foreach my $node (@nodes) {
			my @fields = map {my $value = $dataset->{$node}->{$_}*$scale; sprintf("%.4f\t%.4f", $value, $GRAPH_WIDTH{$_}-$value);} @labels;
			$this->{DATALINE}->[$dataset_index]->{$node} = join("\t", @fields);
		}

		print {$this->config_fh} <<EOF;
dataset${dataset_index}Color = \\#FF0000
dataset${dataset_index}File = $dataset_file
dataset${dataset_index}Type = multibar
dataset${dataset_index}BarSizeMax = $picture_size
EOF
		$this->{FINAL_GRAPH_WIDTH}->[$dataset_index] = \%GRAPH_WIDTH;
	} elsif ($graph_type eq "simplebar") {

		warn "NOTE: simplebar can only display one dataset!\n";

		########################
		# set the labels and colors
		########################

		$this->{LABELLINE}->[$dataset_index] = $labels[0];
		my $color = $this->get_color_for_label($labels[0]);

		########################
		# set the dataset lines
		########################

		# print the data. for each value to be a bar chart, it needs to be padded upto the graphwidth.

		foreach my $node (@nodes) {
			$this->{DATALINE}->[$dataset_index]->{$node} = $dataset->{$node}->{$labels[0]}*$scale;
		}

		print {$this->config_fh} <<EOF;
dataset${dataset_index}Color = \\$color
dataset${dataset_index}File = $dataset_file
dataset${dataset_index}Type = simplebar
dataset${dataset_index}BarSizeMax = $picture_size
EOF
		$this->{FINAL_GRAPH_WIDTH}->[$dataset_index] = \%GRAPH_WIDTH;

	} elsif ($graph_type eq "heatmap") {
		########################
		########################

		# set the labels and colors

		$this->{LABELLINE}->[$dataset_index] = join("\t", @labels);

		my ($min, $mid, $max);
		my ($min_color, $mid_color, $max_color);

		# set the dataset

		$min = 10000000;
		$max = 0;
		my @values = ();
		foreach my $node (@nodes) {
			my @fields = map {my $value = $dataset->{$node}->{$_}*$scale; push(@values, $value); sprintf "%.4f", $value} @labels;
			$this->{DATALINE}->[$dataset_index]->{$node} = join("\t", @fields);
		}

		# keep track of min and max

		@values = sort {$a <=> $b} @values;
		($min, $max) = (0, $values[$#values]);
		$mid = get_median(\@values);
		my $l = 10**int(log($max)/log(10));
		$max  = $l*POSIX::ceil($max/$l);

		# set min,mid,max colors
		# if the colors are passed, use them

		($min_color, $mid_color, $max_color) = ("#0000FF", "#FFFFFF", "#FF0000");
		if ($params{HEATMAP_COLORS}) {
			my $colors = $params{HEATMAP_COLORS};
			($min, $mid, $max) = sort {$a <=> $b} keys %$colors;
			($min_color, $mid_color, $max_color) = map {$colors->{$_}} ($min, $mid, $max);
		}

		print {$this->config_fh} <<EOF;
dataset${dataset_index}Color = \\#FF0000
dataset${dataset_index}File = $dataset_file
dataset${dataset_index}Type = heatmap
dataset${dataset_index}HeatmapBoxWidth = $graph_width
dataset${dataset_index}MinPointValue   = $min
dataset${dataset_index}MidPointValue   = $mid
dataset${dataset_index}MaxPointValue   = $max
dataset${dataset_index}MinPointColor   = \\$min_color
dataset${dataset_index}MidPointColor   = \\$mid_color
dataset${dataset_index}MaxPointColor   = \\$max_color
EOF
	} elsif ($graph_type eq "piechart") {

		# set the labels and colors

		$this->{LABELLINE}->[$dataset_index] = join("\t", @labels);
		$this->{COLORLINE}->[$dataset_index] = join("\t", map {$this->get_color_for_label($_)} @labels);

		foreach my $node (@nodes) {
			my $total = 0;
			foreach my $label (@labels) {
				my $value = $dataset->{$node}->{$label}*$scale;
				$total += $value;
			}
			my @fields = map {sprintf "%d", $dataset->{$node}->{$_}*$scale} @labels;
			my $radius = sprintf("R%d", sqrt($total));
			$this->{DATALINE}->[$dataset_index]->{$node} = join("\t", $radius, @fields);
		}
		print {$this->config_fh} <<EOF;
dataset${dataset_index}Color = \\#FF0000
dataset${dataset_index}File = $dataset_file
dataset${dataset_index}Type = multibar
dataset${dataset_index}PieRadiusMin = 7
dataset${dataset_index}PieRadiusMax = 13
EOF
	} elsif ($graph_type eq "colorstrip") {

		# set the labels and colors

		#$this->{LABELLINE}->[$dataset_index] = join("\t", @labels);
		#$this->{COLORLINE}->[$dataset_index] = join("\t", map {$this->get_color_for_label($_)} @labels);

		foreach my $node (@nodes) {
			foreach my $label (@labels) {
				my $value = $dataset->{$node}->{$label};
				if ($value) {
					$this->{DATALINE}->[$dataset_index]->{$node} = "#318246";
				}
			}
		}
		print {$this->config_fh} <<EOF;
dataset${dataset_index}Color = \\#FF0000
dataset${dataset_index}File = $dataset_file
dataset${dataset_index}Type = colorstrip
dataset${dataset_index}BranchColoringType = branch
EOF
	} elsif ($graph_type eq "boxplot") {

		foreach my $node (@nodes) {
			my @fields = ();
			foreach my $label (@labels) {
				my $val = $dataset->{$node}->{$label}*$scale;
				if ($val > 0) {
					push(@fields, sprintf "%.4f", $val);
				}
			}
			$this->{DATALINE}->[$dataset_index]->{$node} = join("\t", @fields);
		}

		print {$this->config_fh} <<EOF;
dataset${dataset_index}Color = \\#FF0000
dataset${dataset_index}File = $dataset_file
dataset${dataset_index}Type = boxplot
dataset${dataset_index}BoxplotSize = $picture_size
EOF

		if ($boxplot_raw) {
			print {$this->config_fh} <<EOF;
dataset${dataset_index}Color = \\#FF0000
dataset${dataset_index}BoxplotRawData = 1
dataset${dataset_index}BoxplotLowerPercentile = 25
dataset${dataset_index}BoxplotUpperPercentile = 75
dataset${dataset_index}BoxplotStepConstant = 1.5
dataset${dataset_index}BoxplotCalculateExtremes = 1
EOF
		} else {
			print {$this->config_fh} <<EOF;
dataset${dataset_index}BoxplotRawData = 0
EOF
		}
		$this->{MIN_VALUE}->[$dataset_index] = min_hash($dataset);
		$this->{MAX_VALUE}->[$dataset_index] = max_hash($dataset);
	} elsif ($graph_type eq "connections") {
		my $Alpha = $params{ALPHA} || {};
		my $Color = $params{COLOR} || {};
		foreach my $node1 (@nodes) {
			foreach my $node2 (keys %{$dataset->{$node1}}) {
				my $value = $dataset->{$node1}->{$node2};
				my ($text, $type, $thickness) = split(SEPARATOR, $value);
				$type ||= "";
				$thickness ||= 25;
				my $alpha = $Alpha->{$type} || 200 ;
				my $color = $Color->{$type} || "#ff0000";
				$this->{DATALINE}->[$dataset_index]->{$node1}->{$node2} = join("\t", $thickness, $alpha, $color, $text);
			}
		}
		print {$this->config_fh} <<EOF;
dataset${dataset_index}Color = \\#FF0000
dataset${dataset_index}File = $dataset_file
dataset1Type = connections
EOF
	}

	$this->{GRAPH_TYPE}->[$dataset_index] = $graph_type;
	$this->{SCALE}->[$dataset_index] = $scale;
	$this->{UNIT}->[$dataset_index] = $params{UNIT};
	$this->{DATASET_COUNT}++;
}

1;
