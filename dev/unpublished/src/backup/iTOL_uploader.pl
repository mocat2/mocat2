#!/usr/bin/perl
use strict;
use HTTP::Request::Common;
use LWP::UserAgent;
use Pod::Usage;
use Getopt::Regex;
use Config::File::Simple;


=pod

=head1 NAME

iTOL_uploader  -  batch upload trees and data to iTOL

=head1 SYNOPSIS

B<iTOL_uploader.pl> I<options>

=head1 DESCRIPTION

Use B<iTOL_uploader.pl> to upload your phylogenetic trees and datasets to iTOL through its batch upload interface. Your trees and datasets must be in plain text format. If you want to upload into your iTOL account, you must first enable the batch upload through a link in 'Account settings' on your personal account page. This will create a unique 'Upload ID' which you must use when uploading (through the 'uploadID' option).


=head1 GENERAL OPTIONS

Separate options and their values with a space.

=over 4

=item B<--help>

display this message

=item B<--config>

read upload options from a file (format described below); example file available at http://itol.embl.de/help/upload.cfg

=item B<--treeFile>

file containing the tree

=item B<--treeFormat>

tree format (newick, nexus or phyloxml)

=item B<--treeName>

tree name

=item B<--treeDescription>

longer description of the tree, used only if uploading into an iTOL account

=item B<--uploadID>

Your upload ID, which is generated when you enable batch uploading in your account (in the B<'Account settings'>).

=item B<--projectName>

Required if B<uploadID> is specified. Project name is case sensitive, and should be unique in your account. We recommend creating a separate project for batch uploaded trees.

=item B<--showInternalIDs>

displays internal node IDs in popups when you mouseover them

=item B<--fontStylesFile>

defines font styles for leaf labels in exported trees

=item B<--internalNamesFile>

defines custom names for internal nodes in the tree

=item B<--HGTFile>

defines horizontal gene transfers

=item B<--colorDefinitionFile>

defines color ranges, clade and leaf label colors

=item B<--preCollapsedFile>

defines which clades should be collapsed by default; useful for large trees

=item B<--branchLabelsFile>

defines text labels for nodes in the tree

=item B<--popupInfoFile>

provides custom HTML info for mouse-over popups

=item B<--assignTaxonomy>

if you tree uses NCBI tax IDs, the names for leaves and internal nodes will be automatically assigned during upload

=item B<--midpointRoot>

tree will be midpoint re-rooted during upload

=item B<--ncbiFile>

instead of uploading a tree, iTOL can automatically generate one from a file containing a list of NCBI tax IDs. NCBI taxonomy will be pruned based on your IDs and a Newick tree generated.

=item B<--ncbiFormat>

format of the tree generated using NCBI tax IDs:

=over 5

=item B<'namesFull'>

generated tree will contain scientific names and internal nodes will not be collapsed

=item B<'namesCollapsed'>

scientific names will be used, and internal nodes with only one child removed

=item B<'idsFull'>

tree will contain NCBI taxonomy IDs and internal nodes will not be collapsed

=item B<'idsCollapsed'>

NCBI taxonomy IDs will be used, and internal nodes with only one child removed

=back

=back

=head1 DATASET UPLOAD OPTIONS 

Up to 10 datasets can be uploaded, replace 'dataset1' as needed for others (ie. dataset2 for the second dataset and so on). Dataset file, label, separator and type are required, other parameters are optional.

=over 4

=item B<--dataset1File>

file containing the dataset

=item B<--dataset1Label>

label for the dataset, used in the legends

=item B<--dataset1Separator>

field separator used in the dataset file: 'space', 'tab' or 'comma'

=item B<--dataset1Type>

dataset type

=over 6

=item B<'binary'>

Binary data

=item B<'simplebar'>

Simple bar chart

=item B<'multibar'>

Multi-value bar chart or Pie chart

=item B<'gradient'>

Color gradient

=item B<'colorstrip'>

Color strips

=item B<'time'>

Time series

=item B<'domains'>

Protein domain architecture

=item B<'heatmap'>

Heatmap

=item B<'boxplot'>

Boxplot

=item B<'connections'>

Connections

=item B<'circles'>

Circles

=back

=item B<--dataset1Color>

used in the legends and for datasets where the color is not specified in the dataset file

=item B<--dataset1StripWidth>

used in B<'gradient'> and B<'colorstrip'>; strip will have the specified width in pixels (valid range: 10 - 100)

=item B<--dataset1PreventOverlap>

used in B<'gradient'> and B<'colorstrip'>; if set to 1, all other dataset types will start after the end of the last strip dataset

=item B<--dataset1MultiBarAlign>

used in B<'multibar'>; if set to 1, each field's start will be aligned

=item B<--dataset1PieTransparent>

used in B<'multibar'>; if transprent, you will be able to see overlapping pie charts, but the colors will be slightly changed

=item B<--dataset1PieRadiusMax>

used in B<'multibar'>; pie chart radii will be scaled to match the selected values (valid range 10-500)

=item B<--dataset1PieRadiusMin>

used in B<'multibar'>; pie chart radii will be scaled to match the selected values (valid range 10-500)

=item B<--dataset1BarSizeMax>

used in B<'multibar'>, B<'simplebar'> and B<'time'>; maximum value in the dataset will have this pixel size (valid range 1-5000)

=item B<--dataset1ProtSizeMax>

used in B<'domains'>; longest protein in the dataset will have this pixel size (valid range 1-5000)

=item B<--dataset1BoxplotSize>

used in B<'boxplot'>; largest value in the dataset will be at this pixel (valid range 1-5000)

=item B<--dataset1BoxplotRawData>

used in B<'boxplot'>; if you provide the raw data, iTOL will automatically calculate the values needed to display the boxplots

=item B<--dataset1BoxplotUpperPercentile>

used in B<'boxplot'>; values above the set percentile will be included in the box (valid range 0-99)

=item B<--dataset1BoxplotLowerPercentile>

used in B<'boxplot'>; values below the set percentile will be included in the box (valid range 1-100)

=item B<--dataset1BoxplotStepConstant>

used in B<'boxplot'>; this constant multiplied with the box size will define the related whisker size (valid value >0)

=item B<--dataset1BoxplotCalculateExtremes>

used in B<'boxplot'>; if set to 0, only the boxes and whiskers will be calculated 

=item B<--dataset1BranchColoringType>

used in B<'colorstrip'>; specifies the effect of the dataset on branch coloring:

=over 6

=item B<'none'>

no coloring, show only leaf boxes (Default)

=item B<'branch'>

color branches only, don't display boxes

=item B<'both'>

color branches and show boxes

=back

=item B<--dataset1HeatmapBoxWidth>

used in B<'heatmap'>; pixel width of an individual value's box in the heatmap. Total width for the dataset cannot exceed 5000 px.

=item B<--dataset1MinPointValue>

used in B<'heatmap'>; minimum value

=item B<--dataset1MidPointValue>

used in B<'heatmap'>; mid-point value

=item B<--dataset1MaxPointValue>

used in B<'heatmap'>; maximum value

=item B<--dataset1MinPointColor>

used in B<'heatmap'>; color of the minimum value

=item B<--dataset1MidPointColor>

used in B<'heatmap'>; color of the mid-point value

=item B<--dataset1MaxPointColor>

used in B<'heatmap'>; color of the maximum value

=item B<--dataset1CirclesMinR>

used in B<'circles'>; Maximum radius will depend on the available space between leaves, and minimum radius will be calculated based on the percentage specified here (valid range 0-100)

=back

=item B<--dataset1CirclesHorizontalGrid>

used in B<'circles'>; Display horizontal grid connecting the circles

=back

=item B<--dataset1CirclesVerticalGrid>

used in B<'circles'>; Display vertical grid connecting the circles

=back

=item B<--dataset1CirclesSpacing>

used in B<'circles'>; horizontal spacing between circles (in pixels, valid range 0-500)

=back

=item B<--dataset1CirclesFill>

used in B<'circles'>; Should the circles be filled or empty? Filled circles will be fully opaque by default, unless 'dataset1CirclesTransparent' is defined.

=back

=item B<--dataset1CirclesTransparent>

used in B<'circles'>; filled circles will be transparent; ignored if 'dataset1CirclesFill' is not defined

=back

=back

=head1 USING A CONFIGURATION FILE

If you are uploading several datasets, the command line becomes hard to handle and its length could reach your shell limit. Therefore, we provided an option of using a separate plain text config file to specify the upload parameters. To use a configuration file, simply use the B<--config> option to specify its location.

The configuration file should contain one parameter per line. Use the B<=> sign to assign values to the parameters. Please note that you should backslash the hash sign in color definitions. Here is an example configuration file:

 treeFile = tree_of_life.tree
 treeFormat = newick
 treeName = upload test
 uploadID = aNrrwThf6vrDtWdZeaA1Ya
 projectName = Batch test
 dataset1File = tree_of_life.heatmap200
 dataset1Label = heatmap_batch
 dataset1Separator = tab
 dataset1Type = heatmap
 dataset1Color = \#FF00FF
 dataset1HeatmapBoxWidth = 10
 dataset1MinPointValue = -1000
 dataset1MinPointColor = \#ff0000
 dataset1MidPointValue = 0
 dataset1MidPointColor = \#000000
 dataset1MaxPointValue = 1000
 dataset1MaxPointColor = \#0000ff
 dataset2File = tree_of_life.multi4
 dataset2Label = multibar_batch
 dataset2Separator = tab
 dataset2Type = multibar
 dataset2Color = \#0000FF
 dataset2BarSizeMax = 2000

=head1 SEE ALSO

 iTOL Home page : http://itol.embl.de
 iTOL Help      : http://itol.embl.de/help/help.shtml

=head1 AUTHOR

 Ivica Letunic <ivica@letunic.com>
 Contact me if you have any questions or comments.

=cut

my $upload_url = "http://itol.embl.de/batch_uploader.cgi";
my ($showHelp, $configFile, $treeFile, $treeFormat, $treeName, $treeDescription, $uploadID, $showInternalIDs, $projectName, $fontStylesFile, $internalNamesFile, $HGTFile, $colorDefinitionFile, $preCollapsedFile, $branchLabelsFile, $popupInfoFile, $assignTaxonomy, $midpointRoot, $ncbiFile, $ncbiFormat);
my %datasets;

Getopt::Regex::GetOptions(\@ARGV,
						  ['--?h(elp)?',\$showHelp, 0],
						  ['--?config',\$configFile, 1],
						  ['--?treeFile',\$treeFile, 1],
						  ['--?treeFormat',\$treeFormat, 1],
						  ['--?ncbiFile',\$ncbiFile, 1],
						  ['--?ncbiFormat',\$ncbiFormat, 1],
						  ['--?treeName',\$treeName, 1],
						  ['--?treeDescription',\$treeDescription, 1],
						  ['--?showInternalIDs',\$showInternalIDs, 1],
						  ['--?uploadID',\$uploadID, 1],
						  ['--?projectName',\$projectName, 1],
						  ['--?fontStylesFile',\$fontStylesFile, 1],
						  ['--?internalNamesFile',\$internalNamesFile, 1],
						  ['--?HGTFile',\$HGTFile, 1],
						  ['--?colorDefinitionFile',\$colorDefinitionFile, 1],
						  ['--?preCollapsedFile',\$preCollapsedFile, 1],
						  ['--?branchLabelsFile',\$branchLabelsFile, 1],
						  ['--?popupInfoFile',\$popupInfoFile, 1],
						  ['--?assignTaxonomy',\$assignTaxonomy, 1],
						  ['--?midpointRoot',\$midpointRoot, 1],
						  ['--?dataset(\d+)File', sub { $datasets{"d$1"}{"file"} = $_[0]; }, 1],
						  ['--?dataset(\d+)Label', sub { $datasets{"d$1"}{"label"} = $_[0]; }, 1],
						  ['--?dataset(\d+)Separator', sub { $datasets{"d$1"}{"separator"} = $_[0]; }, 1],					  
						  ['--?dataset(\d+)Type', sub { $datasets{"d$1"}{"type"} = $_[0]; }, 1],
						  ['--?dataset(\d+)Color', sub { $datasets{"d$1"}{"color"} = $_[0]; }, 1],
						  ['--?dataset(\d+)StripWidth', sub { $datasets{"d$1"}{"stripwid"} = $_[0]; }, 1],
						  ['--?dataset(\d+)PreventOverlap', sub { $datasets{"d$1"}{"preventoverlap"} = $_[0]; }, 1],
						  ['--?dataset(\d+)MultiBarAlign', sub { $datasets{"d$1"}{"multialign"} = $_[0]; }, 1],
						  ['--?dataset(\d+)PieTransparent', sub { $datasets{"d$1"}{"pietrans"} = $_[0]; }, 1],
						  ['--?dataset(\d+)PieRadiusMax', sub { $datasets{"d$1"}{"piemax"} = $_[0]; }, 1],
						  ['--?dataset(\d+)PieRadiusMin', sub { $datasets{"d$1"}{"piemin"} = $_[0]; }, 1],
						  ['--?dataset(\d+)BarSizeMax', sub { $datasets{"d$1"}{"barsize"} = $_[0]; }, 1],
						  ['--?dataset(\d+)ProtSizeMax', sub { $datasets{"d$1"}{"protsize"} = $_[0]; }, 1],
						  ['--?dataset(\d+)BoxplotSize', sub { $datasets{"d$1"}{"boxplotsize"} = $_[0]; }, 1],
						  ['--?dataset(\d+)BoxplotRawData', sub { $datasets{"d$1"}{"boxplotraw"} = $_[0]; }, 1],
						  ['--?dataset(\d+)BoxplotUpperPercentile', sub { $datasets{"d$1"}{"boxplotupper"} = $_[0]; }, 1],
						  ['--?dataset(\d+)BoxplotLowerPercentile', sub { $datasets{"d$1"}{"boxplotlower"} = $_[0]; }, 1],
						  ['--?dataset(\d+)BoxplotCalculateExtremes', sub { $datasets{"d$1"}{"boxplotextr"} = $_[0]; }, 1],
						  ['--?dataset(\d+)BoxplotStepConstant', sub { $datasets{"d$1"}{"boxplotstep"} = $_[0]; }, 1],
						  ['--?dataset(\d+)BranchColoringType', sub { $datasets{"d$1"}{"coloringtype"} = $_[0]; }, 1],
						  ['--?dataset(\d+)HeatmapBoxWidth', sub { $datasets{"d$1"}{"hbox_width"} = $_[0]; }, 1],
						  ['--?dataset(\d+)MinPointValue', sub { $datasets{"d$1"}{"min_value"} = $_[0]; }, 1],
						  ['--?dataset(\d+)MidPointValue', sub { $datasets{"d$1"}{"mid_value"} = $_[0]; }, 1],
						  ['--?dataset(\d+)MaxPointValue', sub { $datasets{"d$1"}{"max_value"} = $_[0]; }, 1],
						  ['--?dataset(\d+)MinPointColor', sub { $datasets{"d$1"}{"min_color"} = $_[0]; }, 1],
						  ['--?dataset(\d+)MidPointColor', sub { $datasets{"d$1"}{"mid_color"} = $_[0]; }, 1],
						  ['--?dataset(\d+)MaxPointColor', sub { $datasets{"d$1"}{"max_color"} = $_[0]; }, 1],
						  ['--?dataset(\d+)CirclesMinR', sub { $datasets{"d$1"}{"circ_min_r"} = $_[0]; }, 1],
						  ['--?dataset(\d+)CirclesHorizontalGrid', sub { $datasets{"d$1"}{"circ_h_grid"} = $_[0]; }, 1],
						  ['--?dataset(\d+)CirclesVerticalGrid', sub { $datasets{"d$1"}{"circ_v_grid"} = $_[0]; }, 1],
						  ['--?dataset(\d+)CirclesSpacing', sub { $datasets{"d$1"}{"circ_spacing"} = $_[0]; }, 1],
						  ['--?dataset(\d+)CirclesFill', sub { $datasets{"d$1"}{"circ_fill"} = $_[0]; }, 1],
						  ['--?dataset(\d+)CirclesTransparent', sub { $datasets{"d$1"}{"circ_trans"} = $_[0]; }, 1]


	);

pod2usage(VERBOSE => 2) if ( $showHelp );

if (-e $configFile) {
	my $cfg = new Config::File::Simple($configFile);
	$treeFile = $cfg->read('treeFile');
	$treeFormat = $cfg->read('treeFormat');
	$treeName = $cfg->read('treeName');
	$treeDescription = $cfg->read('treeDescription');
	$ncbiFile = $cfg->read('ncbiFile');
	$ncbiFormat = $cfg->read('ncbiFormat');
	$showInternalIDs = $cfg->read('showInternalIDs');
	$uploadID = $cfg->read('uploadID');
	$projectName = $cfg->read('projectName');
	$preCollapsedFile = $cfg->read('preCollapsedFile');
	$branchLabelsFile = $cfg->read('branchLabelsFile');
	$popupInfoFile = $cfg->read('popupInfoFile');
	$fontStylesFile = $cfg->read('fontStylesFile');
	$internalNamesFile = $cfg->read('internalNamesFile');
	$HGTFile = $cfg->read('HGTFile');
	$colorDefinitionFile = $cfg->read('colorDefinitionFile');
	$assignTaxonomy = $cfg->read('assignTaxonomy');
	$midpointRoot = $cfg->read('midpointRoot');

	for (1..10) {
	  next unless (-e $cfg->read("dataset$_" . "File"));
	  $datasets{"d$_"}{"file"} = $cfg->read("dataset$_" . "File") if ($cfg->read("dataset$_" . "File"));
	  $datasets{"d$_"}{"label"} = $cfg->read("dataset$_" . "Label") if ($cfg->read("dataset$_" . "Label"));
	  $datasets{"d$_"}{"separator"} = $cfg->read("dataset$_" . "Separator") if ($cfg->read("dataset$_" . "Separator"));
	  $datasets{"d$_"}{"type"} = $cfg->read("dataset$_" . "Type") if ($cfg->read("dataset$_" . "Type"));
	  $datasets{"d$_"}{"color"} = $cfg->read("dataset$_" . "Color") if ($cfg->read("dataset$_" . "Color"));
	  $datasets{"d$_"}{"stripwid"} = $cfg->read("dataset$_" . "StripWidth") if ($cfg->read("dataset$_" . "StripWidth"));
	  $datasets{"d$_"}{"preventoverlap"} = $cfg->read("dataset$_" . "PreventOverlap") if ($cfg->read("dataset$_" . "PreventOverlap"));
	  $datasets{"d$_"}{"multialign"} = $cfg->read("dataset$_" . "MultiBarAlign") if ($cfg->read("dataset$_" . "MultiBarAlign"));
	  $datasets{"d$_"}{"pietrans"} = $cfg->read("dataset$_" . "PieTransparent") if ($cfg->read("dataset$_" . "PieTransparent"));
	  $datasets{"d$_"}{"piemax"} = $cfg->read("dataset$_" . "PieRadiusMax") if ($cfg->read("dataset$_" . "PieRadiusMax"));
	  $datasets{"d$_"}{"piemin"} = $cfg->read("dataset$_" . "PieRadiusMin") if ($cfg->read("dataset$_" . "PieRadiusMin"));
	  $datasets{"d$_"}{"barsize"} = $cfg->read("dataset$_" . "BarSizeMax") if ($cfg->read("dataset$_" . "BarSizeMax"));
	  $datasets{"d$_"}{"protsize"} = $cfg->read("dataset$_" . "ProtSizeMax") if ($cfg->read("dataset$_" . "ProtSizeMax"));
	  $datasets{"d$_"}{"boxplotsize"} = $cfg->read("dataset$_" . "BoxplotSize") if ($cfg->read("dataset$_" . "BoxplotSize"));
	  $datasets{"d$_"}{"boxplotraw"} = $cfg->read("dataset$_" . "BoxplotRawData") if ($cfg->read("dataset$_" . "BoxplotRawData"));
	  $datasets{"d$_"}{"boxplotupper"} = $cfg->read("dataset$_" . "BoxplotUpperPercentile") if ($cfg->read("dataset$_" . "BoxplotUpperPercentile"));
	  $datasets{"d$_"}{"boxplotlower"} = $cfg->read("dataset$_" . "BoxplotLowerPercentile") if ($cfg->read("dataset$_" . "BoxplotLowerPercentile"));
	  $datasets{"d$_"}{"boxplotextr"} = $cfg->read("dataset$_" . "BoxplotCalculateExtremes") if ($cfg->read("dataset$_" . "BoxplotCalculateExtremes"));
	  $datasets{"d$_"}{"boxplotstep"} = $cfg->read("dataset$_" . "BoxplotStepConstant") if ($cfg->read("dataset$_" . "BoxplotStepConstant"));
	  $datasets{"d$_"}{"coloringtype"} = $cfg->read("dataset$_" . "BranchColoringType") if ($cfg->read("dataset$_" . "BranchColoringType"));
	  $datasets{"d$_"}{"hbox_width"} = $cfg->read("dataset$_" . "HeatmapBoxWidth") if ($cfg->read("dataset$_" . "HeatmapBoxWidth"));
	  $datasets{"d$_"}{"min_value"} = $cfg->read("dataset$_" . "MinPointValue") if ($cfg->read("dataset$_" . "MinPointValue"));
	  $datasets{"d$_"}{"mid_value"} = $cfg->read("dataset$_" . "MidPointValue") if ($cfg->read("dataset$_" . "MidPointValue"));
	  $datasets{"d$_"}{"max_value"} = $cfg->read("dataset$_" . "MaxPointValue") if ($cfg->read("dataset$_" . "MaxPointValue"));
	  $datasets{"d$_"}{"min_color"} = $cfg->read("dataset$_" . "MinPointColor") if ($cfg->read("dataset$_" . "MinPointColor"));
	  $datasets{"d$_"}{"mid_color"} = $cfg->read("dataset$_" . "MidPointColor") if ($cfg->read("dataset$_" . "MidPointColor"));
	  $datasets{"d$_"}{"max_color"} = $cfg->read("dataset$_" . "MaxPointColor") if ($cfg->read("dataset$_" . "MaxPointColor"));
	  $datasets{"d$_"}{"circ_min_r"} = $cfg->read("dataset$_" . "CirclesMinR") if ($cfg->read("dataset$_" . "CirclesMinR"));
    $datasets{"d$_"}{"circ_h_grid"} = $cfg->read("dataset$_" . "CirclesHorizontalGrid") if ($cfg->read("dataset$_" . "CirclesHorizontalGrid"));
	  $datasets{"d$_"}{"circ_v_grid"} = $cfg->read("dataset$_" . "CirclesVerticalGrid") if ($cfg->read("dataset$_" . "CirclesVerticalGrid"));
	  $datasets{"d$_"}{"circ_spacing"} = $cfg->read("dataset$_" . "CirclesSpacing") if ($cfg->read("dataset$_" . "CirclesSpacing"));
	  $datasets{"d$_"}{"circ_fill"} = $cfg->read("dataset$_" . "CirclesFill") if ($cfg->read("dataset$_" . "CirclesFill"));
	  $datasets{"d$_"}{"circ_trans"} = $cfg->read("dataset$_" . "CirclesTransparent") if ($cfg->read("dataset$_" . "CirclesTransparent"));
	}
}

print "\niTOL batch uploader\n===================\n";
#need a tree
unless ( -e $treeFile or -e $ncbiFile) {
 print STDERR "Missing requried parameters. At least 'treeFile' or 'ncbiFile' must be specified. Use the --help option to display full help.\n";
 exit;
}

#sanity checks
if (defined $uploadID and not defined $projectName) {
	print STDERR "ERROR: When uploading to an account, you must specify the project name.\n\nUpload aborted.";
	exit;
}
if (defined $treeDescription and not defined $uploadID) {
	print STDERR "ERROR: Cannot use tree description without uploadID. You need an iTOL account to use tree descriptions.\n";
	exit;
}
if (defined $treeFormat and not ($treeFormat ne 'newick' or $treeFormat ne 'nexus')) {
  print STDERR "ERROR: Invalid tree format. Supported formats: 'newick' or 'nexus'\n";
  exit;
}
if (-e $treeFile and -e $ncbiFile) {
  print STDERR "ERROR: You cannot upload a tree and generate a new one from NCBI IDs at the same time.\n";
  exit;
}
if ($ncbiFormat and not $ncbiFormat =~ /^namesFull$|^namesCollapsed$|^idsFull$|^idsCollapsed$/) {
  print STDERR "ERROR: Invalid NCBI tree format.\nSupported formats:\n\t\t'namesFull' \t\t: generated tree will contain scientific names and internal nodes will not be collapsed\n\t\t'namesCollapsed' \t: scientific names will be used, and internal nodes with only one child removed\n\t\t'idsFull' \t\t: tree will contain NCBI taxonomy IDs and internal nodes will not be collapsed\n\t\t'idsCollapsed' \t\t: NCBI taxonomy IDs will be used, and internal nodes with only one child removed\n";
	exit;
}

#prepare the basic POST data
my %post_content;
if (-e $treeFile) { $post_content{'treeFile'} = [ $treeFile ]; }
if (-e $ncbiFile) { $post_content{'ncbiFile'} = [ $ncbiFile ]; }
if ($ncbiFormat) { $post_content{'ncbiFormat'} = $ncbiFormat; }

if ($uploadID) {
	$post_content{'uploadID'} = $uploadID;
	$post_content{'projectName'} = $projectName;
}
if ($treeFormat) { $post_content{'treeFormat'} = $treeFormat; }
if ($treeName) { $post_content{'treeName'} = $treeName; }
if ($treeDescription) { $post_content{'treeDescription'} = $treeDescription; }
if ($showInternalIDs) {$post_content{'showInternalIDs'} = 1; }
if ($assignTaxonomy) {$post_content{'assignTaxonomy'} = 1; }
if ($midpointRoot) {$post_content{'midpointRoot'} = 1; }
if (-e $fontStylesFile) { $post_content{'fontStylesFile'} = [ $fontStylesFile ]; }
if (-e $internalNamesFile) { $post_content{'internalNamesFile'} = [ $internalNamesFile ]; }
if (-e $HGTFile) { $post_content{'HGTFile'} = [ $HGTFile ]; }
if (-e $colorDefinitionFile) { $post_content{'colorDefinitionFile'} = [ $colorDefinitionFile ]; }
if (-e $preCollapsedFile) { $post_content{'preCollapsedFile'} = [ $preCollapsedFile ]; }
if (-e $branchLabelsFile) { $post_content{'branchLabelsFile'} = [ $branchLabelsFile ]; }
if (-e $popupInfoFile) { $post_content{'popupInfoFile'} = [ $popupInfoFile ]; }

#check the dataset info and fill the POST data
foreach my $dID (keys %datasets) {
	$dID =~ /d(\d+)/;
	my $pID = "dataset$1";
	#required parameters
	unless (-e $datasets{$dID}{'file'}) {
		print STDERR "ERROR: $pID file does not exist.\n"; exit;
	} else {
		$post_content{$pID . "File"} = [ $datasets{$dID}{'file'} ];
	}
	unless (length $datasets{$dID}{'label'}) {
		print STDERR "ERROR: $pID label not defined.\n"; exit;
	} else {
		$post_content{$pID . "Label"} = $datasets{$dID}{'label'};
	}
	unless ($datasets{$dID}{'separator'} =~ /^tab$|^comma$|^space$/) {
		print STDERR "ERROR: $pID has an invalid field separator.\n Supported separators: 'comma', 'space' or 'tab'\n"; exit;
	} else {
		$post_content{$pID . "Separator"} = $datasets{$dID}{'separator'};
	}
	unless ($datasets{$dID}{'type'} =~ /^binary$|^simplebar$|^multibar$|^gradient$|^colorstrip$|^time$|^domains$|^heatmap$|^boxplot$|^connections$|^circles$/) {
	  print STDERR "ERROR: $pID has an invalid type.\nSupported types:\n\t\t'binary' \t: Binary data\n\t\t'simplebar' \t: Simple bar chart\n\t\t'multibar' \t: Multi-value bar chart or Pie chart\n\t\t'gradient' \t: Color gradient\n\t\t'colorstrip' \t: Color strips\n\t\t'time' \t\t: Time series\n\t\t'domains' \t: Protein domain architecture\n\t\t'heatmap' \t: Heatmap\n\t\t'boxplot' \t: Boxplot\n\t\t'connections' \t: Connections\n\t\t'circles' \t: Circles\n"; exit;
	} else {
	  $post_content{$pID . "Type"} = $datasets{$dID}{'type'};
	}
	unless (length $datasets{$dID}{'label'}) {
	  print STDERR "ERROR: $pID label not defined.\n"; exit;
	}  else {
	  $post_content{$pID . "Label"} = $datasets{$dID}{'label'};
	}

	#optional parameters
	if (length $datasets{$dID}{'color'}) {
	  unless ($datasets{$dID}{'color'} =~ /^#......$/) {
		print STDERR "ERROR: $pID has an invalid color definition (" . $datasets{$dID}{'color'} . " ). Use the standard hexadecimal color notation (for example #ff0000 for red)\n"; exit;
	  }  else {
		$post_content{$pID . "Color"} = $datasets{$dID}{'color'};
	  }
	}
	if (length $datasets{$dID}{'stripwid'}) {
	  $post_content{$pID . "StripWidth"} = $datasets{$dID}{'stripwid'};
	}
	if ($datasets{$dID}{'preventoverlap'}) {
	  $post_content{$pID . "PreventOverlap"} = 1;
	}
	if (length $datasets{$dID}{'multialign'}) {
	  $post_content{$pID . "MultiBarAlign"} = $datasets{$dID}{'multialign'};
	}
	if (length $datasets{$dID}{'pietrans'}) {
	  $post_content{$pID . "PieTransparent"} = $datasets{$dID}{'pietrans'};
	}
	if (length $datasets{$dID}{'piemax'}) {
	  $post_content{$pID . "PieRadiusMax"} = $datasets{$dID}{'piemax'};
	}
	if (length $datasets{$dID}{'piemin'}) {
	  $post_content{$pID . "PieRadiusMin"} = $datasets{$dID}{'piemin'};
	}
	if (length $datasets{$dID}{'barsize'}) {
	  $post_content{$pID . "BarSizeMax"} = $datasets{$dID}{'barsize'};
	}
	if (length $datasets{$dID}{'protsize'}) {
	  $post_content{$pID . "ProtSizeMax"} = $datasets{$dID}{'protsize'};
	}
	if (length $datasets{$dID}{'boxplotsize'}) {
	  $post_content{$pID . "BoxplotSize"} = $datasets{$dID}{'boxplotsize'};
	}
	if (length $datasets{$dID}{'boxplotraw'}) {
	  $post_content{$pID . "BoxplotRawData"} = $datasets{$dID}{'boxplotraw'};
	}
	if (length $datasets{$dID}{'boxplotlower'}) {
	  $post_content{$pID . "BoxplotLowerPercentile"} = $datasets{$dID}{'boxplotlower'};
	}
	if (length $datasets{$dID}{'boxplotupper'}) {
	  $post_content{$pID . "BoxplotUpperPercentile"} = $datasets{$dID}{'boxplotupper'};
	}
	if (length $datasets{$dID}{'boxplotextr'}) {
	  $post_content{$pID . "BoxplotCalculateExtremes"} = $datasets{$dID}{'boxplotextr'};
	}
	if (length $datasets{$dID}{'boxplotstep'}) {
	  $post_content{$pID . "BoxplotStepConstant"} = $datasets{$dID}{'boxplotstep'};
	}
	if (length $datasets{$dID}{'coloringtype'}) {
	  unless ($datasets{$dID}{'coloringtype'} =~ /^none$|^both$|^branch$/) {
      print STDERR "ERROR: $pID has an invalid branch coloring type.\nSupported branch coloring types:\n\t\t'none' \t: no coloring, show only leaf boxes (Default)\n\t\t'branch': color branches only, don't display boxes\n\t\t'both' \t: color branches and show boxes\n"; exit;
    }  else {
      $post_content{$pID . "BranchColoringType"} = $datasets{$dID}{'coloringtype'};
	  }
	}
	if (length $datasets{$dID}{'hbox_width'}) {
	  $post_content{$pID . "HeatmapBoxWidth"} = $datasets{$dID}{'hbox_width'};
	}
	if (length $datasets{$dID}{'min_value'}) {
	  $post_content{$pID . "MinPointValue"} = $datasets{$dID}{'min_value'};
	}
	if (length $datasets{$dID}{'mid_value'}) {
	  $post_content{$pID . "MidPointValue"} = $datasets{$dID}{'mid_value'};
	}
	if (length $datasets{$dID}{'max_value'}) {
	  $post_content{$pID . "MaxPointValue"} = $datasets{$dID}{'max_value'};
	}
	if (length $datasets{$dID}{'min_color'}) {
	  $post_content{$pID . "MinPointColor"} = $datasets{$dID}{'min_color'};
	}
	if (length $datasets{$dID}{'mid_color'}) {
	  $post_content{$pID . "MidPointColor"} = $datasets{$dID}{'mid_color'};
	}
	if (length $datasets{$dID}{'max_color'}) {
	  $post_content{$pID . "MaxPointColor"} = $datasets{$dID}{'max_color'};
	}
	if (length $datasets{$dID}{'circ_min_r'}) {
	  $post_content{$pID . "CirclesMinR"} = $datasets{$dID}{'circ_min_r'};
	}
	if (length $datasets{$dID}{'circ_h_grid'}) {
	  $post_content{$pID . "CirclesHorizontalGrid"} = $datasets{$dID}{'circ_h_grid'};
	}
	if (length $datasets{$dID}{'circ_min_r'}) {
	  $post_content{$pID . "CirclesVerticalGrid"} = $datasets{$dID}{'circ_v_grid'};
	}
	if (length $datasets{$dID}{'circ_min_r'}) {
	  $post_content{$pID . "CirclesSpacing"} = $datasets{$dID}{'circ_spacing'};
	}
	if (length $datasets{$dID}{'circ_min_r'}) {
	  $post_content{$pID . "CirclesFill"} = $datasets{$dID}{'circ_fill'};
	}
	if (length $datasets{$dID}{'circ_min_r'}) {
	  $post_content{$pID . "CirclesTransparent"} = $datasets{$dID}{'circ_trans'};
	}

}

#submit the data
my $ua  = LWP::UserAgent->new();
$ua->agent("iTOLbatchUploader1.0");

my $req = POST $upload_url, Content_Type => 'form-data', Content => [ %post_content ];

my $response = $ua->request($req);

if ($response->is_success()) {
    my @res = split(/\n/, $response->content);
	#check for an upload error
	if ($res[$#res] =~ /^ERR/) {
		print "Upload failed. iTOL returned the following error message:\n\n$res[0]\n\n";
		exit;
	}
	#upload without warnings, ID on first line
	
	if ($res[0] =~ /^SUCCESS: (\S+)/) {
		print "Upload successful. Your tree is accessible using the following iTOL tree ID:\n\n$1\n\n";
		exit;
	}
	print "Upload successful. Warnings occured during upload:\n\n";
	for (0..$#res-1) {
		print $res[$_] . "\n";
	}
	if ($res[$#res] =~ /^SUCCESS: (\S+)/) {
		print "\n\nYour tree is accessible using the following iTOL tree ID:\n\n$1\n\n";
		exit;
	} else {
		print "This shouldn't happen. iTOL did not return an error, but there is no tree ID. Please email the full dump below to letunic\@embl.de:\n\n===DEBUG DUMP==\n";
		print join("\n", @res);
		print  "===END DUMP==\n";
	}
} else {
	print "iTOL returned a web server error. Full message follows:\n\n";
	print $response->as_string;
}

exit;
