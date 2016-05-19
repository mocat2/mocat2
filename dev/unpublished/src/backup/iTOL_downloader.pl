#!/usr/bin/perl
use strict;
use LWP::UserAgent;
use HTTP::Request::Common;
use Pod::Usage;
use Getopt::Regex;
use Config::File::Simple;


=pod

=head1 NAME

iTOL_downloader  -  batch download trees from iTOL in different formats

=head1 SYNOPSIS

B<iTOL_downloader.pl> I<options>

=head1 DESCRIPTION

Use B<iTOL_downloader.pl> to download your trees from iTOL. You can specify the format and other export parameters using the options described below.


=head1 OPTIONS

Separate options and their values with a space.

=over 4

=item B<--help>

display this message

=item B<--config>

read upload options from a file (format described below); example file available at http://itol.embl.de/help/export.cfg

=item B<--tree>

iTOL tree ID

=item B<--outputFile>

save the resulting tree into the specified file

=item B<--format>

export format

=over 6

=item B<'png'>

Portable Network Graphics; bitmap format

=item B<'svg'>

Scalable Vector Graphics; vector format

=item B<'eps'>

Encapsulated Postscript; vector format

=item B<'ps'>

Postscript; vector format

=item B<'pdf'>

Portable Document Format; vector format

=item B<'nexus'>

Nexus tree file; plain text

=item B<'newick'>

Newick tree file; plain text

=back

=item B<--reRoot>

Reroot the tree on the specified node. If your tree does not have internal node IDs, you can use two leaf IDs separated with a vertical line. Last common ancestor of these nodes will be used as the new root (for example --reRoot 'Homo_sapiens|Gallus_gallus')

=back

=head2 OPTIONS FOR GRAPHICAL FORMATS

=over 4

=item B<--displayMode>

tree display mode; normal or circular

=item B<--alignLabels>

should the leaf labels be aligned? (0 or 1)

=item B<--omitLegend>

if set to 1, dataset and other legends will not be displayed

=item B<--ignoreBRL>

branch lengths will be ignored if set to 1

=item B<--showBRL>

branch lenght values will be displayed on the tree. Overrides --showBS if both are specified

=item B<--showBS>

bootstraps will be displayed on the tree

=item B<--BSdisplayValue>

specify which bootstrap will be displayed. You must specify the numeric value and choose the direction using the letters 'M' (for more than) or 'L' (for less than). For example:
 display bootstraps above 60: --BSdisplayValue M60  
 display bootstraps below 0.3: --BSdisplayValue L0.3

=item B<--BSdisplayType>

'text' or 'symbol'

=item B<--BSsymbolMaxSize>

if 'symbol' is used to display bootstraps, use this option to specify the maximum circle size (in pixels)

=item B<--hideRanges>

if the tree has colored ranges defined, you can omit them from the exported tree by setting this to 1 

=item B<--rangesCover>

'clades' or 'leaves'

if the tree has colored ranges defined, specify their extent using this option

=item B<--colorBranches>

if the tree has colored branches, use this option to use the colors in the exported tree

=item B<--rotation>

for circular display mode, specify the rotation angle in degrees

=item B<--arc>

for circular display mode, specify the arc in which the tree should be displayed (in degrees, default is 350)

=item B<--inverted>

if set to 1, the circular tree will be inverted

=item B<--resolution>

image resolution in DPI (dots per inch), used only when exporting to PNG

=item B<--hideLabels>

do not display leaf labels

=item B<--fontSize>

font size to be used for leaf labels

=item B<--lineWidth>

line width in pixels (default is 3)

=item B<--scaleFactor>

the default horizontal tree scale will be multiplied with the specified value

=item B<--showInternalLabels>

internal node IDs will be displayed if this option is set to 1

=item B<--pruneList>

comma separated list of leaf IDs which should be included in the tree. If not specified, complete tree will be exported.

=item B<--collapseList>

comma separated list of internal node IDs which will be collapsed. If the tree does not have internal nodes defined, you can specify them using two leafs IDs instead. Separate leaf IDs with a vertical line ('|'). The last common ancestor of each pair will be collapsed. For example:

 --collapseList Homo_sapiens|Gallus_gallus,Escherichia_coli|Mycoplasma_pneumoniae


=item B<--datasetList>

comma separated list of dataset IDs to include in the exported tree. For example:

 --datasetList dataset1,dataset3,dataset4

=item B<--omitDashedLines>

used only when datasets are included in the exported tree. When set, dataset values will not be connected to tree's leaves using dashed lines.

=back

=head2 OPTIONS FOR PLAIN TEXT FORMATS

=over 4

=item B<--newickExportFormat>

output options for Newick format

=over 6

=item B<'withIDs'>

Include internal node IDs in the tree. Example output: ((A:0.1, B:0.1)INTid1:0.2[95]);

=item B<'withoutIDs'>

Only bootstraps and branch lengths will be included. Example output: ((A:0.1, B:0.1)95:0.2);

=back

=back

=head1 USING A CONFIGURATION FILE

We provided an option of using a separate plain text config file to specify the export parameters. To use a configuration file, simply use the B<--config> option to specify its location.

The configuration file should contain one parameter per line. Use the B<=> sign to assign values to the parameters. Here is an example configuration file:

 outputFile = example.pdf
 format = pdf
 tree = 123456789
 displayMode = circular
 inverted = 1
 arc = 360
 rotate = 70

=head1 SEE ALSO

 iTOL Home page : http://itol.embl.de
 iTOL Help      : http://itol.embl.de/help/help.shtml

=head1 AUTHOR

 Ivica Letunic <ivica@letunic.com>
 Contact me if you have any questions or comments.

=cut

my $download_url = "http://itol.embl.de/batch_downloader.cgi";
my ($configFile, $outFile, $showHelp, $tree, $format,$pruneList,$collapseList,$omitLegend,$ignoreBRL,$newickExportFormat, $colorBranches, $displayMode, $alignLabels,  $showBS, $showBRL, $reRoot, $hideRanges, $rangesCover, $BSdisplayValue, $BSdisplayType, $BSsymbolMaxSize, $rotation, $arc, $inverted, $resolution, $hideLabels, $fontSize, $lineWidth, $scaleFactor, $showInternalLabels, $datasetList, $omitDashedLines);


Getopt::Regex::GetOptions(\@ARGV,
						  ['--?h(elp)?',\$showHelp, 0],
						  ['--?config=?',\$configFile, 1],
						  ['--?tree=?',\$tree, 1],
						  ['--?outputFile=?',\$outFile, 1],
						  ['--?format=?',\$format, 1],
						  ['--?displayMode=?',\$displayMode, 1],
						  ['--?alignLabels=?',\$alignLabels, 1],
						  ['--?pruneList=?',\$pruneList, 1],
						  ['--?collapseList=?',\$collapseList, 1],
						  ['--?omitLegend=?',\$omitLegend, 1],
						  ['--?ignoreBRL=?',\$ignoreBRL, 1],
						  ['--?showBS=?',\$showBS, 1],
						  ['--?showBRL=?',\$showBRL, 1],
						  ['--?reRoot=?',\$reRoot, 1],
						  ['--?hideRanges=?',\$hideRanges, 1],
						  ['--?rangesCover=?',\$rangesCover, 1],
						  ['--?colorBranches=?',\$colorBranches, 1],
						  ['--?BSdisplayValue=?',\$BSdisplayValue, 1],
						  ['--?BSdisplayType=?',\$BSdisplayType, 1],
						  ['--?BSsymbolMaxSize=?',\$BSsymbolMaxSize, 1],
						  ['--?rotation=?',\$rotation, 1],
						  ['--?arc=?',\$arc, 1],
						  ['--?inverted=?',\$inverted, 1],
						  ['--?resolution=?',\$resolution, 1],
						  ['--?hideLabels=?',\$hideLabels, 1],
						  ['--?fontSize=?',\$fontSize, 1],
						  ['--?lineWidth=?',\$lineWidth, 1],
						  ['--?scaleFactor=?',\$scaleFactor, 1],
						  ['--?showInternalLabels=?',\$showInternalLabels, 1],
						  ['--?newickExportFormat=?',\$newickExportFormat, 1],
						  ['--?datasetList=?',\$datasetList, 1],
              ['--?omitDashedLines=?',\$omitDashedLines, 1]					  
	);

pod2usage(VERBOSE => 2) if ( $showHelp );

if (-e $configFile) {
	my $cfg = new Config::File::Simple($configFile);
	$tree = $cfg->read('tree') if ($cfg->read('tree'));
	$outFile = $cfg->read('outputFile')  if ($cfg->read('outputFile') );
	$format = $cfg->read('format')  if ($cfg->read('format') );
	$omitLegend = $cfg->read('omitLegend')  if ($cfg->read('omitLegend') );
	$ignoreBRL = $cfg->read('ignoreBRL')  if ($cfg->read('ignoreBRL') );
	$pruneList = $cfg->read('pruneList')  if ($cfg->read('pruneList') );
	$collapseList = $cfg->read('collapseList')  if ($cfg->read('collapseList') );
	$newickExportFormat = $cfg->read('newickExportFormat')  if ($cfg->read('newickExportFormat') );
	$scaleFactor = $cfg->read('scaleFactor')  if ($cfg->read('scaleFactor') );
	$colorBranches = $cfg->read('colorBranches')  if ($cfg->read('colorBranches') );
	$displayMode = $cfg->read('displayMode')  if ($cfg->read('displayMode') );
	$hideLabels = $cfg->read('hideLabels')  if ($cfg->read('hideLabels') );
	$fontSize = $cfg->read('fontSize')  if ($cfg->read('fontSize') );
	$rotation = $cfg->read('rotation')  if ($cfg->read('rotation') );
	$arc = $cfg->read('arc')  if ($cfg->read('arc') );
	$inverted = $cfg->read('inverted')  if ($cfg->read('inverted') );
	$lineWidth = $cfg->read('lineWidth')  if ($cfg->read('lineWidth') );
	$showInternalLabels = $cfg->read('showInternalLabels')  if ($cfg->read('showInternalLabels') );
	$alignLabels = $cfg->read('alignLabels')  if ($cfg->read('alignLabels') );
	$showBS = $cfg->read('showBS')  if ($cfg->read('showBS') );
	$showBRL = $cfg->read('showBRL')  if ($cfg->read('showBRL') );
	$hideRanges = $cfg->read('hideRanges')  if ($cfg->read('hideRanges') );
	$reRoot = $cfg->read('reRoot')  if ($cfg->read('reRoot') );
	$rangesCover = $cfg->read('rangesCover')  if ($cfg->read('rangesCover') );
	$BSdisplayValue = $cfg->read('BSdisplayValue')  if ($cfg->read('BSdisplayValue') );
	$BSdisplayType = $cfg->read('BSdisplayType')  if ($cfg->read('BSdisplayType') );
	$BSsymbolMaxSize = $cfg->read('BSsymbolMaxSize')  if ($cfg->read('BSsymbolMaxSize') );
	$resolution = $cfg->read('resolution')  if ($cfg->read('resolution') );
	$datasetList = $cfg->read('datasetList') if ($cfg->read('datasetList'));
	$omitDashedLines = $cfg->read('omitDashedLines') if ($cfg->read('omitDashedLines'));
}

print "\niTOL batch downloader\n=====================\n";

#need a tree ID, format and outfile
unless ( length($tree) and length($outFile) and length($format) ) {
  print STDERR "Missing requried parameters. At least 'tree', 'format' and 'outputFile' must be specified. Use the --help option to display full help.\n";
 exit;
}

unless ($format eq 'svg' or $format eq 'eps' or $format eq 'ps' or $format eq 'pdf' or $format eq 'png' or $format eq 'newick' or $format eq 'nexus') {
  print STDERR "ERROR: Invalid output format. Supported formats: 'eps', 'svg', 'ps', 'pdf', 'png', 'newick' or 'nexus'\n";
  exit;
}


#prepare the basic POST data
my %post_content;
$post_content{'format'} = $format;
$post_content{'tree'} = $tree;

if (defined $displayMode) { $post_content{'displayMode'} = $displayMode; } 
if (defined $hideLabels) { $post_content{'hideLabels'} = $hideLabels; } 
if (defined $fontSize) { $post_content{'fontSize'} = $fontSize; } 
if (defined $omitLegend) { $post_content{'omitLegend'} = $omitLegend; }
if (defined $ignoreBRL) { $post_content{'ignoreBRL'} = $ignoreBRL; }
if (defined $pruneList) { $post_content{'pruneList'} = $pruneList; }
if (defined $collapseList) { $post_content{'collapseList'} = $collapseList; }
if (defined $scaleFactor) { $post_content{'scaleFactor'} = $scaleFactor; } 
if (defined $colorBranches) { $post_content{'colorBranches'} = $colorBranches; } 
if (defined $newickExportFormat) { $post_content{'newickExportFormat'} = $newickExportFormat; }
if (defined $displayMode) { $post_content{'displayMode'} = $displayMode; }
if (defined $fontSize) { $post_content{'fontSize'} = $fontSize; }
if (defined $rotation) { $post_content{'rotation'} = $rotation; }
if (defined $arc) { $post_content{'arc'} = $arc; }
if (defined $inverted) { $post_content{'inverted'} = $inverted; }
if (defined $lineWidth) { $post_content{'lineWidth'} = $lineWidth; }
if (defined $showInternalLabels) { $post_content{'showInternalLabels'} = $showInternalLabels; }
if (defined $alignLabels) { $post_content{'alignLabels'} = $alignLabels; }
if (defined $showBS) { $post_content{'showBS'} = $showBS; }
if (defined $showBRL) { $post_content{'showBRL'} = $showBRL; }
if (defined $reRoot) { $post_content{'reRoot'} = $reRoot; }
if (defined $hideRanges) { $post_content{'hideRanges'} = $hideRanges; }
if (defined $rangesCover) { $post_content{'rangesCover'} = $rangesCover; }
if (defined $BSdisplayType) { $post_content{'BSdisplayType'} = $BSdisplayType; }
if (defined $BSsymbolMaxSize) { $post_content{'BSsymbolMaxSize'} = $BSsymbolMaxSize; }
if (defined $resolution) { $post_content{'resolution'} = $resolution; }
if (defined $datasetList) { $post_content{'datasetList'} = $datasetList; }
if (defined $omitDashedLines) { $post_content{'omitDashedLines'} = $omitDashedLines; }

#submit the data
my $ua  = LWP::UserAgent->new();
$ua->agent("iTOLbatchDownloader1.0");
my $req = POST $download_url, Content_Type => 'form-data', Content => [ %post_content ];
my $response = $ua->request($req);
if ($response->is_success()) {
  #check for the content type
  #if text/html, there was an error
  #otherwise dump the results into the outfile
  if ($response->header("Content-type") =~ /text\/html/) {
	my @res = split(/\n/, $response->content);
	print "Export failed. iTOL returned the following error message:\n\n$res[0]\n\n";
	exit;
  } else {
	open (OUT, ">$outFile") or die "Cannot write to $outFile";
	binmode OUT;
	print OUT $response->content;
	print "Exported tree saved to $outFile\n";

	# print join("\n", @res);
  }
} else {
	print "iTOL returned a web server error. Full message follows:\n\n";
	print $response->as_string;
}
