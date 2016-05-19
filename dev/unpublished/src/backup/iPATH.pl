#!/usr/bin/perl
use strict;
use LWP::UserAgent;
use HTTP::Request::Common;
use Pod::Usage;
use Getopt::Regex;
use Config::File::Simple;

=pod

=head1 NAME

iPATH interface

=head1 SEE ALSO

iPATH Home page : http://pathways.embl.de

=head1 AUTHOR

Jens Roat Kultima <kultima@embl.de>

=cut

my $download_url = "http://pathways.embl.de/mapping.cgi";
my ( $showHelp, $outFile, $selection, $default_opacity, $keep_colors, $map_type, $tax_filter, $default_width, $default_radius, $background_color, $default_color, $whole_pathways, $query_reactions, $include_metabolic, $include_regulatory, $include_secondary, $png_dpi );

# Defaults
$default_opacity    = 0.1;
$keep_colors        = 0;
$map_type           = "PNG";
$include_metabolic  = 1;
$include_regulatory = 1;
$include_secondary  = 1;
$png_dpi            = 100;

Getopt::Regex::GetOptions(
	\@ARGV,
	[ '--?h(elp)?',              \$showHelp,           0 ],
	[ '--?outfile=?',            \$outFile,            1 ],
	[ '--?infile=?',             \$selection,          1 ],
	[ '--?default_opacity=?',    \$default_opacity,    1 ],
	[ '--?default_width=?',      \$default_width,      1 ],
	[ '--?default_radius=?',     \$default_radius,     1 ],
	[ '--?keep_colors=?',        \$keep_colors,        1 ],
	[ '--?default_color=?',      \$default_color,      1 ],
	[ '--?background_color=?',   \$background_color,   1 ],
	[ '--?whole_pathways=?',     \$whole_pathways,     1 ],
	[ '--?query_reactions=?',    \$query_reactions,    1 ],
	[ '--?tax_filter=?',         \$tax_filter,         1 ],
	[ '--?map_type=?',           \$map_type,           1 ],
	[ '--?include_metabolic=?',  \$include_metabolic,  1 ],
	[ '--?include_regulatory=?', \$include_regulatory, 1 ],
	[ '--?include_secondary=?',  \$include_secondary,  1 ],
	[ '--?png_dpi=?',            \$png_dpi,            1 ],
);

pod2usage( VERBOSE => 2 ) if ($showHelp);

print "\niPATH interface\n===============\n";
print "Loading...\n";

#prepare the basic POST data
my %post_content;
if ( defined $selection )          {
	open FILE, "$selection" or die "ERROR & EXIT: Could not open $selection\n"; 
	my $string = join("", <FILE>); 
	close FILE;	 
	$post_content{'selection'} = $string;
}
unless ($outFile) {
	$outFile = "iPATH.pathway.$map_type";
}
if ( defined $default_opacity )    { $post_content{'default_opacity'}    = $default_opacity; }
if ( defined $keep_colors )        { $post_content{'keep_colors'}        = $keep_colors; }
if ( defined $map_type )           { $post_content{'map_type'}           = $map_type; }
if ( defined $tax_filter )         { $post_content{'tax_filter'}         = $tax_filter; }
if ( defined $default_width )      { $post_content{'default_width'}      = $default_width; }
if ( defined $default_radius )     { $post_content{'default_radius'}     = $default_radius; }
if ( defined $background_color )   { $post_content{'background_color'}   = $background_color; }
if ( defined $default_color )      { $post_content{'default_color'}      = $default_color; }
if ( defined $whole_pathways )     { $post_content{'whole_pathways'}     = $whole_pathways; }
if ( defined $query_reactions )    { $post_content{'query_reactions'}    = $query_reactions; }
if ( defined $include_metabolic )  { $post_content{'include_metabolic'}  = $include_metabolic; }
if ( defined $include_regulatory ) { $post_content{'include_regulatory'} = $include_regulatory; }
if ( defined $include_secondary )  { $post_content{'include_secondary'}  = $include_secondary; }
if ( defined $png_dpi )            { $post_content{'png_dpi'}            = $png_dpi; }

#submit the data
my $ua = LWP::UserAgent->new();
$ua->agent("iPATHinterface");
my $req = POST $download_url, Content_Type => 'form-data', Content => [%post_content];
print "Sending...\n";
my $response = $ua->request($req);
if ( $response->is_success() && $map_type eq 'interactive') {
	
	my @res = split( /\n/, $response->content );
	open OUT, ">$outFile";
	foreach my $res (@res){
		print OUT "$res\n";
	}
	close OUT;
	
	
	if ( $response->header("Content-type") =~ /text\/html/ ) {
		my @res = split( /\n/, $response->content );
		die "Export failed. iPATH interface returned the following error message:\n\n$res[0]\n\n";
	}
	else {
		my @res = split( /\n/, $response->content );
		die "Export failed. iPATH interface returned the following error message:\n\n$res[0]\n\n";
	}
}
else {
	print "Checking...\n";
	if ($response->as_string =~ m/Location: (http.*)\n/) {
	my $address = $1;
	print "Downloading...\n";
	system "wget -O $outFile \"$address\"  ";
	print "===============\n";
	print "Success! File saved in $outFile\n";
	} else {
		print "An error occured.\nResponse:\n";
		print "$response->as_string\n";
		die;
	}
}

exit 0;
