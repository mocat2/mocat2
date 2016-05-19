package Smash::Utils::HTML;
use strict;
use warnings;
use LWP::UserAgent;

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(replace_macros trim_string strip_html get_attribute get_html_page get_html_page_persistent post_html_page retrieve_remote_file retrieve_remote_file_persistent parameterize_form_html get_tag get_tag_contents remove_andhash_tags);
our %EXPORT_TAGS = ('all' => [@EXPORT_OK]);

my $agent;

sub get_agent {
	use LWP::ConnCache;
	$agent = LWP::UserAgent->new(conn_cache => LWP::ConnCache->new());
}

sub replace_macros {
	my $string = shift;
	#$string =~ s///sg;
	$string =~ s/&quot;/'/sg;
	$string =~ s/&lt;/</sg;
	$string =~ s/&gt;/>/sg;
	$string =~ s/&nbsp;/ /sg;
	$string =~ s/&alpha;/alpha/sg;
	$string =~ s/&beta;/beta/sg;
	$string =~ s/&#40;/(/g;
	$string =~ s/&#41;/)/g;
	return $string;
}

sub trim_string {
	my $string = shift;
	$string =~ s/^\s+//s;
	$string =~ s/\s+$//s;
	return $string;
}

sub strip_html {
	my $string = shift;
	$string =~ s/<[^>]*?>//sg;
	return $string;
}

sub remove_andhash_tags {
	my $string = shift;
	return $string;
}

sub get_tag_contents {
	my ($string, $tag) = @_;
	my ($value) = $string =~ m|<\s*${tag}\b[^>]*>([^<]+)</${tag}>|si;
	return $value;
}

sub get_tag {
	my ($string, $tag) = @_;
	my ($value) = $string =~ m|(<\s*${tag}\b[^>]*>[^<]+</${tag}>)|si;
	if (!$value) {
		($value) = $string =~ m|(<\s*${tag}\b[^>]*/>)|si;
	}
	return $value;
}

sub get_attribute {
	my ($string, $attribute) = @_;
	my (undef, $value) = $string =~ m/\b${attribute}=(["'])([^\1]*?)\1/si;
	return $value;
}

sub process_remote_url_persistent {
	my ($url, $func) = @_;
	my $response;
	my $attempt = 0;
	get_agent() unless $agent;
	do {
		$response = $agent->get($url, ":content_cb" => $func);
		if (!$response->is_success) {
			#print $response->content;
			warn "ATTEMPT $attempt: ".$response->status_line."\n";
		}
		$attempt++;
		if ($attempt % 20 == 0) {
			die "Giving up after $attempt attempts!\n" if ($attempt == 200);
			$agent = undef;
			get_agent();
		}
	} while (!$response->is_success && sleep 60);
	return 1;
}

sub retrieve_remote_file_persistent {
	my ($url, $file) = @_;
	my $response;
	my $attempt = 0;
	get_agent() unless $agent;
	do {
		$response = $agent->get($url, ":content_file" => $file);
		if (!$response->is_success) {
			#print $response->content;
			warn "ATTEMPT $attempt: ".$response->status_line."\n";
		}
		$attempt++;
		if ($attempt % 20 == 0) {
			die "Giving up after $attempt attempts!\n" if ($attempt == 200);
			$agent = undef;
			get_agent();
		}
	} while (!$response->is_success && sleep 60);
	return 1;
}

sub retrieve_remote_file {
	my ($url, $file) = @_;
	my $response;
	get_agent() unless $agent;
	$response = $agent->get($url, ":content_file" => $file);
	if (!$response->is_success) {
		print $response->content;
		die $response->status_line;
	}
	return 1;
}

sub get_html_page {
	my $response;
	get_agent() unless $agent;
	push (@{ $agent->requests_redirectable }, 'GET');
	$response = $agent->get(@_);
	if (!$response->is_success) {
		die $response->status_line;
	}
	return $response->content;
}

sub get_html_page_persistent {
	my $response;
	my $attempt = 0;
	get_agent() unless $agent;
	do {
		$response = $agent->get(@_);
		if (!$response->is_success) {
			#print $response->content;
			warn "ATTEMPT $attempt: ".$response->status_line."\n";
		}
		$attempt++;
		if ($attempt % 20 == 0) {
			die "Giving up after $attempt attempts!\n" if ($attempt == 200);
			$agent = undef;
			get_agent();
		}
	} while (!$response->is_success && sleep 30);
	return $response->content;
}

sub post_html_page {
	my $response;
	get_agent() unless $agent;
	push (@{ $agent->requests_redirectable }, 'POST');
	$response = $agent->post(@_);
	if (!$response->is_success) {
		die $response->status_line;
	}
	return $response->content;
}

sub parameterize_form_html {
	my $html = shift;
	my %params;
	my @fields = ($html =~ m@<input .*?/>@sg);
	foreach my $field (@fields) {
		my $name = get_attribute($field, "name");
		my $value = get_attribute($field, "value");
		$params{$name} = $value;
	}
	return %params;
}

1;
