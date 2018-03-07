package Smash::CommandLineParser;

use strict;
use warnings;
use Smash::Global qw(:all);

use base 'Exporter';
our @EXPORT_OK = qw(parse_options check_required_options print_options);
our %EXPORT_TAGS = ('all' => [@EXPORT_OK], 'standard' => [qw(parse_options check_required_options print_options)]);

=head1 NAME

Smash::CommandLineParser - Specialized command line parser for Smash

=head1 SYNOPSIS

	use Smash::CommandLineParser qw(parse_options check_required_options);

	my @expected = qw(name=s number=f verbose);
	my @required = qw(name number);

	my ($status, %options) = parse_options(\@expected);
	if ($status != 1) {
		die "Error parsing command-line!";
	}

	my ($cstatus, $missing) = check_required_options(\%options, \@required);
	if ($cstatus != 1) {
		die "Missing argument: $missing";
	}

	printf "name=%s\n", $options{name};
	printf "number=%f\n", $options{number};

=head1 DESCRIPTION

Smash::CommandLineParser is an attempt to avoid several lines of code for command-line parsing
using C<Getopt::Long> and its robust function C<GetOptions()>. Every script in Smash needs to parse
command line and I wanted to make it a simpler process. The design is simple:

=over 4

=item 1

command line arguments should populate a hash

=item 2

option name-value pairs should be converted to hash key-value pairs

=item 3

a list of expected options should be sent to the parser, so it can ignore others

=item 4

for this list of expected options, they should exist in the hash even if no value was 
specified in the command line
(this is to avoid an C<exists()> check on every option to be clean)

=item 5

(optionally) a list of required options can be tested for their presence in the hash

=item 6

the caller function should know if every call went ok

=back

Smash::CommandLineParser implements it in a series of functions described below.

=head1 FUNCTIONS

=over 4

=item B<add_required_args(@options)>

=cut

sub add_required_args {
	my $this = shift;
	my @options = @_;
}

=item B<add_optional_args(@options)>

=cut

=item B<parse_options>

This function takes a (reference to a) list of allowed options and parses the command line 
and creates a hash where the option name-value pairs are converted into key-value pairs. 
It is analogous to a similar option in C<Getopt::Long> where you pass
a reference to a hash to C<GetOptions>, but slightly better since it brings to existence
an entry for every option in the list. If the option is not set, then the value for that key
is undefined, but it still exists. This makes sure that every option can be tested without 
worrying about testing it with C<exists()> or C<defined()>.

Returns a status code (1 for success, and others for failure) and the hash itself. The return
values should be handled like so:

	($status, %options) = parse_options([qw(name=s file=s)]);

=cut

sub parse_options {
	use Getopt::Long;
	use File::Basename;
	use Config;
	my $allowed  = shift;

	($SMASH_SCRIPT_NAME, $SMASH_SCRIPT_LOCATION) = fileparse($0);
	$SMASH_PERL = $Config{perlpath};

	my %options;

	my $status = 1;
	my %types = map {my ($key, $value) = split('='); $key => $value} @$allowed;

	# construct the specs for Getopt::Long.
	# if just a flag, keep it as "option" => \$options{option}
	# else, make it "option=s" => \$options{option}
	my %specs  = map {
				my $key;
				if (defined($types{$_})) { # with values
					$key = "$_=".$types{$_};
					my $type = $types{$_};
					if ($type =~ /[siof]{\d?,?\d?}/ || $type =~ /[siof]@/) {
						$options{$_} = [];
					} elsif ($type =~ /[siof]%/) {
						$options{$_} = {};
					} else {
						$options{$_} = undef;
					}
				} else { # without values
					$key = $_;
					$key =~ s/!$//;
					$key =~ s/\+$//;
					$options{$key} = undef;
				}
				$key => \$options{$_};
			} keys %types;

	# Get the arguments

	GetOptions(\%options, @$allowed) || ($status = 0);

	# Remove options that werent set

#	map {
#		if (!defined($options{$_})) {
#			delete $options{$_};
#		} elsif (ref($options{$_}) eq "ARRAY" && scalar(@{$options{$_}}) == 0) {
#			delete $options{$_};
#		} elsif (ref($options{$_}) eq "HASH" && scalar(keys %{$options{$_}}) == 0) {
#			delete $options{$_};
#		}
#	    } keys %options;

	# wipeout implicitly assumes unload as well!
	if ($options{wipeout}) {
		$options{unload} = 1;
	}

	return ($status, %options);
}

=item B<check_required_options>

This function takes a reference to a list of required options and a reference to a hash. It checks if each of the 
options is defined in the hash and quits on the first that does not. 

It returns a status code (1 for success, 0 if one of the options is missing). The missing options is also returned
as the second entry in the return array. For example,

	($status, $missing) = check_required_options([qw(name=s file=s)], \%options);

will return C<(0, 'file')> if C<--file> was not set in the command line.

=cut

sub check_required_options {
	my ($required, %options) = @_;
	# Check for required args
	foreach my $opt (@$required) {
		my $value = $options{$opt};
		if (!ref($value)) { # just a scalar
			if (!defined($value)) {
				return (0, "$opt");
			}
		} elsif (ref($value) eq "ARRAY") {
			if (scalar(@$value) == 0) {
				return (0, "$opt");
			}
		} elsif (ref($value) eq "HASH") {
			if (scalar(keys %$value) == 0) {
				return (0, "$opt");
			}
		}
	}
	return 1;
}

=item B<print_options>

This is a utility function that prints option name-value pairs. Useful for debugging purposes.

=back

=cut

sub print_options {
	my %options = @_;
	foreach my $key (sort {$a cmp $b} keys %options) {
		my $value = $options{$key};
		if (ref($value) eq "ARRAY") {
			$value = join(",", @$value);
		} elsif (ref($value) eq "HASH") {
			$value = "{".join(",", map { "$_=".$value->{$_}} keys %$value)."}";
		}
		print "$key=$value\n";
	}
}

1;
