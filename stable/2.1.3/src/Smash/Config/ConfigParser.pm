package Smash::Config::ConfigParser;
use strict;
use warnings;

=head1 NAME

Smash::Config::Parser - Generic config parser.

=head1 SYNOPSIS

	use Smash::Config::Parser;
	my $parser = new Smash::Config::Parser("smash.conf");
	my $conf1  = $parser->parse();    # defaults to : as delimiter
	my $conf2  = $parser->parse("="); # explicitly sets = as delimiter
	print $conf1->{"Taxonomy"}->{"data_repos"};


=head1 DESCRIPTION

C<Smash::Config::Parser> is a generic parser with specific application to parse
the Smash config file.  The format of this file is:

	# Comments can go anywhere in a line that starts with a '#'.
	# Internal comments like perl, where a '#' symbol followed 
	#   by comments at the end of the line, are not recognized.
	# And empty lines don't matter.

	[Section1]
	key1	: value1
	# This comment is fine
	key2	: value2

	[Section2]
	key1	: value1
	key2	: value2 # This comment is NOT fine.
	key3	: value3 # This is part of value.

=head1 FUNCTIONS

=over 4

=item B<parse($delimiter)>

Once the object is made with the FILE attribute, this function parses
the file using the given delimiter, or by default colon (:). If there 
are multiple instances of delimiter, the first one is used as the real
delimiter and what follows (including the other instances of the 
delimiter) will be stored as the value of the key. All keys and 
values will be trimmed for whiltespace in the ends. It then returns a
pointer to a multilevel hash that can be queried as 
C<$conf->{section}->{key}>.

=back

=cut

# new

sub new {
	my $class  = shift;
	my $file   = shift;
	my %object = (FILE => $file);
	bless {%object}, $class;
}

sub file {shift->{FILE}}

sub parse {
	my $this    = shift;
	my $delim   = shift || ":";
	my $file    = $this->file;
	my $config  = {};
	my $section = "global";
	open(FILE, "<$file") || die "Cannot open config file $file: $!";
	while(<FILE>) {

		# Remove whitespace

		chomp();
		s/^\s+//;
		s/\s+$//;

		if (m/^\s*#/ || m/^\s*$/) {                      # skip empty lines and comment lines
			next;
		} elsif (m/^\s*\[([^\]]*)\]\s*$/) {                 # get the section name
			$section = $1;
			$config->{$section} = ();
		} elsif (m/${delim}/) {  # get name value pairs with delimiter
			my @fields = split($delim);
			$fields[0] =~ s/^\s+//;
			$fields[0] =~ s/\s+$//;
			if (@fields == 2) {
				$fields[1] =~ s/^\s+//;
				$fields[1] =~ s/\s+$//;
			}
			my ($key, $value);
			$key = shift @fields;
			$value = join($delim, @fields);
			$config->{$section}->{$key} = $value;
			if (@fields > 2 && $fields[0] ne "remote_repository") {
				warn "WARNING: Multiple ${delim}'s in following line in config file $file:\n";
				warn "\t$_\n";
				warn "First instance of $delim is used.\n";
			}
		}
	}
	close(FILE);
	return $config;
}

1;
