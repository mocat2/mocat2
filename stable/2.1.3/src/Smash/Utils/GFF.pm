package Smash::Utils::GFF;

use strict;
use warnings;

use base 'Exporter';
our @EXPORT_OK = qw(parse_gff text2flag);

sub line2feature {
	my @field_names = qw(SEQNAME SOURCE FEATURE START END SCORE STRAND FRAME ATTRIB);
	my $line = shift;
	my @words = split(/\t/, $line);
	if (@words < 8) { # attribute field could be empty
		die "Incorrect GFF feature line: $line";
	}
	my %hash  = map {$field_names[$_] => $words[$_]} 0..7; # attribs is yet to be parsed
	$hash{LOCATION} = [bless {START => $hash{START}, END => $hash{END}}, "Smash::Utils::GFF::Location"];

	########
	# parse the attribute field in GFF
	########

	if ($words[8]) {
		my @attrib_words = split(/\s*;\s*/, $words[8]);
		my $attributes = {};
		map {
				########
				# this definition comes from GFF specification on Sanger's website:
				# http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml (as on 13.05.2009)
				# http://www.sanger.ac.uk/resources/software/gff/spec.html    (as on 30.07.2010)
				########
				if (m/^([A-Za-z][A-Za-z0-9_]*)\s+([\"]?)(.*)\2/) {  # .*? is probably not necessary since it is split already, but just in case!
					my ($key, $value) = ($1, $3); # this also works as the '$key => $value' assignment
					$attributes->{$key} = $value if ($value);
				}
			   } @attrib_words;
		$hash{ATTRIBUTES} = $attributes;
	}
	my $feature = bless {%hash}, "Smash::Utils::GFF::Feature";
	return $feature;
}

sub raw_parse_gff {
	my $gff_file = shift;
	my @features;
	open(GFF, "<$gff_file") || die "Cannot open $gff_file: $!";
	while(<GFF>) {
		next if (m/^#/);
		next if (m/^$/);
		chomp();
		my $feature = line2feature($_);
		push(@features, $feature);
	}
	close(GFF);
	return \@features;
}

=head2 parse_gff

Parses a GFF file and returns a reference to an array containing L<Smash::Utils::GFF::Feature|Utils::GFF> object. Normally called
as

	my $gff = parse_gff($file);
	foreach my $feature (@$gff) {
		# do something on $feature
	}

but can also be called with an additional group field. Passing group field will make this function group the locations
of multiple features with the same value for group field. For example,

my $gff = parse_gff($file, FIELD => "attribute", ATTRIBUTE => "gene_id");

will group all CDS features belonging to the same gene. Beware that it only groups the locations by making the locations
attribute into an array of all the locations from the ungrouped features. It does not check if the other fields are the same,
that responsibility lies with the caller. For example, if there are two unrelated genes with the same value for the "gene"
attribute, they will be grouped as though they were the same entity.

=cut

sub parse_gff {
	my $gff_file = shift;
	my %options  = @_;

	my $group_attr  = $options{ATTRIBUTE};
	my $group_field = $options{FIELD};

	# Parse GFF file 

	my $gff_features = raw_parse_gff($gff_file);

	# If there is no grouping specified, return the raw gff
	if (!$group_field) {
		return $gff_features;
	}

	my @grouped; # Final array where grouped features are stored
	my @buffer;
	my %Keys; # Store all available keys to group on
	my $must_group = 0;

	# And see if any of them is split and if we really need to do the grouping procedure

	if ($group_field eq "attribute") {
		if ($group_attr) {
			my @groupable = grep {defined($_->get_attribute($group_attr))} @$gff_features;
			map {my $key = $_->get_attribute($group_attr); $Keys{$key}++; if ($Keys{$key} == 2) {$must_group = 1;}} @groupable;
		}

		if ($must_group == 0) {
			#return $gff_features;
		}

		# Anything without $group_attr, or just one occurrence is unchanged
		# Others are put in a buffer to process later
		foreach my $feature (@$gff_features) {
			my $key = $feature->get_attribute($group_attr);
			if (!defined($key) || $Keys{$key} == 1) {
				push(@grouped, $feature);
			} else {
				push(@buffer, $feature);
			}
		}

		# Process the buffer
		foreach my $key (grep {$Keys{$_} > 1} keys %Keys) {
			my @features = grep {$_->get_attribute($group_attr) eq $key} @buffer;
			my @locations = ();
			foreach my $feature (@features) {
				push(@locations, $feature->location->[0]);
			}
			$features[0]->{LOCATION} = \@locations;
			push(@grouped, $features[0]);
		}

		@grouped = sort {$a->seqname cmp $b->seqname || $a->start <=> $b->start || $a->end <=> $b->end || $a->get_attribute($group_attr) cmp $b->get_attribute($group_attr)} @grouped;

	} elsif ($group_field eq "seqname") {
		map {$Keys{$_->seqname}=1;} @$gff_features;

		# Process the gff
		foreach my $key (keys %Keys) {
			my @features = grep {$_->seqname eq $key} @$gff_features;
			my @locations = ();
			foreach my $feature (@features) {
				push(@locations, $feature->location->[0]);
			}
			$features[0]->{LOCATION} = \@locations;
			push(@grouped, $features[0]);
		}

		# Sorting like in attribute based grouping is not necessary, since this is not a usable GFF file anyway.
		# This is a quick hack to group features from a seqname, without worrying about the details
	}

	return \@grouped;
}

sub text2flag {
	my $text = shift;
	if (defined($text) && ($text eq "yes" || $text eq "1")) {
		return 1;
	} else {
		return 0;
	}
}

1;

package Smash::Utils::GFF::Feature;

use strict;
use warnings;

sub new {
	my $class = shift;
	my %param = @_;
	bless {%param}, $class;
}

sub seqname {shift->{SEQNAME}}
sub source  {shift->{SOURCE}}
sub feature {shift->{FEATURE}}
sub type    {shift->{TYPE}}
sub strand  {shift->{STRAND}}
sub score   {shift->{SCORE} || "."}
sub frame   {shift->{FRAME} || "."}
sub name    {shift->{NAME}}
sub gene_name {shift->{NAME}}

sub location  {shift->{LOCATION}}
sub start     {shift->{LOCATION}->[0]->start;}
sub begin     {shift->{LOCATION}->[0]->start;}
sub end {
	my $location = shift->{LOCATION};
	return $location->[$#$location]->end;
}

sub set_start {
	my $this  = shift;
	my $value = shift;
	$this->{LOCATION}->[0]->set_property("START", $value);
}
sub set_begin {shift->set_start(shift)}
sub set_end {
	my $this  = shift;
	my $value = shift;
	my $location = $this->{LOCATION};
	$location->[$#$location]->set_property("END", $value);
}

sub set_property {
	my $this = shift;
	my $key  = shift;
	my $value= shift;
	$this->{$key} = $value;
}
sub get_property {
	my $this = shift;
	my $key  = shift;
	return $this->{$key};
}

sub attributes {shift->{ATTRIBUTES}}
sub get_attribute {
	my $this = shift;
	my $key  = shift;
	return $this->{ATTRIBUTES}->{$key};
}
sub set_attribute {
	my $this = shift;
	my $key  = shift;
	my $value= shift;
	$this->{ATTRIBUTES}->{$key} = $value;
}
sub flattened_attributes {
	my $this   = shift;
	my $attrib = $this->attributes;
	return join("; ", map {sprintf("%s \"%s\"", $_, $attrib->{$_})} sort {$a cmp $b} keys %{$attrib}).";";
}

sub print_feature_gff {
	my $this   = shift;
	my $STREAM = shift || \*STDOUT;
	foreach my $location (@{$this->location}) {
		print $STREAM join("\t", $this->seqname, $this->source, $this->feature, $location->start, $location->end, $this->score, $this->strand, $this->frame, $this->flattened_attributes)."\n";
	}
}

sub copy_feature {
	my $this = shift;
	my $params = {};
	map {$params->{$_} = $this->{$_}} qw(SEQNAME SOURCE FEATURE STRAND SCORE FRAME NAME);
	map {$params->{ATTRIBUTES}->{$_} = $this->{ATTRIBUTES}->{$_}} keys %{$this->{ATTRIBUTES}};
	my $location = $this->{LOCATION};
	map {$params->{LOCATION}->[$_] = $location->[$_]->copy_location()} (0..$#$location);
	my $copy = new Smash::Utils::GFF::Feature(%$params);
	return $copy;
}

sub reverse_feature {
	my $this = shift;
	foreach my $l (@{$this->location}) {
		($l->{START}, $l->{END}) = ($l->{END}, $l->{START});
	}
	my @locations = reverse(@{$this->location});
	$this->{LOCATION} = \@locations;
	$this->{STRAND} = ($this->strand eq "+")?"-":"+";
}

sub offset_feature {
	my $this   = shift;
	my $offset = shift;
	foreach my $i (0..scalar(@{$this->location})-1) {
		$this->location->[$i]->offset($offset);
	}
}

1;

package Smash::Utils::GFF::Location;
use strict;
use warnings;

sub new {
	my $class = shift;
	my %params = @_;
	bless {%params}, $class;
}
sub start {shift->{START}};
sub begin {shift->{START}};
sub end   {shift->{END}};

sub copy_location {
	my $this = shift;
	my $params = {};
	map {$params->{$_} = $this->{$_}} qw(START END);
	my $copy = new Smash::Utils::GFF::Location(%$params);
	return $copy;
}

sub offset {
	my $this   = shift;
	my $offset = shift;
	$this->{START} += $offset;
	$this->{END} += $offset;
}

sub set_property {
	my $this = shift;
	my $key  = shift;
	my $value= shift;
	$this->{$key} = $value;
}

1;
