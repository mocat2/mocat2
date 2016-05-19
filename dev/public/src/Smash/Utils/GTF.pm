package Smash::Utils::GTF;

use strict;
use warnings;
use Smash::Utils::GFF;

use base 'Exporter';
our @EXPORT_OK = qw(parse_gtf);

sub parse_gtf {
	my $gtf_file = shift;

	my $group_attr  = "gene_id";
	my $group_field = "attribute";

	my %Feature2Type = ("CDS" => "CDS", "cds" => "CDS", "start_codon" => "CDS", "stop_codon" => "CDS", "tRNA" => "tRNA", "rRNA" => "rRNA", "misc_RNA" => "miscRNA");

	# Parse GFF file 

	my $gff_features = Smash::Utils::GFF::raw_parse_gff($gtf_file);

	my @genes; # Final array where grouped features are stored
	my @buffer;
	my %Keys; # Store all available keys to group on

	# Find who needs to be grouped and who doesnt

	my @groupable = grep {defined($_->get_attribute($group_attr))} @$gff_features;
	map {my $key = $_->get_attribute($group_attr); $Keys{$key}++;} @groupable;

	# Here's the hash of all the genes that must be grouped
	my %ToGroup = map {$_ => 1} grep {$Keys{$_} > 1} keys %Keys;

	# Anything without $group_attr, or just one occurrence is unchanged
	# Others are put in a buffer to process later

	my $gene_class_name = __PACKAGE__."::Gene";

	foreach my $feature (@$gff_features) {
		my $key = $feature->get_attribute($group_attr);
		die "$group_attr not defined!" unless defined $key;
		my $type = $Feature2Type{$feature->feature} || "unknown";
		$feature->{TYPE} = $type;
		if (!defined($ToGroup{$key})) {
			my $gene = $gene_class_name->new(	NAME          => $key, 
								START         => $feature->start, 
								END           => $feature->end, 
								TYPE          => $type,
								FEATURE_COUNT => 1,
								FEATURES      => [$feature]);
			push(@genes, $gene);
		} else {
			push(@buffer, $feature);
		}
	}

	# Process the buffer
	foreach my $key (keys %ToGroup) {
		my %params = (NAME => $key);
		my $features = [grep {$_->get_attribute($group_attr) eq $key} @buffer];

		my $type = "unknown";
		my @locations = ();
		foreach my $feature (@$features) {
			push(@locations, $feature->location->[0]);
			$type = $Feature2Type{$feature->feature} if $Feature2Type{$feature->feature};
		}

		$params{TYPE} = $type;

		# get gene start and end

		@locations = sort {$a->start <=> $b->start} @locations;
		$params{START} = $locations[0]->start;

		@locations = sort {$b->end <=> $a->end} @locations; # watch out - reverse sort
		$params{END}   = $locations[0]->end;

		# assign features

		$params{FEATURES}      = $features;
		$params{FEATURE_COUNT} = scalar(@$features);

		push(@genes, $gene_class_name->new(%params));
	}

	@genes = sort {$a->seqname cmp $b->seqname || $a->start <=> $b->start || $a->end <=> $b->end || $a->get_attribute($group_attr) cmp $b->get_attribute($group_attr)} @genes;

	return \@genes;
}

1;

package Smash::Utils::GTF::Gene;
use strict;
use warnings;

use overload '""' => 'name';

sub new {
	my $class = shift;
	my %param = @_;
	bless {%param}, $class;
}

sub features {shift->{FEATURES}}
sub feature_count {shift->{FEATURE_COUNT}}

sub type  {shift->{TYPE}}
sub start {shift->{START}}
sub begin {shift->{START}}
sub end   {shift->{END}}
sub name  {shift->{NAME}}

sub _get_value {
	my $this = shift;
	my $key  = shift;
	my $f    = $this->{FEATURES}->[0];
	return $f->{$key};
}

sub seqname {shift->_get_value("SEQNAME")}
sub source  {shift->_get_value("SOURCE")}
sub strand  {shift->_get_value("STRAND")}

sub print_gene_gff {
	my $this   = shift;
	my $STREAM = shift;
	foreach my $feature (@{$this->features}) {
		$feature->print_feature_gff($STREAM);
	}
}

sub set_attribute {
	my $this = shift;
	my $key  = shift;
	my $value= shift;
	foreach my $f (@{$this->features}) {
		$f->set_attribute($key, $value);
	}
}

sub set_property {
	my $this = shift;
	my $key  = shift;
	my $value= shift;
	$this->{$key} = $value;
	foreach my $f (@{$this->features}) {
		$f->set_property($key, $value);
	}
}

sub copy_gene {
	my $this = shift;
	my $params = {};
	map {$params->{$_} = $this->{$_}} qw(SEQNAME SOURCE TYPE STRAND NAME FEATURE_COUNT);
	foreach my $i (0..($this->feature_count-1)) {
		$params->{FEATURES}->[$i] = $this->features->[$i]->copy_feature();
	}
	my $copy = new Smash::Utils::GTF::Gene(%$params);
	return $copy;
}

sub offset_gene {
	my $this   = shift;
	my $offset = shift;
	foreach my $i (0..($this->feature_count-1)) {
		$this->features->[$i]->offset_feature($offset);
	}
}

1;
