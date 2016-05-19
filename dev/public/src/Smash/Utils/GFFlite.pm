package Smash::Utils::GFFlite;
use strict;
use warnings;
use Smash::Utils::GFF;

sub new {
	my $class = shift;
	my $fh    = shift;
	if (ref $fh !~ /GLOB/) {
		die "GFFlite needs a GLOB reference!";
	}
	my $this = bless {FH => $fh}, $class;
	$this->{BUFFER} = $this->read_next_line();
	return $this;
}

sub nextFeature {
	my $this = shift;

	return 0 unless $this->{BUFFER};
	my $feature = Smash::Utils::GFF::line2feature($this->{BUFFER});
	$this->{BUFFER} = $this->read_next_line();
	return $feature;

}

sub nextSequence {
	my $this = shift;

	if (!$this->{LASTFEATURE}) {
		$this->{LASTFEATURE} = $this->nextFeature();
	}
	return 0 unless $this->{LASTFEATURE};
	my $last_feature  = $this->{LASTFEATURE};
	my $features = [$last_feature];
	my $f;
	while ($f = $this->nextFeature) {
		if ($f->seqname ne $last_feature->seqname) {
			$this->{LASTFEATURE} = $f;
			return $features;
		}
		push(@$features, $f);
	}

	# if I reach here, the $this->nextFeature failed, because of EOF

	$this->{LASTFEATURE} = undef;
	return $features;
}

sub read_next_line {
	my $this = shift;
	my $fh   = $this->{FH};
	my $line;

	return 0 if eof $fh;

	# get the next non-comment line

	do {
		$line = <$fh>;
	} while (!eof($fh) && $line =~ /^#/);

	return 0 unless $line;
	return 0 if ($line =~ /^#/);
	chomp($line);
	return $line;
}

1;
