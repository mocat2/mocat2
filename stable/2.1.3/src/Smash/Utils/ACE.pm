package Smash::Utils::ACE;

use strict;
use warnings;

=head1 NAME

Smash::Utils::ACE - ACE assembly format file parser.

=head1 SYNOPSIS

	# create the object

	my $ace = new Smash::Utils::ACE("454Contigs.ace");

	print $ace->contigs."\n";
	print $ace->reads."\n";

	# process each contig

	while (my $contig = $ace->nextContig) {
		printf ">%s\n", $contig->name;
		printf Smash::Core->pretty_fasta($contig->sequence);
	}

=cut

sub new {
	my $class    = shift;
	my %param    = @_;
	my $this     = {%param};
	bless $this, $class;

	my $ace_file = $this->file;
	open(FH, "<$ace_file") || die "Cannot open $ace_file: $!";
	$this->{FH}   = \*FH;
	$this->_parse_header() || die "Incorrect header in $ace_file\n";
	$this->getNextLine();
	return $this;
}

sub test {
	my $file = shift;
	my $ace = new Smash::Utils::ACE(FILE => $file || "454Contigs.ace", SOURCE => "Newbler");

	print $ace->contigs."\n";
	print $ace->reads."\n";

	while (my $contig = $ace->nextContig) {
		printf ">%s\n", $contig->name;
		printf Smash::Core->pretty_fasta($contig->sequence);
		my @features = $contig->getReadMappingGFF();
		foreach my $f (@features) {
			$f->print_feature_gff(\*STDOUT);
		}
	}
}

=head1 DESCRIPTION

=head2 Member variables

=over 4

=item B<C<file>>

=item B<C<source>>

=item B<C<contigs>>

=item B<C<reads>>

=back

=cut

sub _fh       {shift->{FH}}         # file-handle
sub _lastline {shift->{LAST_LINE}}  # last-read line
sub file      {shift->{FILE}}       # ace file
sub contigs   {shift->{CONTIGS}}    # number of contigs in this assembly
sub reads     {shift->{READS}}      # number of reads in this assembly
sub source    {shift->{SOURCE} || "ACEparser"}     # Who made this ace-file?

=head2 Member functions

=over 4

=item B<C<_parse_header()>>

parses the ACE header

=cut

sub _parse_header {
	my $this = shift;
	my $header = $this->getNextLine();
	if ($header =~ /^AS\s+(\d+)\s+(\d+)\s*$/) {
		$this->{CONTIGS} = $1;
		$this->{READS}   = $2;
		return 1;
	} else {
		return 0;
	}
}

=item B<C<nextContig()>>

returns the next contig in the ACE file as a C<Smash::Utils::ACE::Contig> object.

=cut

sub nextContig {
	my $this = shift;

	if (!$this->_lastline || $this->_lastline !~ /^CO\s+/) {
		return 0;
	}

	# parse contig header
	# CO contig00001 491 82 5 U

	my ($name, $length, $reads, $dummy, $dir);
	if ($this->_lastline =~ /^CO\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([UC]\s*$)/) {
		($name, $length, $reads, $dummy, $dir) = ($1, $2, $3, $4, $5);
	} else {
		die "Incorrect format: ".$this->_lastline."\n";
	}

	my $contig = Smash::Utils::ACE::Contig->new(NAME => $name, SOURCE => $this->source, PADDED_LENGTH => $length, READS => $reads, DIR => $dir);

	# parse the contig sequence
	# after CO

	my $buffer = $this->getSequence();

	die "Length mismatch: CO=>$length, Sequence=>".length($buffer)."\n" unless length($buffer) == $length;
	$contig->{PADDED_SEQUENCE} = $buffer;
	$buffer =~ s/\*//g;
	$contig->{SEQUENCE} = $buffer;
	$length = length($buffer);
	$contig->{LENGTH} = $length;

	$this->getNextLine();

	# parse contig quality values, if available
	# BQ

	my @buffer = ();
	if ($this->_lastline =~ /^BQ$/) {
		BQ:while (my $line = $this->getNextLine()) {
			last BQ if ($line =~ /^(AF|BS|RD)\s+/);
			push(@buffer, split(/\s+/, $line));
		}
		die "Length mismatch: CO=>$length, Quality=>".scalar(@buffer)."\n" unless scalar(@buffer) == $length;
	}
	$contig->{QUAL} = \@buffer;

	# parse read mappings
	# AF GDAUG0J08JHBOS U 1

	if ($this->_lastline =~ /^AF\s+/) {
		my $ReadMap = {};
		my $ReadDir = {};
		my @AF = ($this->_lastline);

		AF:while (my $line = $this->getNextLine()) {
			last AF unless ($line =~ /^AF\s+/);
			push(@AF, $line);
		}

		foreach my $AF (@AF) {
			my @w = split(/\s+/, $AF);
			my $read          = $w[1];
			$ReadDir->{$read} = ($w[2] eq "C")?"-":"+";
			$ReadMap->{$read} = $w[3];
		}
		$contig->{READMAPSTART} = $ReadMap;
		$contig->{READDIR} = $ReadDir;
	}

	# parse highqual read mappings
	# BS 1 232 GDAUG0J08JIE0W

	if ($this->_lastline =~ /^BS\s+/) {
		my $BSReadStart = {};
		my $BSReadEnd   = {};

		my @BS = ($this->_lastline);

		BS:while (my $line = $this->getNextLine()) {
			last BS unless ($line =~ /^BS\s+/);
			push(@BS, $line);
		}

		foreach my $BS (@BS) {
			my @w = split(/\s+/, $BS);
			my $read              = $w[3];
			$BSReadStart->{$read} = $w[1];
			$BSReadEnd->{$read}   = $w[2];
		}
		$contig->{BSSTART} = $BSReadStart;
		$contig->{BSEND}   = $BSReadEnd;
	}

	# parse RD tag

	RD:while ($this->_lastline =~ /^RD\s+/) {
		my (undef, $read, $length, $wa, $rt) = split(/\s+/, $this->_lastline);
		my $sequence = $this->getSequence();
		$contig->{$read}->{PADDED_SEQUENCE} = $sequence;
		$sequence =~ s/\*//g;
		$contig->{$read}->{SEQUENCE} = $sequence;

		while (my $line = $this->getNextLine()) {
			# QA and DS are ignored for now
			last unless ($line =~ /^(QA|DS)\s+/);
			if ($line =~ /^QA\s+/) {
				my ($undef, $qb, $qe, $alnBegin, $alnEnd) = split(/\s+/, $line);
				$contig->{READMAPSTART}->{$read} += ($alnBegin - 1);
				$contig->{READMAPEND}->{$read} = $contig->{READMAPSTART}->{$read} + ($alnEnd - $alnBegin);
			}
		}
	}

	# process the contig now!

	if ($this->_lastline =~ /^CO\s+/ || eof $this->_fh) {

		# Map the padded positions to unpadded positions.
		# We only need unpadded positions for the assembly summary.

		my $Padded2Unpadded = [];
		my @bases = split(//, $contig->{PADDED_SEQUENCE});

		# For now, map * to -1. We will fix it later

		for (my ($i, $j) = (0, 0); $i <= $#bases; $i++) {
			if ($bases[$i] eq "*") {
				$Padded2Unpadded->[$i] = -1; # cannot start mapping at *
			} else {
				$Padded2Unpadded->[$i] = $j;
				$j++;
			}
		}

		# Fix all the mappings of * to the mapping of the next non-*

		# I assume that the last padded base is NOT a *, which is reasonable?

		my $nextGoodOne = $#$Padded2Unpadded;
		for (my $i=$#$Padded2Unpadded; $i >= 0; $i--) {
			if ($Padded2Unpadded->[$i] == -1) {
				$Padded2Unpadded->[$i] = $nextGoodOne;
			} else {
				$nextGoodOne = $Padded2Unpadded->[$i];
			}
		}

		$contig->{PADDED2UNPADDED} = $Padded2Unpadded;
		return $contig;
	}
}

=item B<C<getNextLine()>>

reads the next non-empty line and stores it in C<LAST_LINE>.

=cut

sub getNextLine {
	my $this = shift;
	my $fh   = $this->{FH};

	return 0 if eof $fh;

	# Skip empty lines

	my $line;
	LINE:while ($line = <$fh>) {
		last LINE if ($line !~ /^\s*$/);
	}

	return 0 if eof $fh;

	# here's a non empty line

	chomp($line);
	$this->{LAST_LINE} = $line;
	return $line;
}

=item B<C<getSequence()>>

gets the sequence in the ACE file (sequence continues
until the next empty line in ACE format)

=cut

sub getSequence {
	my $this = shift;
	my $fh   = $this->{FH};
	my $buffer = "";
	SEQ:while (my $line = <$fh>) {
		last SEQ if (eof $fh);      # reached end of file
		$line =~ s/\s//g;
		last SEQ if $line =~ /^$/;  # reached end of sequence
		$buffer .= $line;
	}
	return $buffer;
}

=back

=cut

1;

package Smash::Utils::ACE::Contig;

use strict;
use warnings;

=head1 NAME

Smash::Utils::ACE::Contig - Class to hold a contig from the ACE file.

=head1 SYNOPSIS

	my $contig = Smash::Utils::ACE::Contig->new( \
				NAME => $name, \
				PADDED_LENGTH => $length, \
				READS => $reads, \
				DIR => $dir\
			);

	# print the unpadded contig sequence in fasta format

	printf ">%s\n%s", $contig->name, Smash::Core->pretty_fasta($contig->sequence);

=cut

sub new {
	my $class = shift;
	my %param = @_;
	my $this  = bless {%param}, $class;
}

=head2 Member variables

=over 4

=item B<C<name>>

=item B<C<source>>

=item B<C<reads>>

=item B<C<length>>

=item B<C<paddedLength>>

=item B<C<sequence>>

=item B<C<paddedSequence>>

=back

=cut

sub name           {shift->{NAME}}
sub source         {shift->{SOURCE} || "ACEparser"}
sub reads          {shift->{READS}}
sub length         {shift->{LENGTH}}
sub paddedLength   {shift->{PADDED_LENGTH}}
sub sequence       {shift->{SEQUENCE}}
sub paddedSequence {shift->{PADDED_SEQUENCE}}

=head2 Member functions

=over 4

=item B<C<getReadMappingGFF()>>

returns the contig2read mappings for this contig as an L<Utils::GFF> object.

=cut

sub getReadMappingGFF {
	use Smash::Utils::GFF;
	my $this = shift;
	my $Padded2Unpadded = $this->{PADDED2UNPADDED};
	my $ReadMapStart = $this->{READMAPSTART};
	my $ReadMapEnd   = $this->{READMAPEND};
	my $ReadDir = $this->{READDIR};
	my @features;
	my @reads   = sort {$ReadMapStart->{$a} <=> $ReadMapStart->{$b} || 
			    $ReadMapEnd->{$a} <=> $ReadMapEnd->{$b}} keys %$ReadMapStart;
	READ:foreach my $read (@reads) {
		my $p_start  = $ReadMapStart->{$read};
		my $up_start = $Padded2Unpadded->[$p_start-1]; # P2UP is 0-based
		my $p_end    = $ReadMapEnd->{$read};
		my $up_end   = $Padded2Unpadded->[$p_end-1]; # P2UP is 0-based
		if (!defined($up_start) || !defined($up_end)) {
			# Skip this read, since something's wrong
			#print "Oops!\n";
			next READ;
		}
		$up_start++; $up_end++;  # Back to 1-based
		my $attribute = sprintf("read \"%s\";", $read);
		my $line = join("\t", $this->name, $this->source, "read", $up_start, $up_end, ".", $ReadDir->{$read}, ".", $attribute);
		my $f = Smash::Utils::GFF::line2feature($line);
		push(@features, $f);
	}
	return @features;
}

=back

=cut

1;
