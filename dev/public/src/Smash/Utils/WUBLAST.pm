package Smash::Utils::WUBLAST;
use strict;
use warnings;
use base("Smash::Utils");

=head1 Name

	Smash::Utils::WUBLAST - A light-weight parser for WU-BLAST tabular output

=head1 Synopsis

	my $blast = Smash::Utils::WUBLAST->new("file.blastp");
	while (my $report = $blast->nextReport) {
		print "Report for Query: ".$report->queryName."\n";
		while (my $sbjct = $report->nextSbjct) {
			print "Subject: ".$sbjct->sbjctName."\n";
			while (my $group = $sbjct->nextGroup) {
				print "HSP Group: ".$group->groupId."\n";
				while (my $hsp = $group->nextHSP) {
					print $hsp->line."\n";
					$hsp->query;
					$hsp->sbjct;
					$hsp->qb;
					$hsp->qe;
					$hsp->sb;
					$hsp->se;
					$hsp->strand;
					$hsp->percent;
					$hsp->positive;
				}
			}
		}
	}

=head1 Description

Smash::Utils::WUBLAST is a light-weight parser for WU-BLAST tabular output, designed after the
BPlite parser from Ian Korf. It provides most of the functionalities (except the sequences and
alignments) using the same interfaces.

The C<group> functionality is useful when you specify C<topComboN> or C<topComboE> together with
C<links> to WU-BLAST. It will then link groups using a given id, and that can be used to track the
groups. If you do not specify this option, a default group will be created for all HSPs between the
query and a subject.

=cut

sub new {
	my $class = shift;
	my $file  = shift;
	my $this   = bless {}, $class;
	my $FH;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, "<$file") || die "Cannot open $file: $!";
	}
	$this->{FH} = $FH;

	# read one line

	my $line = <$FH>;
	return 0 unless $line; # empty file

	$this->{BUFFER} = new HSP($line);
	return $this;
}

sub nextReport {
	my $this = shift;
	my $FH = $this->{FH};

	# get the hsp for the last line read: from the buffer
	# if it is unset, then last line has been processed

	my $last_hsp = $this->{BUFFER};
	return 0 unless $last_hsp;
	my @hsps = ($last_hsp);

	# if there is no more line left, just process the buffer and 
	# set buffer to undef

	if (eof $FH) {
		delete $this->{BUFFER};
		return Report->new(@hsps) 
	}

	# if there are more lines, process them

	LINE: while (my $hsp = $this->getNextLineAsHSP) {
		if ($hsp->{QUERY} ne $last_hsp->{QUERY}) {
			last LINE;
		}
		push(@hsps, $hsp);
	}
	return Report->new(@hsps);
}

sub getNextLineAsHSP {
	my $this = shift;
	my $fh   = $this->{FH};
	my $line;

	# no more lines to read

	if (eof $fh) {
		delete $this->{BUFFER};
		return 0;
	}

	# get the next non-comment line

	do {
		$line = <$fh>;
	} while (!eof($fh) && $line =~ /^#/);

	# only unusable lines since last read usable line

	if (!$line || $line =~ /^#/) {
		delete $this->{BUFFER};
		return 0;
	}

	chomp($line);
	my $hsp = HSP->new($line);
	$this->{BUFFER} = $hsp;
	return $hsp;
}

sub test {
	my $this = shift;
	while (my $report = $this->nextReport) {
		print "REPORT:".$report->queryName."\n";
		while (my $sbjct = $report->nextSbjct) {
			print "SBJCT:".$sbjct->sbjctName."\n";
			while (my $group = $sbjct->nextGroup) {
				print "GROUP:".$group->groupId."\n";
				while (my $hsp = $group->nextHSP) {
					print $hsp->line."\n";
				}
			}
		}
	}
}

1;

package Report;
use strict;
use warnings;

sub new {
	my $class = shift;
	my @in_hsps = @_;
	my $this = bless {}, $class;
	my $sbjcts = [];
	my $prev = shift(@in_hsps);
	$this->{QUERY_NAME} = $prev->{QUERY};
	my @sbjct_hsps = ($prev);
	while (my $hsp = shift(@in_hsps)) {
		if ($hsp->{SBJCT} ne $prev->{SBJCT}) {
			push(@$sbjcts, new Subject(@sbjct_hsps));
			@sbjct_hsps = ();
		}
		push(@sbjct_hsps, $hsp);
		$prev = $hsp;
	}
	push(@$sbjcts, new Subject(@sbjct_hsps));
	$this->{SBJCTS} = $sbjcts;

	# populate buffer

	my $buffer = [];
	foreach my $sbjct (@$sbjcts) {
		push(@$buffer, $sbjct);
	}
	$this->{SBJCTS_BUFFER} = $buffer;

	return $this;
}

sub queryName {shift->{QUERY_NAME}}

sub nextSbjct {
	my $this = shift;
	my $buffer = $this->{SBJCTS_BUFFER};
	return shift @$buffer;
}

1;

package Subject;
use strict;
use warnings;

sub new {
	my $class = shift;
	my @in_hsps = @_;
	my $this = bless {}, $class; 
	my $groups = [];
	my $prev = shift(@in_hsps);
	$this->{SBJCT_NAME} = $prev->{SBJCT};
	my @group_hsps = ($prev);
	while (my $hsp = shift(@in_hsps)) {
		if ($hsp->{GROUP}) {
			if ($hsp->{GROUP} ne $prev->{GROUP}) {
				push(@$groups, new Group(@group_hsps));
				@group_hsps = ();
			}
		}
		push(@group_hsps, $hsp);
		$prev = $hsp;
	}
	push(@$groups, new Group(@group_hsps));
	$this->{GROUPS} = $groups;

	# populate buffer

	my $buffer = [];
	foreach my $group (@$groups) {
		push(@$buffer, $group);
	}
	$this->{GROUPS_BUFFER} = $buffer;

	return $this;
}

sub sbjctName {shift->{SBJCT_NAME}}

sub nextGroup {
	my $this = shift;
	my $buffer = $this->{GROUPS_BUFFER};
	return shift @$buffer;
}

1;

package Group;
use strict;
use warnings;

sub new {
	my $class = shift;
	my @in_hsps = @_;
	my $this = bless {}, $class;

	# populate hsps and buffer

	my $hsps = [];
	my $buffer = [];
	foreach my $hsp (@in_hsps) {
		push(@$hsps,   $hsp);
		push(@$buffer, $hsp);
	}
	$this->{HSPS}        = $hsps;
	$this->{HSPS_BUFFER} = $buffer;
	$this->{GROUP_ID} = $hsps->[0]->{GROUP};

	return $this;
}

sub groupId {return shift->{GROUP_ID} || "NA";}

sub nextHSP {
	my $this = shift;
	my $buffer = $this->{HSPS_BUFFER};
	return shift @$buffer;
}

1;

package HSP;
use strict;
use warnings;

sub new {
	my $class = shift;
	my $line = shift;
	my $this = bless {}, $class;
	chomp($line);
	my @w = split(/\t/, $line);
	my @labels = qw(QUERY SBJCT E N BITS RAW ALNLEN NIDENT NPOS NMISM PCIDENT PCPOS QGAPS QGAPLEN SGAPS SGAPLEN QFRAME QSTART QEND SFRAME SSTART SEND GROUP LINKS);
	map {$this->{$labels[$_]} = $w[$_]} (0..$#w);
	$this->{LINE} = $line;
	return $this;
}

sub line {shift->{LINE}}
sub query{shift->{QUERY}}
sub sbjct{shift->{SBJCT}}
sub qb   {shift->{QSTART}}
sub qe   {shift->{QEND}}
sub sb   {shift->{SSTART}}
sub se   {shift->{SEND}}
sub percent  {shift->{PCIDENT}}
sub positive {shift->{PCPOS}}
sub strand {
	my $this = shift;
	if ($this->{QFRAME} < 0) {
		return "-";
	} else {
		return "+";
	}
}

1;
