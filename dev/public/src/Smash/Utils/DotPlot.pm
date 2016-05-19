package Smash::Utils::DotPlot;
use strict;
use warnings;
use Smash::Utils::SVG qw(:all);

#our @ISA = qw(Smash::Utils::Plot);

#########################################################
# arguments sent in a hash
# required arguments:
# BLAST_FILE, 
# SIMILARITY, BITS, E, MINALIGNSLEN
# QLENFILE, SLENFILE, 
#########################################################

# QLEN and SLEN are like [["seq1", 1000], ["seq2", 2000"]]

sub new {
	my $class  = shift;
	my %params = @_;
	my $this   = bless {%params}, $class;
	return $this;
}

sub init {
	my $this = shift;
	$this->{SYNTENY}       = {};
	$this->{SYNTENY_COUNT} = {};
	$this->{QORIENTATION}  = {};
	$this->{PROTEIN_SYNTENY}       = {};
	$this->{PROTEIN_SYNTENY_COUNT} = {};
}

sub draw {
	my $this = shift;
	$this->init();
	$this->parse_genome_blast(); # we now have SYNTENY
	$this->parse_lengths(); # we now have SLENGTH, QLENGTH
	foreach my $s (keys %{$this->{SLENGTH}}) {
		foreach my $q (keys %{$this->{QLENGTH}}) {
			$this->{SYNTENY_COUNT}->{$s}->{$q} = 0 unless defined $this->{SYNTENY_COUNT}->{$s}->{$q};
		}
	}
	if ($this->{PROTEIN_BLAST_FILE}) {
		$this->parse_protein_blast();
	}
	map {$this->{$_} = undef} qw(BLAST_FILE PROTEIN_BLAST_FILE SIMILARITY BITS MINALIGNSLEN);
	$this->draw_dot_plot();
}

# parse the sequence lengths

sub parse_lengths {
	my $this          = shift;
	my $QLENGTH_FILE  = $this->{QLENFILE};
	my $SLENGTH_FILE  = $this->{SLENFILE};
	my $QLEN          = $this->{QLEN};
	my $SLEN          = $this->{SLEN};
	my $SyntenyCount  = $this->{SYNTENY_COUNT};

	my $show_hits_only = 1;
	my $KeepQuery = {};
	my $KeepSbjct = {};

	if ($show_hits_only) {
		foreach my $sbjct (keys %$SyntenyCount) {
			$KeepSbjct->{$sbjct} = 1;
			foreach my $query (keys %{$SyntenyCount->{$sbjct}}) {
				if ($SyntenyCount->{$sbjct}->{$query}) {
					$KeepQuery->{$query} = 1;
				}
			}
		}
	}

	# query/sbjct lengths

	my %QLength;
	my %SLength;

	# Get subject length information

	if ($SLEN) {
		my @Length = @$SLEN;
		my $sum = 0;
		for (my $i=0; $i<=$#Length; $i++) {
			my ($seq, $length) = map {$Length[$i][$_]} (0,1);
			$seq = parse_genbank_header($seq);
			$SLength{$seq} = $length;
		}
	} elsif ($SLENGTH_FILE) {
		my $residue_length = 1;
		if ($SLENGTH_FILE =~ /protein/) {
			$residue_length = 3;
		}
		open(LEN, "<$SLENGTH_FILE") || die "Cannot open $SLENGTH_FILE: $!";
		while (<LEN>) {
			chomp();
			my ($seq, $length) = split(/\t/);
			$seq = parse_genbank_header($seq);
			$length *= $residue_length;
			$SLength{$seq} = $length;
		}
		close(LEN);
	} else {
		die "Sbjct lengths not defined";
	}

	# Get query length information

	if ($QLEN) {
		my @Length = @$QLEN;
		for (my $i=0; $i<=$#Length; $i++) {
			my ($seq, $length) = map {$Length[$i][$_]} (0,1);
			$QLength{$seq} = $length;
		}
	} elsif ($QLENGTH_FILE) {
		my $residue_length = 1;
		if ($QLENGTH_FILE =~ /protein/) {
			$residue_length = 3;
		}
		open(LEN, "<$QLENGTH_FILE") || die "Cannot open $QLENGTH_FILE: $!";
		while (<LEN>) {
			chomp();
			my ($seq, $length) = split(/\t/);
			$length *= $residue_length;
			$QLength{$seq} = $length;
		}
		close(LEN);
	} else {
		die "Query lengths not defined";
	}

	if ($show_hits_only) {
		my @queries = keys %QLength;
		foreach my $query (@queries) {
			delete $QLength{$query} unless $KeepQuery->{$query};
		}
		my @sbjcts = keys %SLength;
		foreach my $sbjct (@sbjcts) {
			delete $SLength{$sbjct} unless $KeepSbjct->{$sbjct};
		}
	}

	# clear parameters that are not useful anymore

	map {$this->{$_} = undef} qw(QLENFILE SLENFILE QLEN SLEN);

	$this->{SLENGTH} = \%SLength;
	$this->{QLENGTH} = \%QLength;
}

sub parse_protein_blast {
	my $this                = shift;
	my $BLAST_FILE          = $this->{PROTEIN_BLAST_FILE} || return;
	my $ProteinSynteny      = $this->{PROTEIN_SYNTENY}; 
	my $ProteinSyntenyCount = $this->{PROTEIN_SYNTENY_COUNT}; 

	my ($scds, $qcds) = $this->parse_gff_files();

	my %SProteinStart;
	my %SProteinEnd;
	my %SProteinStrand;
	my %SProteinSeq;
	my %QProteinStart;
	my %QProteinEnd;
	my %QProteinStrand;
	my %QProteinSeq;

	foreach my $feature (@$scds) {
		my $gene = $feature->get_attribute("gene_id");
		$SProteinStart{$gene}  = $feature->begin;
		$SProteinEnd{$gene}    = $feature->end;
		$SProteinStrand{$gene} = $feature->strand;
		$SProteinSeq{$gene}    = $feature->seqname;
	}

	foreach my $feature (@$qcds) {
		my $gene = $feature->get_attribute("gene_id");
		$QProteinStart{$gene}  = $feature->begin;
		$QProteinEnd{$gene}    = $feature->end;
		$QProteinStrand{$gene} = $feature->strand;
		$QProteinSeq{$gene}    = $feature->seqname;
	}

	open(BLAST, "<$BLAST_FILE") || die "Cannot open $BLAST_FILE: $!";
	QRY:while (<BLAST>) {
		chomp();
		my @w = split(/\t/);
		my ($query, $sbjct, $bits, $similarity) = ($w[0], $w[1], $w[4], $w[11], $w[17], $w[18], $w[20], $w[21]);
		$sbjct = parse_genbank_header($sbjct);

		# remap coordinates in genome coordinates
		# TO DO: use EXACT coordinates of protein hits by using locations in gff object to calculate the real begin and end

		my ($sStart, $sEnd) = ($SProteinStart{$sbjct}, $SProteinEnd{$sbjct});
		my ($qStart, $qEnd) = ($QProteinStart{$query}, $QProteinEnd{$query});

		# if the genes are on different strands of genome, reverse the locus
		# rotating the genome itself will be taken care of later by $QOrient

		if ($QProteinStrand{$query} ne $SProteinStrand{$sbjct}) {
			($qStart, $qEnd) = ($qEnd, $qStart);
		}
if (!$qStart || !$qEnd) {
	print "$query uwe!\n";
	next QRY;
}

		$ProteinSynteny->{$SProteinSeq{$sbjct}}->{$QProteinSeq{$query}}->{$ProteinSyntenyCount->{$SProteinSeq{$sbjct}}->{$QProteinSeq{$query}}++} = "$sStart,$sEnd:$qStart,$qEnd:$bits";
	}
	close(BLAST);
}

sub parse_gff_files {
	use Smash::Utils::GFF;
	my $this = shift;
	my $qgff = $this->{QGTFFILE};
	my $sgff = $this->{SGTFFILE};
	my $qcds = Smash::Utils::GFF::parse_gff($qgff, (field => "attribute", ATTRIBUTE => "gene_id"));
	my $scds = Smash::Utils::GFF::parse_gff($sgff, (field => "attribute", ATTRIBUTE => "gene_id"));
	return ($scds, $qcds);
}

sub parse_genome_blast {
	my $this         = shift;
	my $Synteny      = $this->{SYNTENY}; 
	my $SyntenyCount = $this->{SYNTENY_COUNT}; 
	my $BLAST_FILE   = $this->{BLAST_FILE};
	open(BLAST, "<$BLAST_FILE") || die "Cannot open $BLAST_FILE: $!";
	QRY:while (<BLAST>) {
		chomp();
		my @w = split(/\t/);
		my ($query, $sbjct, $bits, $similarity, $qStart, $qEnd, $sStart, $sEnd) = ($w[0], $w[1], $w[4], $w[11], $w[17], $w[18], $w[20], $w[21]);
		if ($similarity < $this->{SIMILARITY} || $bits < $this->{BITS}) {
			next QRY;
		}
		if ($sEnd - $sStart < $this->{MINALIGNSLEN}) {
			next QRY;
		}

		$sbjct = parse_genbank_header($sbjct);

		$Synteny->{$sbjct}->{$query}->{$SyntenyCount->{$sbjct}->{$query}++} = "$sStart,$sEnd:$qStart,$qEnd:$bits";
	}
	close(BLAST);
}

sub parse_genbank_header {
	my $sbjct = shift;
	$sbjct =~ s/gi\|\d+\|//;
	$sbjct =~ s/ref\|//;
	$sbjct =~ s/\|$//;
	$sbjct =~ s/\.\d+$//;
	return $sbjct;
}

# QLEN and SLEN are like [["seq1", 1000], ["seq2", 2000"]]

sub prettify {
	my $length = shift;
	my %Keys = (10**12 => "Tb", 10**9 => "Gb", 10**6 => "Mb", 10**3 => "Kb", 1 => "bp");
	ATTEMPT:foreach my $round (sort {$b <=> $a} keys %Keys) {
		if ($length > $round) {
			return sprintf("%.1f%s", $length/$round, $Keys{$round});
		}
	}
	return "${length}bp";
}

sub draw_dot_plot {
	my $this = shift;
	my %params = @_;

	my $sort_by_hits = 1;

	# running sum of sizes, used for offset

	my $total_query_length  = 0;
	my $total_sbjct_length  = 0;

	# We have to make a hash of length and offset for display
	# Sometimes sequences would be grouped together so we display
	# a different set of sequences than the ones from BLAST report.
	# These hashes are useful to do that, while keeping the original
	# in tact.

	# query/sbjct lengths and offsets in concat string

	my %SOffset;
	my %QOffset;

	# length information

	my $QLength = $this->{QLENGTH};
	my $SLength = $this->{SLENGTH};

	# synteny information

	my $Synteny = $this->{SYNTENY};
	my $SyntenyCount = $this->{SYNTENY_COUNT};

	# orientation of the query against subject

	my $QOrientation = $this->{QORIENTATION};

	# Get concatenated subject and the start positions of each subject in it

	my $subject_total_length = 0;
	foreach my $sbjct (sort {$SLength->{$b} <=> $SLength->{$a}} keys %$SLength) {
		$SOffset{$sbjct} = $total_sbjct_length;
		$total_sbjct_length += $SLength->{$sbjct};
	}

	# Sort queries (X-axis) by position in Y-axis, to make the dotplot look nice

	if ($sort_by_hits == 1) {
		my %QueryStart;
		QUERY:foreach my $query (keys %$QLength) {
			my $max_bits = 0;
			$QOrientation->{$query} = 1;
			# get the beginning of the highest bitscoring hsp in the concatenated string
			SBJCT:foreach my $sbjct (keys %$SLength) {
				for (my $i=0; $i<$SyntenyCount->{$sbjct}->{$query}; $i++) {
					my $map = $Synteny->{$sbjct}->{$query}->{$i};
					my ($s, $q, $bits) = split(/:/, $map);
					my ($start_x, $end_x) = map {$SOffset{$sbjct} + $_} split(/,/, $s);

					# choose the best hit's start point
					if (!defined($QueryStart{$query}) || $bits > $max_bits) {
						($QueryStart{$query}) = sort {$a <=> $b} ($start_x, $end_x);
						$max_bits = $bits;
						my ($qb, $qe) = split(/,/, $q);
						$QOrientation->{$query} = ($qb < $qe)?(1):(-1);
					}
				}
			}
		}

		# Sort by start position

		foreach my $query (sort {($QueryStart{$a} || 10**10) <=> ($QueryStart{$b} || 10**10)} keys %$QLength) {
			my $len = $QLength->{$query};
			$QOffset{$query} = $total_query_length;
			$total_query_length += $len;
		}

	}

	# external units

	my $START_X     = 0;
	my $END_X       = $total_sbjct_length;
	my $XSPAN       = $END_X-$START_X+1;
	my $MAJOR_XTIC  = 10**(int(log($XSPAN-1)/log(10)));
	my $MINOR_XTIC  = $MAJOR_XTIC/10;

	my $START_Y     = 0;
	my $END_Y       = $total_query_length;
	my $YSPAN       = $END_Y-$START_Y+1;
	my $MAJOR_YTIC  = 10**(int(log($YSPAN-1)/log(10)));
	my $MINOR_YTIC  = $MAJOR_YTIC/10;

	my $X_SCALE_DOWN  = $this->{X_SCALE_DOWN} || nearest(10, $XSPAN/1000);
	my $Y_SCALE_DOWN  = $this->{Y_SCALE_DOWN} || $X_SCALE_DOWN;

	# 2px is optimal if you let the script decide automatically using the lines above. Otherwise, you are on your own!

	my $line_width    = 4; 

	my $FH            = $this->{OUT_FILE_HANDLE};
	my $DISPLAY_LABEL = $this->{DISPLAY_LABEL};
	my $DISPLAY_TICS  = $this->{DISPLAY_TICS};

	########
	# draw the labels for the axes
	########

	if ($DISPLAY_TICS == 1) {
		$GLOBAL_X_SCALE_DOWN = $X_SCALE_DOWN;
		print $FH "<g transform=\"translate(0, -15)\">\n";
		my $label = sprintf("<tspan style=\"font-family:sans; font-style:italic;\">%s</tspan> (%dbp)", $this->{SNAME}, $total_sbjct_length);
		$this->draw_dot_plot_labels(OUT_FILE_HANDLE => $FH, LABEL => $label, AXIS => "x", ORIGIN => $START_X, END => $END_X, MAJOR_TIC => $MAJOR_XTIC, MINOR_TIC => $MINOR_XTIC, OFFSETS => \%SOffset, LENGTHS => $SLength, MIN_LABEL_LENGTH => 250000, TICS_ONLY => !$DISPLAY_LABEL);
		print $FH "</g>\n";

		$GLOBAL_X_SCALE_DOWN = $Y_SCALE_DOWN;
		print $FH "<g transform=\"translate(-15, 0) rotate(90 0,0)\">\n";
		$label = sprintf("<tspan style=\"font-family:sans; font-style:italic;\">%s</tspan> (%dbp)", $this->{QNAME}, $total_query_length);
		$this->draw_dot_plot_labels(OUT_FILE_HANDLE => $FH, LABEL => $label, AXIS => "y", ORIGIN => $START_Y, END => $END_Y, MAJOR_TIC => $MAJOR_YTIC, MINOR_TIC => $MINOR_YTIC, OFFSETS => \%QOffset, LENGTHS => $QLength, MIN_LABEL_LENGTH => 10000, TICS_ONLY => !$DISPLAY_LABEL);
		print $FH "</g>\n";
	}

	# from mySVG

	$POINT_PRECISION     = 0.00001;
	$GLOBAL_X_SCALE_DOWN = $X_SCALE_DOWN;
	$GLOBAL_Y_SCALE_DOWN = $Y_SCALE_DOWN;

	########
	# Draw the bounding rectangle
	# The graph goes from (0,0) to ($XSPAN, $YSPAN) since YSPAN = [$START_Y, $END_Y]
	########

	draw_rect_scaled($FH, 0, 0, $XSPAN, $YSPAN, "stroke:black;stroke-width:1px;fill:white;");

	########
	# Draw minor lines for the contig boundaries
	########

	foreach my $sbjct (keys %$SLength) {
		draw_line_scaled($FH, $SOffset{$sbjct}-$START_X, 0, $SOffset{$sbjct}-$START_X, $YSPAN, "stroke:rgb(64,64,64); stroke-width:@{[$line_width/20]}px");
	}
	foreach my $query (keys %$QLength) {
		draw_line_scaled($FH, 0, $QOffset{$query}-$START_Y, $XSPAN, $QOffset{$query}-$START_Y, "stroke:rgb(64,64,64); stroke-width:@{[$line_width/20]}px");
	}

	########
	# Draw genome dot plot
	########

	print $FH "<g id=\"genome_dot_plot\">\n";
	foreach my $sbjct (keys %$Synteny) {
		foreach my $query (keys %{$Synteny->{$sbjct}}) {
			for (my $i=0; $i<$SyntenyCount->{$sbjct}->{$query}; $i++) {
				my $map = $Synteny->{$sbjct}->{$query}->{$i};
				my ($s, $q, $bits) = split(/:/, $map);
				my ($start_x, $end_x) = map {$SOffset{$sbjct} + $_} split(/,/, $s);

				# if the query is reverse, assume reverse complement of the sequence to make it look better
				my ($qb, $qe) = split(/,/, $q);
				if ($QOrientation->{$query} == -1) {
					($qb, $qe) = map {$QLength->{$query} - $_ + 1} ($qb, $qe);
				}

				my ($start_y, $end_y) = map {$QOffset{$query} + $_} ($qb, $qe);

				#my $color = ($qStart < $qEnd)?"0,255,0":"255,0,0";
				my $color = "0,0,0";
				draw_line_scaled($FH, $start_x-$START_X, $start_y-$START_Y, $end_x-$START_X, $end_y-$START_Y, "stroke:rgb($color);stroke-width:${line_width}px");
			}
		}
	}
	print $FH "</g>\n";

	########
	# Draw gene dot plot
	########

	print $FH "<g id=\"gene_dot_plot\">\n";
	my $ProteinSynteny      = $this->{PROTEIN_SYNTENY}; 
	my $ProteinSyntenyCount = $this->{PROTEIN_SYNTENY_COUNT}; 
	foreach my $sbjct (keys %$ProteinSynteny) {
		foreach my $query (keys %{$ProteinSynteny->{$sbjct}}) {
			for (my $i=0; $i<$ProteinSyntenyCount->{$sbjct}->{$query}; $i++) {
				my $map = $ProteinSynteny->{$sbjct}->{$query}->{$i};
				my ($s, $q, $bits) = split(/:/, $map);
				my ($start_x, $end_x) = map {$SOffset{$sbjct} + $_} split(/,/, $s);

				# if the query is reverse, assume reverse complement of the sequence to make it look better
				my ($qb, $qe) = split(/,/, $q);
				if ($QOrientation->{$query} == -1) {
					($qb, $qe) = map {$QLength->{$query} - $_ + 1} ($qb, $qe);
				}

				my ($start_y, $end_y) = map {$QOffset{$query} + $_} ($qb, $qe);

				#my $color = ($qStart < $qEnd)?"0,255,0":"255,0,0";
				my $color = "200,50,100";
				draw_line_scaled($FH, $start_x-$START_X, $start_y-$START_Y, $end_x-$START_X, $end_y-$START_Y, "stroke:rgb($color);stroke-width:${line_width}px");
			}
		}
	}
	print $FH "</g>\n";

}

=begin DEPRECATED

# QLEN and SLEN are like [["seq1", 1000], ["seq2", 2000"]]
sub draw_dot_plot_old {
	my %params = @_;

	my $FH            = $params{OUT_FILE_HANDLE};
	my $QLENGTH_FILE  = $params{QLENFILE};
	my $SLENGTH_FILE  = $params{SLENFILE};
	my $QLEN          = $params{QLEN};
	my $SLEN          = $params{SLEN};
	my $BLAST_FILE    = $params{BLAST_FILE};
	my $DISPLAY_LABEL = $params{DISPLAY_LABEL};
	my $DISPLAY_TICS  = $params{DISPLAY_TICS};
	my $X_SCALE_DOWN  = $params{X_SCALE_DOWN};
	my $Y_SCALE_DOWN  = $params{Y_SCALE_DOWN};

	my ($query_name, $sbjct_name);
	my $sort_by_hits = 1;

	# running sum of sizes, used for offset

	my $query_size  = 0;
	my $sbjct_size  = 0;

	# query/sbjct lengths and offsets in concat string

	my @QueryLength;
	my @SbjctLength;
	my %QueryOffset;
	my %SbjctOffset;

	# make a hash of length and offset for display
	# Sometimes sequences would be grouped together so we display
	# a different set of sequences than the ones from BLAST report.
	# These hashes are useful to do that, while keeping the original
	# in tact.

	my %QueryDisplayLength;
	my %QueryDisplayOffset;
	my %SbjctDisplayLength;
	my %SbjctDisplayOffset;

	# Get subject length information

	if ($SLEN) {
		@SbjctLength = @$SLEN;
		my $sum = 0;
		for (my $i=0; $i<=$#SbjctLength; $i++) {
			my ($seq, $length) = map {$SbjctLength[$i][$_]} (0,1);
			$SbjctOffset{$seq} = $sbjct_size;
			$SbjctDisplayLength{$seq} = $length;
			$sbjct_size += ($length+5*$X_SCALE_DOWN);
			$sbjct_name = $seq;
		}
	} elsif ($SLENGTH_FILE) {
		my $residue_length = 1;
		if ($SLENGTH_FILE =~ /protein/) {
			$residue_length = 3;
		}
		open(LEN, "<$SLENGTH_FILE") || die "Cannot open $SLENGTH_FILE: $!";
		while (<LEN>) {
			chomp();
			my ($seq, $length) = split(/\t/);
			$length *= $residue_length;
			$SbjctOffset{$seq} = $sbjct_size;
			$SbjctDisplayLength{$seq} = $length;
			$sbjct_size += ($length+5*$X_SCALE_DOWN);
		}
		close(LEN);
		($sbjct_name) = map {my ($n) = m/(.*)\.len/; $n =~ s/^(.*)_/$1. /; $n} ($SLENGTH_FILE);
	} else {
		die "Sbjct lengths not defined";
	}
	%SbjctDisplayOffset = %SbjctOffset;

	# Get query length information

	if ($QLEN) {
		@QueryLength = @$QLEN;
		for (my $i=0; $i<=$#QueryLength; $i++) {
			my ($seq, $length) = map {$QueryLength[$i][$_]} (0,1);
			$QueryOffset{$seq} = $query_size;
			$QueryDisplayLength{$seq} = $length;
			$query_name = $seq;
			$query_size += $length;
		}
	} elsif ($QLENGTH_FILE) {
		my $residue_length = 1;
		open(LEN, "<$QLENGTH_FILE") || die "Cannot open $QLENGTH_FILE: $!";
		if ($QLENGTH_FILE =~ /protein/) {
			$residue_length = 3;
		}
		while (<LEN>) {
			chomp();
			my ($seq, $length) = split(/\t/);
			$length *= $residue_length;
			$QueryOffset{$seq} = $query_size;
			$QueryDisplayLength{$seq} = $length;
			$query_size += $length;
		}
		close(LEN);
		($query_name) = map {my ($n) = m/(.*)\.len/; $n =~ s/^(.)_/$1. /; $n} ($QLENGTH_FILE);
	} else {
		die "Query lengths not defined";
	}


	# Sort queries (X-axis) by position in Y-axis, to make the dotplot look nice

	# Get start position for each query

	my %QueryStart;
	if ($sort_by_hits == 1) {
		my %Bits;
		open(BLAST, "<$BLAST_FILE") || die "Cannot open $BLAST_FILE: $!";
		QRY:while (<BLAST>) {
			chomp();
			my @w = split(/\t/);
			my ($query, $sbjct, $bits, $similarity, $qStart, $qEnd, $sStart, $sEnd) = ($w[0], $w[1], $w[4], $w[11], $w[17], $w[18], $w[20], $w[21]);
			if ($similarity < $params{SIMILARITY} || $bits < $params{BITS}) {
				next QRY;
			}
			if ($sEnd - $sStart < $params{MINALIGNSLEN}) {
				next QRY;
			}
			$sbjct =~ s/ref\|//;
			$sbjct =~ s/\|$//;
			$sbjct =~ s/\.\d+$//;
			if (!defined($SbjctOffset{$sbjct})) {
				next QRY;
			}
			my ($start_x, $end_x) = map {$SbjctOffset{$sbjct} + $_} ($sStart, $sEnd);

			# choose the best hit's start point
			if (!defined($QueryStart{$query}) || $bits > $Bits{$query}) {
				($QueryStart{$query}) = sort {$a <=> $b} ($start_x, $end_x);
				$Bits{$query} = $bits;
			}
		}
		close(BLAST);

		# Sort by start position

		foreach my $query (sort {($QueryStart{$a} || 10**10) <=> ($QueryStart{$b} || 10**10)} keys %QueryDisplayLength) {
			my $len = $QueryDisplayLength{$query};
			$QueryOffset{$query} = $query_size;
			$query_size += $len;
		}

	}

	%QueryDisplayOffset = %QueryOffset;

	# external units

	my $START_X     = $params{START_X} || 0;
	my $END_X       = $params{END_X}   || $sbjct_size;
	my $XSPAN       = $END_X-$START_X+1;
	my $MAJOR_XTIC  = 10**(floor(log($XSPAN-1)/log(10)));
	my $MINOR_XTIC  = $MAJOR_XTIC/10;

	my $START_Y     = 0;
	my $END_Y       = $query_size;
	my $YSPAN       = $END_Y-$START_Y+1;
	my $MAJOR_YTIC  = 10**(floor(log($YSPAN-1)/log(10)));
	my $MINOR_YTIC  = $MAJOR_YTIC/10;

	# if fragments have to be grouped together
	if ($params{REMAP_DISPLAY} == 1) {
		for (my $i=0; $i<=$#SbjctLength; $i++) {
			my ($seq, $length) = map {$SbjctLength[$i][$_]} (0,1);
			my ($disp) = $seq =~ /gnl\|GLV\|(.{6})/;
			$SbjctDisplayLength{$disp} += $length;
		}
		my $offset = 0;
		for my $sub (sort {$a cmp $b} keys %SbjctDisplayLength) {
			$SbjctDisplayOffset{$sub} = $offset;
			$offset += $SbjctDisplayLength{$sub};
		}
	}

	# draw the labels for the axes

	if ($DISPLAY_TICS == 1) {
		$GLOBAL_X_SCALE_DOWN = $X_SCALE_DOWN;
		print $FH "<g transform=\"translate(0, -15)\">\n";
		$this->draw_dot_plot_labels(OUT_FILE_HANDLE => $FH, LABEL => "$sbjct_name (${sbjct_size}bp)", AXIS => "x", ORIGIN => $START_X, END => $END_X, MAJOR_TIC => $MAJOR_XTIC, MINOR_TIC => $MINOR_XTIC, OFFSETS => \%SbjctDisplayOffset, LENGTHS => \%SbjctDisplayLength, TICS_ONLY => !$DISPLAY_LABEL);
		print $FH "</g>\n";

		$GLOBAL_X_SCALE_DOWN = $Y_SCALE_DOWN;
		print $FH "<g transform=\"translate(-15, 0) rotate(90 0,0)\">\n";
		$this->draw_dot_plot_labels(OUT_FILE_HANDLE => $FH, LABEL => "<tspan style=\"font-style:italic;\">$query_name</tspan> (${query_size}bp)", AXIS => "y", ORIGIN => $START_Y, END => $END_Y, MAJOR_TIC => $MAJOR_YTIC, MINOR_TIC => $MINOR_YTIC, OFFSETS => \%QueryDisplayOffset, LENGTHS => \%QueryDisplayLength, TICS_ONLY => !$DISPLAY_LABEL);
		print $FH "</g>\n";
	}

	# from mySVG

	$POINT_PRECISION     = 0.00001;
	$GLOBAL_X_SCALE_DOWN = $X_SCALE_DOWN;
	$GLOBAL_Y_SCALE_DOWN = $Y_SCALE_DOWN;

	# The graph goes from (0,0) to ($XSPAN, $YSPAN) since YSPAN = [$START_Y, $END_Y]

	draw_rect_scaled($FH, 0, 0, $XSPAN, $YSPAN, "stroke:black;stroke-width:1px;fill:white;");

	# draw reads

	#my ($line_width) = scale_down_y($opt_salnlen);
	my $line_width = 3;

	open(BLAST, "<$BLAST_FILE") || die "Cannot open $BLAST_FILE: $!";
	QRY:while (<BLAST>) {
		chomp();
		my @w = split(/\t/);
		my ($query, $sbjct, $bits, $similarity, $qStart, $qEnd, $sStart, $sEnd) = ($w[0], $w[1], $w[4], $w[11], $w[17], $w[18], $w[20], $w[21]);
		if ($similarity < $params{SIMILARITY} || $bits < $params{BITS}) {
			next QRY;
		}
		if ($sEnd - $sStart < $params{MINALIGNSLEN}) {
			next QRY;
		}
		$sbjct =~ s/ref\|//;
		$sbjct =~ s/\|$//;
		$sbjct =~ s/\.\d+$//;
		if (!defined($SbjctOffset{$sbjct}) || !defined($QueryOffset{$query})) {
			next QRY;
		}
		my ($start_x, $end_x) = map {$SbjctOffset{$sbjct} + $_} ($sStart, $sEnd);
		my ($start_y, $end_y) = map {$QueryOffset{$query} + $_} ($qStart, $qEnd);
		if ($end_x < $START_X || $start_x > $END_X) {
			next QRY;
		}
		if ($end_y < $START_Y || $start_y > $END_Y) {
			next QRY;
		}
		#my $color = ($qStart < $qEnd)?"0,255,0":"255,0,0";
		my $color = "0,0,0";
		draw_line_scaled($FH, $start_x-$START_X, $start_y-$START_Y, $end_x-$START_X, $end_y-$START_Y, "stroke:rgb($color);stroke-width:${line_width}px");
	}
	close(BLAST);

}

=cut

sub draw_dot_plot_labels {

	use Math::Round qw(nearest);

	my $this           = shift;
	my %params         = @_;
	my $FH             = $params{OUT_FILE_HANDLE};
	my $SEQUENCE_LABEL = $params{LABEL};
	my $AXIS           = $params{AXIS};
	my $ORIGIN         = $params{ORIGIN};
	my $END            = $params{END};
	my $MAJOR_TIC      = $params{MAJOR_TIC};
	my $MINOR_TIC      = $params{MINOR_TIC};
	my $TICS_ONLY      = $params{TICS_ONLY};
	my $SequenceLength = $params{LENGTHS};
	my $SequenceOffset = $params{OFFSETS};
	my $MIN_LABEL_LENGTH = $params{MIN_LABEL_LENGTH};

	my $span   = $END-$ORIGIN+1;
	my $order  = ($AXIS eq "x")?1:-1;
	my $font_size;
	my $level;
	
	######
	# sequence scale line
	######

	$level = 0;
	draw_line_scaled_x($FH, 0, $level, $span, $level, "stroke:rgb(0,0,0);stroke-width:2px");
	for (my $i=$ORIGIN; $i <= $END; $i+=$MAJOR_TIC) {
		draw_line_scaled_x($FH, $i-$ORIGIN, $level-9, $i-$ORIGIN, $level+9, "stroke:rgb(0,0,0);stroke-width=1px");
	}
	for (my $i=$ORIGIN; $i <= $END; $i+=$MINOR_TIC) {
		draw_line_scaled_x($FH, $i-$ORIGIN, $level-5, $i-$ORIGIN, $level+5, "stroke:rgb(0,0,0);stroke-width=0.25px");
	}
	$level -= $order*10; # move the width of what was drawn now

	######
	# sequence scale tic labels
	######

	$font_size = 14;
	if ($order == -1) { # if drawing reverse, move the font size before drawing (0.8 since font size is not centered vertically)
		$level -= $order*int($font_size*0.8);
	}
	for (my $i=$ORIGIN; $i <= $END; $i+=$MAJOR_TIC) {
		draw_text_scaled_x($FH, $i-$ORIGIN, $level, nearest(0.01, $i/1000000)."Mb", ("font-size" => "${font_size}px"));
	}

	if ($order == 1) { # if drawing not reverse, move the font size after drawing
		$level -= $order*($font_size+5); # label for the seq scale
	}

	if ($TICS_ONLY == 1) {
		return;
	}

	######
	# sequence label
	######

	$font_size = 25;
	if ($order == -1) { # if drawing reverse, move the font size before drawing (0.8 since font size is not centered vertically)
		$level -= $order*int($font_size*0.8+5);
	}
	draw_text($FH, 0, $level, "$SEQUENCE_LABEL", ("font-size" => "${font_size}px"));
	$level -= $order*10;

	if ($order == 1) {
		$level -= $order*($font_size);
	}

	######
	# sequence fragments split lines
	######

	my $sign  = 1;
	for my $seq (sort {$SequenceOffset->{$a} <=> $SequenceOffset->{$b}} keys %$SequenceOffset) {
		draw_line_scaled_x($FH, $SequenceOffset->{$seq}-$ORIGIN, $level+$sign*2, $SequenceOffset->{$seq}+$SequenceLength->{$seq}-$ORIGIN, $level+$sign*2, "stroke:rgb(0,0,0);stroke-width:2px");
		$sign *= -1;
	}
	$level -= $order*10; # gap between labels and the line(s)

	######
	# sequence fragments split labels
	######

	$font_size = 20;
	my $rotate = ($AXIS eq "x")?90:270;
	for my $seq (sort {$SequenceOffset->{$a} <=> $SequenceOffset->{$b}} keys %$SequenceOffset) {

		# sometimes sequences are so short that you cannot write them, even rotated
		# so we should reduce the font size for them, so that we can zoom in the svg file to see who these are

		my $local_font_size = $font_size;
		if ($SequenceLength->{$seq}/$GLOBAL_X_SCALE_DOWN < $font_size/2) {
			$local_font_size = 2+int($SequenceLength->{$seq}/$GLOBAL_X_SCALE_DOWN);
		}
		my $offset = ($AXIS eq "x")?0:$local_font_size;
		my $display = $seq;
		$display =~ s/ref\|//;
		$display =~ s/\|$//;

		#beg hacks
		$display =~ s/scaffold0+//;
		$display =~ s/scf71800000118/scaffold_/;
		$display =~ s/supercontig_/scaffold_/;
		#end hacks

		$display .= sprintf(" (%s)", prettify($SequenceLength->{$seq}));
		if (defined($this->{QORIENTATION}->{$seq}) && $this->{QORIENTATION}->{$seq} == -1) {
			$display .= " RC";
		}
		draw_text_scaled_x($FH, $SequenceOffset->{$seq}-$ORIGIN+$offset*$GLOBAL_X_SCALE_DOWN, $level, $display, ("font-size" => "${local_font_size}px", "style" => "text-anchor:end;", rotate => $rotate));
	}
	$level -= $order*($font_size+5); # for the query seq label
}

1;
