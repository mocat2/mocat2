package Smash::Utils::Plot;

use strict;
use warnings;
use POSIX;
use Math::Round;
use Smash::Utils::SVG qw(:all);
use base "Exporter";

our @EXPORT_OK = qw(draw_contig_assembly draw_fragment_recruitment draw_genome_line);
our %EXPORT_TAGS = (all => [@EXPORT_OK]);

#sub draw_y_labels {
#	my ($sequence_name, $sequence_size);
#	my %SequenceLength;
#	my %SequenceOffset;
#	my @levels = (-90, -50, -20);
#	######
#	# y axis guide
#	# use internal (absolute) coordinates in x axis, and external coordinates in y axis
#	######
#
#	######
#	# sequence label
#	######
#
#	$level = $levels[0];
#	draw_text_scaled_y($FH, $level, 0, "$sequence_name (${sequence_size}bp)", 90);
#
#	######
#	# sequence split
#	######
#
#	$level = $levels[1];
#	$sign = 1;
#	for my $seq (sort {$SequenceOffset{$a} <=> $SequenceOffset{$b}} keys %SequenceOffset) {
#		draw_line_scaled_y($FH, $level+$sign*2, $SequenceOffset{$seq}-$START_Y, $level+$sign*2, $SequenceOffset{$seq}+$SequenceLength{$seq}-$START_Y, "stroke:rgb(0,0,0);stroke-width:2px");
#		if ($SequenceLength{$seq} > 250000) {
#			my $display = $seq;
#			$display =~ s/scaffold0+//;
#			draw_text_scaled_y($FH, $level-20, $SequenceOffset{$seq}-$START_Y, $display, 90);
#		}
#		$sign *= -1;
#	}
#
#	$level = $levels[2];
#
#	######
#	# query line
#	######
#
#	draw_line_scaled_y($FH, $level, 0, $level, $YSPAN, "stroke:rgb(0,0,0);stroke-width:2px");
#
#	######
#	# ytics and ytic labels
#	######
#
#	for (my $i=$START_Y; $i <= $END_Y; $i+=$MAJOR_YTIC) {
#		draw_line_scaled_y($FH, $level-9, $i-$START_Y, $level+9, $i-$START_Y, "stroke:rgb(0,0,0);stroke-width=1px");
#		draw_text_scaled_y($FH, $level-20, $i-$START_Y, nearest(0.01, $i/1000000)."Mb", 90);
#	}
#	for (my $i=$START_Y; $i <= $END_Y; $i+=$MINOR_YTIC) {
#		draw_line_scaled_y($FH, $level-5, $i-$START_Y, $level+5, $i-$START_Y, "stroke:rgb(0,0,0);stroke-width=0.25px");
#	}
#}

sub draw_contig_assembly {
	use Smash::Core;
	use Smash::Utils::SVG qw(:all);

	my ($FH, $collection, $contig_id) = @_;
	my ($smash, $dbh, $sth);
	my ($contig_name, $contig_length);

	$smash = new Smash::Core(COLLECTION => $collection);
	$smash->init();
	$dbh = $smash->get_db_handle;

	# get contig name and length

	$sth = $dbh->prepare("SELECT external_id, length FROM contig WHERE contig_id=?");
	$sth->execute($contig_id);
	($contig_name, $contig_length) = $sth->fetchrow_array();
	$sth->finish();

	$sth = $dbh->prepare("SELECT read_id, start, end FROM contig2read WHERE contig_id=?");

	# level tracking
	my $MIN_GAP = 100;
	my $max_level=1;
	my @Used;
	my $width=2;

	# init level tracking

	map {$Used[$max_level][$_] = 0} 0..$contig_length;

	# draw the contig

	draw_line_scaled_x($FH, 0, 0, $contig_length, 0);

	# draw the reads

	$sth->execute($contig_id);
	while (my ($read, $read_start, $read_end) = $sth->fetchrow_array()) {
		my $level;
		LEVEL:for ($level=1; $level <= $max_level; $level++) {
			if (level_available(\@Used, $level, max($read_start+1-$MIN_GAP, 0), min($read_end+1+$MIN_GAP, $contig_length))) {
				last LEVEL;
			}
		}

		# make new level if necessary

		if ($level > $max_level) {
			$max_level = $level;
			map {$Used[$level][$_] = 0} 0..$contig_length;
		}

		# place the read in $level

		draw_line_scaled_x($FH, $read_start+1, $level*$width, $read_end+1, $level*$width);
		map {$Used[$level][$_] = 1} ($read_start+1)..($read_end+1);
	}
	$sth->finish();
	$smash->finish();
}

sub level_available {
	my ($UsedRef, $level, $begin, $end) = @_;
	for (my $i=$begin; $i<=$end; $i++) {
		if ($UsedRef->[$level][$i] == 1) {
			return 0;
		}
	}
	return 1;
}

sub draw_genome_line {
	my ($FH, $y, $genome_label, $genome_size, $start, $end, $major_tic, $minor_tic) = @_;
	# genome label

	my $level;
	my $font_size;
	my ($DISPLAY_LABEL, $DISPLAY_TICS) = (1,1);

	print $FH "<g id=\"genome_scaffold\">\n";
	if ($DISPLAY_LABEL == 1) {
		$level = $y - 25;
		$font_size = 15;
		draw_text_scaled_x($FH, 0, $level, "$genome_label (${genome_size}bp)", ('font-size' => $font_size));
	}

	# genome line and xtics

	$level = $y;
	$font_size = 10;
	draw_line_scaled_x($FH, 1, $level, $end-$start+1, $level, "stroke:rgb(0,0,0);stroke-width:2px");
	print $FH "</g>\n";

	if ($DISPLAY_TICS == 1) {
		my ($unit, $factor);
		if ($major_tic >= 100000000) {
			$unit = "Gb";
			$factor = 1000000000;
		} elsif ($major_tic >= 100000) {
			$unit = "Mb";
			$factor = 1000000;
		} else {
			$unit = "Kb";
			$factor = 1000;
		}
		print $FH "<g id=\"ytics\">\n";
		for (my $i=$start; $i <= $end; $i+=$minor_tic) {
			if ($i % $major_tic == 0) {
				draw_line_scaled_x($FH, $i-$start+1, $level-9, $i-$start+1, $level+9, "stroke:rgb(0,0,0);stroke-width=1px");
				draw_text_scaled_x($FH, $i-$start+1, $level-10, nearest(0.01, $i/$factor).$unit, ('font-size' => $font_size));
			} else {
				draw_line_scaled_x($FH, $i-$start+1, $level-5, $i-$start+1, $level+5, "stroke:rgb(0,0,0);stroke-width=0.25px");
			}
		}
		print $FH "</g>\n";
	}
}

sub draw_fragment_recruitment {
	my %params = @_;

	my $FH            = $params{OUT_FILE_HANDLE};
	my $ALIGN_FILE    = $params{ALIGN_FILE};
	my $ALIGN_FORMAT  = $params{ALIGN_FORMAT};
	my $DISPLAY_LABEL = $params{DISPLAY_LABEL};
	my $DISPLAY_TICS  = $params{DISPLAY_TICS};
	my $X_SCALE_DOWN  = $params{X_SCALE_DOWN};
	my $Y_SCALE_DOWN  = $params{Y_SCALE_DOWN};
	my $GENOME_SIZE   = $params{GENOME_SIZE};
	my $GENOME_NAME   = $params{GENOME_NAME};
	my $GENOME_LABEL  = $params{GENOME_LABEL};
	my $COLORS        = $params{COLORS};

	my $START_POS     = $params{START_X} || 0;
	my $END_POS       = $params{END_X}   || $GENOME_SIZE;
	my $XSPAN         = $END_POS-$START_POS+1;
	my $MAJOR_XTIC    = 10**floor(log($XSPAN-1)/log(10));
	my $MINOR_XTIC    = $MAJOR_XTIC/10;


	# external units - usually percent identity

	my $END_Y         = $params{SIMILARITY} || 70;
	my $START_Y       = 100;
	my $YSPAN         = abs($END_Y-$START_Y);
	my $MAJOR_YTIC    =  5;
	my $MINOR_YTIC    =  $MAJOR_YTIC/5;

	# The graph goes from (0,0) to ($XSPAN, $YSPAN) since YSPAN = [$START_Y, $END_Y]

	$POINT_PRECISION     = 0.1;
	$GLOBAL_X_SCALE_DOWN = $X_SCALE_DOWN;
	$GLOBAL_Y_SCALE_DOWN = $Y_SCALE_DOWN;

	# Colors

	my %Colors  = %$COLORS;
	my $percent_precision = 1;                    # pi values will be rounded with this precision
	my $jitter_level      = 0.1*$percent_precision;        # Add jitter that is within +- 10% of precision
	my $line_stroke_width = 1;

	my $level = -$jitter_level*$YSPAN;
	my $font_size;

	# use internal (absolute) coordinates in y axis, and external coordinates in x axis

	# draw genome sequence line with tics

	draw_genome_line($FH, -$jitter_level/$Y_SCALE_DOWN - 15, $GENOME_LABEL, $GENOME_SIZE, $START_POS, $END_POS, $MAJOR_XTIC, $MINOR_XTIC);

	# draw bounding box for frag-plot

	draw_rect_scaled_y($FH, -5, 0-1.1*$jitter_level, $XSPAN/$X_SCALE_DOWN+10, $YSPAN+2.2*$jitter_level, "stroke:black;stroke-width:1px;fill:white;");

	# ytics
	# mandatory, since this is important
	# use internal (absolute) coordinates in x axis, and external coordinates in y axis

	my $width;
	print $FH "<g id=\"similarity_scale-y-axis\">\n";

	# minor tics

	$width = -7;
	$font_size = 5;
	for (my $i=$START_Y; $i >= $END_Y; $i-=$MINOR_YTIC) {
		my $ylevel = nearest($percent_precision, ($START_Y-$i));
		draw_rect_scaled_y($FH, $width-5, $ylevel-$jitter_level, abs($width), 2*$jitter_level, "fill:gray;");
		#draw_text_scaled_y($FH, $width-35, $ylevel, "$i%", ('font-size' => $font_size));
	}

	# major tics

	$width = -15;
	$font_size = 10;
	for (my $i=$START_Y; $i >= $END_Y; $i-=$MAJOR_YTIC) {
		my $ylevel = nearest($percent_precision, ($START_Y-$i));
		draw_rect_scaled_y($FH, $width-5, $ylevel-$jitter_level, abs($width), 2*$jitter_level, "fill:black;");
		draw_text_scaled_y($FH, $width-35, $ylevel, "$i%", ('font-size' => $font_size));
	}
	print $FH "</g>\n";

	# draw reads

	my @wu_indices   = (0,1,2,4,5,10,11,17,18,20,21);
	my @ncbi_indices = (0,1,10,11,11,2,2,6,7,8,9);

	print $FH "<g id=\"fragments\">\n";
	my $minslen    = $params{MINALIGNSLEN};
	my $similarity = $params{SIMILARITY};
	open(ALIGN, "<$ALIGN_FILE") || die "Cannot open $ALIGN_FILE: $!";
	my $Count = {};
	my $Drawn = {};
	if ($ALIGN_FORMAT eq "SAM") { 
		QRY:while (<ALIGN>) {
			chomp();
			my @words = split(/\t/);
			#my ($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send);
			my ($query, $qlen, $subject, $sstart, $send, $match, $edit) = @words;
			die "Error: SAM SUMMARY requires exactly 7 columns. Check the following line:\n$_\n" unless @words == 7;

			if ($subject ne $GENOME_NAME) {
				next QRY;
			}
			my $pcid = 100*(1 - $edit/($edit+$match));
			if ($pcid < $similarity) {
				next QRY;
			}
			#if ($qlen < $minslen || $qlen > $minslen+1) {}
			if ($qlen < $minslen) {
				next QRY;
			}
			if ($send < $START_POS || $sstart > $END_POS) {
				next QRY;
			}
			#my $sample = join('.', (split(/\./, $query))[0,1]);
			#my $color  = $Colors{$sample} || "#000000";
			my $color = "#000000";
			my $ylevel = nearest($percent_precision, ($START_Y-$pcid));
			$Count->{$ylevel}++;
			my $jitter = rand()*$jitter_level*((-1)**(int(rand(2)))); 
			$ylevel += $jitter;
			$ylevel = nearest($POINT_PRECISION, $ylevel/$GLOBAL_Y_SCALE_DOWN);
			if (!$Drawn->{$sstart}->{$ylevel}->{$send}) { # No point drawing it again
				draw_line_scaled_x($FH, $sstart-$START_POS+1, $ylevel, $send-$START_POS+1, $ylevel, "class:f");
				$Drawn->{$sstart}->{$ylevel}->{$send} = 1;
			}
		}
		close(ALIGN);
		print $FH "</g>\n";
use Smash::Utils::MatrixIO qw(:all);
write_two_column_hash_sorted(\*STDOUT, $Count);
	} else {
		QRY:while (<ALIGN>) {
			chomp();
			my @words = split(/\t/);
			my ($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send);
			if ($ALIGN_FORMAT eq "WU") {
				die "Error: WU BLAST requires at least 22 columns. Check the following line:\n$_\n" unless @words > 21;
				($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send) = @words[@wu_indices];
			} elsif ($ALIGN_FORMAT eq "NCBI") {
				die "Error: NCBI BLAST requires exactly 12 columns. Check the following line:\n$_\n" unless @words == 12;
				($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send) = @words[@ncbi_indices];
			}

			if ($subject ne $GENOME_NAME) {
				next QRY;
			}
			if ($pcid < $params{SIMILARITY} || $bits < $params{BITS} || $E > $params{E}) {
				next QRY;
			}
			if (abs($send - $sstart + 1) < $params{MINALIGNSLEN}) {
				next QRY;
			}
			if ($send < $START_POS || $sstart > $END_POS) {
				next QRY;
			}
			my $sample = join('.', (split(/\./, $query))[0,1]);
			my $color  = $Colors{$sample} || "#000000";
			my $ylevel = nearest($percent_precision, ($START_Y-$pcid));
			my $jitter = $percent_precision*0.1*((-1)**(int(rand(2))))*rand(); # Add jitter that is +- 10% of precision
			$ylevel += $jitter;
			if (!$Drawn->{$sstart}->{$ylevel}->{$send}) {
				draw_line_scaled($FH, $sstart-$START_POS+1, $ylevel, $send-$START_POS+1, $ylevel, "stroke:$color;stroke-width:${line_stroke_width}px");
				$Drawn->{$sstart}->{$ylevel}->{$send} = 1;
			}
		}
		close(ALIGN);
		print $FH "</g>\n";
	}
}

sub min {
	my ($a, $b) = @_;
	if ($a < $b) {
		return $a;
	} else {
		return $b;
	}
}

sub max {
	my ($a, $b) = @_;
	if ($a > $b) {
		return $a;
	} else {
		return $b;
	}
}

1;
