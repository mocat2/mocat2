package Smash::Utils::SVG;

use strict;
use warnings;
use POSIX;
use Math::Round;

use base 'Exporter';
our @EXPORT_OK = qw(open_svg close_svg open_style close_style add_class draw_polygon draw_rect draw_rect_scaled_x draw_rect_scaled_y draw_rect_scaled draw_line draw_line_scaled_x draw_line_scaled_y draw_line_scaled draw_text draw_text_scaled_x draw_text_scaled_y draw_text_scaled draw_pie_chart $POINT_PRECISION $GLOBAL_X_SCALE_DOWN $GLOBAL_Y_SCALE_DOWN $GLOBAL_X_OFFSET);
our %EXPORT_TAGS = ('all' => [@EXPORT_OK]);

our $POINT_PRECISION     = 0.00001;
our $GLOBAL_X_SCALE_DOWN = 1;
our $GLOBAL_Y_SCALE_DOWN = 1;
our $GLOBAL_X_OFFSET = 0;

sub open_svg {
	my ($FH, $width, $height) = @_;
	$width || ($width = "100%");
	$height || ($height = "100%");
	print $FH <<EOF;
<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg width="$width" height="$height" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">
EOF
	return;
}

sub close_svg {
	my $FH = shift;
	print $FH <<EOF;
</svg>
EOF
	return;
}

sub open_style {
	my $FH = shift;
	my $type = shift;
	print $FH <<EOF;
<defs>
<style type="$type"><![CDATA[
EOF
	return;
}

sub close_style {
	my $FH = shift;
	print $FH <<EOF;
]]></style>
</defs>
EOF
	return;
}

sub draw_polygon {
	my ($FH, $points, $style) = @_;
	$style || ($style = "stroke:rgb(99,99,99);stroke-width:1px");
	print $FH <<EOF;
<polygon points="$points" style="$style"/>
EOF
	return;
}

sub draw_rect {
	my ($FH, $x, $y, $w, $h, $style) = @_;
	$x += scale_down_x($GLOBAL_X_OFFSET);
	($x, $y, $w, $h) = map { nearest($POINT_PRECISION, $_) } ($x, $y, $w, $h);
	$style || ($style = "stroke:rgb(99,99,99);stroke-width:1px");
	print $FH <<EOF;
<rect x="$x" y="$y" width="$w" height="$h" style="$style"/>
EOF
	return;
}

sub draw_rect_scaled_x {
	my ($FH, $x, $y, $w, $h, $style) = @_;
	($x, $w) = scale_down_x($x, $w);
	draw_rect($FH, $x, $y, $w, $h, $style);
	return;
}

sub draw_rect_scaled_y {
	my ($FH, $x, $y, $w, $h, $style) = @_;
	($y, $h) = scale_down_y($y, $h);
	draw_rect($FH, $x, $y, $w, $h, $style);
	return;
}

sub draw_rect_scaled {
	my ($FH, $x, $y, $w, $h, $style) = @_;
	($x, $w) = scale_down_x($x, $w);
	($y, $h) = scale_down_y($y, $h);
	draw_rect($FH, $x, $y, $w, $h, $style);
	return;
}

sub draw_line {
	my ($FH, $x1, $y1, $x2, $y2, $style) = @_;
	$x1 += scale_down_x($GLOBAL_X_OFFSET);
	$x2 += scale_down_x($GLOBAL_X_OFFSET);
	($x1, $y1, $x2, $y2) = map { my $a = nearest($POINT_PRECISION, $_); $a =~ s/0+$// if ($a =~ /\.[0-9]+0/); $a } ($x1, $y1, $x2, $y2);
	$style || ($style = "stroke:rgb(99,99,99);stroke-width:1px");
	if ($style =~ /^class:/) {
		$style =~ s/^class://;
		$style = "class=\"$style\"";
	} else {
		$style = "style=\"$style\"";
	}
	print $FH <<EOF;
<line x1="$x1" y1="$y1" x2="$x2" y2="$y2" $style/>
EOF
	return;
}

sub draw_line_scaled_x {
	my ($FH, $x1, $y1, $x2, $y2, $style) = @_;
	($x1, $x2) = scale_down_x($x1, $x2);
	draw_line($FH, $x1, $y1, $x2, $y2, $style);
	return;
}

sub draw_line_scaled_y {
	my ($FH, $x1, $y1, $x2, $y2, $style) = @_;
	($y1, $y2) = scale_down_y($y1, $y2);
	draw_line($FH, $x1, $y1, $x2, $y2, $style);
	return;
}

sub draw_line_scaled {
	my ($FH, $x1, $y1, $x2, $y2, $style) = @_;
	($x1, $x2) = scale_down_x($x1, $x2);
	($y1, $y2) = scale_down_y($y1, $y2);
	draw_line($FH, $x1, $y1, $x2, $y2, $style);
	return;
}

sub draw_text {
	my ($FH, $x, $y, $text, %options) = @_;
	$x += scale_down_x($GLOBAL_X_OFFSET);
	$options{"font-size"} || die "draw_text needs font-size to be set";
	my $attribute = "";
	my $rotate    = $options{rotate};
	if ($rotate) {
		$attribute = "transform=\"rotate($rotate $x,$y)\" ";
		delete $options{rotate};
	}
	while (my ($key, $value) = each %options) {
		$attribute .= "$key=\"$value\" ";
	}
	($x, $y) = map { nearest($POINT_PRECISION, $_) } ($x, $y);
	print $FH <<EOF;
<text x="$x" y="$y" fill="navy" $attribute>
$text
</text>
EOF
	return;
}

sub draw_text_scaled_x {
	my ($FH, $x, $y, $text, %options) = @_;
	($x) = scale_down_x($x);
	draw_text($FH, $x, $y, $text, %options);
	return;
}

sub draw_text_scaled_y {
	my ($FH, $x, $y, $text, %options) = @_;
	($y) = scale_down_y($y);
	draw_text($FH, $x, $y, $text, %options);
	return;
}

sub draw_text_scaled {
	my ($FH, $x, $y, $text, %options) = @_;
	($x) = scale_down_x($x);
	($y) = scale_down_y($y);
	draw_text($FH, $x, $y, $text, %options);
	return;
}

sub draw_pie_chart {
	my $FH  = shift;
	my $r   = shift;
	my $pie = shift;
	my $colors = shift;

	print $FH <<EOF;
<g transform="translate($r, $r)">
EOF

	my $sum = 0;
	foreach my $piece (@$pie) {
		my $theta;
		my $sweep_flag = 0;
		my $large_arc_flag = 0;
		if ($piece > 0.5) {
			$large_arc_flag=1;
		}
		if ($piece == 1) {
			$sweep_flag = 1;
		}

		my $color = shift(@$colors);
		push(@$colors, $color);

		$theta    = $sum*2*3.1416;
		my $start_x = $r*cos($theta);
		my $start_y = -$r*sin($theta);

		$sum += $piece;
		$theta    = $sum*2*3.1416;
		my $end_x   = $r*cos($theta);
		my $end_y   = -$r*sin($theta);
		

		print $FH <<EOF;
<path d="M 0,0 L $start_x,$start_y A $r,$r 0 $large_arc_flag,$sweep_flag $end_x,$end_y z" fill="$color" />
EOF
	}

	print $FH <<EOF;
</g>
EOF
}

sub scale_down_x {
	return map { $_/$GLOBAL_X_SCALE_DOWN } @_;
}

sub scale_down_y {
	return map { $_/$GLOBAL_Y_SCALE_DOWN } @_;
}

1;
