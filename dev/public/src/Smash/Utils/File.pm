#!/usr/bin/env perl

use strict;
use warnings;

package Smash::Utils::File;

use overload '""' => 'filename';

sub new {
	my $class  = shift;
	my %params = @_;
	%params = map {uc($_) => $params{$_}} keys %params;

	my $name = $params{NAME};
	my $this = open_file($name, $params{GZIP});

	# scalar context

	${*$this} = $name;

	# hash context, from constructor

	%{*$this} = %params;

	bless $this, $class;

	return $this;
}

sub gzip {shift->{GZIP}}

sub filename {
	my $this = shift;
	return ${*$this};
}

sub open_file {
	my $name = shift;
	my $gzip = shift;
	my $FH;
	if ($gzip) {
		open($FH, "zcat $name | ") || die "Cannot pipe uncompressed $name: $!";
	} else {
		open($FH, "<$name") || die "Cannot open $name: $!";
	}
	return $FH;
}

sub close {
	my $this = shift;
	close($this);
	return 1;
}

sub reopen {
	my $this = shift;
	if ($this->fh) {
		$this->close();
	}
	return $this->open();
}

sub read_file_as_string {
	my $this = shift;
	my $FH = $this->fh;
	my $string = "";
	while (<$FH>) {
		$string .= $_;
		$string .= "\n";
	}
	return $string;
}

sub write_string_to_file {
	my ($file, $string) = @_;
	my $FH;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, ">$file") || die "Cannot open $file: $!";
	}
	print {$FH} $string;
	if (ref($file) !~ /GLOB/) {
		CORE::close($FH);
	}
}

1;
