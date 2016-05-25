#!/usr/bin/env perl
use warnings; 
use strict;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

unless(@ARGV==3) {
	&usage;
	exit;
}
my($file,$filter,$seq,$len,@new);
open IN,$ARGV[0] or die "$ARGV[0] $!\n";
$filter = $ARGV[1];
open OT,">$ARGV[2].more$ARGV[1]" or die "$ARGV[2].more$ARGV[1] $!\n";
$/ = ">";
#<IN>;
while(<IN>) {
	chomp;
	@new = split /\n/;
	shift @new;
	$seq = join "",@new;
	$len = rindex $seq."\$","\$";
	if($len >= $filter) {
		print OT ">",$_;
	}
}
$/ = "\n";
close IN;
close OT;

sub usage{
	print STDERR<<USAGE;
	
	perl $0 [fa.file] [cut by length] [out name]
	
USAGE
}
