#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

die
"Usage: perl $0 <eindel.xls> <input.fsa> <output.fsa> <cutoff> <reformat file> <quality>\n"
  unless ( scalar @ARGV == 6 and -f $ARGV[0] and -f $ARGV[1] );

our $reformat_pl = $ARGV[4];
my ( $cutoff, $quality ) = ( $ARGV[3], $ARGV[5] );
if ( $ARGV[0] =~ m/\.gz/ ) {
	open IN, "gzip -dc $ARGV[0] |";
}
else {
	open IN, "$ARGV[0]";
}
open OUT, ">$ARGV[2].log";

my ( @l, %hash, %h );
while ( @l = split( /\s+/, <IN> ) ) {
	next unless ( $#l == 9 or $#l == 14 );
	next if ( $l[7] <= 10 or $l[4] < $quality or $l[5] < $quality );
	if ( $#l == 9 ) {
		my ( $ctg, $pos, $aln ) = @l[ 0, 1, 8 ];
		if ( $aln =~ m/[\+-](\d+)[atcgATCGnN]+/ ) {
			$aln =~ s/[\+-]\d+[atcgATCGnN]{$1}//g;
		}
		$aln =~ s/[\$\^\?]//g;

		my $ref = $aln =~ tr/.,/.,/;
		next unless ( $ref <= 10 );

		my $a = $aln =~ tr/aA/aA/;
		my $t = $aln =~ tr/tT/tT/;
		my $c = $aln =~ tr/cC/cC/;
		my $g = $aln =~ tr/gG/gG/;

		my ( $s, $ans );
		$s   = ( $a > $t ) ? $a   : $t;
		$ans = ( $a > $t ) ? "A"  : "T";
		$s   = ( $s > $c ) ? $s   : $c;
		$ans = ( $s > $c ) ? $ans : "C";
		$s   = ( $s > $g ) ? $s   : $g;
		$ans = ( $s > $g ) ? $ans : "G";
		next unless ( $s / $l[7] >= 0.7 );

		print OUT join( "\t", @l ) . "\n";
		$hash{$ctg}{$pos} = $ans;

	}
	elsif ( $#l == 14 ) {
		my ( $ctg, $pos, $dep, $a1, $a2, $qual, $qual2 ) =
		  @l[ 0, 1, 7, 8, 9, 10, 11 ];
		next
		  unless ( ( $qual / $dep > $cutoff and $a1 =~ m/^[\+\-]/ )
			or
			( $qual2 / $dep > $cutoff and $a2 =~ m/^[\+\-]/ )
		  );
		next if ( defined $hash{$ctg}{$pos} );
		print OUT join( "\t", @l ) . "\n";
		if ( $qual / $dep > $cutoff and $a1 =~ m/^[\+\-]/ ) {
			$hash{$ctg}{$pos} = $a1;
		}
		elsif ( $qual2 / $dep > $cutoff and $a2 =~ m/^[\+\-]/ ) {
			$hash{$ctg}{$pos} = $a2;
		}
	}
}

close IN;
close OUT;

##
# eindel revision
open IN,  "$ARGV[1]";
open OUT, ">$ARGV[2].tmp";
local $/ = '>';

while (<IN>) {
	chomp;
	next unless $_;
	my ( $ids, $seq ) = split( /\n/, $_, 2 );
	$seq =~ s/\n//g;

	my ( $id, $dummy ) = split( /\s+/, $ids, 2 );
	my $h = $hash{$id};
	my @lst = sort { $a <=> $b } keys %$h;

	print OUT ">$id\n";
	my ( $prev, $len, $x ) = ( 1, 0, 0 );
	foreach my $l (@lst) {
		$len = $l - $prev;
		print OUT substr( $seq, $prev - 1, $len );
		my $c = $h->{$l};

		#print "l=$l, prev=$prev, len=$len, c=$c\n";

		if ( $c =~ m/^-/ ) {
			print OUT substr( $seq, $l - 1, 1 );

			#print substr($seq, $l -1, 1), "1";
			$x = length( substr( $c, 1 ) );
			my $lx = $l + $x;
			if ( exists( $h->{$lx} ) ) {
				$prev = $l + 1;
			}
			else {
				$prev = $l + $x + 1;
			}
		}
		elsif ( $c =~ m/^\+/ ) {
			print OUT substr( $seq, $l - 1, 1 );
			print OUT substr( $c, 1 );

			#print substr($seq, $l -1, 1), "2";
			#print substr($c, 1), "3";
			$prev = $l + 1;
		}
		else {
			print OUT $c;

			#print $c, "4";
			$prev = $l + 1;
		}

		#print "---------\n";
	}
	print OUT substr( $seq, $prev - 1 );
	print OUT "\n";
}

system("$reformat_pl $ARGV[2].tmp $ARGV[2] >/dev/null; rm $ARGV[2].tmp");

