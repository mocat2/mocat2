#!/usr/bin/perl -w

###############################################################
#This is a commented version of a perl script that will       #
#convert files from fasta to tab delimited format. 06/28/07   #
###############################################################

use strict;
use warnings;

#Die, if script is executed without argument (fastafile) 
unless ($ARGV[0]) {die "Usage: fasta2tabdel.pl fastafile\n"}

#Assign (first) argument to variable $fasta
my $fasta = shift @ARGV;

#Add .tabdel to outputfile
my $outfile = "$fasta.tabdel";

#Die if fastafile does not exist
unless (open(FASTA, $fasta)) {
	die "Cannot open file \"$fasta\"\n"
};

#Open (create) outputfile
open OUT, ">$outfile";

#Initiate variables $count and $len
my $count=0;
my $len=0;

#Read fastafile, delete all return or newline; all tabs are replaced by space
while(<FASTA>) {
s/\r?\n//g; s/\t/ /g;

#if ">" is found at beginning of line  
	if (m/^>/){

#if ">" is found in first line of fastafile, skip this loop, else add newline
		if ($. != 1) {
			print OUT "\n"
		}	

#count occurrences of ">" and add tab to end of line
		$count++; $_ .= "\t"
	}
#for all other lines, remove all space characters and sum up the length of line
	else {
		s/ //g; $len += length($_)
	} 

#print out the line	 
print OUT $_;
}

#print  "\n"; warn "\nConverted $count FASTA records in $. lines to tabular format\nTotal sequence length: $len\n\n"; 
