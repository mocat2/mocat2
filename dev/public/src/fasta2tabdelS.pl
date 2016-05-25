#!/usr/bin/perl -w

###############################################################
#This is a commented version of a perl script that will       #
#convert files from fasta to tab delimited format. 06/28/07   #
###############################################################

use strict;
use warnings;

#Initiate variables $count and $len
my $count=0;
my $len=0;

#Read fastafile, delete all return or newline; all tabs are replaced by space
while(<STDIN>) {
s/\r?\n//g; s/\t/ /g;

#if ">" is found at beginning of line  
	if (m/^>/){

#if ">" is found in first line of fastafile, skip this loop, else add newline
		if ($. != 1) {
			print STDOUT "\n"
		}	

#count occurrences of ">" and add tab to end of line
		$count++; $_ .= "\t"
	}
#for all other lines, remove all space characters and sum up the length of line
	else {
		s/ //g; $len += length($_)
	} 

#print out the line	 
print STDOUT $_;
}


