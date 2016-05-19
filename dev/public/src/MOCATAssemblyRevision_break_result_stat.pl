#!/usr/bin/env perl
use warnings;
use strict;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

my $usage =<<USAGE;#******* Instruction of this program *********#

Usage: perl $0 <break.num.list> <out>

USAGE
die $usage unless(@ARGV==2);

############################################################ VARIABLES ############################################################
my ($break_num, $lack, $i, $j, $name, @temp);
my ($list, $out) = @ARGV;
$break_num = 0;
$lack = 0;
############################################################ CONSTANTS ############################################################
#use constant SAMPLESIZE => 73;
############################################################ MAIN ############################################################
my $date;

open L, $list or die "$list $! \n";
open O, ">$out" or die "$out $! \n";
print O "name\tbreak_num\tlack_base_num\n";
while( <L> ){
    chomp;
    @temp = split /\//;
    $temp[-1] =~ m/(.*)\.scaftig/;
    $name = $1;
    open I, $_ or die "$_ $! \n";
    <I>;
    while( <I> ){
        chomp;
        @temp = split;
        my(@a1, @a2, @a3);
        @a1 = split /,/, $temp[-3];
        @a2 = split /[,-]/, $temp[-2];
        @a3 = split /,/, $temp[-1];
        push @a2, $temp[1]+1;
        unshift @a2, 0;
        for($i=0; $i<@a1; $i++){
            if($a1[$i] <= 50 and $a3[$i] >= 20){
                next;
            }
            $break_num++;
            if($a2[2*$i+1]-$a2[2*$i]-1<500){
                $lack += $a2[2*$i+1]-$a2[2*$i]-1;
            }
            $lack += $a1[$i];
            if($i == $#a1){
                if($a2[-1]-$a2[-2]-1<500){
                    $lack += $a2[-1]-$a2[-2]-1;
                }
            }
        }
    }
    print O "$name\t$break_num\t$lack\n";
    $break_num = 0;
    $lack = 0;
    close I;
}
close L;
close O;


print STDERR "break result stat program finished!\t",$date = localtime,"\n";
############################################################ SUB FUNCTION ############################################################

sub function{
	return  ;
}
