#!/usr/bin/env perl
use warnings;
use strict;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

my $usage =<<USAGE;#******* Instruction of this program *********#

Usage: perl $0 <contig.fa> <depth.out> <pe.list> <filter.contig.fa> <cut off length> <sample>

USAGE
die $usage unless(@ARGV==6);

############################################################ VARIABLES ############################################################
my (%contig, @lack, %length, %block_0, %pair, $num_0, %num_0, @temp, @t, $temp, $i, $j, $k);
my ($contig, $depth, $pe_list, $out, $lengthcutoff, $sample) = @ARGV;
############################################################ CONSTANTS ############################################################
#use constant SAMPLESIZE => 73;
############################################################ MAIN ############################################################
my $date;

open I, $contig or die "$contig $! \n";
$/ = ">";
while(<I>){
    chomp;
    next if($_ eq "");
    @temp = split /\n/;
    @t = split /\s+/, $temp[0];
    shift @temp;
    $contig{$t[0]} = join "", @temp;
    $length{$t[0]} = length $contig{$t[0]};
}
close I;

open I, $depth or die "$depth $! \n";
while(<I>){
    chomp;
    next if($_ eq "");
    @temp = split /\n/;
    @t = split /\s+/, $temp[0];
    shift @temp;
    $temp = join " ", @temp;
    @temp = split /\s+/, $temp;
    @lack = ();
    push @lack, 1;
    $num_0 = 0;
    for($i=100; $i<@temp-100; $i++){
        if($temp[$i] eq "0"){
            push @lack, $i;
            while($i<@temp-100 and $temp[$i] eq "0"){
                $i ++;
                $num_0 ++;
            }
            push @lack, $i+1;
        }
    }
    push @lack, scalar(@temp);
    $num_0{$t[0]} = $num_0;
    next if($num_0 == 0);
    for($i=1; $i<@lack-1; $i+=2){
        push @{$block_0{$t[0]}}, ($lack[$i]+1)."-".($lack[$i+1]-1);
    }
}
close I;

$num_0 = 0;
$/ = "\n";
open L, $pe_list or die "$pe_list $! \n";
while( <L> ){
    chomp;
    if(m/\.pe\.gz$/){
        open I, "gzip -dc $_ |" or die "$_ $! \n";
    }elsif(m/\.pe$/){
        open I, $_ or die "$_ $! \n";
    }else{
        next;
    }
    my $name = "";
    my %h = ();
    chomp($temp = <I>);
    unless ($temp) {
    	die "ERROR & EXIT: Something is wrong. Program died in AssemblyRevision_filter_contigs.pl on line 82. Run revision again? Sorry, don't know what's wrong...";
    }
    @temp = split /\s+/, $temp;
    $temp[0] =~ s/\/\d$//;
    $name = $temp[0];
    push @{$h{$temp[7]}}, $temp[8];
    while( <I> ){
      chomp;
        @temp = split;
        $temp[0] =~ s/\/\d$//;
        if($name ne $temp[0]){
	  foreach $i (keys %h){
	    if(defined $block_0{$i}){
	      foreach $j (@{$block_0{$i}}){
		@t = split /-/, $j;
		if(($h{$i}[0] <= $t[0] and $h{$i}[1] >= $t[1]) or
		   ($h{$i}[1] <= $t[0] and $h{$i}[0] >= $t[1])){
		  $pair{$i}{$t[0]} ++;
		  $num_0 ++;
		}
	      }
	    }
	  }
	  $name = $temp[0];
	  %h = ();
	  #            next;
        }
        push @{$h{$temp[7]}}, $temp[8];
      }
    close I;
  }
close L;
print $num_0, "\n";

open O, ">$out.num" or die "$out.num $! \n";
print O "name\tlength\t0_num\t0_block_num\t0_block_length\tposition\tpair_num\n";
foreach $i (keys %block_0){
    print O "$i\t$length{$i}\t$num_0{$i}\t", scalar(@{$block_0{$i}}), "\t";
    if(@{$block_0{$i}} == 1){
        @t = split /\-/, $block_0{$i}[0];
        print O $t[1]-$t[0]+1;
    }else{
    foreach $j (@{$block_0{$i}}){
        @t = split /\-/, $j;
        print O $t[1]-$t[0]+1,",";
    }
    }
    print O "\t";
    if(@{$block_0{$i}} == 1){
        print O $block_0{$i}[0];
    }else{
        foreach $j (@{$block_0{$i}}){
            print O "$j,";
        }
    }
    print O "\t";
    if(defined $pair{$i}){
            foreach $j (@{$block_0{$i}}){
                @t = split /-/, $j;
                if(defined $pair{$i}{$t[0]}){
                    print O "$pair{$i}{$t[0]},";
                }else{
                    print O "0,";
                    $pair{$i}{$t[0]} = 0;
                }
            }
        print O "\n";
    }else{
        foreach $j (@{$block_0{$i}}){
            @t = split /-/, $j;
            print O "0,";
            $pair{$i}{$t[0]} = 0;
        }
        print O "\n";
    }
}
close O;

open O, ">$out" or die "$out $! \n";
foreach $i (keys %contig){
    @temp = split //, $contig{$i};
    next if(length($contig{$i}) < $lengthcutoff);
    unless(defined $block_0{$i}){
      my $PRINT = "$sample\_revised_$i";
      $PRINT =~ s/_$sample//;
        print O ">$PRINT  length=", length($contig{$i});
        for($j=0; $j<@temp; $j++){
            print O "\n" if($j%60 == 0);
            print O $temp[$j];
        }
        print O "\n";
        next;
    }
    @lack = ();
    for($j=0; $j<@{$block_0{$i}}; $j++){
        @t = split /-/, $block_0{$i}[$j];
        if($t[1]-$t[0]<=50 and $pair{$i}{$t[0]}>=20){
            next;
        }
        if(@lack < 2){
            push @lack, 0;
        }
        push @lack, $t[0]-2;
        push @lack, $t[1];
    }
    push @lack, scalar(@temp)-1;
    if(@lack < 2){
      my $PRINT = "$sample\_revised_$i";
      $PRINT =~ s/_$sample//;
        print O ">$PRINT  length=", length($contig{$i});
        for($k=0; $k<@temp; $k++){
            print O "\n" if($k%60 == 0);
            print O $temp[$k];
        }
        print O "\n";
        next;
    }
    $temp = 0;
    for($j=0; $j<@lack; $j+=2){
        if($lack[$j+1]-$lack[$j]+1 >= $lengthcutoff){
            $temp ++;
	    my $PRINT = "$sample\_revised_$i";
	    $PRINT =~ s/_$sample//;
            print O ">$PRINT\_$temp  length=", $lack[$j+1]-$lack[$j]+1;
            for($k=$lack[$j];$k<=$lack[$j+1];$k++){
                if(($k-$lack[$j]) % 60 == 0){
                    print O "\n";
                }
                print O $temp[$k];
            }
            print O "\n";
        }
    }
}

print STDERR "program finished!\t",$date = localtime,"\n";
############################################################ SUB FUNCTION ############################################################

