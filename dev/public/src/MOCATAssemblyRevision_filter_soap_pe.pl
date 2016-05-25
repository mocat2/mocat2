#!/usr/bin/env perl
use warnings;
use strict;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

my($pair_pe,$identity);
my $opt;
while($opt = shift){
    if($opt eq "-p"){$pair_pe = shift;}
    elsif($opt eq "-i"){$identity = shift;}
    elsif($opt eq "-h"){&usage();exit;}
}

unless($pair_pe){
    &usage();
    exit;
}
$identity = 0.9 unless(defined($identity));

my $filter_pe = "$pair_pe.filter.pe";
my $filter_se = "$pair_pe.filter.se";
my $filter_stat= "$pair_pe.filter.stat";

open OUPP,">$filter_pe" or die "can\'t open filter of pe result: $filter_pe\n";
open OUPS,">$filter_se" or die "can\'t open filter of se result: $filter_se\n";
open IPP,$pair_pe or die "can\'t open pe result: $pair_pe\n";
my $pe_num = 0;
my $se_num = 0;
my @pe=();
my $old_name = "";
while(<IPP>){
    chomp;
    my @temp = split;
    my $reads=$temp[0];    
    $reads=~s/\/[12]$//;
    if($old_name eq $reads){
        push @pe,$_;
    }else{
        &filter_pe(@pe) if(@pe);
        @pe = ();
        push @pe,$_;
        $old_name = $reads;
    }
}
 &filter_pe(@pe) if(@pe);

close IPP;
close OUPP;
close OUPS;

open OUS,">$filter_stat" or die "can\'t open filter stat: $filter_stat\n";
print OUS "Paired:\t$pe_num\tPE\n";
print OUS "Singled:\t$se_num\tSE\n";
close OUS;
#-------------------------------------------------------------------------------
sub filter_pe{
    my @array = @_;
    my $num= @array/2;
    my (@pe_reads1,@pe_reads2,@se_reads1,@se_reads2);
    my $pe_mark = 0;
    my $flag1 = 0;
    my $flag2 = 0;
    my @se = ();
    for(my $i=0;$i<$num;$i++){
        my $reads1=$array[$i];
        my $reads2=$array[$i+$num];
        my @temp1=split /\s+/,$reads1;
        my @temp2=split /\s+/,$reads2;
        my $mismatch1 = $temp1[-1] =~ tr/[a-zA-Z]/[a-zA-Z]/;
        my $mismatch2 = $temp2[-1] =~ tr/[a-zA-Z]/[a-zA-Z]/;
        $flag1=(($temp1[5] - $mismatch1)/$temp1[5] >= $identity)?1:0;
        $flag2=(($temp2[5] - $mismatch2)/$temp2[5] >= $identity)?1:0;
        if($flag1 == 1 and $flag2 == 1){
            push @pe_reads1,$reads1;
            push @pe_reads2,$reads2;
            $pe_mark = 1;
        }elsif($flag1 == 1){
            push @se_reads1,$reads1;
        }elsif($flag2 == 1){
            push @se_reads2,$reads2;
        }
    }
    if($pe_mark){
        $pe_num++;
        my $comp_num = @pe_reads1;
        foreach my $comparison1(@pe_reads1){
            my @info1 = split /\s+/,$comparison1;
            $info1[3] = $comp_num;
            print OUPP join "\t",@info1;
            print OUPP "\n";
        }
        foreach my $comparison2(@pe_reads2){
            my @info2 = split /\s+/,$comparison2;
            $info2[3] = $comp_num;
            print OUPP join "\t",@info2;
            print OUPP "\n";
        }
    }else{
        $se_num++ if(@se_reads1);
        $se_num++ if(@se_reads2);
        my $comp_num1 = @se_reads1;
        foreach my $comparison1(@se_reads1){
            my @info1 = split /\s+/,$comparison1;
            $info1[3] = $comp_num1;
            print OUPS join "\t",@info1;
            print OUPS "\n";
        }
        my $comp_num2 = @se_reads2;
        foreach my $comparison2(@se_reads2){
            my @info2 = split /\s+/,$comparison2;
            $info2[3] = $comp_num2;
            print OUPS join "\t",@info2;
            print OUPS "\n";
        }
    }
}

sub usage{
    print <<EOD
    usage: perl $0 -p soap.pe -i identity
        -p the pair-end soap result of pair-end reads, required
        -i identity, default 0.9
EOD
}
