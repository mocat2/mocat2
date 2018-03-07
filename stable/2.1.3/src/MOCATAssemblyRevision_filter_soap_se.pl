#!/usr/bin/env perl
use warnings;
use strict;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

my($se,$identity);
my $opt;
while($opt = shift){
    if($opt eq "-s"){$se = shift;}
    elsif($opt eq "-i"){$identity = shift;}
    elsif($opt eq "-h"){&usage();exit;}
}

unless($se){
    &usage();
    exit;
}
$identity = 0.9 unless(defined($identity));

my $filter_se = "$se.filter.se";
my $filter_stat = "$se.filter.stat";
my $se_num = 0;
open OUPS,">$filter_se" or die "can\'t open filter of se result: $filter_se\n";
open IPS,$se or die "can\'t open se result: $se\n";
my $last_name = "";
my @info;
while(<IPS>){
    chomp;
    my @temp = split;
    if($temp[0] eq $last_name){
        push @info,$_;
    }else{
        &trim_se(@info) if(@info);
        @info = ();
        push @info,$_;
        $last_name = $temp[0];
    }
}
&trim_se(@info) if(@info);
close IPS;
close OUPS;

#my $se_num = keys %se;
open OUS,">$filter_stat" or die "can\'t open filter stat: $filter_stat\n";
print OUS "Singled:\t$se_num\n";
close OUS;

sub trim_se{
    my @note = @_;
    my %se = ();
    my @result = ();
    for(0..$#note){
        my @temp = split /\s+/,$note[$_];
        my $mismatch = $temp[-1] =~ tr/[a-zA-Z]/[a-zA-Z]/;
        if(($temp[5] - $mismatch)/$temp[5] >= $identity){
            push @result,$note[$_];
            $se{$temp[0]}++;
        }
    }
    for(0..$#result){
        my @temp = split /\s+/,$result[$_];
        $temp[3] = $se{$temp[0]};
        print OUPS join "\t",@temp;
        print OUPS "\n";
    }
    $se_num++ if(@result);
}

sub usage{
    print <<EOD
    usage: perl $0 -s soap.se -i identity
        -s soap se result, required
        -i identity, default 0.9
EOD
}
