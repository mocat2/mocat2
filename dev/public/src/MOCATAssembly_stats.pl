#!/usr/bin/env perl
use strict;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

die "usage: perl $0 contig.fa scafSeq.fa scaftig.fa output MOCATbindir\n" unless(@ARGV == 5);
my($contig,$scafSeq,$scaftig,$output,$bin_dir) = @ARGV;

`echo "contig" >$output`;
`echo " " >>$output`;
`perl $bin_dir/MOCATAssembly_contig_stats.pl $contig 0:200:500:1000:2000:4000 >>$output`;

`echo " " >>$output`;
`echo " " >>$output`;
`echo "scafold" >>$output`;
`echo " " >>$output`;
`perl $bin_dir/MOCATAssembly_contig_stats.pl $scafSeq 0:200:500:1000:2000:4000 >>$output`;

`echo " " >>$output`;
`echo " " >>$output`;
`echo "scaftig" >>$output`;
`echo " " >>$output`;
`perl $bin_dir/MOCATAssembly_contig_stats.pl $scaftig 0:200:500:1000:2000:4000 >>$output`;
