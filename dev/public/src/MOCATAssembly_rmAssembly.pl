#!/usr/bin/env perl
use warnings;
use strict;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

die "usage: perl $0 assemlbyDir\n" unless(@ARGV == 1);

print "rm -f $ARGV[0]/*edge* $ARGV[0]/*Arc* $ARGV[0]/*ContigIndex* $ARGV[0]/*vertex* $ARGV[0]/*preGraphBasic* $ARGV[0]/*ContigIndex* $ARGV[0]/*updated.edge* $ARGV[0]/*readOnContig* $ARGV[0]/*peGrads* $ARGV[0]/*readInGap* $ARGV[0]/*newContigIndex* $ARGV[0]/*links $ARGV[0]/*scaf_gap* $ARGV[0]/*contigPosInscaff* $ARGV[0]/*contig4asm* $ARGV[0]/*asm* $ARGV[0]/*kmerFreq* $ARGV[0]/*scaf* $ARGV[0]/*gapSeq* $ARGV[0]/*shortreadInGap* $ARGV[0]/*PEreadOnContig* $ARGV[0]/*lanes $ARGV[0]/*inserts* $ARGV[0]/*lengths*\n";
`rm -f $ARGV[0]/*edge* $ARGV[0]/*Arc* $ARGV[0]/*ContigIndex* $ARGV[0]/*vertex* $ARGV[0]/*preGraphBasic* $ARGV[0]/*ContigIndex* $ARGV[0]/*updated.edge* $ARGV[0]/*readOnContig* $ARGV[0]/*peGrads* $ARGV[0]/*readInGap* $ARGV[0]/*newContigIndex* $ARGV[0]/*links $ARGV[0]/*scaf_gap* $ARGV[0]/*contigPosInscaff* $ARGV[0]/*contig4asm* $ARGV[0]/*asm* $ARGV[0]/*kmerFreq* $ARGV[0]/*scaf* $ARGV[0]/*gapSeq* $ARGV[0]/*shortreadInGap* $ARGV[0]/*PEreadOnContig* $ARGV[0]/*lanes $ARGV[0]/*inserts* $ARGV[0]/*lengths*`;
