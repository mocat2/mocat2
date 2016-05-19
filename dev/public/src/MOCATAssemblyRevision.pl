#!/usr/bin/env perl
use warnings;
use strict;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

my $usage = <<USAGE;    #******* Instruction of this program *********#

Usage: perl $0 <contig> <fq.list> <insertsize.list> <length.info> <outdir> <bindir> <scr folder> <processors> <scaftig min length> <sample> <StatsOutFile>

Version 2.0
	For each sample, 3 fq file required, 1.fq, 2.fq and single.fq.
	Speed up in BWA step.
	
USAGE
die $usage unless ( @ARGV == 11 );

############################################################ VARIABLES ############################################################
my ( %fq, @temp, $tmp, $sample_name, %insert, %len );
my ( $scaftig, $fq_list, $insert, $len_info, $outdir, $bin_folder, $scr_folder, $processors, $scaftigminlength, $sample, $tatsOut ) = @ARGV;
############################################################ CONSTANTS ############################################################
my $cutmore       = "$scr_folder/MOCATAssemblyRevision_cut_more.pl";
my $reformat      = "$scr_folder/MOCATAssemblyRevision_reformat.pl";
my $reads_mapping = "$scr_folder/MOCATAssemblyRevision_reads_mapping.pl";
my $correct_it    = "$scr_folder/MOCATAssemblyRevision_correct_it.pl";
my $samtools      = "$bin_folder/samtools";
my $soap          = "$bin_folder/soap2.21";
my $soap_coverage = "$bin_folder/soap.coverage";
my $bwt           = "$bin_folder/2bwt-builder";
my $filter_pe     = "$scr_folder/MOCATAssemblyRevision_filter_soap_pe.pl";
my $filter_se     = "$scr_folder/MOCATAssemblyRevision_filter_soap_se.pl";
my $filter        = "$scr_folder/MOCATAssemblyRevision_filter_contig.pl";
my $correct       = "$scr_folder/MOCATAssemblyRevision_count.sh";
my $stat          = "$scr_folder/MOCATAssemblyRevision_stat.sh";
my $cutoff        = 0.5;
my $quality       = 20;
############################################################ MAIN ############################################################
my $date;

print STDERR "program begin!\t", $date = localtime, "\n";

#print "Bin folder is $bin_folder\n";

system_("mkdir -p $outdir")                unless ( -e $outdir );
system_("mkdir -p $outdir/index")          unless ( -e "$outdir/index" );
system_("mkdir -p $outdir/bwa")            unless ( -e "$outdir/bwa" );
system_("mkdir -p $outdir/correct_Result") unless ( -e "$outdir/correct_Result" );
system_("mkdir -p $outdir/soap")           unless ( -e "$outdir/soap" );
system_("mkdir -p $outdir/soap_coverage")  unless ( -e "$outdir/soap_coverage" );
system_("mkdir -p $outdir/break_contig")   unless ( -e "$outdir/break_contig`unless" );

chomp( $tmp = `head -1 $fq_list` );
@temp = split /\//, $tmp;
$sample_name = $temp[-3];
system_("mkdir $outdir/soap/$sample_name") unless ( -e "$outdir/soap/$sample_name" );
open IN, $fq_list or die "$fq_list $! \n";
while (<IN>) {
	chomp;
	@temp = split /\//, $_;
	if ( $temp[-1] =~ /1.fq/ ) {
		my $x = $temp[-1];
		$x =~ s/.pair.1.*//;
		$fq{$x}[0] = $_;
	}
	elsif ( $temp[-1] =~ /2.fq/ ) {
		my $x = $temp[-1];
		$x =~ s/.pair.2.*//;
		$fq{$x}[1] = $_;
	}
	else {
		my $x = $temp[-1];
		$x =~ s/.single.*//;
		$fq{$x}[2] = $_;
	}
}
close IN;

my $contig = ( split /\//, $scaftig )[-1];
$contig = "$outdir/$contig";
system_("perl $cutmore $scaftig $scaftigminlength $contig");

#print "perl $cutmore $scaftig $scaftigminlength $contig\n";
$contig = "$contig.more$scaftigminlength";

@temp = split /\//, $contig;

#print "perl $reformat $contig $outdir/index/$temp[-1].newformat\n";
system_("perl $reformat $contig $outdir/index/$temp[-1].newformat");

print STDERR "\nreformat complete!  ", $date = localtime, " \n";
$contig = "$outdir/index/$temp[-1].newformat";
foreach my $key ( keys %fq ) {
	die "Error in fq.list, 3 fastq file for each sample required: 1.fq, 2.fq, single.fq" unless ( @{ $fq{$key} } == 3 );
	my $t_q1 = $fq{$key}[0];
	my $t_q2 = $fq{$key}[1];
	my $t_qs = $fq{$key}[2];

	#print "perl $reads_mapping --ref $contig --q1 $t_q1 --q2 $t_q2 --qs $t_qs --outdir $outdir/bwa --bindir $bin_folder --processors $processors";
	system_("perl $reads_mapping --ref $contig --q1 $t_q1 --q2 $t_q2 --qs $t_qs --outdir $outdir/bwa --bindir $bin_folder --processors $processors");
}
print STDERR "reads_mapping complete!  ", $date = localtime, "\n";
if ( ( keys %fq ) >= 1 ) {
	$tmp = " $outdir/bwa/$sample_name\_allmerge.bam";
	foreach my $k ( keys %fq ) {
		foreach my $n ( @{ $fq{$k} } ) {
			@temp = split /\//, $n;
			if ( $temp[-1] =~ /^(.*pair.*)\.1\.fq/ ) {
				$tmp .= " $outdir/bwa/$1.sampe.bam";
			}
			elsif ( $temp[-1] =~ /^(.*single.*)\.fq/ ) {
				$tmp .= " $outdir/bwa/$1.samse.bam";
			}

			#$temp[-1] =~ m/(.*)\.fq/;
			#           $tmp .= " $outdir/bwa/$1.samse.bam";
		}
	}
	system_("$samtools merge $tmp");

	#print "$samtools merge $tmp\n";

	print STDERR "merge complete!  ", $date = localtime, "\n";

	#print "$samtools pileup -c -f $contig $outdir/bwa/$sample_name\_allmerge.bam >$outdir/bwa/$sample_name.xls\n";
	system_("$samtools pileup -c -f $contig $outdir/bwa/$sample_name\_allmerge.bam >$outdir/bwa/$sample_name.xls");

	print STDERR "pileup complete!  ", $date = localtime, "\n";
	@temp = split /\//, $contig;

	#print "perl $correct_it $outdir/bwa/$sample_name.xls $contig $outdir/correct_Result/$temp[-1].revised $cutoff $reformat $quality\n";
	system_("perl $correct_it $outdir/bwa/$sample_name.xls $contig $outdir/correct_Result/$temp[-1].revised $cutoff $reformat $quality");
	$contig = "$outdir/correct_Result/$temp[-1].revised";
	print STDERR "correct complete!  ", $date = localtime, "\n";
}
else { die "error 2 \n"; }

system_("rm $outdir/bwa/*sai -r");
system_("rm $outdir/bwa/*.bam -f");
system_("rm $outdir/bwa/*xls -f");

print "CONTIG: $contig\n";
unless ( -s $contig ) { die "ERROR & EXIT: Missing internal file $contig.\nAre there any scaftigs longer than $scaftigminlength bp in the assembly?\n" }
system_("cp -f $contig $outdir/index/");
print "cp -f $contig $outdir/index/\n";
$contig = "$outdir/index/$temp[-1].revised";
print "CONTIG: $contig\n";
print "$bwt $contig\n";
system_("$bwt $contig");
print "INSERT: $insert\n";
open I, $insert or die "$insert $! \n";

while (<I>) {
	@temp = split;
	$insert{ $temp[0] } = $temp[1];
}
close I;
open I, "$len_info" or die "$insert $!\n";
<I>;
while (<I>) {
	my @temp = split;
	$len{ $temp[0] } = $temp[1];
}
close I;
my ( $min, $max );
foreach my $key ( keys %fq ) {
	print STDERR"\n\n processing $key\n";
	die "Error in fq.list, 3 fastq file for each sample required: 1.fq, 2.fq, single.fq" unless ( @{ $fq{$key} } == 3 );
	my $t_q1 = $fq{$key}[0];
	my $t_q2 = $fq{$key}[1];
	my $t_qs = $fq{$key}[2];
	my @temp = split /\//, $t_q1;
	my $lane;
	if ( $temp[-1] =~ /1.fq/ ) {
		$lane = $temp[-1];
		$lane =~ s/.pair.1.*//;
	}
	elsif ( $temp[-1] =~ /2.fq/ ) {
		$lane = $temp[-1];
		$lane =~ s/.pair.2.*//;
	}
	else {
		$lane = $temp[-1];
		$lane =~ s/.single.*//;
	}
	system_("mkdir -p $outdir/soap/$sample_name/$lane") unless ( -e "$outdir/soap/$sample_name/$lane" );
	if ( defined $insert{"$sample_name\_$lane"} ) {
		$max = $insert{"$sample_name\_$lane"} + 50;
		$min = $insert{"$sample_name\_$lane"} - 50;
	}
	else {
		die "ERROR & EXIT: Insert size missing for $sample_name : $lane";
	}
	$temp[-3] = $sample_name;
	$temp[-2] = $lane;

	#my $v = int($len{$sample_name}*0.1)-2
	my $v = 8;    # SET BY JRK. 10% or read length - 2. This wasn't correctly set before anyway, I guess.
	              #print "$soap -a $t_q1 -b $t_q2 -D $contig.index -m $min -x $max -o $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.pe -2 $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.se -r 2 -l 30 -v $v -M 4 -p $processors 2>$outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.log\n";
	system_("$soap -a $t_q1 -b $t_q2 -D $contig.index -m $min -x $max -o $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.pe -2 $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.se -r 2 -l 30 -v $v -M 4 -p $processors 2>$outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.log");

	#print "$soap -a $t_qs -D $contig.index -o $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].single.soap.se -r 2 -l 30 -v $v -M 4 -p $processors 2>$outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].single.soap.log\n";
	system_("$soap -a $t_qs -D $contig.index -o $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].single.soap.se -r 2 -l 30 -v $v -M 4 -p $processors 2>$outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].single.soap.log");

	#print "perl $filter_pe -p $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.pe -i 0.9;rm $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.pe\n";
	system_("perl $filter_pe -p $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.pe -i 0.9");
	system_("rm $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.pe -f");

	#print "perl $filter_se -s $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.se -i 0.9;rm $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.se\n";
	system_("perl $filter_se -s $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.se -i 0.9");
	system_("rm $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.se -f");

	#print "perl $filter_se -s $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].single.soap.se -i 0.9;rm $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].single.soap.se\n";
	system_("perl $filter_se -s $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].single.soap.se -i 0.9");
	system_("rm $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].single.soap.se -f");

	#print "rm $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.pe -f\nrm $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].pair.soap.se -f\nrm $outdir/soap/$sample_name/$temp[-2]/$temp[-3]_$temp[-2].single.soap.se -f\n";
}

print STDERR "soap complete! ", $date = localtime, "\n";
system_("mkdir $outdir/soap/soap_list") unless(-e "$outdir / soap / soap_list ");
system_("find $outdir/soap/$sample_name/*/*filter*e >$outdir/soap/soap_list/$sample_name.all.trimsoap.list");

#print "$soap_coverage -cvg -refsingle $contig -p $processors -il $outdir/soap/soap_list/$sample_name.all.trimsoap.list -o $outdir/soap_coverage/$sample_name.cvg.out -depthsingle $outdir/soap_coverage/$sample_name.depth.out\n";
system_("$soap_coverage -cvg -refsingle $contig -p $processors -il $outdir/soap/soap_list/$sample_name.all.trimsoap.list -o $outdir/soap_coverage/$sample_name.cvg.out -depthsingle $outdir/soap_coverage/$sample_name.depth.out");
print STDERR "soap coverage complete! ", $date = localtime, "\n";

#print "perl $filter $contig $outdir/soap_coverage/$sample_name.depth.out $outdir/soap/soap_list/$sample_name.all.trimsoap.list $outdir/break_contig/$sample_name.scaftig.more$scaftigminlength.revised.break.fa $scaftigminlength $sample\n";
system_("perl $filter $contig $outdir/soap_coverage/$sample_name.depth.out $outdir/soap/soap_list/$sample_name.all.trimsoap.list $outdir/break_contig/$sample_name.scaftig.more$scaftigminlength.revised.break.fa $scaftigminlength $sample");

#print "bash $correct $scr_folder $outdir\n";
system_("bash $correct $scr_folder $outdir");

#print "bash $stat $scr_folder $outdir/break_contig/$sample_name.scaftig.more$scaftigminlength.revised.break.fa $sample_name $tatsOut\n";
system_("bash $stat $scr_folder $outdir/break_contig/$sample_name.scaftig.more$scaftigminlength.revised.break.fa $sample_name $tatsOut");

#`rm $outdir/soap/$name/$temp[-2]/$temp[-3]_$temp[-2].soap.pe -f`;
#`rm $outdir/soap/$name/$temp[-2]/$temp[-3]_$temp[-2].soap.se -f`;

print STDERR "All programs finished!\t", $date = localtime, "\n";

exit 0;

sub system_ {
	my $cmd = shift;
	print STDERR "SYSTEM CALL: $cmd";
	system($cmd) == 0 or die("\nERROR & EXIT: system($cmd) failed: $!\n");
}

