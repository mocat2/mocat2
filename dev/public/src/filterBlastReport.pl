#! /usr/bin/env perl
use strict;
use warnings;
use Smash::Utils::BLAST;
use Getopt::Long;

my $opt_flavor   = "WU";
my $opt_bits     = 40;
my $opt_exp      = 10;
my $opt_qalnlen  = 10;
my $opt_salnlen  = 10;
my $opt_qalnfrac = 0.0;
my $opt_salnfrac = 0.0;
my $opt_sim      = 10;
my $opt_id       = 10;
my $opt_top_bits;
my $opt_sbjct;
my $opt_qlength_file;
my $opt_slength_file;
my $opt_check_complete = 0;
my $opt_nofilter = 0;
my $opt_fasta;
my $opt_help;

my $usage = "
$0 - Filter WU-BLAST output on combinations 
of bit score, E-value, aligned query length, percent similarity
and max subjects (cf -B option in WU-BLAST). Subjects are sorted by 
individual HSP scores and the best <n> subjects are chosen.

Usage: $0 [options] <blast-file> 

Program Options:
 --flavor=<string>     flavor of BLAST used - one of WU or NCBI (default: $opt_flavor)
 --bits=<num>          bit threshold for results to be included (default: $opt_bits)
 --topbits=<frac>      bit threshold in terms of percent distance from best hit (default: none)
 --exp=<num>           E-value threshold for results to be included (default: $opt_exp)
 --sim=<num>           percent similarity threshold for hsp to be included (default: $opt_sim)
 --id=<num>            percent identity threshold for hsp to be included (default: $opt_id)
 --sbjct=<num>         max number of subjects whose hsp hits should be reported (default: none)
 --minqlen=<num>       minimum number of query bases for hsp to be included (default: $opt_qalnlen)
 --minqfrac=<frac>     fraction of the query that must be aligned (default: $opt_qalnfrac)
 --qlen=<file>         file containing lengths of the query sequences (required for --minqfrac)
 --minslen=<num>       minimum number of subject bases for hsp to be included (default: $opt_salnlen)
 --minsfrac=<frac>     fraction of the subject that must be aligned (default: $opt_salnfrac)
 --slen=<file>         file containing lengths of the subject sequences (required for --minsfrac)
 --check_complete      check if the blast file is complete (default: false)
 --nofilter            do not perform filter (default: false)
		        can be combined with --check_complete to only check completeness 
 --fasta=<file>        fasta file containing query sequences to check completeness (required for --check_complete)
";

GetOptions(
	"flavor=s"   => \$opt_flavor,
	"bits=n"     => \$opt_bits,
	"exp=f"      => \$opt_exp,
	"topbits=f"  => \$opt_top_bits,
	"qlen=s"     => \$opt_qlength_file,
	"minqlen=n"  => \$opt_qalnlen,
	"minqfrac=f" => \$opt_qalnfrac,
	"slen=s"     => \$opt_slength_file,
	"minslen=n"  => \$opt_salnlen,
	"minsfrac=f" => \$opt_salnfrac,
	"sbjct=n"    => \$opt_sbjct,
	"sim=n"      => \$opt_sim,
	"id=n"       => \$opt_id,
	"check_complete" => \$opt_check_complete,
	"nofilter"   => \$opt_nofilter,
	"fasta=s"    => \$opt_fasta,
	"help"       => \$opt_help);

my ($blast_file) = @ARGV;

if ($opt_help || !defined($blast_file)) {
	die $usage;
}

if (! ($blast_file eq "-" || -f $blast_file)) {
	warn "File $blast_file not found";
	die $usage;
}

if ($opt_qalnfrac && !$opt_qlength_file) {
	warn "--qlen is required with --minqfrac";
	die $usage;
}

if ($opt_salnfrac && !$opt_slength_file) {
	warn "--slen is required with --minsfrac";
	die $usage;
}

my $complete = 1;
if ($opt_check_complete) {
	$complete = Smash::Utils::BLAST::check_complete_blast($blast_file, $opt_fasta, $opt_flavor);
}
#HACK
#$opt_sbjct = undef;
#END_HACK
if ($complete == 1 && !$opt_nofilter) {
	Smash::Utils::BLAST::write_sorted_blast_report(
						"BLAST_FILE"  => $blast_file, 
						"OUT_FH"      => \*STDOUT, 
						"FLAVOR"      => $opt_flavor,
						"BITSCORE"    => $opt_bits, 
						"TOPBITS"     => $opt_top_bits, 
						"EXP"         => $opt_exp,
						"SIM"         => $opt_sim,
						"IDENTITY"    => $opt_id,
						"MAXSUBJECTS" => $opt_sbjct,
						"MINQLEN"     => $opt_qalnlen,
						"MINQFRAC"    => $opt_qalnfrac,
						"QLENFILE"    => $opt_qlength_file,
						"MINSLEN"     => $opt_salnlen,
						"MINSFRAC"    => $opt_salnfrac,
						"SLENFILE"    => $opt_slength_file
						);
}
