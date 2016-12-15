#!/usr/bin/env perl

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

#####################################################################################################################################
# UPDATES
#####################################################################################################################################
#
# 20161215: Fixed a bug in Funcitonal.pm around line 51
#
# 20160330: I changed how SAM input lengths are calculated
#
# 20151005: I noted that the old length file was still used, despite this info in map file. So changed that.
#           I also noted the .hashlevels were not needed when running NCBI or mOTU mode, so removed those.
#           Specifically ensuring that the length for NCBI in the map file are used, and also removed the hashlevel files,
#           as these ewre not actually used.
#
# 20151002: Added the date variable, so that rowname files are printed to the.date and then moved
#
# 20150815: Fixed this line (x2) in Main.pm @{ $multipleMapper{$i}{"$prev_insert"} }[$position2]       += ( 1 / $count ) * $value
#           Fixed minor bug, now not printing norm gene -1 in NCBI mode
# 20150814: v1.5.2 t16: fixed a bug that gave the wrong values for mm.dist.among.unique.norm for NCBI taxa.
#
# In version 1.4.6:t15 I fixed the bug that normalized length counts were incorrect for NCBI mode. In all previous modes the
# -1 wasn't correctly calculated. The number used instead of the correct one was the full -1, but it should've been divided by the
# average summarized taxa id length
#
# 050615:
# Added support for horizontal gene coverage. Currently stores a 100 char long string for each gene.
# This is stored for genes in HASH[INIT_NUM+1]
#
# 010715:
# Fixed a bug that calculated the wrong unannotated abundances for functional profiles and mm.dist.among.unique

#####################################################################################################################################
# INTRODUCTION
#####################################################################################################################################
# This is a copilation of 11 scripts that does the main calculation of MOCAT. Base and inert counts are calculated from a BAM/SAM/SOAP file.
# This version of the script (MOCAT v 1.4.0 and above) is faster than v 1.3 and below, and has been completely rewritten from scratch.
# Some features have been removed from within these scripts, e.g. the automatic making of a coord file. This should be done elsewhere.
# A number of bugs have been fixed compared to previous version, but mainly bases are now counted from the portion of the gene that is covered.
# v 1.3.8:
# GENE: -----------XXXXXXXXXXXXXX--------------
# READ:        ========
# Then the entire portion above was used for calculaing base coverages. In this evrsion, the portion below is used:
# v 1.4.0+:
# GENE: -----------XXXXXXXXXXXXXX--------------
# READ:            ====
# This is more correct, but also means that if a read has multiple mappers, mapping to multiple targets, it may be that 100% of the read isn't actually counted for.
#
# Before short genes, less than 500 bp were over estimates from 2.1x to 1.2x, now the effect is an under estimate of 0.8x max (min)
#
# These are the bug fixes that have been made since v 1.3.8
# 1)  The scaling factor was incorrect for mm dist files, affecting both base and insert counts. The effect of this bug is hard to
#    tell, but in general all scaled numbers would be somewhat slightly affected.
# 2) The total mapped in functional headers referred to the numbers in the abundance files without PE filter, this affected the raw,
#    norm, scaled files for both abse and inserts. This means the total number of bases written was the bases before paired end filtering,
#    but it should've been the bases after PE filtering. The effect is probably minimal.
# 3) Previously the scaling factor for cog, ko, module, pathway, taxonomic levels were based on total cog/ko/... abundances + the unannotated
#    fraction. I think now the scaling factors should be based on only the total cog/ko/... abundance. This affects only the scaled files for
#    both insert and bases for functional and taxonomic profiles. The effect ranges from minimal to somewhat higher.
# 4) A line specific bug (line 1625) meant that the raw mm dist for bases sum not annotated header line was incorrect
# 5) A minor bug in the code because some variables were not reset, when distributing to multiple mappers to 2 taxa or more.
# 6) If the same read matched the same reference more than once, the value of the read was only counted as the last match in the file,
#    now this has been fixed so that if a read matches the same reference multiple of times, the highest scoring match is counted.
#    We don't want to count the same read twice, or take an average.
#
# Combining bugs 1-4, the effect is the following on correltions between v 1.3.8 and 1.4.0 (tested on one random sample):
# Correlation is > .99, but the actual values, with a certain number of signif digits is not 100%
# base norm mm.dist.among.unique. gene With 3 significant digits 90.372 %
# base scaled only.unique. gene With 3 significant digits 17.365 % (98.443 % with 1 digit)
# And similar effects for inserts. Generally raw and norm files have identical values, and mm dist files have more different values, and scaled
# files have further different values. Also, the higher levels (cog, ko, mod, pathway, kingdom, phylum, class, ...) are affected more because
# they are calculated form a number of genes

#####################################################################################################################################
# MORE RANDOM NEWS
#####################################################################################################################################
# It's possible to specifu coord file like: file.coord, but only file.coord.gz exists, and the it will load the gz version.
# The same applies for the map file

#####################################################################################################################################
# GENERAL PROGRAM LAYOUT
#####################################################################################################################################
# 1. Run Initalize::initialize to load stuff and check files
# 2. Run Load::loadStatsFile to load stats file. This file may be missing
# 3. Open input file for reading
# 3. Create the threads for processing the input file, MOCATProfiling::Main2::main. It's important this is done before loading the coordinate files
# 4. Load coord files
# 5. If in NCBI mode, load map files
#
# Part 1
# 6. Launch MainQueuer::mainQueuer($command)
# 6.1 Each line is read, and then after X number of lines these are packaged, and the needed coord files are also saved and then these are enqueued for the Main2::main queuer to process this batch
# 6.2 While running, here is also a check to see if some files have finished, and if so, the results are processed using MOCATProfiling::Threading::merge($file);
#     The merge subroutine adds data to the main HASH and HASH_LEVELS
# 6.3 Once all threads have finished, cadd up total sums using MOCATProfiling::Threading::parseThreadResult( $returned, $i );
# 6.4 Run MOCATProfiling::Threading::merge($file); for the jobs that haven't been merged
#
# Part 2
# 7. Print some stats and die if nothing was mapped
#
# Part 3
# 8. Distribute multiple mappers, this is done for each thread and each level. For gene and funcitonal mode the level is only gene. For mOTU mode it's gene and mOTU,
#    and for NCBI mode it's gene, kingdom, phylum, ... levels
#
# Part 4
# 9. Print gene abundance files
# 10. If mode is funcitonal, we print the functional levels
#
# Part 5
# 11. Print taxonomic levels, e.g. mOTU level for mOTU mode, and taxonomic levels for NCBI mode
#
# Part 6
# 12. Print stats files
#
# Part 7
# 13. Move files securely
# 14. Done with program

#####################################################################################################################################
# HASH AND HASH_LEVELS LAYOUT
#####################################################################################################################################
#The HASH and HASH_LEVELS are the two main hashes, in which data is stored. The layout of these differ slightly between the different
#modes gene, funcitonal, motu and NCBI. The information below may not be 100% correct, but will hopefully point to a general idea of the layout.
#The HASH contains an entry for each gene, and this will be true for all modes. HASH_LEVELS contains a hash for each of the other levels,
#e.g. cog, ko, mod and pathway for funcitonal and each of the taxonomic levels for NCBI mode, and mOTU, as a single level for mOTU mode.
#
# An array entry in the HASH, for normal gene mode, looks like this:
# 0: start position
# 1: stop position
# 2: base raw counts
# 3: insert raw counts
# 4: base unique counts
# 5: insert unique counts
# 6: base mm dist raw, these are added in the DistMM function
# 7: insert mm dist raw, these are added in the DistMM function
# 8: base mm dist norm, this is not added in the Main loop, but in the DistMM function
# 9: insert mm dist norm, this is not added in the Main loop, but in the DistMM function
# 10: for mOTU, this is the ID the gene map to :: @{ $HASH{"$line[0].0"} }[$INIT_NUM] = "$line[2].$line[3]"
# 11: horizontal gene coverage, a string 100 char long
# The norm values are not stored, except for mm dist, because these come from multiple sources
#
# For HASH_LEVELS, the array layout is:
# 0: start position
# 1: stop position
# 2: base raw counts
# 3: insert raw counts
# 4: base unique counts
# 5: insert unique counts
# 6: base mm dist raw, these are added in the DistMM function
# 7: insert mm dist raw, these are added in the DistMM function
# 8: base mm dist norm, these are added in the DistMM function
# 9: insert mm dist norm, these are added in the DistMM function
# 10: base norm
# 11: insert norm
# 12: base unique only norm
# 13: insert unique only norm
# 14: ??? mm dist among unique base norm ???
# 15: ??? mm dist among unique insert norm ???
# 16: horizontal coverage (total)
# 17: horizontal coverage (number of COGs matched)
# For taxonomic levels, the start and stop are both = 1, because the length is stored in the %MAP hash.

#####################################################################################################################################
# VARIABLE DESFRIPTIONS
#####################################################################################################################################
# @fields : contains the index field to load for a specific taxa in the MAP file
# $input_file_format : this can be either, BAM och SOAP, or SAM, but SAM is mostly used while developing
# $INIT_NUM : this is defined as the first position in the vector of the hash, where maps are started,
#             this value is unknown when starting scripting this script, but will become apparent
# %HASH : this is the main hash where EVERYTHING is stored, from coordinates for each gene, to the values to the mappings
#         currently in position $INIT_NUM and forward are the mapping saved. Note that the mapping are for KEGG comma
#         seperated KOs, MODs and PATHs
# @STATS_FILE_DATA : contains the data in the stats file
# $reference_for_read_exists_in_coord_file : keeps track how many of the reads have a reference
# $reference_for_read_do_not_exists_in_coord_file : keeps track of how many dont have a ref
# $TOTAL_LENGTH_MAPPED : the total length mapped
# $seen_types{$i}{$ref} : this is a hash that is reset for each new insert, it stores which e.g. taxa has been seen at the different levels
# %intersection : a temporary hash that contains the taxa at individual levels that passes the paired end filtering
# %PEaffectedInserts : this hash contains stats on how many genes on each level were affected by the paired end filtering
# $TOTAL_INSERTS_MAPPED2, $TOTAL_INSERTS_MAPPED3 : a development variable to see if this gives same result as $TOTAL_INSERTS_MAPPED
# $INSERT_COUNTER : this is a counter that keeps in increasing throughout the script, it's used in the multipleMapper and multipleMapperValue hashes,
#                  for exactly what, I can't remember...
# %LENGTHS : if the mode is RefMG or moTU the lengths are stored here, otherwise tey are retrieves from the HASH[1]-HASH[0]+1
# $ALL_TAXA_MEAN_LENGTH : Used in RefMG mode
# $countThisTaxa : internal vairbale used to check whether we should count this taxa, in use if PE_filter enabled
# %HASH_LEVELS : this is the main HASH where speicies up to family, and mOTU and final cog and ko abundances are stored
# %multipleMapperValue : this is used for storing the total value of count that should be distributed to multiple mappers
# %multipleMapper : this stores hashes with the taxa that should have the value distributed
# $FRACTION_MAPPED_BASES : well...
# %multipleMapperMap : this is a hash that contains a hash reference for each level, the idea is to avoid having hashes of hashes for multipleMapper hash
#                      and instead use hash references, we'll see if it works. hopefully this uses a lot less memory
# $temp_file : if this is specified we save MM data in this file rather than in arrays in memory, and then load this file line by line (it's actually one file per level)
# %open_temp_files : a list of filehandles for open temp files, one for each level
# %did_PE_filtering : i had to introduce this variable to keep track of whether PE filtering for a specific level was performed or not, just checking if intersection>0 isn't enough, because it could be wrong
# $temp_folder : if specified as input, we create this folder, useful for when jobs are submitted to another node that initialized on
# @LEVELS : only created for funciotnal profiles, contains the levels except gene, that is going to be summarized into, this is then saved as @levels
#####################################################################################################################################

#####################################################################################################################################
# LOAD USE
#####################################################################################################################################
use strict;
use warnings;
use threads;
use Storable;
use threads::shared;
use MOCATProfiling::QueueJ;
use Getopt::Long;
use MOCATProfiling::Variables;
use MOCATProfiling::Threading;
use MOCATProfiling::DistMM;
use MOCATProfiling::Misc;
use MOCATProfiling::Main;
use MOCATProfiling::Load;
use MOCATProfiling::Functional;
use MOCATProfiling::Print;
use MOCATProfiling::Initialize;
use MOCATProfiling::MainQueuer;
use File::Sync qw(fsync sync);
#####################################################################################################################################

#####################################################################################################################################
# MAKE STDOUT AND STDERR NON CHACHED
#####################################################################################################################################
select(STDERR);
$| = 1;
select(STDOUT);    # default
$| = 1;
#####################################################################################################################################

#####################################################################################################################################
# DEFINE VARIABLES
#####################################################################################################################################
$reference_for_read_exists_in_coord_file        = 0;
$reference_for_read_do_not_exists_in_coord_file = 0;

#$gene_output_folder                             = "";
#$profiles_output_folder                         = "";
$TOTAL_LENGTH_MAPPED = 0;
$INSERT_COUNTER      = 0;
$INIT_NUM            = 10;
$prev_insert         = "";
$seen_this_insert    = 0;
$THREADS             = 1;
$temp_file           = "MP.tmp";
$LOADED_MAP          = 0;          # this is set to 1 once the map file has been loaded, and then processing the files can start
$BLOCKSIZE           = 2000000;

# DEFINE SHARED VARIABLES
share($TERMINATE);
share($LOADED_MAP);
share($FINISHED_THREADS);
$TERMINATE        = 0;
$FINISHED_THREADS = 0;
#####################################################################################################################################

#####################################################################################################################################
# INITIAL SETUP
#####################################################################################################################################
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : ============ MOCAT MAGICAL MULTI THREADING PROFILING SCRIPT FOR GETTING YOUR PROFILES. MOSTLY BUG FREE.  ============\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : ============================== PART 0 : initialize and load stuff, always a good idea! ==============================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";

# SAVE ARGS, USED FOR PRINTING LATER
@argv = @ARGV;

# LOAD OPTIONS
GetOptions(
	'coord_files=s{,}'       => \@coord_files,
	'input_file=s{,}'        => \@input_files,
	'stats_file=s'           => \$stats_file,
	'temp_file=s'            => \$temp_file,
	'map_file=s'             => \$map_file,
	'output_file=s'          => \$output_file,
	'rownames_file=s'        => \$rownames_file,
	'mode=s'                 => \$mode,
	'samtools_executable=s'  => \$samtools_executable,
	'PE_filter=s'            => \$PE_filter,
	'sample_name=s'          => \$sample_name,
	'zip=s'                  => \$zip_execute,
	'threads=s'              => \$THREADS,
	'blocksize=s'            => \$BLOCKSIZE,
	'old_functional_headers' => \$USE_OLD_HEADERS,
	'out_stats_file=s'       => \$out_stats_file,
	'out_PE_stats_file=s'    => \$out_PE_stats_file,
	'temp_folder=s'          => \$temp_folder,
	'print_rownames_file'    => \$print_rownames_file,
	'sam'                    => \$sam,
	'bin:s'                  => \$bin_dir,
	'version:s'              => \$version,
	'verbose'                => \$VERBOSE,
	'firstColumnName:s'      => \$firstColumnName,
	'levels:s{,}'            => \@CONFIG_LEVELS,
	'horizon'                => \$CALCULATE_HORIZONTAL_COVERAGE,
	'date=s'                 => \$date
);

# INITIALIZE AND LOAD STATS FILE
MOCATProfiling::Initialize::initialize();
MOCATProfiling::Load::loadStatsFile();

####################################################################################################################################
# PART 1 get gene abundances and, if applicable, abundances for different taxonomic levels
#####################################################################################################################################
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : ==================================== PART 1 : calculate abundances in each thread ===================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";

# OPEN INPUT
my $i = 0;    # this is 0 because, we only load the first entry in the input_files array, because right now only support one input file
my $command;
if ( $input_file_format eq 'BAM' ) {
	$command = "$samtools_executable view $input_files[$i] | ";
}
elsif ( $input_file_format eq 'SAM' ) {
	if ( $input_files[$i] =~ m/\.gz$/ ) {
		$command = "gunzip -c $input_files[$i] | ";
	}
	else {
		$command = "<$input_files[$i]";
	}
}
elsif ( $input_file_format eq 'SOAP' ) {
	$command = "gunzip -c $input_files[$i] | ";
}
open $READ, "$command" or die "ERROR & EXIT: Could not execute $command";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [MAIN THREAD] OPEN $command\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [MAIN THREAD] Submitting inserts to threads\n";

######### THREADING and RUNNING #########
# CREATE THREAD HANDLES
for my $thread ( 1 .. $THREADS ) {
	my $OUT;
	open( $OUT, $READ );
	$threadHandles[$thread] = $OUT;
}

# This is one queue per THREAD
$request_q = MOCATProfiling::QueueJ->new();
$return_q  = MOCATProfiling::QueueJ->new();
my @threads;
for my $i ( 1 .. $THREADS ) {
	push @threads, threads->create( \&MOCATProfiling::Main::main );
}
######### THREADING and RUNNING #########

######### LOAD MAPS #########
# this loads the levels into @LEVELS;
if ( $mode eq 'functional' ) {
	MOCATProfiling::Load::loadFunctLevels();
}

# we have to load these AFTER we've dispatched the threads, coord file needs to be loaded before map file
MOCATProfiling::Load::loadCoordFile(-1);    # we just set -1 because of semantic for not printing text THREAD X later on, but also because we want to loop over all entries, see Loader script.

# Load map files, but not if it's functional mode, then we load them later on the fly and only for those data points with values
if ( $mode ne 'functional' && $mode ne 'gene' ) {
	MOCATProfiling::Load::loadMapFiles();
}
######## LOAD MAPS #########

######### WAITING WHILE ALL THREADS FINISHED #########
# this is the thread that sends jobs to other threads, we expect and want it to finish first
# currently this thread has a copy of the large HASH, hopefully this memory is cleared before we continue downstream

# create jobs
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Processing input file\n";
MOCATProfiling::MainQueuer::mainQueuer($command);

# This tells the main program to keep running until all threads have finished.
foreach my $i ( 0 .. $THREADS - 1 ) {
	my $returned = $threads[$i]->join();
	unless ($returned) {
		die "ERROR & EXIT: Main routine shutdown because thread " . ( $i + 1 ) . " failed";
	}
	MOCATProfiling::Threading::parseThreadResult( $returned, $i );
}

# this takes care of the summaries, that have been saved while the threads have been running
$return_q->enqueue(undef);
while ( defined( my $file = $return_q->dequeue() ) ) {
	MOCATProfiling::Threading::merge($file);
}

# CLOSE PIPES
close $READ or die("ERROR & EXIT: Closing INPUT PIPE ($command) failed, this means something is wrong. Corrput files?");
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [MAIN THREAD] CLOSE $command\n";
@threads = ();    # EMPTY THREADS
######### WAITING WHILE ALL THREADS FINISHED #########
####################################################################################################################################

####################################################################################################################################
# PART 2 print some simple stats, and if nothing matched, die, and close some files
####################################################################################################################################
unless ( $TOTAL_MAPPED{'base'}{'gene'} ) {
	$TOTAL_MAPPED{'base'}{'gene'}   = 0;
	$TOTAL_MAPPED{'insert'}{'gene'} = 0;
}

print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =============================== PART 2 : some basic stats summarized for all threads ================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
print STDERR <<"EOF";
                                       Uncounted reads                 : $reference_for_read_do_not_exists_in_coord_file
    ("`-''-/").___..--''"`-._          Counted   reads                 : $reference_for_read_exists_in_coord_file
     `6_ 6  )   `-.  (     ).`-.__.`)  TOTAL BASES   MAPPED (gene)     : $TOTAL_MAPPED{'base'}{'gene'}
     (_Y_.)'  ._   )  `._ `. ``-..-'   TOTAL BASES   (from stats file) : $TOTAL{'base'}
   _..`--'_..-_/  /--'_.' ,'           TOTAL INSERTS MAPPED (gene)     : $TOTAL_MAPPED{'insert'}{'gene'}
  (il),-''  (li),'  ((!.-'             TOTAL INSERTS (from stats file) : $TOTAL{'insert'}
EOF

# Maybe die?
if ( $TOTAL_MAPPED{'base'}{'gene'} > 0 ) {
}
else {
	warn "WARNING: Total bases mapped was 0. Most likely no reads were mapped to the databases because the dataset was small or the number of reads to start with were not that many.";
}

# if we are using temp_file, which we ARE, because, I only support that, we need to close the files now
if ($temp_file) {
	foreach my $level (@levels) {
		for my $i ( 1 .. $THREADS ) {
			close $open_temp_files{"$level.$i"} or die "ERROR & EXIT: Closing pipe failed, this means something is wrong. Corrput files?";
		}
	}
}
####################################################################################################################################

####################################################################################################################################
# PART 3 distribute multiple mappers among genes or taxonomic levels
####################################################################################################################################
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : ================================ PART 3 : distribute multiple mappers for each thread ===============================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";

# DIST MMs without threads
for my $thread ( 1 .. $THREADS ) {
	MOCATProfiling::DistMM::distributeMultipleMappersThreadWrapper($thread);    # here we distribute the multiple mappers
}
####################################################################################################################################

####################################################################################################################################
# PART 4 if it's funcitonal we need to print the gene files and then summarize it at the different functional levels
####################################################################################################################################
# if mode is functional we have to summarize at the different levels other than gene
# before printing we have to reload all the entries for the coord files
if ( $mode eq 'functional' ) {
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =========================== PART 4 : print gene files and also summarize functional levels ==========================\n";
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";

	# we need to print the gene files before we do the summary step, because the summary step removes them from the main HASH, this has to be done in order to gte the stats for sum_not_annotated
	MOCATProfiling::Print::printFiles( \%HASH, 'gene', 1 );    # this prints the raw, norm and scaled files for the abundances for specified level level (gene level only for functional)

	# dist to HASH_LEVELS
	#@levels = ( 'cog', 'ko', 'module', 'pathway' );                                                     # before we had levels set to only genes for functional because it's on the gene level we distribute multiple mappers
	@levels = @LEVELS;
	MOCATProfiling::Functional::summarizeFunctionalLevels();
}
else {
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : ============================================= PART 4 : print gene files  ============================================\n";
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
	MOCATProfiling::Print::printFiles( \%HASH, 'gene', 1 );    # this prints the raw, norm and scaled files for the abundances for specified level level (gene level only for functional)
	@levels = @levels[ 1 .. scalar @levels - 1 ];           # lets remove gene from the levels

}
####################################################################################################################################

####################################################################################################################################
# PART 5 print taxonomic files
####################################################################################################################################
if ( scalar @levels >= 1 ) {
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : ========================================== PART 5 : print additional files  =========================================\n";
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
	foreach my $level (@levels) {
		MOCATProfiling::Print::printFiles( $HASH_LEVELS{$level}, $level, 6 );    # the gene abundances are stored in HASH, whereas eg species and cog abundances and other levels are stored in HASH_LEVELS
	}
}
####################################################################################################################################

####################################################################################################################################
# PART 6 print stats files
####################################################################################################################################
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =========================================== PART 6 : print statistics files  ========================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
MOCATProfiling::Print::printStats();
####################################################################################################################################

####################################################################################################################################
# PART 7 move files
####################################################################################################################################
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =========================================== PART 7 : move all created files  ========================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
chomp( my $file_size_all = `du -hsLc $temp_file* | tail -n 1 | cut -f 1` );
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : output files and other temp files uses $file_size_all in $temp_file*\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : tar gz gene files\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : started zipping files\n";
MOCATProfiling::Misc::zipFiles("$output_file.zip");
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : completed zipping files\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : remove temp files $temp_file*\n";
system "rm -fr $temp_file*";
chomp( $file_size_all = `du -hsLc $temp_file* 2>/dev/null | tail -n 1 | cut -f 1` );
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : output files and other temp files uses $file_size_all in $temp_file*\n";
####################################################################################################################################

####################################################################################################################################
# DONE
####################################################################################################################################
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : ========================= THE END. Seriously. It's over now. Enjoy the sun! Have an ice cream! ======================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
exit 0;
####################################################################################################################################

