use strict;
use warnings;
use Getopt::Long;
use File::Sync qw(fsync sync);
use MOCATProfiling::Variables;
use MOCATProfiling::Threading;
use MOCATProfiling::DistMM;
use MOCATProfiling::Misc;
use MOCATProfiling::Main2test;
use MOCATProfiling::Load;
use MOCATProfiling::Functional;
use MOCATProfiling::Print;
use MOCATProfiling::Initialize;
use MOCATProfiling::MainQueuer;
use Storable;

select(STDERR);
$| = 1;
select(STDOUT);    # default
$| = 1;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

#####################################################################################################################################
# Define variables

$reference_for_read_exists_in_coord_file        = 0;
$reference_for_read_do_not_exists_in_coord_file = 0;
$base_output_folder                             = "";
$insert_output_folder                           = "";
$profiles_output_folder                         = "";
$TOTAL_LENGTH_MAPPED                            = 0;
$INSERT_COUNTER                                 = 0;
$INIT_NUM                                       = 10;
$prev_insert                                    = "";
$seen_this_insert                               = 0;
$THREADS                                        = 1;
$temp_file                                      = "MP.tmp";
$LOADED_MAP                                     = 0;          # this is set to 1 once the map file has been loaded, and then processing the files can start
$BLOCKSIZE                                      = 2000000;
$TERMINATE        = 0;
$FINISHED_THREADS = 0;

#$SIG{INT} = $SIG{TERM} = \&MOCATProfiling::Threading::signal_handler;
#####################################################################################################################################
# Variable descriptions
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
#####################################################################################################################################

#####################################################################################################################################
# Initial setup

print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : ============ MOCAT MAGICAL MULTI THREADING PROFILING SCRIPT FOR GETTING YOUR PROFILES. MOSTLY BUG FREE.  ============\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : ============================== PART 0 : initialize and load stuff, always a good idea! ==============================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";

@argv = @ARGV;

# Load options
GetOptions(
	'coord_files=s{,}'         => \@coord_files,
	'input_file=s{,}'          => \@input_files,
	'stats_file=s'             => \$stats_file,
	'temp_file=s'              => \$temp_file,
	'map_files=s'              => \$map_file,
	'output_file=s'            => \$output_file,
	'rownames_file=s'          => \$rownames_file,
	'refmg_length_file'        => \$refmg_length_file,
	'mode=s'                   => \$mode,
	'samtools_executable=s'    => \$samtools_executable,
	'PE_filter=s'              => \$PE_filter,
	'sample_name=s'            => \$sample_name,
	'zip=s'                    => \$zip_execute,
	'base_output_folder=s'     => \$base_output_folder,
	'insert_output_folder=s'   => \$insert_output_folder,
	'profiles_output_folder=s' => \$profiles_output_folder,
	'threads=s'                => \$THREADS,
	'blocksize=s'              => \$BLOCKSIZE,
	'old_functional_headers'   => \$USE_OLD_HEADERS,
	'NCBI_length_file=s'       => \$NCBI_length_file,
	'out_stats_file=s'         => \$out_stats_file,
	'out_PE_stats_file=s'      => \$out_PE_stats_file
	
);
MOCATProfiling::Initialize::initialize();

# End initial setup
####################################################################################################################################

####################################################################################################################################
# PART 1 get gene abundances and, if applicable, abundances for different taxonomic levels

print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : ================================= PART 1 : calculate gene abundances in each thread =================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";

#$LOADED_MAP = 1;

# open input
#!# $i=0 hard coded for now, need to change to support only one input file perhaps
my $i = 0;
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

#for my $thread ( 1 .. $THREADS ) {
#	my $OUT;
#	open( $OUT, $READ );
#	$threadHandles[$thread] = $OUT;
#}
#
########## THREADING and RUNNING #########
#
## THREADING
## THREADING
## THREADING
## This is one queue per THREAD
#$request_q = Thread::Queue->new();
#$return_q  = Thread::Queue->new();
#my @threads;
#for my $i ( 1 .. $THREADS ) {
#	push @threads, threads->create( \&MOCATProfiling::Main::main );
#}
## THREADING
## THREADING
## THREADING

#NON THREADING
#my @threads = ();
#my $returned = MOCATProfiling::Main::main;
#MOCATProfiling::Threading::parseThreadResult( $returned, $i, 'abund' );
#NON THREADING



# This is the main queue that handles submitted jobs to threads
# currently we load the coord file in here, and then it's released from within here and then reloaded in the
#$main_q = Thread::Queue->new();
#my @mainThread;
#push @mainThread, threads->create( \&MOCATProfiling::MainQueuer::mainQueuer, $command );

######### THREADING and RUNNING #########

######### LOAD MAPS #########
# we have to load these AFTER we\ve dispatched the threads
# coord file needs to be loaded before map file
MOCATProfiling::Load::loadCoordFile(-1);    # we just set -1 because of semantic for not printing text THREAD X later on, but also because we want to loop over all entries, see Loader script.
# Load map files, but not if it's functional mode, then we load them later on the fly and only for those data points with values
if ( $mode ne 'functional' && $mode ne 'gene') {
	MOCATProfiling::Load::loadMapFiles();
}
######## LOAD MAPS #########

######### WAITING WHILE ALL THREADS FINISHED #########
# this is the thread that sends jobs to other threads, we expect and want it to finish first
# currently this thread has a copy of the large HASH, hopefully this memory is cleared
# before we continue downstream
#my $returned = $mainThread[0]->join();    # join mainQueue, once the others have finished
#unless ($returned) {
#	die "ERROR & EXIT: Main routine shutdown because main queue failed";
#}
#@mainThread = ();


## THREADING
## THREADING
## THREADING
## create jobs
#MOCATProfiling::MainQueuer::mainQueuer($command);
#
## This tells the main program to keep running until all threads have finished.
#foreach my $i ( 0 .. $THREADS - 1 ) {
#	my $returned = $threads[$i]->join();
#	unless ($returned) {
#		die "ERROR & EXIT: Main routine shutdown because thread " . ( $i + 1 ) . " failed";
#	}
#	MOCATProfiling::Threading::parseThreadResult( $returned, $i, 'abund' );
#}
#
## this takes care of the summaries, that have been saved while the threads have been running
#$return_q->enqueue(undef);
#while ( defined( my $file = $return_q->dequeue() ) ) {
#	MOCATProfiling::Threading::merge($file);
#}
## THREADING
## THREADING
## THREADING

######### WAITING WHILE ALL THREADS FINISHED #########

MOCATProfiling::Main2test::main();
close $READ or die("ERROR & EXIT: Closing INPUT PIPE ($command) failed, this means something is wrong. Corrput files?");
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [MAIN THREAD] CLOSE $command\n";

####### SUMMARIZE HASH local files into HASH ######
#foreach my $file (@RETURNED) {
#	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [MAIN THREAD] LOADING $file\n";
#	my $hashref = retrieve($file);
#	foreach my $key1 ( keys %{$hashref} ) {
#		my @array = @{ $$hashref{$key1} };
#		for my $i ( 0 .. scalar @array - 1 ) {
#			if ( $array[$i] ) {
#				@{ $HASH{$key1} }[$i] += $array[$i];
#			}
#		}
#	}
#}
#
####### SUMMARIZE HASH local files into HASH ######

# Add the start and stop and fill in the HASH
# Load map
#@threads = ();

#not needed if we run mainQueuer without a thread
#MOCATProfiling::Load::loadCoordFile(0);

####################################################################################################################################

####################################################################################################################################
# PART 2 print some simple stats, and if nothing matched, die, and close some files
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
	die "ERROR & EXIT: Total bases mapped was 0. Most likely no reads were mapped to the databases because the dataset was small.";
}

# if we are using temp_file, which we ARE, because, I only support that, we need to close the files now
if ($temp_file) {
	foreach my $level (@levels) {
		for my $i ( 1 .. $THREADS ) {
			close $open_temp_files{"$level.$i"} or die "ERROR & EXIT: Closing pipe failed, this means something is wrong. Corrput files?";
		}
	}
}

#		# TEMP
#		foreach my $key1 ( keys %HASH ) {
#				print "HASH $key1 : ", join(" " ,@{$HASH{$key1}} ) . "\n";
#		}

####################################################################################################################################

####################################################################################################################################
# PART 3 distribute multiple mappers among genes or taxonomic levels
#foreach my $key1 ( keys %HASH_LEVELS ) {
#	foreach my $key2 ( keys %{$HASH_LEVELS{$key1}} ) {
#				print "HASH_LEVELS $key1:$key2 : ", join("\t" ,@{$HASH_LEVELS{$key1}{$key2}} ) . "\n";
#	}
#		}

print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : ================================ PART 3 : distribute multiple mappers for each thread ===============================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
###############################################
# ALT 1
##### this is for running DistMM using multiple threads, which we don't want to do anymore because it uses a lot of memory
#### NOTE that this implementation will not work properly for mm dist values unless the HASH is shared, which it isn't anymore
# # Create and start threads for main sub
# @threads = MOCATProfiling::Threading::initThreads();
# for my $i ( 0 .. $THREADS - 1 ) {
# 	my $thread = $i + 1;
# 	$threads[$i] = threads->create( \&MOCATProfiling::DistMM::distributeMultipleMappersThreadWrapper, $thread );
# }
# # This tells the main program to keep running until all threads have finished.
# foreach my $i ( 0 .. $THREADS - 1 ) {
# 	my $returned = $threads[$i]->join();
# 	unless ($returned) {
# 		die "ERROR & EXIT: Main routine shutdown because thread " . ( $i + 1 ) . " failed";
# 	}
# 	MOCATProfiling::Threading::parseThreadResult( $returned, $i, 'mm dist' );
# }
###############################################

###############################################
# ALT 2, dist MMs without threads
for my $thread ( 1 .. $THREADS ) {
	MOCATProfiling::DistMM::distributeMultipleMappersThreadWrapper($thread);    # here we distribute the multiple mappers
}
###############################################

#foreach my $key1 ( keys %HASH_LEVELS ) {
#	foreach my $key2 ( keys %{$HASH_LEVELS{$key1}} ) {
#				print "HASH_LEVELS $key1:$key2 : ", join("\t" ,@{$HASH_LEVELS{$key1}{$key2}} ) . "\n";
#	}
#		}

####################################################################################################################################

####################################################################################################################################
# PART 4 if it's funcitonal we need to print the gene files and then summarize it at the different functional levels
# if mode is functional we have to summarize at the different levels other than gene
# before printing we have to reload all the entries for the coord files

if ( $mode eq 'functional' ) {
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =========================== PART 4 : print gene files and also summarize functional levels ==========================\n";
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";

	# we need to print the gene files before we do the summary step, because the summary step removes them from the main HASH, this has to be done in order to gte the stats for sum_not_annotated
	MOCATProfiling::Print::printFiles( \%HASH, 'gene', $base_output_folder, $insert_output_folder );    # this prints the raw, norm and scaled files for the abundances for specified level level (gene level only for functional)

	# dist to HASH_LEVELS
	@levels = ( 'cog', 'ko', 'module', 'pathway' );                                                     # before we had levels set to only genes for functional because it's on the gene level we distribute multiple mappers
	MOCATProfiling::Functional::summarizeFunctionalLevels();
}
else {
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : ============================================= PART 4 : print gene files  ============================================\n";
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
	MOCATProfiling::Print::printFiles( \%HASH, 'gene', $base_output_folder, $insert_output_folder );    # this prints the raw, norm and scaled files for the abundances for specified level level (gene level only for functional)
	@levels = @levels[ 1 .. scalar @levels - 1 ];                                                       # lets remove gene from the levels

}
####################################################################################################################################

####################################################################################################################################
# PART 5 print and exit - yay! done!
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =========================== PART 5 : print all the files and finally have that coffee break =========================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";

foreach my $level (@levels) {
	MOCATProfiling::Print::printFiles( $HASH_LEVELS{$level}, $level, $profiles_output_folder, $profiles_output_folder );    # the gene abundances are stored in HASH, whereas eg species and cog abundances and other levels are stored in HASH_LEVELS
}
####################################################################################################################################

####################################################################################################################################
# DONE
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : ========================= THE END. Seriously. It's over now. Enjoy the sun! Have an ice cream! ======================\n";
print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : =====================================================================================================================\n";
exit 0;
####################################################################################################################################

