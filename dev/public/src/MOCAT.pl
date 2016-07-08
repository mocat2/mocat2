#!/usr/bin/env perl

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

use strict;
use warnings;
use MOCATbwa;
use MOCATCore;
use Pod::Usage;
use MOCATModule;
use MOCATScreen;
use MOCATBowtie2;
use Getopt::Long qw(:config pass_through);
use MOCATAssembly;
use MOCATFetchMGs;
use MOCATFilterV2;
use MOCATVariables;
use MOCATProfiling;
use Term::ANSIColor;
use MOCATStatsFastqc;
use MOCATSampleStatus;
use MOCATProfilingPaste;
use MOCATReadTrimFilter;
use MOCATSampleStatus_dev;
use MOCATAssemblyRevision;
use MOCATSampleStatusExcel2;
use MOCATGenePredictionProdigal;
use MOCATGenePredictionMetaGeneMark;
use MOCATUnpublished;    #DEV VER#

$MOCAT_DESCRIPTION_FILE = "~/MOCAT.versions";
$MOCAT_DESCRIPTION      = "Version 2.1";
$INTERNAL_MOCAT_VERSION = "v2.1.1";
$version                = $INTERNAL_MOCAT_VERSION;
my $space = " ";
$INTERNAL_NCBI_VERSION       = "4";        #20150928: updated ncbi map to v3, #20151002: updated ncbi map to v4 because v3 had bugs and this version should have correct names.dmp names (except 5)
$INTERNAL_TAXONOMIC_VERSION  = "t19";      #t14->t15: fixed important bug that correctly calculates -1 for norm files in NCBI mode, t16: fixed bug in profiling affected mm.dist.among.unique norm for NCBI; t17: 20150928 updated ncbi map to v3; t17: 20151002 updated ncbi map to v4
$INTERNAL_FUNCTIONAL_VERSION = "f16";      #14->15: added horizontal coverage; 15->16: Fixed a bug that calculated the wrong unannotated abundances for functional profiles and mm.dist.among.unique
$INTERNAL_MOTU_VERSION       = "m14";
$INTERNAL_GENE_VERSION       = "g15";      #14->15: added horizontal coverage
$INTERNAL_RESISTANCE_VERSION = "r14";
$INTERNAL_VIRULENCE_VERSION  = "v14";
$INTERNAL_GENERAL_VERSION    = "gen20";    # gen15: uses FilterV2 instead of V1. gen16: has fixed the bug in v2, gen17: fixed minor bug in FilterV2 where len file wasn't sorted, gen18: added horizontal coverage, gen91: fixed bug in profiling for mm.dist.among.unique and NCBI
$INTERNAL_A_VERSION          = "A4";       # this is the assembly version A3: supported up to SOAPdenovo 1.06, A4: SOAPdenovo 2.04
$INTERNAL_B_VERSION          = "B4";       # from version 1.4.2 and onwards this version will incerase if there are changes to the /PROFILES and /SUMMARIES folders
$INTERNAL_C_VERSION          = "C4";       # this is filtering V2. In v3 there was still a minor bug, which was fixed in v4.
$INTERNAL_X_VERSION          = "X3";
$INTERNAL_Y_VERSION          = "Y3";
$INTERNAL_Z_VERSION          = "Z3";

# v1.5.3: added support for using pair and single files in raw data input folders
# v1.5.4: fixed major bugs in new RTF using 3 files
# v1.5.6: t18->t19; removed unnecessary code in the NCBI profiling steps
# v1.5.7: Support for pair.1 pair.2 single in screen fasta file
# v2.1.0: Added support for SLURM queue and BWA for mapping
# v2.1.1: Updated version of fastq_trim_filter_v5_EMBL, fixing FastX bug

##DEV VER#
$MOCAT_DESCRIPTION      = "Major updates";                            #DEV VER#

### DEFINE VARIABLES AND LOG ###
our @args = @ARGV;
our @all_commands;
$print_final_text = 1;
$sep              = "===============================================================================";
$sample_file      = "";
$no_temp          = "";
my $USAGE;
our $cpu = "0";
our @screen;
our $date = `date +\%Y\%m\%d_\%H\%M\%S`;
chomp($date);

$|++;
### AVAILABLE ROUTINES ###
my $do_read_trim_filter;
my @do_screen;
my @do_screen_fasta_file;
my $do_assembly          = "NOTSET";
my $do_assembly_revision = "NOTSET";
my $do_gene_prediction;
my $do_stats_fastqc;
my $do_sample_status;
my $BOWTIE2;
my $BWA;
my @do_filter;
my $help;
my $manual;
my $do_version;

### ADDITIONAL VARIABLES ###
my $SETWD;
our @databases;
our $mod_single_job    = "";
our $screen_fasta_file = "";
our $reads             = "";
our $assembly          = "";
our $use_mem           = "";
our $realtime_status   = "";
my $LENGTH_CUTOFF  = "";
my $PERCENT_CUTOFF = "";
our $run_on_host = "";
our $pcf_uniq;
our $taxo_profiling_mode = "";
my $memory;
my $TC;
my $no_queue;
my @CONFIG_OVERRIDE;
our $PUBLIC = 1;
my $cite;
my $tiger;
my $NO_CALCULATE_HORIZONTAL_COVERAGE;
$CALCULATE_HORIZONTAL_COVERAGE = 1;
$RUN_MODULE                    = "0";

$PUBLIC = 0;                                #DEV VER#
our $do_cluster_genes;                      #DEV VER#
our $do_import_old;                         #DEV VER#
our $do_redo_scaf;                          #DEV VER#
our $do_extract_genes;                      #DEV VER#
our $do_import;                             #DEV VER#
our @do_calculate_multiplemappers;          #DEV VER#
our $do_calculate_taxonomic_composition;    #DEV VER#
our $do_run_carma_annotation;               #DEV VER#
our $analyze;                               #DEV VER#
our $resistance_screen_reads2         = "";           #DEV VER#
our $importOldScafLength              = "NO";         #DEV VER#
our $coverage_version_1               = "";           #DEV VER#
our $carma_file                       = "";           #DEV VER#
our $analyze_project_name             = "";           #DEV VER#
our $run_carma_annotate_host1         = "epsilon";    #DEV VER#
our $run_carma_annotate_host2         = "eta";        #DEV VER#
our $multiple_mappers_summarize_level = "";           #DEV VER#
my $filter_psort_buffer = "";                         #DEV VER#
our $identity2                        = "";                            #DEV VER#
our $resistance_screen_refgenecatalog = '263RefGeneCatalog.padded';    #DEV VER#
our $resistance_screen_suffix         = "";                            #DEV VER#
our $OUTPUT_FOLDER                    = "OUTPUT";
our $MOCAT_ID                         = "identifier unset";

chomp( $username = `whoami` );
chomp( $hostname = `hostname` );

### LOAD OPTIONS ###
GetOptions(
  'tiger'                     => \$tiger,
  'cite'                      => \$cite,
  'h|help|?'                  => \$help,
  'man|manual'                => \$manual,
  'rtf|read_trim_filter'      => \$do_read_trim_filter,
  'sfq|stats_fastqc'          => \$do_stats_fastqc,
  'ss|sample_status'          => \$do_sample_status,
  'x|no_execute'              => \$only_print,
  'cfg:s'                     => \$config_file,
  'sf|sample_file:s'          => \$sample_file,
  'cpus:i'                    => \$cpu,
  'nt|no_temp'                => \$no_temp,
  'soap:s{,}'             => \@do_screen,
  'sff|screen_fastafile:s{,}' => \@do_screen_fasta_file,
  'r|reads:s{,}'              => \@reads,
  'a|assembly'                => \$do_assembly,
  'ar|assembly_revision'      => \$do_assembly_revision,
  'gp|gene_prediction:s'      => \$do_gene_prediction,
  'fmg|fetch_mg:s'            => \$do_fetch_mg,
  'fsoap:s{,}'             => \@do_filter,
  'e|extracted'               => \$use_extracted_reads,
  'old'                       => \$use_old_version,
  'use_mem'                   => \$use_mem,
  'length:s'                  => \$LENGTH_CUTOFF,
  'id|identity:s'             => \$PERCENT_CUTOFF,
  'host:s'                    => \$run_on_host,
  'SETWD:s'                   => \$SETWD,
  'version'                   => \$do_version,
  'unique'                    => \$pcf_uniq,
  'm|mode:s'                  => \$profiling_mode,
  'screened_files'            => \$screened_files,
  'extracted_files'           => \$extracted_files,
  'temp'                      => \$TEMP,
  'shm'                       => \$SHM,
  'bowtie2:s{,}'              => \@BOWTIE2,
  'bwa:s{,}'              => \@BWA,

  'v1_pcf|v1_paste_coverage_files=s{,}' => \@do_paste_coverage_files,                #DEV VER#
  'v1_t|v1_taxonomic_profiling=s{,}'    => \@do_taxo_profiling,                      #DEV VER#
  'v1_c|v1_coverage=s{,}'               => \@do_calculate_coverage,                  #DEV VER#
  'v1_functional_profiling=s{,}'        => \@do_functional_profiling,                #DEV VER#
  'cluster_genes_by_sample:s'           => \$do_cluster_genes,                       #DEV VER#
  'create_carma_annotation_file=s'      => \$do_run_carma_annotation,                #DEV VER#
  'host1:s'                             => \$run_carma_annotate_host1,               #DEV VER#
  'host2:s'                             => \$run_carma_annotate_host2,               #DEV VER#
  'calculate_multiple_mappers=s{,}'     => \@do_calculate_multiplemappers,           #DEV VER#
  'import=s'                            => \$import_path,                            #DEV VER#
  'import_sl=s'                         => \$import_max_sample_name,                 #DEV VER#
  'scaf=s'                              => \$importOldScafLength,                    #DEV VER#
  'scafa=s'                             => \$assembly,                               #DEV VER#
  'eg_genes_file=s'                     => \$extract_genes,                          #DEV VER#
  'annotate_using_carma_file=s{,}'      => \@do_annotate_using_carma_file,           #DEV VER#
  'file=s'                              => \$carma_file,                             #DEV VER#
  'summarize_carma_coverages:s{,}'      => \@do_summarize_carma_coverages,           #DEV VER#
  'map:s'                               => \$gbc_map,                                #DEV VER#
  'analyze=s{,}'                        => \@do_analyze,                             #DEV VER#
  'analyze_function:s'                  => \$analyze_function,                       #DEV VER#
  'analyze_downsample:s'                => \$analyze_downsample,                     #DEV VER#
  'analyze_type:s'                      => \$analyze_type,                           #DEV VER#
  'analyze_norm:s'                      => \$analyze_norm,                           #DEV VER#
  'analyze_taxa:s'                      => \$analyze_taxa,                           #DEV VER#
  'analyze_desc:s'                      => \$analyze_desc,                           #DEV VER#
  'analyze_file:s'                      => \$analyze_file,                           #DEV VER#
  'analyze_output:s{,}'                 => \@analyze_output,                         #DEV VER#
  'analyze_metadata:s{,}'               => \@analyze_metadata,                       #DEV VER#
  'analyze_group:s'                     => \$analyze_group,                          #DEV VER#
  'analyze_condA:s'                     => \$analyze_condA,                          #DEV VER#
  'analyze_condB:s'                     => \$analyze_condB,                          #DEV VER#
  'analyze_settings:s{,}'               => \@analyze_settings,                       #DEV VER#
  'analyze_project_name:s'              => \$analyze_project_name,                   #DEV VER#
  'analyze_feature_annotations:s{,}'    => \@analyze_feature_annotations,            #DEV VER#
  'analyze_running_galaxy'              => \$analyze_running_galaxy,                 #DEV VER#
  'eggnog_map:s'                        => \$functional_profiling_eggnog_map,        #DEV VER#
  'kegg_map:s'                          => \$functional_profiling_kegg_map,          #DEV VER#
  'resistance_screen:s{,}'              => \@do_resistance_screen,                   #DEV VER#
  'r2:s'                                => \$resistance_screen_reads2,               #DEV VER#
  'e2'                                  => \$use_extracted_reads2,                   #DEV VER#
  'analyze_downsample_size:s'           => \$analyze_downsample_size,                #DEV VER#
  'gbc=s{,}'                            => \@do_group_by_column,                     #DEV VER#
  'extract_genes_genefile:s'            => \$do_extract_genes,                       #DEV VER#
  'extract_genes_scaftigfile:s'         => \$extract_genes_scaftig,                  #DEV VER#
  'status_for_assembly'                 => \$SampleStatusExcel2_useAssembly,         #DEV VER#
  'only_regenerate_reads'               => \$only_regenerate_reads,                  #DEV VER#
  'summarize_level:s'                   => \$multiple_mappers_summarize_level,       #DEV VER#
  'manual_stats_file:s'                 => \$calculateTaxonomy_manual_stats_file,    #DEV VER#
  'virulence_screen:s{,}'               => \@do_virulence_screen,                    #DEV VER#
  'identity2:s'                         => \$identity2,                              #DEV VER#
  'filter_psort_buffer:s'               => \$filter_psort_buffer,                    #DEV VER#
  'refgenecatalog:s'                    => \$resistance_screen_refgenecatalog,       #DEV VER#
  'SOLiD'                               => \$preprocessSOLiD,                        #DEV VER#
  'sam'                                 => \$profiling_SAM,                          #DEV VER#
  'v2_t|v2_calculate_taxonomy:s{,}'     => \@do_calculate_taxonomy,                  #DEV VER#
  'database_name:s'                     => \$database_name,                          #DEV VER#

  'previous_db_calc_tax_stats_file' => \$calculateTaxonomy_previous_calc_coverage_stats_file,
  'no_paste'                        => \$no_paste,
  'only_paste'                      => \$only_paste,
  'config:s{,}'                     => \@CONFIG_OVERRIDE,
  'memory:s'                        => \$memory,
  'tc:s'                            => \$TC,
  'noqueue'                         => \$no_queue,
  'psoap:s{,}'                => \@do_profiling,
  'o|output:s'                      => \$OUTPUT_FOLDER,
  'verbose'                         => \$VERBOSE,
  'no_horizontal'                   => \$NO_CALCULATE_HORIZONTAL_COVERAGE,
  'pbwa:s{,}'                => \@do_profiling_bwa,
  
);

# module support
if ( scalar @ARGV > 0 )
{
  $MODULE_OPTION = $ARGV[0];
  $MODULE_OPTION =~ s/^-*//;
  @MODULE_PARAM = @ARGV[ 1 .. scalar @ARGV - 1 ];
}

if ($NO_CALCULATE_HORIZONTAL_COVERAGE)
{
  $CALCULATE_HORIZONTAL_COVERAGE = 0;
}

if ($tiger)
{
  tiger();
  exit 0;
}

# PRINT VERSION AND EXIT
if ($do_version)
{
  print "MOCAT $version\n";
  exit(0);
}

# SET WD AND CFG #
if ($SETWD)
{
  $cwd = $SETWD;
} else
{
  chomp( $cwd = `pwd` );
}
unless ($config_file)
{
  $config_file = "$cwd/MOCAT.cfg.$username";
  unless ( -e $config_file )
  {
    $config_file = "$cwd/MOCAT.cfg";
  }
}

$jobdir = "$cwd";
system "mkdir -p $jobdir";

### CHECK AND ORGANIZE CPUS ###
my $cpu_read_trim_filter     = 1;
my $cpu_screen               = 8;
my $cpu_bowtie2              = 8;
my $cpu_bwa                  = 8;
my $cpu_screen_fastafile     = 1;
my $cpu_assembly             = 8;
my $cpu_assembly_revision    = 8;
my $cpu_gene_prediction      = 1;
my $cpu_filter               = 8;
my $cpu_paste_coverage_files = 1;
my $cpu_fetch_mg             = 4;
my $cpu_profiling            = 4;
$cpu_module = 8;
$cpus       = 1;

if ( $cpu > 0 )
{
  $cpu_read_trim_filter     = $cpu;
  $cpu_screen               = $cpu;
  $cpu_bowtie2              = $cpu;
  $cpu_bwa                  = $cpu;
  $cpu_screen_fastafile     = $cpu;
  $cpu_assembly             = $cpu;
  $cpu_assembly_revision    = $cpu;
  $cpu_gene_prediction      = $cpu;
  $cpu_filter               = $cpu;
  $cpu_paste_coverage_files = $cpu;
  $cpu_fetch_mg             = $cpu;
  $cpu_profiling            = $cpu;
  $cpus                     = $cpu;
  $cpu_module               = $cpu;
}

### CHECK IF SOME OPTION HAS BEEN STARTED ###
if ( $do_assembly eq "NOTSET" )
{
  undef $do_assembly;
} elsif ( $do_assembly eq "" )
{
  $do_assembly = "1";
}
if ( $do_assembly_revision eq "NOTSET" )
{
  undef $do_assembly_revision;
} elsif ( $do_assembly_revision eq "" )
{
  $do_assembly_revision = "1";
}

if ($do_profiling_bwa[0]) {
 @do_profiling = @do_profiling_bwa;
}

### CHECK IF WE WANT TO PRINT HELP ###
### CHECK IF ANY ROUTINES ARE BEING RUN OR HELP ###
my $routine_sum = defined($MODULE_OPTION) + defined($do_read_trim_filter) + defined( $BWA[0] ) + defined( $BOWTIE2[0] ) + defined( $do_screen[0] ) + defined( $do_screen_fasta_file[0] ) + defined($do_assembly) + defined($do_assembly_revision) + defined($do_gene_prediction) + defined($do_sample_status) + defined($do_stats_fastqc) + defined( $do_filter[0] ) + defined( $do_paste_coverage_files[0] ) + defined($do_fetch_mg) + defined( $do_profiling[0] );

$routine_sum +=                                    #DEV VER#
  defined($do_import) +                            #DEV VER#
  defined($do_redo_scaf) +                         #DEV VER#
  defined($do_cluster_genes) +                     #DEV VER#
  defined($do_extract_genes) +                     #DEV VER#
  defined( $do_calculate_multiplemappers[0] ) +    #DEV VER#
  defined($do_run_carma_annotation) +              #DEV VER#
  defined( $do_annotate_using_carma_file[0] ) +    #DEV VER#
  defined( $do_summarize_carma_coverages[0] ) +    #DEV VER#
  defined( $do_analyze[0] ) +                      #DEV VER#
  defined($analyze_file) +                         #DEV VER#
  defined( $do_functional_profiling[0] ) +         #DEV VER#
  defined( $do_resistance_screen[0] ) +            #DEV VER#
  defined( $do_virulence_screen[0] ) +             #DEV VER#
  defined( $do_group_by_column[0] ) +              #DEV VER#
  defined( $do_calculate_coverage[0] ) +           #DEV VER#
  defined( $do_calculate_taxonomy[0] ) +           #DEV VER#
  defined( $do_taxo_profiling[0] );                #DEV VER#

if ($cite)
{
  Cite();
}
if ($help)
{
  Usage();
} elsif ($manual)
{
  pod2usage( -exitstatus => 0, -verbose => 2 );
} elsif ( $routine_sum == 0 || $routine_sum > 1 )
{
  Usage();
}

### PRINT WELCOME ###
sub welcome
{
  print "$sep\n";
  print color 'bold';
  print "                  MOCAT";
  print color 'reset';
  print " - Metagenomics Analysis Toolkit                $space";
  print color 'bold';
  print "$version\n";
  print color 'reset';
  print " by Jens Roat Kultima, Luis Pedro Coelho, Shinichi Sunagawa @ Bork Group, EMBL
$sep\n";
}
welcome();
### LOAD SAMPLE AND CONFIG FILE ###
unless ( $do_run_carma_annotation || $analyze_file || $do_extract_genes )
{    #DEV VER#
  @samples = MOCATCore::read_samples($sample_file);
  ### CHECK NO OF SAMPLES ###
  if ( scalar @samples == 0 )
  {
    die "ERROR & EXIT: Number of samples is 0. Correct sample file?";
  }
  chomp( my $bnsf = `basename $sample_file` );
  $date = "$bnsf\_$date";
}    #DEV VER#

%conf         = MOCATCore::read_config( $config_file, "dont_print_config" );
$mod_dir      = $conf{MOCAT_dir} . "/mod";
$bin_dir      = $conf{MOCAT_dir} . "/bin";
$scr_dir      = $conf{MOCAT_dir} . "/src";
$src_dir      = $conf{MOCAT_dir} . "/src";
$ext_dir      = $conf{MOCAT_dir} . "/ext";
$data_dir     = $conf{MOCAT_data_dir};
$qsub         = $conf{MOCAT_qsub_system};
$ziplevel     = $conf{MOCAT_zip_level};
$rownames_dir = $conf{MOCAT_rownames_dir};

print "EXECUTE COMMAND   : " . join( " ", @args ) . "\n";
print "CURRENT USER      : $username\n";
print "CURRENT HOST      : $hostname\n";
print "WORKING DIRECTORY : $cwd\n";
print "CONFIG FILE       : $config_file\n";
print "NUMBER OF SAMPLES : " . scalar @samples . "\n";
print "SAMPLE NAMES      : " . join( ", ", @samples ) . "\n";
print "SAMPLE FILE       : $sample_file\n";
print "$sep\n";

# Manual override of settings
if ( $LENGTH_CUTOFF ne "" )
{
  $conf{readtrimfilter_length_cutoff} = $LENGTH_CUTOFF;
  $conf{screen_length_cutoff}         = $LENGTH_CUTOFF;
  $conf{filter_length_cutoff}         = $LENGTH_CUTOFF;
}
if ( $PERCENT_CUTOFF ne "" )
{
  $conf{screen_percent_cutoff} = $PERCENT_CUTOFF;
  $conf{filter_percent_cutoff} = $PERCENT_CUTOFF;
  $conf{bwa_percent_cutoff} = $PERCENT_CUTOFF;
  $conf{bowtie2_percent_cutoff} = $PERCENT_CUTOFF;
}

#DEV OVERRIDES# #DEV VER#
foreach my $override (@CONFIG_OVERRIDE)
{
  my @overrideX = $override =~ /^([a-zA-Z_]*)=(.*)/; 
  if ( scalar @overrideX == 2 )
  {
    $conf{ $overrideX[0] } = $overrideX[1];
  }
}
if ( $hostname =~ /submaster/ )
{    #DEV VER#
  $conf{MOCAT_qsub_system} = 'LSF';    #DEV VER#
  $qsub = $conf{MOCAT_qsub_system};    #DEV VER#
}    #DEV VER#
if ($no_queue)
{
  $conf{MOCAT_qsub_system} = 'none';
  $qsub = $conf{MOCAT_qsub_system};
}
if ($memory)
{
  $conf{MOCAT_SGE_qsub_add_param} = $conf{MOCAT_SGE_qsub_add_param} . " -l h_vmem=$memory ";
  $conf{MOCAT_LSF_memory_limit}   = $memory;
  $conf{MOCAT_LSF_memory_limit} =~ s/G/000/;
  my $mem2 = $memory;
  $mem2 =~ s/G/000/;
  $conf{MOCAT_SLURM_qsub_add_param} = $conf{MOCAT_SLURM_qsub_add_param} . " --mem=$mem2 ";
  my $tmem = $memory;
  $tmem =~ s/G//;
  $conf{filter_psort_buffer} = ( $tmem - 5 ) . "G";
}
if ( $filter_psort_buffer ne "" )
{    #DEV VER#
  $conf{filter_psort_buffer} = $filter_psort_buffer;    #DEV VER#
}    #DEV VER#
if ($TC)
{
  $conf{MOCAT_SGE_qsub_add_param} = $conf{MOCAT_SGE_qsub_add_param} . " -tc $TC ";
}

#DEV OVERRIDES# #DEV VER#

MOCATCore::print_settings("MOCAT");
### Get reads, if not specified on command line, but in config ###
if ( $reads[0] )
{
  $reads = $reads[0];
} elsif ( $do_assembly || $do_assembly_revision )
{
  die "ERROR & EXIT: For the assembly and assembly revision steps, please specify -r [READS] explicitly (to minimize risk of using wrong reads for assembly)";
} elsif ( $conf{MOCAT_default_reads} ne '-define if desired-' && $reads eq '' )
{
  $reads = $conf{MOCAT_default_reads};
  $reads[0] = $reads;

}
### CHECK BIN ###
my @bin = qw(
  2bwt-builder
  bwa
  cd-hit-est
  cdbfasta
  cdbyank
  diamond
  faSomeRecords
  fastq_quality_filter
  fastq_quality_trimmer
  fastq_trim_filter_v5_EMBL
  fastqc
  fastx_quality_stats
  fastx_trimmer
  hmmsearch
  msamtools
  nt2aa
  pigz
  prodigal
  ps
  psort
  psort_OSX
  pstree
  samtools
  soap.coverage
  soap2.21
  sqlite3_linux
  sqlite3_osx
  unpigz
  unzip
  zip
  SOAPdenovo2/SOAPdenovo-127mer
  SOAPdenovo2/SOAPdenovo-127mer_OSX
  SOAPdenovo2/SOAPdenovo-63mer
  SOAPdenovo2/SOAPdenovo-63mer_OSX
  soap1.05/SOAPdenovo
  soap1.05/SOAPdenovo127mer
  soap1.05/SOAPdenovo31mer
  soap1.05/SOAPdenovo63mer
  soap1.06/SOAPdenovo-127mer-static
  soap1.06/SOAPdenovo-31mer-nobam-noAIO
  soap1.06/SOAPdenovo-31mer-static
  soap1.06/SOAPdenovo-63mer-nobam-noAIO
  soap1.06/SOAPdenovo-63mer-static
  soap1.06OSX/SOAPdenovo-127mer
  soap1.06OSX/SOAPdenovo-31mer
  soap1.06OSX/SOAPdenovo-63mer
  soap1.06OSXnobamaio/SOAPdenovo-127mer
  soap1.06OSXnobamaio/SOAPdenovo-31mer
  soap1.06OSXnobamaio/SOAPdenovo-63mer
);

foreach my $bin (@bin)
{
  unless ( -s "$bin_dir/$bin" )
  {
    die "ERROR & EXIT: Missing executable $bin_dir/$bin.\nThis file should be shipped with the MOCAT zip.";
  }
}
### CHECK CONFIG ###
my @configs = qw(
  MOCAT_pre_execute
  MOCAT_umask
  MOCAT_qsub_system
  MOCAT_SGE_qsub_add_param
  MOCAT_SGE_qsub_add_param
  MOCAT_SGE_parallell_env
  MOCAT_dir
  MOCAT_data_type
  MOCAT_paired_end
  MOCAT_zip_program
  MOCAT_zip_level
  MOCAT_default_reads
  MOCAT_mapping_mode
  MOCAT_prompt_before_run
  readtrimfilter_length_cutoff
  readtrimfilter_qual_cutoff
  readtrimfilter_use_sanger_scale
  readtrimfilter_trim_5_prime
  readtrimfilter_use_precalc_5prime_trimming
  screen_length_cutoff
  screen_percent_cutoff
  screen_soap_seed_length
  screen_soap_max_mm
  screen_soap_cmd
  screen_save_format
  screen_fasta_file_blast_e_value
  screen_fasta_file_blast_read_min_length
  screen_fasta_file_additional_usearch_cmd
  assembly_soap_version
  assembly_calculate_insert_size_using
  assembly_db_for_calc_insertsize
  assembly_scaftig_min_length
  assembly_revision_scaftig_min_length
  gene_prediction_software
  gene_prediction_input
  filter_psort_buffer
  filter_length_cutoff
  filter_percent_cutoff
  filter_paired_end_filtering
  filter_remove_mapping_files
  filter_samtools_memory
  filter_make_unique_sorted
  realtime_status_use
  realtime_status_timer
  realtime_status_log
  screen_fasta_file_usearch_version_5_exe
  screen_fasta_file_usearch_version_6_exe
  screen_fasta_file_usearch_version
  profiling_paired_end_filtering
);

foreach my $c (@configs)
{
  unless ( exists $conf{$c} )
  {
    die "ERROR & EXIT: Missing setting '$c' in config file. Is the config file correct?";
  }
}
### OTHER CHECKS ###
our $ZCAT = "gunzip -c";
our $ZIP  = "gzip";
chomp( our $systemType = `uname -s` );
if ( $systemType =~ m/Darwin/ )
{
  $ZCAT = "gunzip -c";
}
if ( $conf{MOCAT_zip_program} eq 'pigz' )
{
  $ZIP = "$bin_dir/pigz variableswillbesetelsewhere";
}
### CHECK SGE ###
#chomp( my $SGE_ROOT = `echo \$SGE_ROOT` );
#if ( $SGE_ROOT eq "" & $conf{MOCAT_qsub_system} eq "SGE" ) {
#	die "ERROR & EXIT: MOCAT is set to run under the SGE queuing system, but UNIX environmental variable SGE_ROOT isn't set.
#              This usually indicates that there is no SGE queuing system available.
#              Please change 'qsub_system' to 'none' in the config file.\n";
#}
### ALL CHECKS PASSED, print run command :) ###
chomp( my $cd = `date +'\%F \%T'` );
chomp( my $h  = `hostname` );
chomp( my $un = `whoami` );
system "echo 'MOCAT.pl " . join( " ", @args ) . " @ $cd ($un on $h)' >> MOCAT.executed";
if ( $sample_file ne '' )
{
  chomp( $sample_file_basename = `basename $sample_file` );
} else
{
  $sample_file_basename = 'samples';
}

# If executing on a different host
if ( $run_on_host ne "" )
{
  print "$sep\nEXECUTING COMMANDS ON A DIFFERENT HOST: $run_on_host\n";
}

# Set umask
umask $conf{MOCAT_umask};

# Record script versions
if ( $conf{MOCAT_record_script_versions} eq 'yes' )
{
  MOCATCore::record_script_version();
}

#LOAD AND RUN DEV              #DEV VER#
MOCATUnpublished::Launch();    #DEV VER#

### READ TRIM AND FILTER ###
if ( defined $do_read_trim_filter )
{
  MOCATCore::print_settings("readtrimfilter");
  print localtime() . ": PERFORMING READ TRIM AND FILTER...\n";
  MOCATReadTrimFilter::pre_check_files();
  MOCATReadTrimFilter::create_job( "readtrimfilter", $cpu_read_trim_filter );
  MOCATCore::execute_job( "readtrimfilter", $cpu_read_trim_filter, scalar @samples, "15gb" );
  unless ($only_print)
  {
    MOCATReadTrimFilter::post_check_files();
    print localtime() . ": COMPLETED READ TRIM AND FILTER.\n";
    print "Files stored in $cwd/<SAMPLE>/reads.processed.$conf{MOCAT_data_type}/*[pair|single]*.fq.gz\n";
  }
}

### SCREEN ###
if ( defined $do_screen[0] )
{
  if ( $do_screen[0] eq "" )
  {
    die "ERROR & EXIT: No database defined with the -s option. Please specify -s DATABASE (either as single file name, or full path)";
  }
  @screen = @do_screen;
  MOCATCore::print_settings("screen");
  print localtime() . ": PERFORMING SCREEN USING DATABASE...\n";
  MOCATScreen::pre_check_files();
  if ($use_mem)
  {
    unless ($only_print)
    {
      foreach my $screen (@screen)
      {
        MOCATCore::fake_data_init("$screen.index");
      }
    }
  }
  MOCATScreen::create_job( "screen", $cpu_screen );
  MOCATCore::execute_job( "screen", $cpu_screen, scalar @samples, "15gb" );
  if ($use_mem)
  {
    unless ($only_print)
    {
      foreach my $screen (@screen)
      {
        MOCATCore::fake_data_done("$screen.index");
      }
    }
  }
  unless ($only_print)
  {
    MOCATScreen::post_check_files();
    print localtime() . ": COMPLETED SCREEN USING DATABASE.\n";
    if (    $screen[0] eq 's'
         || $screen[0] eq 'c'
         || $screen[0] eq 'f'
         || $screen[0] eq 'r' )
    {
      print "Files stored in $cwd/<SAMPLE/reads.mapped.<assembly|revised.assembly>.K<kmer>.$reads.$conf{MOCAT_data_type}/\n";
      if ($screened_files)
      {
        print "Files stored in $cwd/<SAMPLE/reads.screened.<assembly|revised.assembly>.K<kmer>.$reads.$conf{MOCAT_data_type}/\n";
      }
      if ($extracted_files)
      {
        print "Files stored in $cwd/<SAMPLE/reads.extracted.<assembly|revised.assembly>.K<kmer>.$reads.$conf{MOCAT_data_type}/\n";
      }
    } else
    {
      foreach my $screen (@screen)
      {
        if ($screened_files)
        {
          print "Files stored in $cwd/<SAMPLE>/reads.screened.$screen.$conf{MOCAT_data_type}/*[pair|single]*.fq.gz\n";
        }
        if ($extracted_files)
        {
          print "Files stored in $cwd/<SAMPLE>/reads.extracted.$screen.$conf{MOCAT_data_type}/*[pair|single]*.fq.gz\n";
        }
        print "Files stored in $cwd/<SAMPLE>/reads.mapped.$screen.$conf{MOCAT_data_type}/\n";
      }
    }
  }
}

### BOWTIE2 ###
if ( defined $BOWTIE2[0] )
{
 die "The support for Bowtie2 hasn not been fully implemented and probably will not ever be. Sorry.";
  if ( $BOWTIE2[0] eq "" )
  {
    die "ERROR & EXIT: No database defined with the -bowtie2 option. Please specify -bowtie2 DATABASE (either as single file name, or full path)";
  }
  @screen = @BOWTIE2;
  MOCATCore::print_settings("bowtie2");
  print localtime() . ": PERFORMING BOWTIE2 MAPPING+FILTERING...\n";
  MOCATBowtie2::pre_check_files();
  MOCATBowtie2::create_job( "bowtie2", $cpu_screen );
  MOCATCore::execute_job( "bowtie2", $cpu_screen, scalar @samples, "15gb" );
  unless ($only_print)
  {
    MOCATBowtie2::post_check_files();
    print localtime() . ": COMPLETED SCREEN+FILTERING USING BOWTIE2.\n";
    if (    $screen[0] eq 's'
         || $screen[0] eq 'c'
         || $screen[0] eq 'f'
         || $screen[0] eq 'r' )
    {
      print "Files stored in $cwd/<SAMPLE/reads.filtered.<assembly|revised.assembly>.K<kmer>.$reads.$conf{MOCAT_data_type}/\n";
      if ($screened_files)
      {
        print "Files stored in $cwd/<SAMPLE/reads.screened.<assembly|revised.assembly>.K<kmer>.$reads.$conf{MOCAT_data_type}/\n";
      }
      if ($extracted_files)
      {
        print "Files stored in $cwd/<SAMPLE/reads.extracted.<assembly|revised.assembly>.K<kmer>.$reads.$conf{MOCAT_data_type}/\n";
      }
    } else
    {
      foreach my $screen (@screen)
      {
        if ($screened_files)
        {
          print "Files stored in $cwd/<SAMPLE>/reads.screened.$screen.$conf{MOCAT_data_type}/*[pair|single]*.fq.gz\n";
        }
        if ($extracted_files)
        {
          print "Files stored in $cwd/<SAMPLE>/reads.extracted.$screen.$conf{MOCAT_data_type}/*[pair|single]*.fq.gz\n";
        }
        print "Files stored in $cwd/<SAMPLE>/reads.filtered.$screen.$conf{MOCAT_data_type}/\n";
      }
    }
  }
}


### BWA ###
if ( defined $BWA[0] )
{
  if ( $BWA[0] eq "" )
  {
    die "ERROR & EXIT: No database defined with the -bwa option. Please specify -bwa DATABASE (either as single file name, or full path)";
  }
  @screen = @BWA;
  MOCATCore::print_settings("bwa");
  print localtime() . ": PERFORMING BWA MAPPING+FILTERING...\n";
  MOCATbwa::pre_check_files();
  MOCATbwa::create_job( "bwa", $cpu_screen );
  MOCATCore::execute_job( "bwa", $cpu_screen, scalar @samples, "15gb" );
  unless ($only_print)
  {
    MOCATbwa::post_check_files();
    print localtime() . ": COMPLETED SCREEN+FILTERING USING BWA.\n";
    if (    $screen[0] eq 's'
         || $screen[0] eq 'c'
         || $screen[0] eq 'f'
         || $screen[0] eq 'r' )
    {
      print "Files stored in $cwd/<SAMPLE/reads.filtered.<assembly|revised.assembly>.K<kmer>.$reads.$conf{MOCAT_data_type}/\n";
      if ($screened_files)
      {
        print "Files stored in $cwd/<SAMPLE/reads.screened.<assembly|revised.assembly>.K<kmer>.$reads.$conf{MOCAT_data_type}/\n";
      }
      if ($extracted_files)
      {
        print "Files stored in $cwd/<SAMPLE/reads.extracted.<assembly|revised.assembly>.K<kmer>.$reads.$conf{MOCAT_data_type}/\n";
      }
    } else
    {
      foreach my $screen (@screen)
      {
        if ($screened_files)
        {
          print "Files stored in $cwd/<SAMPLE>/reads.screened.$screen.$conf{MOCAT_data_type}/*[pair|single]*.fq.gz\n";
        }
        if ($extracted_files)
        {
          print "Files stored in $cwd/<SAMPLE>/reads.extracted.$screen.$conf{MOCAT_data_type}/*[pair|single]*.fq.gz\n";
        }
        print "Files stored in $cwd/<SAMPLE>/reads.filtered.$screen.$conf{MOCAT_data_type}/\n";
      }
    }
  }
}


### FASTA FILE SCREEN ###
if ( defined $do_screen_fasta_file[0] )
{
  if ( $do_screen_fasta_file[0] eq "" )
  {
    die "ERROR & EXIT: No fasta file defined with the -sff option. Please specify -sff FASTA FILE (either as single file name, or full path)";
  }
  if ( scalar @do_screen_fasta_file > 1 )
  {
    die "ERROR & EXIT: Only screening against 1 fasta file at a time is supported. Please concatenate the fasta files.";
  }
  ### Print first settings and check fasta file ###
  MOCATCore::print_settings("screen_fasta_file");
  $screen[0] = $do_screen_fasta_file[0];
  $SCREEN_FASTA_FILE = 1;
  MOCATScreen::pre_check_files();
  MOCATScreen::create_job( "screen_fasta_file", $cpu_screen_fastafile );
  MOCATCore::execute_job( "screen_fasta_file", $cpu_screen_fastafile, scalar @samples, "15gb" );
  ### Done ##
  unless ($only_print)
  {
    MOCATScreen::post_check_files();
    print localtime() . ": COMPLETED SCREEN USING FASTA FILE.\n";
    print "Files stored in $cwd/<SAMPLE>/reads.mapped.'FASTA FILE'.$conf{MOCAT_data_type}/*[pair|single]*.fq.gz\n";
    if ($screened_files)
    {
      print "Files stored in $cwd/<SAMPLE>/reads.screened.'FASTA FILE'.$conf{MOCAT_data_type}/*[pair|single]*.fq.gz\n";
    }
    if ($extracted_files)
    {
      print "Files stored in $cwd/<SAMPLE>/reads.extracted.'FASTA FILE'.$conf{MOCAT_data_type}/*[pair|single]*.fq.gz\n";
    }
  }
}
### ASSEMBLY ###
if ($do_assembly)
{
  ### Variables ###
  push @databases, $conf{assembly_db_for_calc_insertsize};
  ### Assembly ###
  MOCATCore::print_settings("assembly");
  print localtime() . ": PERFORMING ASSEMBLY...\n";
  MOCATAssembly::pre_check_files();
  MOCATCore::check_databases( \@databases );
  MOCATAssembly::create_job( "assembly", $cpu_assembly );
  MOCATCore::execute_job( "assembly", $cpu_assembly, scalar @samples, "15gb" );
  unless ($only_print)
  {
    MOCATAssembly::post_check_files();
    print localtime() . ": COMPLETED ASSEMBLY.\n";
    my $reads2 = join( "_AND_", @reads );
    print "Files stored in $cwd/<SAMPLE>/assembly.$reads2.$conf{MOCAT_data_type}.KXX/<SAMPLE>.assembly.$reads2.$conf{MOCAT_data_type}.KXX.*\n";
  }
}
### ASSEMBLY CORRECTION ###
if ($do_assembly_revision)
{
  MOCATCore::print_settings("assembly_revision");
  print localtime() . ": PERFORMING ASSEMBLY REVISION...\n";
  if ( $conf{MOCAT_paired_end} ne 'yes' )
  {
    die "ERROR & EXIT: Assmelby revision can only be performed on paired-end data.";
  }
  MOCATAssemblyRevision::pre_check();
  MOCATAssemblyRevision::create_job( "assembly_revision", $cpu_assembly_revision );
  MOCATCore::execute_job( "assembly_revision", $cpu_assembly_revision, scalar @samples, "15gb" );
  unless ($only_print)
  {
    MOCATAssemblyRevision::post_check_files();
    print localtime() . ": COMPLETED ASSEMBLY REVISION.\n";
    print "Files stored in $cwd/<SAMPLE>/assembly.revised.$reads.$conf{MOCAT_data_type}.KXX/<SAMPLE>.assebmly.revised.$reads.$conf{MOCAT_data_type}.KXX.scaftig\n";
  }
}
### GENE PREDICTION ###
if ( defined $do_gene_prediction )
{
  if ( $do_gene_prediction eq "" || !( $do_gene_prediction eq "assembly" || $do_gene_prediction eq "assembly.revised" ) )
  {
    die "ERROR & EXIT: Please specify as -gp assembly or -gp assembly.revised";
  }
  if ($only_print)
  {
    die "EXIT: Gene prediction step does not support option -x. Sorry.\n";
  } else
  {
    $assembly = $do_gene_prediction;
    MOCATCore::print_settings("gene_prediction");
    if ( $conf{gene_prediction_software} eq 'MetaGeneMark' )
    {
      print localtime() . ": PERFORMING GENE PREDICTION USING METAGENEMARK...\n";
      MOCATGenePredictionMetaGeneMark::pre_check();
      MOCATGenePredictionMetaGeneMark::create_job( "gene_prediction", $cpu_gene_prediction );
    } elsif ( $conf{gene_prediction_software} eq 'Prodigal' )
    {
      print localtime() . ": PERFORMING GENE PREDICTION USING PRODIGAL...\n";
      MOCATGenePredictionProdigal::pre_check();
      MOCATGenePredictionProdigal::create_job( "gene_prediction", $cpu_gene_prediction );
    } else
    {
      die "ERROR & EXIT: 'gene_prediction_software' in the MOCAT.cfg file must be either 'MetaGeneMark' or 'Prodigal'";
    }
    MOCATCore::execute_job( "gene_prediction", $cpu_gene_prediction, scalar @samples, "15gb" );
    if ( $conf{gene_prediction_software} eq 'MetaGeneMark' )
    {
      MOCATGenePredictionMetaGeneMark::post_check_files();
    } elsif ( $conf{gene_prediction_software} eq 'Prodigal' )
    {
      MOCATGenePredictionProdigal::post_check_files();
    }
    print localtime() . ": COMPLETED GENE PREDICTION.\n";
    print "Files stored in $cwd/<SAMPLE>/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.KXX/<SAMPLE>.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.KXX.*\n";
  }
}
### FILTER ###
if ( $do_filter[0] )
{
  if ( $systemType =~ m/Darwin/ && $conf{MOCAT_mapping_mode} eq 'allbest' )
  {
    die "ERROR & EXIT: Filtering in mode 'allbest' is not supported under OSX. Simply because there are no OSX executables of psort and msamtools. Not even the authos could make this work. Sorry.";
  }
  @databases = @do_filter;
  MOCATCore::print_settings("filter");
  print localtime() . ": PERFORMING FILTERING...\n";
  MOCATCore::pre_check( "", $reads, "f" );
  MOCATFilterV2::create_job( 'filter', $cpu_filter );
  MOCATCore::execute_job( "filter", $cpu_filter, scalar @samples, "15gb" );
  print localtime() . ": COMPLETED FILTERING.\n";
}

### STATS FASTQC ###
if ($do_stats_fastqc)
{
  MOCATCore::print_settings("stats_fastqc");
  print localtime() . ": EXECUTING FASTQC...\n";
  MOCATStatsFastqc::create_job( "stats_fastqc", 1 );
  my $fqs = MOCATCore::count_fq();
  MOCATCore::execute_job( "stats_fastqc", 1, $fqs, "15gb" );
  unless ($only_print)
  {
    MOCATStatsFastqc::post_check_files();
    print "Files stored in $cwd/<SAMPLE>/<LANE>_fastqc/*\n";
    print localtime() . ": COMPLETED FASTQC.\n";
  }
}
### SAMPLE STATUS ###
if ($do_sample_status)
{
  MOCATCore::print_settings("sample_status");
  print localtime() . ": CALCULATING SAMPLE STATUS...\n";
  MOCATSampleStatus_dev::run();
  MOCATSampleStatus::run();
  MOCATSampleStatusExcel2::run();
}

### FETCH MG ###
if ( defined $do_fetch_mg )
{
  $assembly = $do_fetch_mg;
  MOCATCore::print_settings("fetch_mg");
  print localtime() . ": FETCHING MARKER GENES...\n";
  MOCATFetchMGs::pre_check();
  MOCATFetchMGs::create_job( "fetch_mg", $cpu_fetch_mg );
  MOCATCore::execute_job( "fetch_mg", $cpu_fetch_mg, scalar @samples, "15gb" );
  print localtime() . ": COMPLETED FETCHING MARKER GENES.\n";
  print "Files stored in $cwd/<SAMPLE>/gene.fetched.MGs.$assembly.$reads.$conf{MOCAT_data_type}.KXX/<SAMPLE>.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.KXX.*\n";
}

### PROFILING ###
if ( $do_profiling[0] )
{
  @databases = @do_profiling;    # This line is used in both CalculateTaxonomy and PasteTaxonomyCoverageFiles
  unless ($only_paste)
  {
    MOCATCore::print_settings("profiling");
    print localtime() . ": PERFORMING PROFILING...\n";
    MOCATCore::pre_check( "", $reads, "c" );
    MOCATProfiling::create_job( 'profiling', $cpu_profiling );
    MOCATCore::execute_job( "profiling", $cpu_profiling, scalar @samples, "15gb" );
    print "$sep\n";
  }
  unless ( $no_paste || $only_print )
  {
    MOCATProfilingPaste::paste2( 'paste_profiling_files', 1 );
  }
  print localtime() . ": COMPLETED PROFILING.\n";
  print "$sep\nSOME PROFILES ARE SAVED IN (note, depending on settings in config file, these profiles may not have been created): $OUTPUT_FOLDER\n";

}

### MODULE SUPPORT ###
if ($MODULE_OPTION)
{
  print "MOCAT MODULE MODE ACTIVATED\n";
  print "      MODULE = $MODULE_OPTION\n";
  unless ( -e "$mod_dir/$MODULE_OPTION/$MODULE_OPTION.sh" )
  {
    die "ERROR & EXIT: Tried to launch module '$MODULE_OPTION', but could not find expected file $mod_dir/$MODULE_OPTION/$MODULE_OPTION.sh";
  }

  $modcomment = "";
  $MODCFG     = "$mod_dir/$MODULE_OPTION/$MODULE_OPTION.cfg";
  foreach my $n ( 0 .. scalar @MODULE_PARAM - 1 )
  {
    if ( $MODULE_PARAM[$n] =~ m/^-*modcfg/ )
    {
      $MODCFG = $MODULE_PARAM[ $n + 1 ];
      print "      MOD CFG FILE = $MODCFG\n";
    }
  }

  foreach my $n ( 0 .. scalar @MODULE_PARAM - 1 )
  {
    if ( $MODULE_PARAM[$n] =~ m/^-*name/ )
    {
      $modcomment = $MODULE_PARAM[ $n + 1 ];
      print "      ADD NAME = $modcomment\n";
    }
  }

  unless ( -e $MODCFG )
  {
    die "ERROR & EXIT: Tried to launch module '$MODULE_OPTION', but could not find expected cfg file $MODCFG";
  }
  chomp( $mod_single_job = `grep -v '^#' $MODCFG | grep 'MODULE_RUNS_SINGLE_JOB' | cut -f 2 -d'=' | sed 's/ //g'` );
  unless ( $mod_single_job eq 'FALSE' || $mod_single_job eq 'TRUE' )
  {
    die "ERROR & EXIT: MODULE_RUNS_SINGLE_JOB=$mod_single_job in $MODCFG. Expected it to be either TRUE or FALSE";
  }
  chomp( @mod_requested_all = `grep -v '^#' $MODCFG | grep '^REQUEST' ` );
  chomp( @mod_requested     = `grep -v '^#' $MODCFG | grep '^REQUEST' | cut -f 2 -d' ' ` );
  chomp( my @mod_requested2 = `grep -v '^#' $MODCFG | grep '^REQUEST' | sed 's%.*\\[%[%' ` );
  chomp( @mod_optional_all  = `grep -v '^#' $MODCFG | grep '^OPTIONAL' ` );
  chomp( @mod_optional      = `grep -v '^#' $MODCFG | grep '^OPTIONAL' | cut -f 2 -d' ' ` );
  chomp( @mod_overrides_all = `grep -v '^#' $MODCFG | grep '^OVERRIDE' | cut -f 3 -d' ' ` );
  chomp( @mod_overrides     = `grep -v '^#' $MODCFG | grep '^OVERRIDE' | cut -f 2 -d' ' ` );

  if ( scalar @mod_optional_all > 0 )
  {
    print "      " . ( join "\n      ", @mod_optional_all ) . "\n";
  }

  my $paste = 1;

  my %req_desc;
  for my $i ( 0 .. scalar @mod_requested - 1 )
  {
    $req_desc{ $mod_requested[$i] } = $mod_requested2[$i];
  }

  foreach my $n ( 0 .. scalar @MODULE_PARAM - 1 )
  {
    if ( $MODULE_PARAM[$n] =~ m/^-/ )
    {
      $paste = 0;
    }
    if ( $paste == 1 )
    {
      push @MODULE_DB, $MODULE_PARAM[$n];
    }
    foreach my $m (@mod_requested)
    {
      if ( $MODULE_PARAM[$n] =~ m/^-*$m$/ )
      {
        $mod_requests{$m} = $MODULE_PARAM[ $n + 1 ];
      }
    }
    foreach my $m (@mod_optional)
    {
      if ( $MODULE_PARAM[$n] =~ m/^-*$m$/ )
      {
        $mod_optionals{$m} = $MODULE_PARAM[ $n + 1 ];
      }
    }
  }
  foreach my $m ( 0 .. scalar @mod_overrides - 1 )
  {
    $mod_overrides{ $mod_overrides[$m] } = $mod_overrides_all[$m];
  }

  foreach my $k ( sort keys %mod_requests )
  {
    unless ( $mod_requests{$k} )
    {
      $mod_requests{$k} = "";
    }
    print "      SET      $k = $mod_requests{$k}\n";
  }
  foreach my $k ( sort keys %mod_optionals )
  {
    print "      SET      $k = $mod_optionals{$k}\n";
  }

  foreach my $k (@mod_requested)
  {
    unless ( $mod_requests{$k} )
    {
      die "ERROR & EXIT: Missing (at least) option -$k, please specify it; $req_desc{$k}";
    }
  }

  print "$sep\n";
  print localtime() . ": EXECUTING $MODULE_OPTION...\n";
  $RUN_MODULE = 1;
  MOCATModule::create_job( $MODULE_OPTION, $mod_single_job, $cpu_module );
  if ( $mod_single_job eq "TRUE" )
  {
    MOCATCore::execute_job( $MODULE_OPTION, $cpu_module, 1, "15gb" );
  } else
  {
    MOCATCore::execute_job( $MODULE_OPTION, $cpu_module, scalar @samples, "15gb" );
  }

  print "$sep\n";
  print "Output can generally be found in $cwd/$mod_overrides{outfolder}\n";
  print "Run log file: $MODLOG\n";
  print "$sep\n";
}

### EXIT ###
unless ($only_print)
{
  if ($print_final_text)
  {
    print "$sep\nCOMPLETED ALL TASKS SUCCESSFULLY!\n$sep\n";
  }
}
close(STDOUT);
exit(0);

sub Cite
{
  print "Latest information about MOCAT: http://mocat.embl.de

Citing MOCAT

   If you have used MOCAT2 in your work, please cite:
          Kultima, J. R., Coelho, L. P., Forslund, K., Huerta-Cepas, J., Li, S. S., Driessen, Bork, P. (2016).
          MOCAT2: a metagenomic assembly, annotation and profiling framework. Bioinformatics . http://doi.org/10.1093/bioinformatics/btw183
   
   
   If you have used MOCAT in your work, please cite:
           Kultima JR, Sunagawa S, Li J, Chen W, Chen H, et al. (2012)
               MOCAT: A Metagenomics Assembly and Gene Prediction Toolkit. PLoS ONE 7(10): e47656. doi:10.1371/journal.pone.0047656

   
   MOCAT2 is a wrapper for 3rd party software. Therefore we strongly suggest you also cite the following papers if you use MOCAT2:

   Initial read trimming and quality control
      - Cox MP, Peterson DA, Biggs PJ (2010) SolexaQA: At-a-glance quality assessment of Illumina second-generation sequencing data. BMC bioinformatics 11: 485 doi:10.1186/1471-2105-11-485. FastX program

   Mapping reads
     - Li R, Yu C, Li Y, Lam T-W, Yiu S-M, et al. (2009) SOAP2: an improved ultrafast tool for short read alignment. Bioinformatics (Oxford, England) 25: 1966–1967 doi:10.1093/bioinformatics/btp336. doi: 10.1093/bioinformatics/btp336
     - Edgar RC (2010) Search and clustering orders of magnitude faster than BLAST. Bioinformatics. 2010 Oct 1;26(19):2460-1. doi: 10.1093/bioinformatics/btq461

   Assembly
     - Li R, Zhu H, Ruan J, Qian W, Fang X, et al. (2010) De novo assembly of human genomes with massively parallel short read sequencing. Genome research 20: 265–272 doi:10.1101/gr.097261.109.

   Assembly revision
     - Li H, Durbin R (2009) Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics (Oxford, England) 25: 1754–1760 doi:10.1093/bioinformatics/btp324.

   Gene Prediciton
     - Hyatt D, Chen G-L, Locascio PF, Land ML, Larimer FW, et al. (2010) Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC bioinformatics 11: 119 doi:10.1186/1471-2105-11-119. 
     - Zhu W, Lomsadze A, Borodovsky M (2010) Ab initio gene identification in metagenomic sequences. Nucleic acids research 38: 1–15 doi:10.1093/nar/gkq275.

    Retrieving Marker Genes
     - Sunagawa S., et al. (2013) Metagenomic species profiling using universal phylogenetic marker genes. Nature Methods 10, 1196–1199 doi:10.1038/nmeth.2693

   Clustering gene catalogs
     - Limin Fu, et al. (2012) CD-HIT: accelerated for clustering the next generation sequencing data. Bioinformatics. doi: 10.1093/bioinformatics/bts565

   Annotating gene catalogs
     - Buchfink, B., et al. (2014) Fast and sensitive protein alignment using DIAMOND. Nat. Methods, 12, 59–60.
     - Arumugam M., et al. (2010) SmashCommunity: a metagenomic annotation and analysis tool. Bioinformatics. doi: 10.1093/bioinformatics/btq536.

   Taxnomic profiles
     - mOTU-LG: Sunagawa,S. et al. (2013) Metagenomic species profiling using universal phylogenetic marker genes. Nat. Methods, 10, 1196–9
     - specI & NCBI: Mende,D.R. et al. (2013) Accurate and universal delineation of prokaryotic species. Nat. Methods, 10, 881–4.

   Pre-compiled reference gene catalogs
     - IGC (human gut): Li,J. et al. (2014) An integrated catalog of reference genes in the human gut microbiome. Nat. Biotechnol., 32, 834–41.
     - CRC-RGC (human gut): Zeller,G. et al. (2014) Potential of fecal microbiota for early-stage detection of colorectal cancer. Mol. Syst. Biol., 10, 766.
     - skin-RGC (human skin): Oh,J. et al. (2014) Biogeography and individuality shape function in the human skin metagenome. Nature, 514, 59–64.
     - mouse-RGC (human skin): Xiao,L. et al. (2015) A catalog of the mouse gut metagenome. Nat Biotech, 33, 1103–1108.
     - OM-RGC (ocean): Sunagawa,S. et al. (2015) Structure and function of the global ocean microbiome. Science, 348 (6237), 1:10

   Functional profiles
     - eggNOG: Huerta-Cepas, J., et al. eggNOG 4.5: a hierarchical orthology framework with improved functional annotations for eukaryotic, prokaryotic and viral sequences. Nucleic Acids Res 2015. 10.1093/nar/gkv1248.
     - ARDB: Liu,B. and Pop,M. (2009) ARDB—Antibiotic Resistance Genes Database. Nucleic Acids Res. , 37 , D443–D447.
     - CARD: McArthur, A.G., et al. The comprehensive antibiotic resistance database. Antimicrob Agents Chemother 2013;57(7):3348-3357. 10.1128/AAC.00419-13.
     - DBETH: Chakraborty, A., et al. DBETH: a Database of Bacterial Exotoxins for Human. Nucleic Acid Res. 2012;40 Database issue):D615-20 0.1093/nar/gkr942
     - dbCAN: Yin, Y., et al. dbCAN: a web resource for automated carbohydrate-active enzyme annotation. Nucleic Acids Res 2012;40(Web Server issue):W445-451. 10.1093/nar/gks479.
     - DrugBank: Knox C., et al. DrugBank 3.0: a comprehensive resource for 'omics' research on drugs. Nucleic Acids Res 2011;39(Database issue):D1035-41. 10.1093/nar/gkq1126
     - ICEberg: Bi D., et al. ICEberg: a web-based resource for integrative and conjugative elements found in Bacteria. Nucleic Acids Res. 2012 Jan;40(Database issue):D621-6. 10.1093/nar/gkr846.
     - KEGG: Kanehisa, M., et al. Data, information, knowledge and principle: back to metabolism in KEGG. Nucleic Acids Res 2014;42(Database issue):D199-205. 10.1093/nar/gkt1076.
     - MetaCyc: Caspi, R., et al. The MetaCyc database of metabolic pathways and enzymes and the BioCyc collection of pathway/genome databases. Nucleic Acids Res 2015. 10.1093/nar/gkv1164.
     - MvirDB: Zhou, C.E., et al. MvirDB--a microbial database of protein toxins, virulence factors and antibiotic resistance genes for bio-defence applications. Nucleic Acids Res 2007;35(Database issue):D391-394. 10.1093/nar/gkl791.
     - PATRIC: Mao, C., et al. Curation, integration and visualization of bacterial virulence factors in PATRIC. Bioinformatics 2015;31(2):252-258. 10.1093/bioinformatics/btu631.
     - Pfam: Finn, R.D., et al. Pfam: the protein families database. Nucleic Acids Res 2014;42(Database issue):D222-230. 10.1093/nar/gkt1223.
     - Prophages: Waller, A.S., et al. Classification and quantification of bacteriophage taxa in human gut metagenomes. ISME J 2014;8(7):1391-1402. 10.1038/ismej.2014.30.
     - Resfams: Gibson, M.K., et al. Improved annotation of antibiotic resistance determinants reveals microbial resistomes cluster by ecology. ISME J 9(1). 10.1038/ismej.2014.106.
     - SEED subsystems: Overbeek, R., et al. The SEED and the Rapid Annotation of microbial genomes using Subsystems Technology (RAST). Nucleic Acids Res 2014;42(Database issue):D206-214. 10.1093/nar/gkt1226.
     - Superfamily: Gough, J., et al. Assignment of homology to genome sequences using a library of hidden Markov models that represent all proteins of known structure. J Mol Biol 2001;313(4):903-919. 10.1006/jmbi.2001.5080.
     - vFam: Skewes-Cox, P., et al. Profile hidden Markov models for the detection of viruses within metagenomic sequence data. PLoS One 2014;9(8):e105067. 10.1371/journal.pone.0105067.
     - VFDB: Chen, L., et al. VFDB 2012 update: toward the genetic diversity and molecular evolution of bacterial virulence factors. Nucleic Acids Res 2012;40(Database issue):D641-645. 10.1093/nar/gkr989.
     - Victors: Mao, C., et al. Curation, integration and visualization of bacterial virulence factors in PATRIC. Bioinformatics 2015;31(2):252-258. 10.1093/bioinformatics/btu631.
     
";
  exit 0;
}

sub Usage
{
  welcome();
  print "\n                    ";
  print color 'bold';
  print color 'underline';
  print "Full manual & FAQ: MOCAT.pl -man\n";
  print color 'reset';
  print "\n                    ";
  print color 'bold';
  print color 'underline';
  print "How to cite MOCAT: MOCAT.pl -cite\n";
  print color 'reset';
  print "\n            ";
  print color 'bold';
  print color 'underline';
  print "Have you tried the wrapper runMOCAT.sh? Try it!\n";
  print color 'reset';
  print "\nUsage: ";
  print color 'bold';
  print "MOCAT.pl";
  print color 'reset';
  print " -sf|sample_file 'FILE' [";
  print color 'bold';
  print "Pipeline";
  print color 'reset';
  print ", Statistics, & Additional Options]\n\n";
  print " 'FILE'
   Contains the list of folder names (sample names), one per line,\n   in which the raw sample data is located\n\n";
  print color 'underline';
  print "Examples\n\n";
  print color 'reset';
  print "Process, Assemble, Revise Assembly, Predict Genes, cluster genes into gene catalog, annotate gene catalog, profile against gene catalog
                            MOCAT.pl -sf my.samples -rtf
                            MOCAT.pl -sf my.samples -a
                            MOCAT.pl -sf my.samples -gp assembly
                            MOCAT.pl -sf my.samples -make_gene_catalog -assembly_type assembly
                            MOCAT.pl -sf my.samples -annotate_gene_catalog
                            MOCAT.pl -sf my.samples -s my.samples.padded -identity 95
                            MOCAT.pl -sf my.samples -f my.samples.padded -identity 95
                            MOCAT.pl -sf my.samples -p my.samples.padded -identity 95 -mode functional\n\n";
  print "Assemble and predict genes: MOCAT.pl -sf my.samples -rtf
  (no screen)               MOCAT.pl -sf my.samples -a
                            MOCAT.pl -sf my.samples -gp assembly
  fetch marker genes:       MOCAT.pl -sf my.samples -fmg assembly
                            MOCAT.pl -sf my.samples -ss\n\n";
  print "Assemble and predict genes: MOCAT.pl -sf my.samples -rtf
  (DB screen)               MOCAT.pl -sf my.samples -s hg19 -screened_files -identity 90
                            MOCAT.pl -sf my.samples -a -r hg19
                            MOCAT.pl -sf my.samples -gp assembly -r hg19
                            MOCAT.pl -sf my.samples -ss\n\n";
  print "Assemble and predict genes: MOCAT.pl -sf my.samples -rtf
  (remove eg. adapters      MOCAT.pl -sf my.samples -sff adapters.fa -screened_files
   and then DB screen)      MOCAT.pl -sf my.samples -bwa hg19 -r adapters.fa  -screened_files
                            MOCAT.pl -sf my.samples -a -r screened.adapters.fa.on.hg19
                            MOCAT.pl -sf my.samples -gp assembly -r screened.adapters.fa.on.hg19
                            MOCAT.pl -sf my.samples -ss\n\n";
  print color 'underline';
  print "Pipeline Options\n\n";
  print color 'reset';

  print " -r|reads ['reads.processed', 'DATABASE' or 'FASTA FILE']
   Required for all pipeline options, except rtf|read_trim_filter
   Specify whether processing trim & filtered, or screened reads.
   A default value to this setting can also be specified in config file\n\n";
  print " -e|extracted
   Optional for all pipeline options, except rtf|read_trim_filter, see full manual\n\n";
  print "\n";
  print color 'bold';
  print " -rtf|read_trim_filter\n";
  print color 'reset';
  print "   performs trimming and filtering of reads\n\n";
  print color 'bold';
  print " -a|assembly\n";
  print color 'reset';
  print "   Performs assembly of reads\n\n";
  print color 'bold';
  print " -ar|assembly_revision\n";
  print color 'reset';
  print "   Further improves assemblies\n\n";
  print color 'bold';
  print " -gp|gene_prediction ['assembly', 'assembly.revised']\n";
  print color 'reset';
  print "   Predicts protein coding genes on assemblies\n\n";
  print color 'bold';
  print " -fmg|fetch_mg ['assembly', 'assembly.revised']\n";
  print color 'reset';
  print "   Extracts marker genes among the predicted genes\n\n";
  print color 'bold';
  print " -soap|bwa ['DB1 DB2 ...',s,c,f,r]\n";
  print color 'reset';
  print "   Screen, extract and map reads against a reference databse (hg19 is provided) or (s)acftigs,
   (c)ontigs, sca(f)folds from an assembly, or scaftigs from a (r)evised assembly.
   This mapping step uses SOAPaligner2 (soap) or BWA (bwa).
   Additional options:
    -screened_files : If set, screened read files are generated, these are reads not matching the DB
    -extracted_files : If set, extracted read files are generated, these are reads matching the DB
    -use_mem  : If set, copies the DB into memory for faster loading\n\n";
  print color 'bold';
  print " -sff|screen_fastafile 'FASTA FILE'\n";
  print color 'reset';
  print "   Same as 's|screen' above, but uses USearch, rather than SOAPaligner2.\n\n";
  print color 'bold';
  print " -fsoap ['DB1 DB2 ...',s,c,f,r]\n";
  print color 'reset';
  print "   Filter screened reads, (s)caftigs, (c)ontigs, sca(f)folds or (r)evised assembly scaftigs
    at higher \%ID and length cutoff. This step has to be run before calculating profiles if the option soap was used\n\n";
  print color 'reset';
  print "   Additional options:
    -shm   : If set, faster, but saves data for the filtering step in /dev/shm/<USER>
	\n";
  print color 'bold';
  print " -psoap|pbwa ['DB1 DB2 ...',s,c,f,r] -m|mode [gene, NCBI, mOTU, functional] -o [OUTPUT FOLDER]\n";
  print color 'reset';
  print "   Generate gene, mOTU, NCBI or functional profiles on filtered reads,
   (s)caftigs, (c)ontigs, sca(f)folds or (r)evised assembly scaftigs. 
   If -mode is set to either NCBI or mOTU, it is expected that the 
   reads have been correctly mapped to the corresponding databases.
   Specify psoap if you used the command 'soap' previously, and 'pbwa' if you used 'bwa'.
   Additional options:
    -no_horizontal : No not calculate horizontal gene & functional coverages
    -verbose       : Prints extra information about status of profiling steps
    -shm           : Faster, but saves 2-5 GB of data for the profiling step in /dev/shm/<USER>
    -uniq          : Specify this flag if you find duplicated row names
                     (e.g. if you have mapped to a DB where the same reference appears multiple times)\n\n";

  my %mods;
  my %mod_text;
  chomp( $mod_dir = `which MOCAT.pl | sed 's|src/MOCAT.pl||' | sed 's/\$/mod/'` );
  chomp( my @tmp = `ls $mod_dir 2>/dev/null` );
  foreach my $tmp (@tmp)
  {
    if ( -e "$mod_dir/$tmp/$tmp.sh" && -e "$mod_dir/$tmp/$tmp.cfg" )
    {
      chomp( @mod_requested_all = `grep -v '^#' $mod_dir/$tmp/$tmp.cfg | grep '^REQUEST' | sed 's/REQUEST //'` );
      @{ $mods{$tmp} } = @mod_requested_all;
      chomp( my $mod_text = `grep -v '^#' $mod_dir/$tmp/$tmp.cfg | grep '^TEXT' | sed 's/TEXT //'` );
      $mod_text{$tmp} = $mod_text;

    }

  }

  print color 'reset';
  print color 'underline';
  print "Available modules";
  print color 'reset';
  print "\n\n These are installed in the folder $mod_dir\n";
  print " Each module requires a NAME.sh and NAME.cfg file inside the NAME folder\n\n";

  foreach my $mod ( sort keys %mods )
  {

    print color 'bold';
    print " -$mod $mod_text{$mod}\n";
    print color 'reset';
    print "   Required options:\n";
    foreach my $mod2 ( @{ $mods{$mod} } )
    {
      print "    -$mod2\n";
    }
    print "\n";
  }

  print "\n";
  print color 'reset';
  print color 'underline';
  print "Statistics Options";
  print color 'reset';
  print "

 -sfq|stats_fastqc
   Produces statistics for each lane with raw reads using the FastQC toolkit
 -ss|sample_status
   Prints a simple view how the processing status of each sample,
   and stores this in <sample_file>.status\n\n";
  print color 'underline';
  print "Additional Options";
  print color 'reset';
  print "

 -cfg|config [file]
   Specify another config file than MOCAT.cfg
 -x|no_execute
   Only create job scripts, but don't execute them
 -nt|no_temp
   Overrides any specified temp folders config file
 -cpus [integer]
   Not recommended, but specifies a fixed number of cores for each job,
   please read the full manual using MOCAT.pl -man
 -host [hostname]
   Runs the jobs on a different host machine
 -identity [integer]
   Overrides any percentage cutoff setting in cfg file
 -length [integer]
   Overrides any length cutoff setting in cfg file
 -memory XGB
   If queuing system is SGE or LSF, it will require XGB of RAM for the job
   This can also be set with the respective memory options by adding these
   to the param fields in the config file
 -config A=b C=d
   Overrides setting A from the config file with b, etc
   \n";

  MOCATUnpublished::Usage();    #DEV VER#

  exit 0;
}

sub tiger
{

  print STDOUT <<"EOF";
                                 ___..........__
           _,...._           _."'_,.++8n.n8898n.`"._        _....._
         .'       `".     _.'_.'" _.98n.68n. `"88n. `'.   ,"       `.
        /        .   `. ,'. "  -'" __.68`""'""=._`+8.  `.'     .     `.
       .       `   .   `.   ,d86+889" 8"""+898n, j8 9 ,"    .          \
      :     '       .,   ,d"'"   _..d88b..__ `"868' .'  . '            :
      :     .      .    _    ,n8""88":8"888."8.  "               '     :
       \     , '  , . .88" ,8P'     ,d8. _   `"8n  `+.      `.   .     '
        `.  .. .     d89' "  _..n689+^'8n88n.._ `+  . `  .  , '      ,'
          `.  . , '  8'    .d88+"    j:""' `886n.    b`.  ' .' .   ."
           '       , .j            ,d'8.         `  ."8.`.   `.  ':
            .    .' n8    ,_      .f A 6.      ,..    `8b, '.   .'_
          .' _    ,88'   :8"8    6'.d`i.`b.   d8"8     688.  ".    `'
        ," .88  .d868  _         ,9:' `8.`8   "'  ` _  8+""      b   `,
      _.  d8P  d'  .d :88.     .8'`j   ;+. "     n888b 8  .     ,88.   .
     `   :68' ,8   88     `.   '   :   l `     .'   `" jb  .`   688b.   ',
    .'  .688  6P   98  =+""`.      `   '       ,-"`+"'+88b 'b.  8689  `   '
   ;  .'"888 .8;  ."+b. : `" ;               .: "' ; ,n  `8 q8, '88:       \
   .   . 898  8:  :    `.`--"8.              d8`--' '   .d'  ;8  898        '
  ,      689  9:  8._       ,68 .        .  :89    ..n88+'   89  689,' `     .
  :     ,88'  88  `+88n  -   . .           .        " _.     6:  `868     '   '
  , '  .68h.  68      `"    . . .        .  . .             ,8'   8P'      .   .
  .      '88  'q.    _.f       .          .  .    '  .._,. .8"   .889        ,
 .'     `898   _8hnd8p'  ,  . ..           . .    ._   `89,8P    j"'  _   `
  \  `   .88, `q9868' ,9      ..           . .  .   8n .8 d8'   +'   n8. ,  '
  ,'    ,+"88n  `"8 .8'     . ..           . .       `8688P"   9'  ,d868'   .  .
  .      . `86b.    " .       .            ..          68'      _.698689;      :
   . '     ,889_.n8. ,  ` .   .___      ___.     .n"  `86n8b._  `8988'b      .,6
    '       q8689'`68.   . `  `:. `.__,' .:'  ,   +   +88 `"688n  `q8 q8.     88
    , .   '  "     `+8 n    .   `:.    .;'   . '    . ,89           "  `q,    `8
   .   .   ,        .    + c  ,   `:.,:"        , "   d8'               d8.    :
    . '  8n           ` , .         ::    . ' "  .  .68h.             .8'`8`.  6
     ,    8b.__. ,  .n8688b., .    .;:._     .___nn898868n.         n868b "`   8
      `.  `6889868n8898886888688898"' "+89n88898868868889'         688898b    .8
       :    q68   `""+8688898P ` " ' . ` '  ' `+688988P"          d8+8P'  `. .d8
       ,     88b.       `+88.     `   ` '     .889"'           ,.88'        .,88
        :    '988b        "88b.._  ,_      . n8p'           .d8"'      '     689
        '.     "888n._,      `"8"+88888n.8,88:`8 .     _ .n88P'   .  `      ;88'
         :8.     "q888.  .            "+888P"  "+888n,8n8'"      .  .     ,d986
         :.`8:     `88986                          `q8"           ,      :688"
         ;.  '8,      "88b .d                        '                  ,889'
         :..   `6n      '8988                                         b.89p
         :. .    '8.      `88b                                        988'
         :. .      8b       `q8.        '                     . '   .d89      '
         . .        `8:       `86n,.       " . ,        , "        ,98P      ,
         .. .         '6n.       +86b.        .      .         _,.n88'     .
           .            `"8b.      'q98n.        ,     .  _..n868688'          .
          ' . .            `"98.     `8868.       .  _.n688868898p"            d
           . .                '88.      "688.       q89888688868"            ,86
            '. .                 88.     `8898        " .889"'              .988
EOF

  print "\n";

}

__END__


=head1 Name

B<MOCAT - Metagenomics Analysis Toolkit>

=head1 Description

B<MOCAT>, is a software package to process Illumina metagenomic data. This software provides a pipeline for processing raw reads, assembly, gene prediction, extracting marker genes and mapping reads to external databases

=head1 Synopsis

> MOCAT.pl [Required Options] [Pipeline Options] [Statistics Options] [Additional Options]

Or, try the wrapper script, which runs multiple MOCAT commands at once:

> runMOCAT.sh

=head1 MOCAT modules

Look in the MOCAT/mod folder for examples how to create modules. Run MOCAT.pl without options for available modules and their options.

=head1 Example - Process, Assemble, Revise Assembly, Predict Genes, cluster genes into gene catalog, annotate gene catalog, profile against gene catalog

=over

=item MOCAT.pl -sf my.samples -rtf

=item MOCAT.pl -sf my.samples -a

=item MOCAT.pl -sf my.samples -gp assembly

=item MOCAT.pl -sf my.samples -make_gene_catalog -assembly_type assembly

=item MOCAT.pl -sf my.samples -annotate_gene_catalog

=item MOCAT.pl -sf my.samples -s my.samples.padded -identity 95

=item MOCAT.pl -sf my.samples -f my.samples.padded -identity 95

=item MOCAT.pl -sf my.samples -p my.samples.padded -identity 95 -mode functional

=back

=head1 Example - Process, Assemble, Predict Genes and fetch marker genes

=over

=item MOCAT.pl -sf my.samples -rtf

=item MOCAT.pl -sf my.samples -a

=item MOCAT.pl -sf my.samples -gp assembly

=item MOCAT.pl -sf my.samples -fmg assembly

=item MOCAT.pl -sf my.samples -ss

=back

=head1 Example - Process, Screen against DB, Assemble, and Predict Genes

=over

=item MOCAT.pl -sf my.samples -rtf

=item MOCAT.pl -sf my.samples -bwa hg19 -screened_files

=item MOCAT.pl -sf my.samples -a -r hg19

=item MOCAT.pl -sf my.samples -gp assembly -r hg19

=item MOCAT.pl -sf my.samples -ss

=back

=head1 Example - Process, screen against fasta file, against DB, Assemble, and Predict Genes

=over 

=item MOCAT.pl -sf my.samples -rtf

=item MOCAT.pl -sf my.samples -sff adapters.fa

=item MOCAT.pl -sf my.samples -bwa hg19 -r adapters.fa -screened_files

=item MOCAT.pl -sf my.samples -a -r screened.adapters.fa.on.hg19

=item MOCAT.pl -sf my.samples -gp assembly -r screened.adapters.fa.on.hg19

=item MOCAT.pl -sf my.samples -ss

=back

=head1 Example - Calculate coverage of reads mapping a custom made database

=over

=item MOCAT.pl -sf my.samples -rtf

=item MOCAT.pl -sf my.samples -s /User/My/Path/REF_DB.fna

=item MOCAT.pl -sf my.samples -f REF_DB.fna

=item MOCAT.pl -sf my.samples -p REF_DB.fna

=item MOCAT.pl -sf my.samples -ss

=back

=head1 Example Sample File

=over

=item MH0001

=item GOS_STATION_2

=item MMMS_0561

=back

=head1 FAQ

=head2 How do I calculate the read coverage of each contigs, scaffold or scaftig of an assembly?

MOCAT can calculate the total read and base coverages (length normalized and on normalized) for each contig, scaftig or scaffold. To do this, first run MOCAT.pl -sf sample -a to assemble your samples. Then execute MOCAT.pl -sf sample -s [s|c|f] to map reads to the assembled scaftigs, contigs or scaffolds. After this run MOCAT.pl -sf sample -f [s|c|f] to filter the mapped reads are selected percentage and length cutoff. Finally run MOCAT.pl -sf sample -p [s|c|f] to calculate the coverages. 

=head1 Options

=head2 Required Options

=over

=over 

=item B<C<-sf|sample_file 'FILE'> (required)>

'FILE' contains the list of folder names (sample names), one per line, in which the raw sample data is located. See section 'Setup MOCAT in a new folder' for more details.

=back

=back

=head2 Pipeline Options

=over

=over

=item B<C<-rtf|read_trim_filter>>

Performs trimming and quality filtering of raw reads and stores them in SAMPLE/reads.proocessed/*[pair|single]*.fq.gz.

B<For supported file formats (old and new Illumina format), see section Supported Formats below.>

NOTE! Both 'fastx' and 'solexaqa' requires the two LANE.1.fq and LANE.2.fq input files to contain reads in the same order (the files may have different number of reads, but their respective order must be the same). This is usually the case when the FastQ file comes directly from the sequencing machine, however if a pre-filtering step has been performed, the order may have changed.

B<Additional note:> If 'readtrimfilter_use_precalc_5prime_trimming' in the config file is set to 'yes', then the file MOCAT.cutoff5prime must exist. For more information, see Config File Section.


=item B<C<-soap|bwa ['DB1 DB2 ...' or 's' or 'c' or 'f'] >>

I<Additionally required options>

=over

=over

=item C<-r|reads ['reads.processed', 'DATABASE NAME']> (required)

This is to specify whether reads should be extracted from reads that have been a) trimmed and filtered, or b) trimmed, filtered and screened against a custom databse. See below.

=back

=back

I<Additional options>

=over

=over

=item C<-e|extract> (optional)

Set this flag if you want to screen the extracted reads (reads that matched the database), rather than the screened (reads tat did not match) from a database.

=item C<-screened_files> (optional)

If set, reads are mapped and the files in the .screened folder are produced, but not the .extracted folder.

=item C<-extracted_files> (optional)

If set, reads are mapped and the files in the .extracted folder are produced, but not the .screened folder.

=item C<-use_mem> (optional)

This option can be set when running the -s option. If set, it copies the selected database into memory and then sets the data folder to the memory location. This will speed up loading the database for each sample, because it\s read from memory instead of disk.

=back

=back

Trimmed and filtered reads will be aligned against a custom database (or multiple databases) located in the MOCAT data folder, using SOAPAligner2 (if command used is soap) or BWA (if command used is bwa). This databse has to be previously generated and the name to specify is the filenames of the database, excluding '.index.xxx'. Provided with MOCAT is the 'hg19' database, which can for example be used to screen for human contamination in procaryotic samples. Note, if a paired end read is filtered, then it's pair mate is also removed. For a detailed description of the database files required and created, see the taxonomic profiling section below. Note that if you wish to use multiple databases (a split database) (only works with option soap) the naming convention must follow DBNAME.1 DBNAME.2 etc.



=item B<C<-sff|screen_fastafile 'FASTA FILE'>>

I<Additionally required options>

=over

=over

=item C<-r|reads ['reads.processed', 'DATABASE NAME']> (required)

This is to specify whether reads should be extracted from reads that have been a) trimmed and filtered, or b) trimmed, filtered and screened against a custom databse. See below.

=back

=back

Reads that have been trimmed and filtered (and possibly screened against a custom database) can be screened against a fasta file with DNA sequences. This screening step is a BLAST search using USearch. 

FASTA FILE is an input file with sequences, against which, the reads will be blasted. a) Trimmed and filter reads OR b) reads that have been both trimmed & filter and additionally screened using a database, can be screened against a fasta file. To screen reads that have been only trimmed and filtered specify 'read_trim_filter'. To screen reads from previsouly screened reads (using -db), specify the database name.

B<IMPORTANT NOTE:> To maximize accuracy, RAW reads are blasted against the specified fasta file. Then, reads that have a blast mastch are removed from either reads that have been a) trimmed and filter only, or b) trimmed and filterd and additionally screened using a database.



=item B<C<-a|assembly>>

I<Additionally required options>

=over

=over

=item C<-r|reads ['reads.processed', 'DATABASE NAME']> (required)

This is to specify whether reads should be extracted from reads that have been a) trimmed and filtered, or b) trimmed, filtered and screened against a custom databse. See below.

=back

=back

The option specifies from which type of reads the assembly should be made. a) 'reads.processed' refers to reads that were only trimmed and filtered. b) 'DATABASE NAME' are reads that were additionally screened using a database. c) 'FASTA FILE' are reads that were screened against a fasta file.

Reads can be assembled using SOAPDenovo. This step assembles trimmed and filtered reads (-ar|assembly_reads 'reads.processed') or reads that were previously screened using a database and/or fasta file (-ar|assembly_reads 'DATABASE NAME' or -ar|assembly_reads 'FASTA FILE') into scaftigs and contigs. The assembly is done in two steps: 1. calculate insert sizes, 2. assembly. Calculating insert sizes can be done by aligning a subset of reads to a database (we have provided '1506MG'), or by pre-assembly. Mapping reads to a database may be more accurate and faster. The choice of aligning reads or assembly, is specified in the MOCAT config file using the tag 'assembly_calculate_insert_size_using' set to 'mapping' or 'assembly'.

B<Additional note:> If 'assembly_calculate_insert_size_using' in the config file is set to 'assembly', the file 'MOCAT.preliminary_insert_sizes' must exist. For more information, see Config File Section.



=item B<C<-ar|assembly_revision>>

I<Additionally required options>

=over

=over

=item C<-r|reads ['reads.processed', 'DATABASE NAME']> (required)

This is to specify whether reads should be extracted from reads that have been a) trimmed and filtered, or b) trimmed, filtered and screened against a custom databse. See below.

=back

=back

The option specifies from which type of reads the assembly was made. a) 'reads.processed' refers to reads that were only trimmed and filtered. b) 'DATABASE NAME' are reads that were additionally screened using a database. c) 'FASTA FILE' are reads that were screened against a fasta file.

Corrects assemblies for base pair and short indel errors. First, reads are aligned to the scaftigs using BWA to detect single base and short indel errors. For regions with low coverage (<10x): a) if there is a majority of one base (>0.7), the base of the scaftig, in that position, is modified accordingly. b) Scaftigs have short indels either inserted or deleted, if the read support is >0.5. Then, reads are aligned to the 'new' scaftigs using SOAPAligner2. The assembly revision can be made on assemblies that were made using reads that were a) trimmed and filtered ('reads.processed'), or b) trimmed, filtered and screened against a custom databse ('DATABASE NAME'), or c) trimmed, filtered, screened against a custom database and against a fasta file ('FASTA FILE')



=item B<C<-gp|gene_prediction ['assembly' or 'assembly.revised']>>

I<Additionally required options>

=over

=over

=item C<-r|reads ['reads.processed', 'DATABASE NAME' or 'FASTA FILE']> (required)

Specifies from which type of reads the assembly was made. Reads that were only trimmed and filtered, or also additionally screened using either a database and/or a fasta file.

=back

=back

Main option specifies whether to use a non revised, or revised assembly.

Predicts genes (ORF and protein sequences) from an assembly using MetaGeneMark or Prodigal (specified in config file). 



=item B<C<-fmg|fetch_mg ['assembly' or 'assembly.revised']>>

I<Additionally required options>

=over

=over

=item C<-r|reads ['reads.processed', 'DATABASE NAME' or 'FASTA FILE']> (required)

Specifies from which type of reads the assembly was made. Reads that were only trimmed and filtered, or also additionally screened using either a database and/or a fasta file.

=back

=back

Main option specifies whether to use a non revised, or revised assembly.

Fetches (extracts) 40 universial marker genes (Ciccarelli et al., Science, 2006 and Sorek et al., Science, 2007) from the set of predicted genes from each sample. In brief, Pre-built HMM models are used to identify protein coding sequencing matching these 40 universal marker genes. It is possible to train this model to retrieve a different set of genes, for instructions how to do this, execute 'MOCATFetchMGs03.pl' found in the MOCAT bin directory. To use these newly configured bit score cutoffs to retrieve sequences, replace the generated 'MG_BitScoreCutoffs.txt', with the original one in the MOCAT lib/fetchMG folder.



=item B<C<-fsoap ['DB1 DB2 ...','s','c','f','r']>>

I<Additionally required options>

=over

=over

=item C<-r|reads ['reads.processed', 'DATABASE NAME' or 'FASTA FILE']> (required)

Specifies from which type of reads the assembly was made. Reads that were only trimmed and filtered, or also additionally screened using either a database and/or a fasta file.

=back

=back

I<Additional options>

=over

=over

=item C<-e|extract>

Set this flag if, in the screen step, this flag was set.

=item C<-shm>

If set, faster, but saves temporary data to /dev/shm/<USER> instead of TEMP directory

=back

=back

Additional filtering of reads that have been mapped using the soap command (required for soap, but not for bwa). Percentage ID and length cut off are defined in the config file. With this option you can map reads at a lower cutoff and then filter them at different (higher) cutoffs. It is required to run this step before running profiling below, even though you do not want to filter at a higher cut off. If you do not want to filter at a higher cut off, set the values in the config file to the same ones as for the screen step.



=item B<C<-psoap|pbwa ['DB1 DB2 ...','s','c','f','r']>>

I<Additionally required options>

=over

=over

=item C<-r|reads ['reads.processed', 'DATABASE NAME' or 'FASTA FILE']> (required)

Specifies from which type of reads the assembly was made. Reads that were only trimmed and filtered, or also additionally screened using either a database and/or a fasta file.

=back

=back

=over

=over

=item C<-m|mode ['gene', 'NCBI', 'mOTU', 'functional']>

Specify psoap if you used soap in the screen step, or pbwa if you used bwa

Basic mode (gene): The (gene) coverages of the database are calculated. Of course, the sequences in the database does not have to be genes, but as MOCAT was designed to be used with gene catalogs, we call this mode 'gene'. :)

Taxonomic modes (NCBI, mOTU): Specify whether the taxonomic profiles to be calculated are based on reference marker genes or reference genomes (NCBI), or metagenomic OTUs (mOTU). These two settings both
have specific requirements, both on the database file and also additionally requireed files.

Functional: If this is specified, and an additional functional.map file exists for the database, the gene abundances are first calculated and then summaried at the different funcitonal levels. By default the levels are cog, ko, module and pathway, but this could be modified for a custom database.

=back

=back

=over

=over

=item C<-o|output [OUTPUT FOLDER]>

Easy to use links to the generated abundances tables are saved in the OUTPUT folder. These files are the main output of the taxonomic and mOTU profiling steps.

=back

=back

I<Additional options>

=over

=over

=item C<-e|extract>

Set this flag if, in the screen step, this flag was set.

=item C<-no_horizontal>

Do not calculate horizontal gene & functional coverages

=item C<-verbose>

Prints additional information about profiling status

=item C<-uniq>

Specify this flag if you find duplicated row names (e.g. if you have mapped to a DB where the same reference appears multiple times)

=item C<-shm>

If set, faster, but saves profiling temporary data (2-5 GB per sample) to /dev/shm/<USER> instead of TEMP directory

=back

=back

If mode is set to 'gene', then calculates the base and insert counts (and normalized by reference ID length) for reads that have been mapped and filtered against a database using first the screen and then filter commands.

Specify 'NCBI' to generates taxonomic profiles by summarizing base and insert coverages into abundances of taxa (kingdom up to species), or specify 'mOTU' to summarize into abundances of mOTUs. The two different modes can not be used on the same type of database. The database need to be constructed in a specific way, in order to run this step. MOCAT ships with to different databases, one for each mode. 'RefMG.v1.padded' is used when you wish to calculate taxonomic profiles with NCBI taxa names, and the 'mOTU.v1.padded' is used for generating mOTU species cluster abundances. If you want to use the provided databases, it is enough to run the screen and filter steps and then this step on either the database RefMG.v1.padded or mOTU.v1.padded, and specifying mode to 'NCBI' and 'mOTU', respectively, when running the this step. 

By specifying mode to be 'functional' is it possible to summarize gene abundances into functional categorizes. Some publicly available database have been annotated to various functional categories. By default MOCAT is expected to be used to summarize into cog, ko, module and pathway abundances. Note that the mapping file DBNAME.functional.map must either be downloaded for the database, or manually generated.

I<NCBI mode: Design and requirements for summarizing gene profiles into kingdom up to species abundances>

The names within () are example names and the convention of the file names need to be used. This means if the database is called 'DB', the required files are called 'DB.xxxxx'.

=over

=over

=item The database file (RefMG.v1.padded)

The fasta headers MUST be on the format '>taxaid.XXXX'. This format is required for the taxonomic profiling step to identify from which taxa the sequence is.

=back

=back


=over

=over

=item The database NCBI map file (RefMG.v1.padded.NCBI.map)

This file must be generated manually. However, MOCAT ships with a version of this file 'RefMG.v1.padded.NCBI.map'. This file can be copied and renamed if you have designed a different databse using the NCBI taxa ids, as described above for the database file. It is a mapping file from taxa ids to different taxonomic levels, and has the format '<TaxID>\t<Kingdom>\t<Phylum>\t<Class>\t<Order>\t<Family>\t<Genus>\t<Species>\t<specI_clusters>\t<length of taxa ID>'

The length field contains the summarized length of all the sequences corresponding to each taxa id in the database.


CuratedSpecies represents the species clusters as described in Mende et al (2013).

=back

=back


=over

=over

=item The database file .coord file (RefMG.v1.padded.coord)

This file contains information of which bases should be counted when calculating base and insert coverages. Normally this file is generated automatically and contains the fasta entry names with start and stop bases, normally start=1 and stop=length of that sequence. Changing this file may be useful if you want to calculate the coverages based on only a part of a gene or genome.

Example line: 1001582.LL3_00013<TAB>101<TAB>1378

=back

=back


=over

=over

=item The database file .len file (RefMG.v1.padded.len)

This file contains the fasta entries and their respective lengths, one per line. It is generated automatically and do not need to be changed.

Example line: 1001582.LL3_00013<TAB>1478

=back

=back


=over

=over

=item The database file .rownames and .rownames.uniq files (RefMG.v1.padded.rownames[.uniq])

These files contains the rownames of the computed base and insert coverages. These files are generated automatically and do not need to be modified. I the .uniq file the _X_X corresponds to the region over which the coverages are calculated. This is normally _1_[length of sequence]. This region can be changed in the .coord file.

=back

=back



I<mOTU mode: Design and requirements for summarizing gene profiles into species cluster (mOTU) abundances>

The names within () are example names and the convention of the file names need to be used. This means if the database is called 'DB', the required files are called 'DB.xxxxx'.

=over

=over

=item The database file (mOTU.v1.padded) (same as for NCBI mode)

The fasta headers must be defined in the map file. As long as the have been defined in this file, they can have any name.

=back

=back


=over

=over

=item The database NCBI map file (mOTU.v1.padded.motu.map) (NOT same as for NCBI mode)

This file must be generated manually. However, MOCAT ships with a version of this file '263MetaRefv9MG.cal.v2.seed.padded.motu.map'. This file can NOT be copied and renamed if you have designed a different databse. It is a mapping file from fasta identifiers in the database to different species clusters, and has the format '<fasta identifier><TAB><sequence length><TAB>COG id><cluster ID>'

Example line: 351605.Gura_3679<TAB>1095<TAB>COG0012<TAB>OTUcal.v2.0

=back

=back


=over

=over

=item The database file .coord file (mOTU.v1.padded.coord) (same as for NCBI mode)

This file contains information of which bases should be counted when calculating base and insert coverages. Normally this file is generated automatically and contains the fasta entry names with start and stop bases, normally start=1 and stop=length of that sequence. Changing this file may be useful if you want to calculate the coverages based on only a part of a gene or genome.

Example line: 1001582.LL3_00013<TAB>101<TAB>1378

=back

=back


=over

=over

=item The database file .len file (mOTU.v1.padded.len) (same as for NCBI mode)

This file contains the fasta entries and their respective lengths, one per line. It is generated automatically and do not need to be changed.

Example line: 1001582.LL3_00013<TAB>1478

=back

=back


=over

=over

=item The database file .rownames and .rownames.uniq files (mOTU.v1.padded.rownames[.uniq]) (same as for NCBI mode)

These files contains the rownames of the computed base and insert coverages. These files are generated automatically and do not need to be modified. I the .uniq file the _X_X corresponds to the region over which the coverages are calculated. This is normally _1_[length of sequence]. This region can be changed in the .coord file.

=back

=back

I<functional mode: Design and requirements for summarizing gene profiles into functional categories>

All that is required to summarize gene profiles into functional profiles, is the additional mapping file: DBNAME.functional.map. Note that this file could be gzipped and saved like: DBNAME.functional.map.gz. MOCAT will first look for the gzipped version. If multiple databases have been mapped against a format example of the name would be: DBNAME.1-3.functional.map.gz.

=over

=over

=item The format of the mapping file should first have a header with a '#' at the beginning of the line (IMPORTANT!):

#gene <TAB> cog <TAB> ko <TAB> module <TAB> pathway

=back

=back

=over

=over

=item Then each line has a gene ID form the database and the corresponding category:

geneID <TAB> cogID <TAB> koID <TAB> moduleID <TAB> pathwayID

=back

=back

=back

=back



=head2 Statistics Options

=over

=over

=item B<C<-sfq|stats_fastqc>>

Runs FastQC on each lane for each sample. The output is described in detail at L<http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/>.

=item B<C<-ss|sample_status>>

This option prints a progress status of each sample into the file <SAMPLE FILE>.status AND summarizes most of the stats files in each sample folder and stores in the files <SAMPLE FILE>.*. The columns represent: if a preset 5' trimming base is defined in the file 'MOCAT.cutoff5prime' (M if set in file, A if not - A means it should be automatically determined), if FastQC has been run on the lanes for the sample, if reads have been: trimmed and filtered, screened against database, screened against fasta file, assembled, assembly revised, genes predicted on an assembly and genes predicted on a revised assembly.

=back

=back

=head2 Additional Options

=over

=over

=item B<C<-x|no_execute> (optional)>

Only create the job files to execute, but don't execute them. Note that it's normal for the program to check for existing files, but not finding them (because no jobs were executed).

=item B<C<-nt|no_temp> (optional)>

If specified, overrides any temp folder settings. All temporary files will be stored in CWD/SAMPLE/temp/.

=item B<C<-cpus> (not recommended)>

MOCAT has been optimized to be run on a cluster of nodes with a larger number of cores (>8). This menas that, without setting the number of available cores for each sample (job), the number of cores listed below will be used in corresponding option. However, if you are using a system with fewer cores that 8, or if you think it will be better to use more cores than 8 when, for example, assembling samples, you can specify this option. Do note, by doing so, other options, such as read_trim_filter or screen, will also occupy the specified number of cores, without actually using them. Number of cores per option:

=over

=item -read_trim_filter     : 3

=item -screen               : 8

=item -screen_fastafile     : 2

=item -assembly             : 8

=item -assembly_revision    : 8

=item -gene_prediction      : 1

=item -filter               : 8

=item -profiling            : 4

=item -fetch_mg             : 4

=item modules               : 8

=back

=back

=back

=head1 Initital setup

=head2 Installing required external software

All steps in MOCAT are dependent on external software. Some of these software you have to download manually, as they require a licence key.

=over

=item Screen Fasta File

If you want to perform a screen using a fasta file (NOT screen against a databse, which uses SOAPaligner2), this step uses Usearch. You can download Usearch from L<http://www.drive5.com/usearch/nonprofit_form.html>. After downloading, extract the Usearch executable into the folder 'MOCAT_PATH/ext/usearch' and rename the file to 'usearch'. Renaming the file is required because you may download a newer version than was initially used in developing MOCAT, and by renaming it to only 'usearch', MOCAT can identify the executable. NOTE: make sure that the file is made 'executable', this can be done by typing (when in the ext/usearch directory) 'chmod u+x usearch'.

=item Gene Prediction

Gene prediction is done by MetaGeneMark, or Prodigal (Default). You can download MetaGeneMark from L<http://exon.gatech.edu/GeneMark/license_download.cgi>. After downloading MetaGeneMark, extract the files into 'MOCAT_PATH/ext/metagenemark'. Normally you also have to copy the 'gm_key' file to your home directory, but MOCAT will do this for you, if needed.

=back

=head2 Required paths in UNIX

Make sure that the MOCAT/src folder is in your UNIX path variable $PATH. This can be done by typing (in the UNIX shell) 'echo $PATH:/usr/me/MOCAT/src > $PATH', assuming you installed MOCAT into the directory /usr/me. By doing this you can execute the MOCAT.pl script from any directory.

Add the MOCAT/src folder to the UNIX variable PERL5LIB. This can be done by typing 'echo $PERL5LIB:/usr/me/MOCAT/src > $PERL5LIB', or if the $PERL5LIB variable is empty: 'echo /usr/me/MOCAT/src > $PERL5LIB', assuming you installed MOCAT to the /usr/me direcotry.

=head2 Setup MOCAT in a new folder

To initiate a new folder structure for MOCAT three steps are required.

=over

=item 1. Copy the config file from the MOCAT/GETTING_STARTED directory into the current direcotry. Make sure all settings in the config are set appropriately. Especially the GLOBAL SETTINGS section.

=item 2. For each sample to be analysed, create a subfolder. Within this subfolder each paired end files should be named: LANE_NAME.1.fq[.gz] and LANE_NAME.2.fq[.gz]. Note that it is possible to have several lanes within one sample.

=item 3. Make a file within the current folder (suggested) where you store the names of the samples to be analysed.

=back

=head2 Example of a folder structure

Current folder is: /usr/me/project. Then the following setup would be appropriate to run MOCAT:

=over

=item ____________________________________

=item Sample files (aoption A):

=item /usr/me/project/SAMPLE_1/lane1.1.fq

=item /usr/me/project/SAMPLE_1/lane1.2.fq

=item /usr/me/project/SAMPLE_1/lane2.1.fq

=item /usr/me/project/SAMPLE_1/lane2.2.fq

=item /usr/me/project/SAMPLE_1/lane3.1.fq

=item /usr/me/project/SAMPLE_1/lane4.2.fq

=item /usr/me/project/SAMPLE_2/laneID.1.fq

=item /usr/me/project/SAMPLE_2/laneID.2.fq

=back

=over


=item ____________________________________

=item Sample files (aoption B):

=item /usr/me/project/SAMPLE_1/lane1.pair.1.fq

=item /usr/me/project/SAMPLE_1/lane1.pair.2.fq

=item /usr/me/project/SAMPLE_1/lane1.single.fq

=item /usr/me/project/SAMPLE_1/lane2.pair.1.fq

=item /usr/me/project/SAMPLE_1/lane2.pair.2.fq

=item /usr/me/project/SAMPLE_1/lane2.single.fq

=back


=over 

=item ____________________________________

=item Required files:

=item /usr/me/project/MOCAT.cfg

=item /usr/me/project/MOCAT.samples

=back

=over

=item ____________________________________

=item Content of MOCAT.samples file:

=item SAMPLE_1

=item SAMPLE_2

=back

=head1 The MOCAT configuration file - MOCAT.cfg

This file contains settings assumed to not be changed equally often as those entered each time when executing the program. The file is divided into different sections, usually one for each possible pipeline option (but not all pipeline options need additional settings).

=head2 Global Settings Section

=over

=item MOCAT_umask        : 0022 [0022]

Sets the UNIX umask options that should be used. If set to 0002, for example, also members of the same group can write to files created by the user. Preferably this is combined with adding 'umask 0002' in the MOCAT_pre_execute option.

=MOCAT_pre_execute       : [umask 0002]

Here you can set a command that is executed before the actual job. Here you can set a specific umask, or perhaps run the command 'newgrp' if you want to run all jobs as a specific group.

=item MOCAT_LSF_qsub_add_param : [-l select=mem=6gb]

Additional parameters for LSF queuing system

=item MOCAT_LSF_queue          : []

Sets the specific LSF '-queue' flag.

=item MOCAT_LSF_memory_limit   : []

Sets the specific LSF '-M' flag

=item qsub_system        : SGE [SGE,PBS,LSF,none]

Currently SGE is supported. If you are using a UNIX cluster without a queuing system set it to 'none'. MOCAT has been tested on SGE queuing system version GE 6.2u4 and GE 6.2u5 and a PBS as well as a LSF queuing system.

=item MOCAT_SGE_qsub_add_param : [-l mem_free=6G -l h_vmem=6G] 

Additional parameters for SGE queuing system

=item MOCAT_PBS_qsub_add_param : [-l select=mem=6gb]

Additional parameters for PBS queuing system

=item MOCAT_SGE_parallell_env : make [make,orte,...]

This options specifies which parallel environment to use. It is set during the installation process but can easily be changed here. If it's incorrect jobs cannot be submitted. What the different options are only the administrator of your system knows.

=item MOCAT_dir          : /bin/MOCAT

the main MOCAT directory, in which the /scr, /bin, /cnf and /data are located.

=item MOCAT_temp_dir_1   : HOSTNAME|/tmp

=item MOCAT_temp_dir_2   : any|/tmp

MOCAT has an advanced system to specify temporary directories. The simplest way, is to not specify any of the tags MOCAT_temp_dir_X, then all temporary files will be stored in CWD/SAMPLE/temp/. The next level is to specify a temp directory for any and all systems you are using, this can be done by specifying 'MOCAT_temp_dir_1 : any_/tmp'. Here, 'any' refers to the hostname of the clusters, by specifying 'any' the temp directory till be the trailing '/tmp' on any system. The '_' is used internally in MOCAT to differentiate between hostname and start of path. If you wish to specify specific temporary directories for specific host computers, this is possible by adding as many MOCAT_temp_dir_1, MOCAT_temp_dir_2 ... MOCAT_temp_dir_N tags as desired. Please note that the tags are evaluated in order. This means if the specified hostname matches the systems's hostname, this temporary directory is used.

=item MOCAT_data_type          : solexaqa [fastx,solexaqa]

Reads can be trimmed and filtered using either FastX 3' trimming or SolexaQA 3' trimming. Very generally, the SolexaQA trimmming is more strict. 

=item MOCAT_paired_end         : yes [yes]

Defines whether sample data is paried end reads or not.

=item MOCAT_zip_program        : pigz [gzip,pigz]

Zip program. PigZ is a parallell version of gzip, provided with MOCAT. pigz is not supported under OSX.

=item MOCAT_default_reads      : reads.processed [-define if desired-,'reads.processed','DATABASE']

By setting this options, you do not have to specify many or the additional 'READS ORIGIN' options. Eg, if you''d do a screen step, you'd normally write MOCAT.pl -sf FILE -sdb DB -sr reads.processed. If this is set to reads.processed, you'd only type: MOCAT.pl -sf FILE -sdb DB.

=item MOCAT_zip_level          : 1 [1-9]

Determines how much the zip files are compressed. 9 is higher but slower. In our tests, fasta and fastq files are compressed very well already using compression level 1 or 2.

=item MOCAT_mapping_mode       : allbest [allbest,random,unique]

Corresponds to SOAPsligner2's 'r' option. Mode 'allbest' is not supported under OSX.

=item MOCAT_prompt_before_run  : no [yes,no]

If set to yes, the settings have to be confirmed each time before running MOCAT

=back

=head2 Read Trim and Filter Section

To remove low quality and short reads, the raw reads are trimmed and filtered. This is done using either the FastX or SolexaQA algorithm and additional 5' trimming.

=over 

=item readtrimfilter_length_cutoff               : 45

Removes reads that are shorter than this cutoff. We have mainly been using 30 or 45.

=item readtrimfilter_qual_cutoff                 : 20

The quality cutoff, at which reads are trimmed. Used internally in the FastX and SolexaQA routines. We recommend 20.

=item readtrimfilter_use_sanger_scale            : auto [yes,no,auto] (if files are in Illumina 1.8+ format, use Sanger scale)

Older Illumina sequences have Illumina quality scale. Later versions have a Sanger quality scale. This would normally be set to 'auto', but can be set to 'yes' or 'no' if you have problems with this for some reason.

=item readtrimfilter_trim_5_prime                : yes [yes,no]

Defines whether the 5' end of each read should be trimmed or not.

=item readtrimfilter_use_precalc_5prime_trimming : no [yes,no]

MOCAT normally does 5' trimming automatically (if 'readtrimfilter_trim_5_prime' set to 'yes'). If you wish to specify at which base the reads should be trimmed for each sample, set this to 'yes'. If set to 'yes', the bases to trim at should be provided in the file 'MOCAT.cutoff5prime' (project folder, same location as MOCAT.cfg). The format of the file is: <SAMPLE> TAB <LANE> TAB <FIRST BASE TO KEEP>. Lane here is the file name up to the '.single', or '.pair'. Note that if this is set to 'no', after running read_trim_filter, MOCAT will produce the file 'MOCAT.cutoff5prime_calculated' in the current folder. If you do not want to perform 5' trimming set 'readtrimfilter_trim_5_prime' to 'no'.

=back

=head2 Screen Section

Custom databases, or provided database 'hg19', can be used to screen reads and exclude these reads from further steps. Here, SOAPAligner2 is used to align reads to the specified database. The reads that match the database are removed.

=over

=item screen_length_cutoff           : 45

After alignment, the SOAPAligner2 output file is filtered using a read length cutoff. If the read is shorter than the cutoff, the alignment will be discarded. 

=item screen_percent_cutoff          : 95

Currently this cutoff should not be below 90. below 93.3%, we cannot guarantee that the criterium is fulfilled. While going below will be unreliable. After alignment, the SOAPAligner2 output file is filtered using an alignment percent identity cutoff. If the % id is lower than the cutoff, the alignment will be discarded.

=item screen_soap_seed_length        : 30

This is the "-l" option in SOAPAligner2. It means that the first 30 bp of a read are taken as a seed. In this seed region, you allow for a maximum of 2 mismatches, which equals 93.3%, see also above. The remainder of the read can have up to x mismatches, specified by the next flag.

=item screen_soap_max_mm             : 10

Maximum number of mismatches allowed on a read, apart from the seed. This will guarantee that reads of 100 bp can be mapped down to 88% identity.

=item screen_soap_cmd                : -M 4 [-M 4]

Set of additional paramters parsed to SOAP. They are provided here, so you can change them - however this is not recommended.

=item screen_save_format             : sam [soap,sam]

Defines the output format of the extracted reads. If set to sam, the soap file is additionally converted into sam format 

=back

=head2 Screen Fasta File Section

Read that have been trimmed and filter (and perhaps screened against a custom database) can be screened against a fasta file or DNA sequences. This is done by BLAST using USearch. 

=over

=item screen_fasta_file_usearch_version       : 6 [5,6]

Specify Usearch version. We recommend Usearch 5, but Usearch v6 is supported.

item screen_fasta_file_usearch_version_5_exe  : usearch (path relative to MOCAT_DIR/ext/usearch/)

Path to executable for Usearch5. These two options have been set if you're running different Usearch versions for different projects, for example.

item screen_fasta_file_usearch_version_6_exe  : usearch (path relative to MOCAT_DIR/ext/usearch/)

Path to executable for Usearch6. These two options have been set if you're running different Usearch versions for different projects, for example.

=item screen_fasta_file_blast_e_value         : 0.00001

If a read\s match has an e-value below this specified value, it is removed.

=item screen_fasta_file_blast_read_min_length : 10

The minimum length of reads to be screened. This applies especially to the sequences in the fasta file. For example, if the fasta file contains very short adapter sequences, this value should be less than, or equal, to the shortest sequence length in the specified fasta file. 

=item screen_fasta_file_additional_usearch_cmd : []

Additional commands that could be sent to Usearch.

=back



=head2 Filter Section

After any screen step, reads are filtered, by removing reads below a certain % ID cut off, a certain length cut off, and removing both paired-end reads if one of the two reads match the database.

=over

=item filter_psort_buffer         : 2G

Maximum memory requierment for psort. If you filter reads mapped against a large database, sorting the reads will take longer time. This time is reduced by increasing the amount of RAM allocated for this step. Change this setting depending on your system, if needed. We sometimes use 50GB of RAM for this step.

=item filter_length_cutoff        : 45

The minimum length of a read

=item filter_percent_cutoff       : 95

The minimum % identity for a read not be removed

=item filter_paired_end_filtering : yes [yes,no]

We recommend setting this to yes. This will remove read A2, if read A1 matches the database. For example if you screen against the hg19 database, and read A1 matches the database. It is likely that also read A2 should be removed, which it will be if this is set to yes.

=item filter_remove_mapping_files : no  [yes.no]

Set to yes, if you want to remove the original mapping files after filtering. This will clear up some disk space.

=back



=head2 Assembly Section

Reads are assembled using any of the supported version of SOAPDenovo. First, if paired end reads, the insert size is calculated using either a pre-assembly, or a pre-mapping (pre-align) step. Then the reads are assembled. 

=over

=item assembly_soap_version                : 1.06 [1.05,1.06]

Specifies which SOAPDenovo version to use. Version 1.06 has several improvements over version 1.05, but version 1.05 is still supported. When running under OSX it is hard coded in the MOCATAssembly.pm script that version 1.06OSXnobamaio is used. This because there are probably library issues when using 1.06 under OSX. The difference between these two versions is that 2 output files are not printed in the nobamaio version. However, these files are not used and deleted when running MOCAT.

=item assembly_calculate_insert_size_using : mapping [assembly,mapping]

If reads are paired end, a more accurate insert size, than the estimate from library preparation, is calculated. This can be done by either an assembly step, or by aligning reads to a database. If 'mapping' is specified a database need to be specified as well (see below). If 'assembly' is specified, the file 'MOCAT.preliminary_insert_sizes' is required. The format of the file is <SAMPLE> TAB <LANE> TAB <insert size>, for each sample and lane, one lane per line. LANE here is the name of the files in the sample folder, up to the .1 and .2, excluding the .1 and .2. This file should be located in the project folder (same location as MOCAT.cfg). Using 'mapping' is faster, but not more accurate.

=item assembly_db_for_calc_insertsize      : 1506MG (used if specified 'mapping' above)

A database, against which reads are aligned to calculate insert sizes, if 'mapping' is specified in the previous setting. The database '1506MG' is provided, which consists of 40 marker genes from 1506 reference genomes.

=item assembly_scaftig_min_length          : 500 [500]

Sets the minimum length of a scaftig to be kept. 

=back

=head2 Assembly Revision Section

Revised assemblies are created by realigning reads to an existing assembly, then correcting for indels and base pair errors.

=over

=item assembly_revision_scaftig_min_length          : 500 [500]

Sets the minimum length of a scaftig to be kept. 

=back

=head2 Gene Prediction

Predicts protein coding genes on assemblies or revised assemblies.

=over

=item gene_prediction_software : Prodigal [MetaGeneMark,Prodigal]

If you want to use Prodigal or MetaGeneMark.

=item gene_prediction_input    : scaftig [scaftig,contig,scafSeq] (if revised assembly, can only be 'scaftig')

Determines whether to predict genes on scaftigs, contigs or scaffolds (scafSeq). If predicting on a revised assembly, only scaftig can be selected.

=item gene_prediction_prodigal_cmd   : -f gff [-f xxx, -none-] (-p, -o, -a, -d, -i already set at runtime in MOCAT)

Additional settings for Prodigal

=back

=head2 Experimental settings

These settings are by default turned off as they are under development and may not function as expected. However turning them on will not affect any result files.

=over

=item realtime_status_use  : no [yes,no]

If set to yes, you will see automatic updates of the status any memory usage of each sample. This will only have affect under the SGE queuing system.

=item realtime_status_timer : 5

A timer in seconds how often the data is updated.

=item realtime_status_log   : no

If sets to yes, print more bug searching information.

=item realtime_status_fix_1 : -hostname-

You can add as many status_fix as you like. Here you would list hostnames that are only on one single machine. It is not needed to set any hosts here, but if you do, the status script will assume that the hosts set here have all jobs submitted to one particular node. Eg, if you have a system with three different nodes and the queuoing system could submit jobs to all these nodes, o not enter the hostname here. But if all jobs from a host will end up on that host, do set that hostname here. It has to do with whether the status script should attempt to SSH into the nodes to get the information about current jobs or not. If a hostname is provided here, no SSH will be attempted.

=item realtime_status_fix_2 : -hostname-

See above.

=back

=head1 Output files

=head2 NOTE: This list of output files may not be exhaustive. We recommend you to investigate the different output files manually to know exactly what youre current version of MOCAT provides.

Each step in the pipeline generates a number of files. Here the main files produced in each step are described. For the example below, we assume that the MOCAT_data_type in the config file is solexaqa, a database named 'hg19' and that we are processing paired end reads. The fasta file screen was performed after a database screen. The output files differ slightly if the data is not paired end (apired end not supported in this evrsion).

Additionally, note, that you can summarize these output files using the option -ss, see above. Then the summary files are written to the files <SAMPLE FILE>.XXX.

=head2 Read Trim Filter

Main folder: reads.processed.solexaqa/

=over

=item LANE.1.fq.gz.qual_stats

=item LANE.2.fq.gz.qual_stats

Internal quality stats files used for trimming.

=item LANE.pair.1.fq.gz

=item LANE.pair.2.fq.gz

=item LANE.single.fq.gz

Contains the processed reads for pair 1, 2, and single reads.

=back

=head2 Screen Database

Main folders: reads.screened.hg19.solexaqa AND reads.extracted.hg19.solexaqa AND reads.mapped.hg19.solexaqa

=over

=item reads.screened.hg19.solexaqa/LANE.screened.hg19.pair.1.fq.gz

=item reads.screened.hg19.solexaqa/LANE.screened.hg19.pair.2.fq.gz

=item reads.screened.hg19.solexaqa/LANE.screened.hg19.single.fq.gz

Contains the screened (kept) reads for pair 1, 2, and single reads.

=item reads.extracted.hg19.solexaqa/LANE.extracted.hg19.pair.1.fq.gz

=item reads.extracted.hg19.solexaqa/LANE.extracted.hg19.pair.2.fq.gz

=item reads.extracted.hg19.solexaqa/LANE.extracted.hg19.single.fq.gz

Contains the extracted reads for pair 1, 2, and single reads.

=item reads.screened.hg19.solexaqa/SAMPLE.all.screened.hg19.aligned.hg19.ids

Contains the identifiers of the reads that matched the database. These reads were removed from the input files and are not in any of the files above.

=item reads.screened.hg19.solexaqa/SAMPLE.all.screened.hg19.aligned.soap

Internal output from SOAPAligner2 containing information about the read that were removed.

=item reads.screened.hg19.solexaqa/SAMPLE.all.screened.hg19.aligned.sam.gz

GZipped SAM formatted file of the SOAP output file mentioned above. These are the reads that were removed/extracted in SAM format.

=back

=head2 Screen Fasta File

Main folder: reads.screened.'fasta file name'.solexaqa AND reads.extracted.'fasta file name'.solexaqa AND reads.mapped.'fasta file name'.solexaqa

=over

=item LANE.screened.'fasta file name'.solexaqa/LANE.screened.'fasta file name'.pair.1.fq.gz

=item LANE.screened.'fasta file name'.solexaqa/LANE.screened.'fasta file name'.pair.1.fq.gz

=item LANE.screened.'fasta file name'.solexaqa/LANE.screened.'fasta file name'.pair.1.fq.gz

Contains the screened (kept) reads for pair 1, 2, and single reads.

=item LANE.extracted.'fasta file name'.solexaqa/LANE.extracted.'fasta file name'.pair.1.fq.gz

=item LANE.extracted.'fasta file name'.solexaqa/LANE.extracted.'fasta file name'.pair.1.fq.gz

=item LANE.extracted.'fasta file name'.solexaqa/LANE.extracted.'fasta file name'.pair.1.fq.gz

Contains the extracted reads for pair 1, 2, and single reads..

=item LANE.1.readsMatchingFastaFile

=item LANE.2.readsMatchingFastaFile

Contains Usearch blast results for each pair (.1, .2). These files are used internally to cerate the .readsToRemove file.

=item SAMPLE.readsToRemove

Contains the read id's of all reads that were removed to create the results files. Note that the reads here are missing the trailing pair identifier (1 or 2). This is because if paired read 1 is filtered, so is also read 2.

=back

=head2 Assembly

Main folder: assembly.hg19.solexaqa.K19/ (19 represents the K-mer size for that assembly)

=over

=item SAMPLE.assembly.hg19.solexaqa.K19.config

This is the config file used in SOAPDenovo to generate the assembly.

=item SAMPLE.assembly.hg19.solexaqa.K19.contig

Contains the contigs produced.

=item SAMPLE.assembly.hg19.solexaqa.K19.scafSeq

Contains the full scafolds.

=item SAMPLE.assembly.hg19.solexaqa.K19.scaftig

Contains all scaftigs longer than specified (normally 500 bp).

=back

=head2 Assembly Revision

Main folder: assembly.revised.hg19.solexaqa.K19/

=over

=item SAMPLE.assembly.revised.hg19.solexaqa.K19.scaftig

Revised scaftigs longer than specified length (normally 500 bp).

=back

=head2 Gene Prediction

Main folder: gene.prediction.assembly.revised.hg19.solexaqa.K19.500.'MetaGeneMark or Prodigal'

Here it is 'assembly.revised' because we used a revised assembly when predicting the genes. The .500 means that genes are predicted on scaftigs that are longer than 500bp, this number is taken from the value in the config file.

=over

=item SAMPLE.gene.prediction.assembly.revised.hg19.solexaqa.K19.faa

Amino acid sequences of predicted proteins.

=item SAMPLE.gene.prediction.assembly.revised.hg19.solexaqa.K19.fna

Nucleotide sequences of predicted genes.

=item SAMPLE.gene.prediction.assembly.revised.hg19.solexaqa.K19.lst

MetaGeneMark output file.

=item SAMPLE.gene.prediction.assembly.revised.hg19.solexaqa.K19.tab

A summary table of all predicted genes. This information is also available in the headers of the predicted genes and proteins.

=back

=head2 Filter

Main folder: reads.filtered.reads.processed.solexaqa

=head2 Profiling

Main folders: [base|insert].coverage.'database'.'solexaqa or fastx'

There are two main folders in each sample folder: base and insert coverages. Inside each of these folders are files with length normalized and raw count data.
The rownames for each database are stored in the data_dir. These files can be concatenated together with the rowname-file in the database directory to generate custom
files with calculated coverages for a combination of samples.

Additionally, if the mode NCBI or mOTU was run, results are saved in the folders taxonomic.profiles.XXX and motu.profiles.XXX

=head2 Fetch Marker Genes

Main folder: gene.fetched.MGs.assembly.revised.hg19.solexaqa.K19.500.'MetaGeneMark or Prodigal'

In this folder the fetched (extracted) marker genes are stored.



=head1 Statistics files

=head2 NOTE: This list of statistics files may not be exhaustive. We recommend you to investigate the different output files manually to know exactly what youre current version of MOCAT provides.

Each step in the pipeline prints different statistics files. These are stoed in the /stats folder in each sample folder.

=over

=item Read Trim Filter

=over

=item LANE.1.fq.gz_raw.reads.stats

=item LANE.2.fq.gz_raw.reads.stats

Number of raw reads and number of raw bases

=item SAMPLE.readtrimfilter.solexaqa.stats

After filtering: number of reads, bases, maximum read length, average read length, estimated K-mer and inserts.

=back

=item Screen Database

=over 

=item SAMPLE.screen.hg19.solexaqa.stats

After screening using database: number of reads, bases, maximum read length, average read length, estimated K-mer and inserts.

=item SAMPLE.extracted.hg19.solexaqa.stats

After screening using database, stats on teh extracted reads: number of reads, bases, maximum read length, average read length, estimated K-mer and inserts.


=back

=item Screen Fasta File

=over

=item SAMPLE.fastafile.hg19.solexaqa.stats

After screening using fasta file: number of reads, bases, maximum read length, average read length and estimated K-mer.

=item LANE.1.screen.fastafile.solexaqa.stats

=item LANE.2.screen.fastafile.solexaqa.stats

Total number of reads screened, total number of hits and total number of unique hits.

=back

=item Assembly

=over

=item SAMPLE.assembly.hg19.solexaqa.K19.assembly.stats

=item SAMPLE.assembly.hg19.solexaqa.K19.scaftig.stats

Specific statistics on each assembly.

=item SAMPLE.assembly.hg19.solexaqa.K19.inserts.stats

Contains the calculated insert sizes for each lane, for that specific sample.

=back

=item Assembly Revision

=over

=item =item SAMPLE.assembly.revised.hg19.solexaqa.K19.baseErrorAndIndelError.stats

Statistics describing how many single base errors, small indels regions were found, and their lengths.

=back

=back


=head1 Additional files used by MOCAT

All these files are stored in the working directory containing the sample folders.

=over

=item MOCAT.cfg

The MOCAT configuration file, described above.

=item MOCAT.samples

Example name of the required sample file used to specify which samples to analyze.

=item MOCAT.cutoff5prime

File required as input if 'readtrimfilter_use_precalc_5prime_trimming' is set to 'yes' on the config file. See config file section.

=item MOCAT.cutoff5prime_calculated

File written by MOCAT, after read trim filter, with the chosen 5' cutoffs for each lane and sample. The number represents the first base to keep.

=item MOCAT.preliminary_insert_sizes

File required if 'assembly_calculate_insert_size_using' is set to 'assembly' in the config file. See config file section.

=back

=head1 The /log directory

The /log directory created contains all temporary output log files created while running MOCAT. These files could be useful for tracking an error. Prior to being loacted in the /log directory, the files are located in the current working directory. All files in the /log directory can safely be removed at any time.

The log folder has a subfolder for each processing step. Eg:

=over

=item assembly

=item assembly_revision

=item profiling

=item filter

=item gene_prediction

=item readtrimfilter

=item screen

=item screen_fasta_file

=item other

Please note that this folder may contain error messages not captured in the other log files. If you run into problems, it could be very fruitful to have a look into this folder!

=item resources

If you run MOCAT using the experimental settings set to 'yes', this folder will contain detailed information about how much memory, disk space and how many processors were used

=back

Each processing folder has in turn these subfolders:

=over

=item commands folder

Here is summarized information including which settings used for each tas performed by MOCAT

=item jobs folder

Here are the specific executed commands for each task performed by MOCAT

=item samples folder

Here are log files for each sample that contains any error messages from any processing step

=item startstop

here is information about when each job started and finished and the exit status

=item arrays

This folder contains the array files used to submit the jobs to the queues.

=back


=head1 Citing MOCAT

=head2 If you have used MOCAT in your work, please cite:

=over

=over

=item Kultima JR, Sunagawa S, Li J, Chen W, Chen H, et al. (2012)

MOCAT: A Metagenomics Assembly and Gene Prediction Toolkit. PLoS ONE 7(10): e47656. doi:10.1371/journal.pone.0047656

=back

=back

=head2 MOCAT is a wrapper for 3rd party software. Therefore we strongly suggest you also cite the following papers if you use MOCAT:

=head3 Initial read trimming and quality control

=over

=over

=item Cox MP, Peterson DA, Biggs PJ (2010)

SolexaQA: At-a-glance quality assessment of Illumina second-generation sequencing data. BMC bioinformatics 11: 485 doi:10.1186/1471-2105-11-485.

=item FastX program:

http://hannonlab.cshl.edu/fastx_toolkit/

=back

=back

=head3 Mapping reads

=over

=over

=item Li R, Yu C, Li Y, Lam T-W, Yiu S-M, et al. (2009)

SOAP2: an improved ultrafast tool for short read alignment. Bioinformatics (Oxford, England) 25: 1966 to 1967 doi:10.1093/bioinformatics/btp336.


=item Edgar RC (2010)

Search and clustering orders of magnitude faster than BLAST. Bioinformatics. 2010 Oct 1;26(19):2460-1. doi: 10.1093/bioinformatics/btq461

=back

=back

=head3 Assembly


=over

=over

=item Li R, Zhu H, Ruan J, Qian W, Fang X, et al. (2010)

De novo assembly of human genomes with massively parallel short read sequencing. Genome research 20: 265 to 272 doi:10.1101/gr.097261.109.

=back

=back

=head3 Assembly revision

=over

=over

=item Li H, Durbin R (2009)

Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics (Oxford, England) 25: 1754 to 1760 doi:10.1093/bioinformatics/btp324.

=back

=back

=head3 Gene Prediciton

=over

=over

=item Hyatt D, Chen G-L, Locascio PF, Land ML, Larimer FW, et al. (2010)

Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC bioinformatics 11: 119 doi:10.1186/1471 to 2105-11-119.

=back

=back

=over

=over

=item Zhu W, Lomsadze A, Borodovsky M (2010)

Ab initio gene identification in metagenomic sequences. Nucleic acids research 38: 1 to 15 doi:10.1093/nar/gkq275.

=back

=back

=head3 Retrieving Marker Genes

=over

=over

=item Sunagawa et al. (2013)

Metagenomic species profiling using universal phylogenetic marker genes. Nature Methods 10, 1196 to 1199 doi:10.1038/nmeth.2693

=back

=back


=head1 Supported Formats

MOCAT supports both older and newer (as of January, 2012) Illumina file formats. The old format has a fasta header descriptor like this: '@HWUSI-EAS100R:6:73:941:1973#0/1'. The new format has a fasta header like this: '@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG'. However, when processing the reads in the first read trim filter step, the headers are converted into the old format. This menas the following information is lost: run id, flowcell id, control bit. Reads that are of low quality, ninth field equals to 'Y', are filtered in the first step. For more information, see L<http://en.wikipedia.org/wiki/FASTQ_format>.

=head1 Author

The B<MOCAT> pipeline was developed by Jens Roat Kultima & Shinichi Sunagawa (Bork Group, EMBL) in collaboration with BGI. External software used by the B<MOCAT> pipeline are copyright respective authors. 

=head1 Copyright

Copyright (c) 2011-2012. Jens Roat Kultima, Shinichi Sunagawa, EMBL and BGI. MOCAT is released under the GNU General Public Licence v3 (http://www.gnu.org/licenses/gpl.html).

=cut
