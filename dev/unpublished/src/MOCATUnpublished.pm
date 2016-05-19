package MOCATUnpublished;
use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor;
#use MOCATCalculateCoverage;
use MOCATTaxoProfiling;
use MOCATCore;
use MOCATImport;
use MOCATCluster;
use MOCATAnalyze;
use MOCATVariables;
use MOCATImportOld;
use MOCATRedoScaftig;
use MOCATExtractGenes;
use MOCATExtractGenes_old;
use MOCATGroupByColumn;
use MOCATVirulenceScreen;
use MOCATResistanceScreen;
use MOCATFunctionalProfiling;
use MOCATCalcMultipleMappers;
use MOCATRunCarmaAnnotation;
use MOCATAnnotateUsingCarmaFile;
use MOCATSummarizeCarmaCoverages;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

sub Usage {
	print "$sep\n";
	print color 'bold';
	print "\n\n                              MOCAT - INTERNAL DEVELOPER VERSION                     \n\n";

	print " _____________________________________________________________________________________________________\n";
	print " | CARMA TAXONOMIC ANNOTATION PIPELINE                                                               |\n";
	print " |                                                                                                   |\n";
	print " |  -create_carma_annotation_file [FILE] -host1 [hostname] -host2 [hostname]                         |\n";
	print " |     host1: a host for running the blast, default is epsilon                                       |\n";
	print " |     host2: a host for running the carma annotations, default is eta                               |\n";
	print " |     STEP 1/3: Runs carma and annotates a blast file.                                              |\n";
	print " |     IMPORTANT! INPUT HAS TO BE PROTEIN SEQUENCE!                                                  |\n";
	print " |     Input file is a set of sequences, example a file with clustered genes.                        |\n";
	print " |     The output is a file with annotations for each input sequence.                                |\n";
	print " |     Note, this step requires a large blast against NR and then running carma.                     |\n";
	print " |     Requires 8 cpus + 24G memory for blast and 2 cpus + 16G memory for the carma                  |\n";
	print " |     Results saved in CWD/annotate.FILE/                                                           |\n";
	print " |                                                                                                   |\n";
	print " |  -annotate_using_carma_file [DB1 DB2 ...] -r [READS ORIGIN] -file [FILE] [-e]                     |\n";
	print " |     STEP 2/3: Uses carma file from step 1 and annotates calculated coverages.                     |\n";
	print " |     On samples which have had insert and base coverages calculates you can                        |\n";
	print " |     summarize these into different taxonomical levels using a carma annotation file.              |\n";
	print " |     This files is created in step 1 above. Please specify same filename is this step.             |\n";
	print " |     Results saved in CWD/SAMPLE/taxonomic.profile/                                                |\n";
	print " |                                                                                                   |\n";
	print " |  -summarize_carma_coverages [DB1 DB2 ...] -r [READS ORIGIN] -file [FILE] [-e]                     |\n";
	print " |     STEP 3/3: Summarizes the .tab files created in step 2 for all the samples                     |\n";
	print " |     This will generate many taxa by sample tables for all the conditions                          |\n";
	print " |     Results saved in CWD/taxonomic.profiles/carma.annotations/FILE/                               |\n";
	print " _____________________________________________________________________________________________________\n\n";

	print " _____________________________________________________________________________________________________\n";
	print " | CALCULATE COVERAGE (additional options not viewed in public version)                              |\n";
	print " |                                                                                                   |\n";
	print " |       Additional options:                                                                         |\n";
	print " |        -mode ['none', 'RefMG', 'mOTU', 'functional', 'identifier']                                |\n";
	print " |           'identifier' is used when we want to summarize on the identifier (gene) level           |\n";
	print " |           'functional' is used to generate functional profiles                                    |\n";
	print " |        -manual_stats_file 'PATH' (PATH should be on format /g/bork/.../something.SAMPLE.stats)    |\n";
	print " |           This overrides the original stats file. Useful when you want to use a manually          |\n";
	print " |           calculated +1                                                                           |\n";
	print " |        -no_paste                                                                                  |\n";
	print " |           If specified, the coverage and taxonomy files are not pasted                            |\n";
	print " |        -only_paste                                                                                |\n";
	print " |           If specified, the coverage and taxonomy files are only pasted. This assumes that the    |\n";
	print " |           calculate_taxonomy step has been run previously for all samples (with or w/o no_paste)  |\n";
	print " _____________________________________________________________________________________________________\n\n";

	#	print " _____________________________________________________________________________________________________\n";
	#	print " | FUNCTIONAL PROFILING PIPELINE                                                                     |\n";
	#	print " |                                                                                                   |\n";
	#	print " |  -functional_profiling [DB1 ...] -r [READS] [-eggnog_map FILE] [-kegg_map FILE] [-e] [-temp]      |\n";
	#	print " |                                                                                                   |\n";
	#	print " |     eggnog_map and kegg_map can be specified in the config file                                   |\n";
	#	print " |     These two files will be project specific and should have been generated by Smash's            |\n";
	#	print " |     doOrthologousMapping script. This Smash script takes as input a blast file, and               |\n";
	#	print " |     output is a mapping file from eg. gene id to EggNOG ID and KEGG IDs                           |\n";
	#	print " |     Note that this step will require at least up to 20-30GB of memory (90 samples and 10M         |\n";
	#	print " |     genes requires 16GB of temporary storage)                                                     |\n";
	#	print " |     Temporary files are stored in /dev/shm/user to speed up processing,                           |\n";
	#	print " |     unless -temp is specified, then large temporary files will be stored in temp folder           |\n";
	#	print " _____________________________________________________________________________________________________\n\n";

	print " _____________________________________________________________________________________________________\n";
	print " | GBC (GROUP BY COLUMN)                                                                             |\n";
	print " |                                                                                                   |\n";
	print " |  This can be used for generating taxonomic or functional profiles, but because it's a more        |\n";
	print " |  general function than taxonomic and functional profiling pipeline, it has its own function.      |\n";
	print " |                                                                                                   |\n";
	print " |  -gbc [DB1 DB2 ...] -r [READS ORIGIN] -map [MAP FILE] [-e] [-temp]                                |\n";
	print " |     The MAP FILE, should be TAB DELIMITED and on the format:                                      |\n";
	print " |     ROW 1: header row. Eg: gene\\tkingdom\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies            |\n";
	print " |     Other rows: The values to be summarized over                                                  |\n";
	print " |     COLUMN 1: These names should correspond to the rownames of the coverage files                 |\n";
	print " |     Temporary files are stored in /dev/shm/user to speed up processing,                           |\n";
	print " |     unless -temp is specified, then large temporary files will be stored in temp folder           |\n";
	print " _____________________________________________________________________________________________________\n\n";

	print " _____________________________________________________________________________________________________\n";
	print " | SCREEN SAMPLES FOR ANTIBIOTICS RESISTANCE or VIRULENCE                                            |\n";
	print " |                                                                                                   |\n";
	print " |  -resistance_screen [DB1 ...] -r [READS] -r2 [READS] -identity2 [%ID] [-e2] -refgenecatalog [DB]  |\n";

	#print " |    -suffix ['', efflux, target, modifier]                                                         |\n"; RUNNING all suffixes at a time currently
	print " |    Calculates the resistance potetional for different antibiotics for each sample                 |\n";
	print " |  -virulence_screen  [DB1 ...] -r [READS] -r2 [READS] -identity2 [%ID] [-e2] -refgenecatalog [DB]  |\n";
	print " |    Calculates virulence factors                                                                   |\n";
	print " |                                                                                                   |\n";
	print " |    -r:  reads used to map against DB (probably RefMGv9) for taxonomical classification            |\n";
	print " |    -r2: reads used to map against 263RefGeneCatalog (NOTE, only 263RefGeneCatalog supported)      |\n";
	print " |    -refgenecatalog: name of reference gene catalog, usually 263RefGeneCatalog.padded or similar   |\n";
	print " |                                                                                                   |\n";
	print " |    Assume you have run RTF, adapter & hg19 screen, and mapped these reads against                 |\n";
	print " |    263RefGeneCatalog, as well as mapped reads first against 263MetaRefv9MG and then mapped the    |\n";
	print " |    extracted reads from this mapping against RefMGv9 (to be able to estimate taxonomical          |\n";
	print " |    abundance), THEN, you can run the following to calculate the resistance potential              |\n";
	print " |    1. Summarize taxonomical profiles for the mappings against RefMGv9: (AN EXAMPLE ONLY)          |\n";
	print " |       mc -sf samples -taxo_profiling RefMGv9.padded                                               |\n";
	print " |          -r screened.screened.adapter.on.hg19.on.263MetaRefv9MG.nr.padded -mode RefMG -e          |\n";
	print " |    2. Summarize the calculated coverages of the mappings against 263RefGeneCatalog:               |\n";
	print " |       mc -sf samples -pcf 263RefGeneCatalog.padded -r screened.adapter.on.hg19                    |\n";
	print " |    3. Calculate the resistance potential / virulence:                                             |\n";
	print " |       mc -sf samples -resistance_screen Ref10MGv9.padded                                          |\n";
	print " |          -r screened.screened.adapter.on.hg19.on.263MetaRef10MGv9.cal.v2.nr.padded -e             |\n";
	print " |          -r2 screened.adapter.on.hg19 -memory 100G -identity 97 -identity2 95                     |\n";
	print " |          -refgenecatalog 263RefGeneCatalog.resistance_virulence_subset.padded                     |\n";
	print " |   NOTE 1: 7 ResistanceScreenData.* files are required in your MOCAT data folder                   |\n";
	print " |   NOTE 2: 7 screen_virulence.* files are required in your MOCAT data folder                       |\n";
	print " _____________________________________________________________________________________________________\n\n";

	print " _____________________________________________________________________________________________________\n";
	print " | ANALYZE                                                                                           |\n";
	print " |                                                                                                   |\n";
	print " |  -analyze [DB1 DB2 ...] -r [READS ORIGIN] [-e]                                                    |\n";
	print " |     Statistics on selected samples using R                                                        |\n";
	print " |     For example you can run DESeq or EdgeR to identify state specific genes/species/...           |\n";
	print " _____________________________________________________________________________________________________\n\n";

	print " _____________________________________________________________________________________________________\n";
	print " | OTHER                                                                                             |\n";
	print " |                                                                                                   |\n";
	print " |  -cluster_genes_by_sample [ASSEMBLY ORIGIN] -r [READS ORIGIN]                                     |\n";
	print " |     Clusters all the genes, sample by sample                                                      |\n";
	print " |                                                                                                   |\n";
	print " |  -calculate_multiple_mappers [DB1 DB2 ...] -r [READS ORIGIN] -summarize_level [LEVEL] -map [FILE] |\n";
	print " |     Calculates the fraction of multiple mappers for mapped reads                                  |\n";
	print " |     summarize_level can be one of 'none', 'taxid', 'kingdom', 'phylum', 'class' 'order',          |\n";
	print " |     'family', 'genus', 'species', 'specI_clusters'                                               |\n";
	print " |     If summarize_level set to anything else but 'none', the identifers mapped to has to be on     |\n";
	print " |     the format TAXID.*                                                                            |\n";
	print " |     Also, if not set to 'none', taxo_profiling_map has to be defined in the config file, or       |\n";
	print " |     specified with -map [FILE]                                                                    |\n";
	print " _____________________________________________________________________________________________________\n\n";

	print " _____________________________________________________________________________________________________\n";
	print " | ADDITIONAL                                                                                        |\n";
	print " |                                                                                                   |\n";
	print " |  -memory [required memory]                                                                        |\n";
	print " |     This will pass the specified memory requirements to the queue on the format -l h_vmem=[X]     |\n";
	print " |     If this is already specified in the SGE_queue paramter in the config file, it will be sent    |\n";
	print " |     twice, probably not a good thing. Also, only implemented for SGE systems.                     |\n";
	print " |  -tc [maximum number of jobs running at the same time]                                            |\n";
	print " |     This will pass the specified memory requirements to the queue on the format -tc [X]           |\n";
	print " |     If this is already specified in the SGE_queue paramter in the config file, it will be sent    |\n";
	print " |     twice, probably not a good thing. Also, only implemented for SGE systems.                     |\n";
	print " |  -noqueue                                                                                         |\n";
	print " |     Sets MOCAT_qsub_system to 'none'                                                              |\n";
	print " |  -filter_percent_cutoff                                                                           |\n";
	print " |     overrides the setting in the config file                                                      |\n";
	print " |  -config A=b AA=bb ...                                                                            |\n";
	print " |     overrides settings form the config file. Eg.: -config MOCAT_qsub_system=none                  |\n";
	print " |  -only_regenerate_reads                                                                           |\n";
	print " |     when running s|screen, no mapping is performed, but the FastQ files are re-generated form the |\n";
	print " |     mapping files directly. This requires the mapping to have been done previously                |\n";
	print " _____________________________________________________________________________________________________\n\n";
	print color 'reset';

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

sub Launch {
	### Check config ###
	my @configs = qw(
	  MOCAT_DEV_DIR
	  cluster_bin
	  cluster_cmd
	  rca_carma_cfg
	  rca_parallel_blast
	  rca_submit_jobs
	  rca_smash_config
	  rca_nr_db_path
	  rca_nr_copy_temp_path
	  analyze_base_and_insert
	  analyze_normalization
	  analyze_taxonomic_level
	  analyze_metadata
	  analyze_feature_annotations
	  analyze_itol_ncbi_names_dmp
	  MOCAT_record_script_versions
	);

	foreach my $c (@configs) {
		unless ( exists $conf{$c} ) {
			die "ERROR & EXIT: Missing setting '$c' in config file. Is the config file correct?";
		}
	}

	### CPU ###
	my $cpu_calculate_coverage        = 1;
	my $cpu_cluster                   = 8;
	my $cpu_calculate_multiplemappers = 2;
	my $cpu_annotate_using_carma_file = 1;
	my $cpu_run_carma_annotation      = 8;
	my $cpu_summarize_carma_coverages = 1;
	my $cpu_analyze                   = 1;
	my $cpu_functional_profiling      = 1;
	my $cpu_resistance_screen         = 1;
	my $cpu_virulence_screen          = 1;
	my $cpu_group_by_column           = 1;
	my $cpu_taxo_profiling            = 1;
	my $cpu_calculate_taxonomy        = 1;
	my $cpu_paste_coverage_files      = 1;
	if ( $cpu > 0 ) {
		$cpu_calculate_coverage        = $cpu;
		$cpu_cluster                   = $cpu;
		$cpu_calculate_multiplemappers = $cpu;
		$cpu_annotate_using_carma_file = $cpu;
		$cpu_run_carma_annotation      = $cpu;
		$cpu_summarize_carma_coverages = $cpu;
		$cpu_analyze                   = $cpu;
		$cpu_functional_profiling      = $cpu;
		$cpu_resistance_screen         = $cpu;
		$cpu_virulence_screen          = $cpu;
		$cpu_group_by_column           = $cpu;
		$cpu_taxo_profiling            = $cpu;
		$cpu_calculate_taxonomy        = $cpu;
		$cpu_paste_coverage_files = $cpu;

	}

	if ( defined $do_cluster_genes ) {
		$assembly = $do_cluster_genes;
		if ( $assembly eq "" || !( $assembly eq 'assembly' || $assembly eq 'assembly.revised' ) ) {
			die "ERROR & EXIT: Input must be either 'assembly' or 'assembly.revised'";
		}
		MOCATCore::print_settings("cluster");
		print localtime() . ": PERFORMING CLUSTERING...\n";
		MOCATCore::pre_check( $assembly, $reads, "cl" );
		MOCATCluster::create_job( 'cluster_genes', $cpu_cluster, $date );
		MOCATCore::execute_job( "cluster_genes", $cpu_cluster, scalar @samples, "15gb" );
		MOCATCluster::post_check_files();
		print localtime() . ": COMPLETED CLUSTERING.\n";
		$print_final_text = 1;
	}
	if ($do_import_old) {
		MOCATCore::print_settings("import_old");
		print localtime() . ": IMPORTING OLD DATA...\n";
		MOCATImportOld::run( 'import_old', 1, $importOldScafLength );
		print localtime() . ": COMPLETED IMPORTING.\n";
		$print_final_text = 1;
	}
	if ($do_import) {
		MOCATCore::print_settings("import");
		print localtime() . ": IMPORTING DATA...\n";
		MOCATImport::run( 'import', 1 );
		print localtime() . ": COMPLETED IMPORTING.\n";
		$print_final_text = 1;
	}
	if ($do_redo_scaf) {
		MOCATCore::print_settings("redo_scaftig");
		print localtime() . ": RE-EXTRACTING SCAFSEQ INTO SCAFTIG...\n";
		MOCATCore::pre_check( $assembly, $reads, "cl" );
		MOCATRedoScaftig::create_job( 'redo_scaftig', 1, $date, $importOldScafLength );
		MOCATCore::execute_job( "redo_scaftig", 1, scalar @samples, "15gb" );
		print localtime() . ": COMPLETED RE-EXTRACTION.\n";
		$print_final_text = 1;
	}
	if ( $do_calculate_multiplemappers[0] ) {
		@databases = @do_calculate_multiplemappers;
		MOCATCore::print_settings("calculate_multiple_mappers");
		print localtime() . ": PERFORMING MULTIPLE MAPPERS CALCULATIONS...\n";
		MOCATCore::pre_check( "", $reads, "calculate_multiple_mappers" );
		MOCATCalcMultipleMappers::create_job( 'calculate_multiple_mappers', $cpu_calculate_multiplemappers, $date );
		MOCATCore::execute_job( "calculate_multiple_mappers", $cpu_calculate_multiplemappers, scalar @samples, "15gb" );
		MOCATCalcMultipleMappers::post_check_files("calculate_multiple_mappers");
		print localtime() . ": COMPLETED CALCULATING MULTIPLE MAPPERS.\n";
		$print_final_text = 1;
	}
	if ($do_run_carma_annotation) {
		MOCATCore::print_settings("rca");
		$carma_file_to_annotate = $do_run_carma_annotation;
		print localtime() . ": PERFORMING BLAST AND CARMA ANNOTATION...\n";
		MOCATRunCarmaAnnotation::create_job( 'run_carma_annotation', $cpu_run_carma_annotation, $date );
		print localtime() . ": NOTE! Running this will take some time as it's first a blast and then carma3.\n";
		MOCATCore::execute_job( "run_carma_annotation", 1, 1, "15gb" );
		print localtime() . ": COMPLETED BLAST AND CARMA ANNOTATION.\n";
		$print_final_text = 1;
	}
	if ( $do_annotate_using_carma_file[0] ) {
		MOCATCore::print_settings("calculate_taxonomic_composition");
		print localtime() . ": PERFORMING TAXONOMIC COMPOSITION CALCULATIONS USING CARMA FILE...\n";
		MOCATCore::pre_check( "", $reads, "r" );
		MOCATAnnotateUsingCarmaFile::create_job( 'annotate_using_carma_file', $cpu_annotate_using_carma_file, $date );
		MOCATCore::execute_job( "annotate_using_carma_file", $cpu_annotate_using_carma_file, scalar @samples, "15gb" );
		print localtime() . ": COMPLETED CALCULATING TAXONOMIC COMPOSITION USING CARMA FILE.\n";
		$print_final_text = 1;
	}
	if ( $do_summarize_carma_coverages[0] ) {
		MOCATCore::print_settings("calculate_taxonomic_composition");
		print localtime() . ": SUMMARIZING CARMA COVERAGES...\n";
		MOCATCore::pre_check( "", $reads, "r" );
		MOCATSummarizeCarmaCoverages::run( 'summarize_carma_coverages', $cpu_summarize_carma_coverages );
		print localtime() . ": COMPLETED SUMMARIZING CARMA COVERAGES.\n";
		$print_final_text = 1;
	}

	if ( $do_functional_profiling[0] ) {
		MOCATCore::print_settings("functional_profiling");
		print localtime() . ": PERFORMING FUNCTIONAL PROFILING...\n";
		MOCATCore::pre_check( "", $reads, "r" );
		MOCATFunctionalProfiling::create_job( 'functional_profiling', $cpu_functional_profiling, $date );
		print localtime() . ": COMPLETED FUNCTIONAL PROFILING.\n";
		$print_final_text = 1;
	}
	if ( $do_resistance_screen[0] ) {
		MOCATCore::print_settings("resistance_screen");
		print localtime() . ": PERFORMING RESISTANCE SCREEN...\n";
		MOCATCore::pre_check( "", $reads, "r" );
		MOCATResistanceScreen::create_job( 'resistance_screen', $cpu_resistance_screen, $date );
		MOCATCore::execute_job( "resistance_screen", $cpu_resistance_screen, 1, "15gb" );
		print "RESULTS SAVED IN $output_file\n";
		print localtime() . ": COMPLETED RESISTANCE SCREEN.\n";
		$print_final_text = 1;
	}
	if ( $do_virulence_screen[0] ) {
		MOCATCore::print_settings("virulence_screen");
		print localtime() . ": PERFORMING VIRULENCE SCREEN...\n";
		MOCATCore::pre_check( "", $reads, "r" );
		MOCATVirulenceScreen::create_job( 'virulence_screen', $cpu_virulence_screen, $date );
		MOCATCore::execute_job( "virulence_screen", $cpu_virulence_screen, 1, "15gb" );
		print "RESULTS SAVED IN $output_file\n";
		print localtime() . ": COMPLETED VIRULENCE SCREEN.\n";
		$print_final_text = 1;
	}
	if ( $do_calculate_coverage[0] ) {
		@databases = @do_calculate_coverage;
		MOCATCore::print_settings("calculate_coverage");
		print localtime() . ": PERFORMING COVERAGE CALCULATIONS...\n";
		MOCATCore::pre_check( "", $reads, "c" );
		MOCATCalculateCoverage::create_job( 'calculate_coverage', $cpu_calculate_coverage );
		MOCATCore::execute_job( "calculate_coverage", $cpu_calculate_coverage, scalar @samples, "15gb" );
		print localtime() . ": COMPLETED CALCULATING COVERAGE.\n";
		$print_final_text = 1;
	}

	if ( $do_taxo_profiling[0] ) {
		MOCATCore::print_settings("taxo_profiling");
		print localtime() . ": PERFORMING TAXONOMIC PROFILING...\n";
		MOCATCore::pre_check( "", $reads, "r" );
		MOCATTaxoProfiling::create_job( 'taxo_profiling', $cpu_taxo_profiling, $date );
		MOCATCore::execute_job( "taxo_profiling", $cpu_taxo_profiling, scalar @samples, "15gb" );
		MOCATTaxoProfiling::summarize("taxo_profiling");
		print localtime() . ": COMPLETED TAXONOMIC PROFILING.\n";
		$print_final_text = 1;
	}

	if ( $do_analyze[0] || $analyze_file ) {
		MOCATCore::print_settings("analyze");
		print localtime() . ": ANALYSING SAMPLES...\n";
		unless ($analyze_file) {
			MOCATCore::pre_check( "", $reads, "r" );
		}
		MOCATAnalyze::run( 'analyze', $cpu_analyze );
		print localtime() . ": COMPLETED ANALYSIS\n";
		$print_final_text = 1;
	}
	if ( $do_group_by_column[0] ) {
		MOCATCore::print_settings("group_by_column");
		print localtime() . ": PERFORMING GROUP BY COLUMN:\n";
		MOCATGroupByColumn::run();
		print localtime() . ": COMPLETED GROUP BY COLUMN.\n";
		$print_final_text = 1;
	}
	if ($do_extract_genes) {
		print localtime() . ": EXTRACTING +- 100bp OF GENES...\n";
		MOCATCore::pre_check( "", $reads, "r" );
		MOCATExtractGenes_old::run();
		print localtime() . ": COMPLETED EXTRACTION.\n";
		$print_final_text = 1;
	}

	### v2 PROFILING (calculate_taxonomy) ###
	if ( $do_calculate_taxonomy[0] ) {
		$taxo_profiling_mode = $profiling_mode;
		@databases = @do_calculate_taxonomy;    # This line is used in both CalculateTaxonomy and PasteTaxonomyCoverageFiles
		unless ($only_paste) {
			MOCATCore::print_settings("profiling");
			print localtime() . ": PERFORMING PROFILING...\n";
			MOCATCore::pre_check( "", $reads, "c" );
			MOCATCalculateTaxonomy::create_job( 'profiling', $cpu_calculate_taxonomy );
			MOCATCore::execute_job( "profiling", $cpu_calculate_taxonomy, scalar @samples, "15gb" );
		}
		unless ( $no_paste || $only_print ) {
			unless ( $taxo_profiling_mode eq 'none' || $taxo_profiling_mode eq '' ) {
				MOCATPasteTaxonomyCoverageFiles::run( 'paste_taxonomy_files', $cpu_calculate_taxonomy );
			}
			MOCATPasteCoverageFiles::paste( 'paste_coverage_files', $cpu_calculate_taxonomy );
		}
		print localtime() . ": COMPLETED PROFILING.\n";
		print "$sep\nGENERATED PROFILES ARE SAVED IN: $OUTPUT_FOLDER\n";
	}

	#### v1 PASTE COVERAGE FILES ###
	if ( $do_paste_coverage_files[0] ) {
		@databases = @do_paste_coverage_files;
		MOCATCore::print_settings("paste_coverage_files");
		print localtime() . ": PASTING COVERAGE FILES:\n";
		MOCATPasteCoverageFiles::paste( 'paste_coverage_files', $cpu_paste_coverage_files );
		print localtime() . ": COMPLETED PASTING COVERAGE FILES.\n";
		$print_final_text = 1;
	}

}
1;
