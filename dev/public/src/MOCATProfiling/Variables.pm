package MOCATProfiling::Variables;
use strict;
use warnings;
use Exporter;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2014
# This code is released under GNU GPL v3.

our @ISA = 'Exporter';
our (
	$date, %NCBI_TAXAID_LENGTH, $CALCULATE_HORIZONTAL_COVERAGE, @CONFIG_LEVELS, $TOTAL_INSERTS_MAPPED, $VERBOSE, $version, @length_files, $NCBI_concated_genelengths, $firstColumnName, $bin_dir, $rownames_file_base, $output_file_base, $sam, $cwd, $HASH, @LEVELS, $print_rownames_file, $temp_folder, $profiling_v2_counting, $db_average_entry_length, $out_MM_stats_file,          $out_PE_stats_file,                       $out_stats_file,      %MAP,             $NCBI_length_file, $request_q,       @input_files,                                    $stats_file,
	$input_file_format,    @coord_files,             $mode,                       @levels,                                  $map_file,            @fields,          @STATS_FILE_DATA,  $output_file,     $rownames_file,                                  $sample_name,
	$samtools_executable,  $PE_filter,               %PEaffectedInserts,          %multipleMapper,                          $zip_execute,         $TOTAL_LENGTH,    %TOTAL_MAPPED2,    %INSERT_COUNTER,  %TOTAL_MAPPED3,                                  %LENGTHS,
	$ALL_TAXA_MEAN_LENGTH, $refmg_length_file,       %TOTAL_MAPPED_NOT_ANNOTATED, %TOTAL_MAPPED_LENGTH_NORM,                $HASH_pointer,        %FRACTION_MAPPED, $temp_file,        %open_temp_files, %TOTAL_MAPPED_ANNOTATED,                         $insert,
	$read,                 $prev_read,               %hash,                       %seen_types,                              $direction,           $ref_id,          $first_base,       $last_base,       $length,                                         %HASH,
	%HASH_LEVELS,          %TOTAL_MAPPED,            %TOTAL,                      $reference_for_read_exists_in_coord_file, @argv,                $QUEUE,           $FINISHED_THREADS, $return_q,        $reference_for_read_do_not_exists_in_coord_file, $THREADS,
	@threadPipes,          @threadHandles,           %gene2other,                              $USE_OLD_HEADERS,     $BLOCKSIZE,       $LOADED_MAP,       $READ,            $main_q,                                         %HASH_local,
	%HASH_LEVELS_local,    $TERMINATE,                                        $TOTAL_LENGTH_MAPPED, $INSERT_COUNTER,  $INIT_NUM,         $prev_insert,     $seen_this_insert,                               $num_of_threads
);
our @EXPORT =
  qw/$date %NCBI_TAXAID_LENGTH $CALCULATE_HORIZONTAL_COVERAGE @CONFIG_LEVELS $VERBOSE $TOTAL_INSERTS_MAPPED @length_files $version $NCBI_concated_genelengths $sam $firstColumnName $bin_dir $rownames_file_base $output_file_base $cwd $HASH @LEVELS $print_rownames_file $temp_folder $profiling_v2_counting $db_average_entry_length $out_MM_stats_file $out_PE_stats_file $out_stats_file $TERMINATE %MAP $NCBI_length_file %HASH_local %gene2other $USE_OLD_HEADERS $BLOCKSIZE $READ $LOADED_MAP $main_q %HASH_LEVELS_local $QUEUE $request_q $THREADS @threadHandles @threadPipes @argv @input_files $stats_file $input_file_format @coord_files $mode @levels $map_file @fields @STATS_FILE_DATA $output_file $rownames_file $sample_name $samtools_executable      $PE_filter      %PEaffectedInserts %multipleMapper $zip_execute $FINISHED_THREADS $return_q         $TOTAL_LENGTH %TOTAL_MAPPED2            %INSERT_COUNTER %TOTAL_MAPPED3     %LENGTHS        $ALL_TAXA_MEAN_LENGTH $refmg_length_file %TOTAL_MAPPED_NOT_ANNOTATED %TOTAL_MAPPED_LENGTH_NORM $HASH_pointer   %FRACTION_MAPPED   $temp_file      %open_temp_files      %TOTAL_MAPPED_ANNOTATED $insert    $read   $prev_read  %hash      %seen_types $direction $ref_id $first_base $last_base $length %HASH %HASH_LEVELS %TOTAL_MAPPED %TOTAL $reference_for_read_exists_in_coord_file $reference_for_read_do_not_exists_in_coord_file $base_output_folder                                                                             $TOTAL_LENGTH_MAPPED                           $INSERT_COUNTER                                $INIT_NUM                                      $prev_insert                                   $seen_this_insert $num_of_threads/;

1;
