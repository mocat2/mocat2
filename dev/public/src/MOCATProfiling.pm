package MOCATProfiling;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

sub create_job
{

  #############################################################################################################################################
  # INITIALIZE
  #############################################################################################################################################
  my $job        = shift;
  my $processors = shift;
  my $jobfile    = "$jobdir/MOCATJob_$job\_$date";
  $ZIP =~ s/pigz.*/pigz -p $processors/;
  open JOB, '>', "$jobfile" or die "ERROR & EXIT: Cannot open $jobfile for writing.";
  print localtime() . ": Creating $job jobs...\n";
  my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
  my $print_rownames;
  my $read_type = 'screened';
  my $temp_folder;    # this will be used for saving to the script file below, similar but not identical to temp_dir

  if ($use_extracted_reads)
  {
    $read_type = 'extracted';
  }
  MOCATCore::mkdir_or_die($rownames_dir);
  my $warnings      = 0;    # THIS IS USED FOR CHECKING IF WE HAVE USED OLD FILES
  my $profiling_map = "";
  my $blocksize;
  my $databases = MOCATCore::checkAndReturnDB( \@do_profiling );
  my $databases_old = join( "_AND_", @do_profiling );

  my $columnName;
  my $profiling_NCBI_map       = "$data_dir/$databases.ncbi.map";
  my $profiling_motu_map       = "$data_dir/$databases.motu.map";
  my $profiling_functional_map = "$data_dir/$databases.functional.map";
  #############################################################################################################################################

  #############################################################################################################################################
  # CHECK MODE AND FILES
  #############################################################################################################################################
  print localtime() . ": Checking mode...";

  if (    $databases[0] eq 's'
       || $databases[0] eq 'c'
       || $databases[0] eq 'f'
       || $databases[0] eq 'r' )
  {

    # set mode
    $profiling_mode = 'gene';

  }

  unless ($profiling_mode)
  {
    die "\nERROR & EXIT: Please specify -mode: available modes are 'gene', 'mOTU', 'NCBI', 'functional'";
  }
  if ( $profiling_mode eq 'RefMG' )
  {
    $profiling_mode = 'NCBI';
  }
  if ( $profiling_mode eq 'NCBI' )
  {
    $profiling_map = $profiling_NCBI_map;
    ( -e "$profiling_map" ) or die "\nERROR & EXIT: Missing map file $profiling_map";
    $blocksize = 200000;
  } elsif ( $profiling_mode eq 'mOTU' )
  {
    $profiling_map = $profiling_motu_map;
    ( -e "$profiling_map" ) or die "\nERROR & EXIT: Missing map file $profiling_map";
    $blocksize = 200000;
  } elsif ( $profiling_mode eq 'gene' )
  {
    $blocksize = 2500000;
  } elsif ( $profiling_mode eq 'functional' )
  {
    $profiling_map = $profiling_functional_map;
    ( -e "$profiling_map" || -e "$profiling_map.gz" ) or die "\nERROR & EXIT: Missing map file $profiling_map (or alternatively $profiling_map.gz)";
    $blocksize = 2500000;
  } else
  {
    die "\nERROR & EXIT: -mode $profiling_mode is an incorrect mode! Available modes are 'gene', 'mOTU', 'NCBI', 'functional'";
  }
  print " OK!\n";
  #############################################################################################################################################

### BWA SUPPORT ###
  my $BWA;
  if ( $do_profiling_bwa[0] )
  {
    unless ( $conf{bwa_options} ) { $conf{bwa_options} = "" }
    my $tmp = $conf{bwa_options};
    $tmp =~ s/ //g;
    $tmp =~ s/\t//g;
    $BWA = "bwa$tmp";
  }
### BWA SUPPORT ###

  my $samples_counter = 0;
  foreach my $sample (@samples)
  {
    #############################################################################################################################################
    # INITIALIZE
    #############################################################################################################################################
    my $LOG = " 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log ";
    my ( $output_gene_folder, $rownames, $coord_files, $output_folder, $temp_file, $fullname, $input_folder, $input_folder_old, $input_file, $input_file_old, $output_gene_file, $output_file, $output_file_rownames, $stats_file, $out_stats_file, $out_PE_stats_file );
    my ( $max, $avg, $kmer );
    MOCATCore::mkdir_or_die("$temp_dir/$sample/temp");
    $samples_counter++;
    #############################################################################################################################################

    #############################################################################################################################################
    # SCAFTIG, CONTIG
    #############################################################################################################################################
    if (    $databases[0] eq 's'
         || $databases[0] eq 'c'
         || $databases[0] eq 'f'
         || $databases[0] eq 'r' )
    {
      my $assembly_type = 'assembly';
      my $end;
      if ( $databases[0] eq 's' )
      {
        $end        = 'scaftig';
        $columnName = 'scaftig';
      }
      if ( $databases[0] eq 'c' )
      {
        $end        = 'contig';
        $columnName = 'contig';
      }
      if ( $databases[0] eq 'f' )
      {
        $end        = 'scafSeq';
        $columnName = 'scafSeq';
      }
      if ( $databases[0] eq 'r' )
      {
        $assembly_type = 'assembly.revised';
        $end           = 'scaftig';
        $columnName    = "revised.scaftig";
      }
      ( $max, $avg, $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-r" );

      # Define input file
      $input_folder = "$cwd/$sample/reads.filtered.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}";
      die "$do_profiling_bwa[0]";
      if ( $do_profiling_bwa[0] )
      {
        $input_file = "$input_folder/$sample.filtered.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$BWA.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.z$conf{bwa_min_perc_query_aligned}.bam";
      } else
      {
        $input_file = "$input_folder/$sample.filtered.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.bam";
      }

      my @parts = split '/', $input_file;

      # Define output files

      if ( $do_profiling_bwa[0] )
      {
        $output_gene_folder   = "coverage.$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end";
        $output_gene_file     = "$sample.filtered.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$BWA.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.z$conf{bwa_min_perc_query_aligned}";
        $output_file          = "$sample.$profiling_mode.profile.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$BWA.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.z$conf{bwa_min_perc_query_aligned}";
        $output_file_rownames = "$sample_file_basename.$profiling_mode.profile.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$BWA.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.z$conf{bwa_min_perc_query_aligned}.rownames";

      } else
      {
        $output_gene_folder   = "coverage.$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end";
        $output_gene_file     = "$sample.filtered.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
        $output_file          = "$sample.$profiling_mode.profile.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
        $output_file_rownames = "$sample_file_basename.$profiling_mode.profile.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.rownames";

      }

      # Define other files
      if ( $do_profiling_bwa[0] )
      {
        $out_stats_file = "$cwd/$sample/stats/$sample.coverage.$reads.on.$end.$assembly_type.$conf{MOCAT_data_type}.$BWA.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.z$conf{bwa_min_perc_query_aligned}.stats";
      } else
      {
        $out_stats_file = "$cwd/$sample/stats/$sample.coverage.$reads.on.$end.$assembly_type.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.stats";
      }
      
      # Check len and coord file
      my $assembly_file = "$cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end";
      ( -e "$assembly_file.gz" ) or die "\nERROR & EXIT: Missing $end file: $assembly_file.gz";
      ( -e $input_file ) or die "\nERROR & EXIT: Missing mapping file: $input_file";
      unless ( -e "$assembly_file.len" )
      {
        print localtime() . ": Creating length file $assembly_file.len";
        system "$scr_dir/MOCATFilter_falen.pl -infile $assembly_file.gz -outfile $assembly_file.len -zip ";
        unless ( -e "$assembly_file.coord" )
        {
          print "\n" . localtime() . ": Creating $assembly_file.coord, and then continue...\n";
          if ( $systemType =~ m/Darwin/ )
          {
            system "sed \'s/[[:space:]]/	1	/\' $assembly_file.len > $assembly_file.coord";
          } else
          {
            system "sed \'s/\\t/\\t1\\t/\' $assembly_file.len > $assembly_file.coord";
          }
          unless ( -e "$assembly_file.coord" )
          {
            print localtime() . ": Error creating .coord file\n";
            die "ERROR & EXIT: Could not create $assembly_file.coord from $assembly_file.len\nDoes the system have write permission?";
          }
        }
      }
      $coord_files = "$assembly_file.coord";

      # Set rownames
      $rownames = "$rownames_dir/$parts[-1].rownames";
      if ( $samples_counter == 1 )
      {
        $print_rownames = " -print_rownames_file";
      } else
      {
        $print_rownames = "";
      }
    }    # End scaftig
    #############################################################################################################################################

    #############################################################################################################################################
    # NORMAL DATABASE
    #############################################################################################################################################
    else
    {

      # Create coord files if needed
      my @coord_files;
      foreach my $databases (@databases)
      {
        push @coord_files, "$data_dir/$databases.coord";
        unless ( -e "$data_dir/$databases.coord" || -e "$data_dir/$databases.coord.gz" )
        {

          die "ERROR & EXIT: Previous version of MOCAT created a missing .coord file automatically, but we found this not to always be desired.
The $data_dir/$databases.coord file is missing. Please manually create it. Each database entry should have a single line entry as such:
ID <TAB> START <TAB> STOP.\nIf there is no padded region of the gene, you can create the coord file by running this:
$scr_dir/MOCATFilter_falen.pl -infile $data_dir/$databases -outfile $data_dir/$databases.len && sed 's/\\t/\\t1\\t/\' $data_dir/$databases.len > $data_dir/$databases.coord";

          #					unless ( -e "$data_dir/$databases.len" ) {
          #						unless ( -e "$data_dir/$databases" ) {
          #							die "\nERROR & EXIT: $data_dir/$databases does not exist. Cannot create .len and .coord files";
          #						}
          #						print "\n" . localtime() . ": Creating length file $data_dir/$databases.len...";
          #						system "$scr_dir/MOCATFilter_falen.pl -infile $data_dir/$databases -outfile $data_dir/$databases.len";
          #						print " OK!";
          #					}
          #					print "\n" . localtime() . ": Creating $data_dir/$databases.coord, and then continue...";
          #					if ( $systemType =~ m/Darwin/ ) {
          #						system "sed \'s/[[:space:]]/	1	/\' $data_dir/$databases.len > $data_dir/$databases.coord";
          #					}
          #					else {
          #						system "sed \'s/\\t/\\t1\\t/\' $data_dir/$databases.len > $data_dir/$databases.coord";
          #					}
          #					unless ( -e "$data_dir/$databases.coord" ) {
          #						print localtime() . ": Error creating .coord file\n";
          #						die "ERROR & EXIT: Could not create $data_dir/$databases.coord from $data_dir/$databases.len\nDoes the system have write permission?";
          #					}
        }
      }

      # Define coord files
      $coord_files = join " ", @coord_files;

      # Set rownames
      $rownames = "$rownames_dir/$databases.rownames";
      if ( $samples_counter == 1 )
      {
        $print_rownames = " -print_rownames_file";
      } else
      {
        $print_rownames = "";
      }

      # Bowtie2 and BWA support
      if ( $conf{MOCAT_mapping_mode} eq "bwa" )
      {
        unless ( $conf{bwa_options} ) { $conf{bwa_options} = "" }
        my $tmp = $conf{bwa_options};
        $tmp =~ s/ //g;
        $tmp =~ s/\t//g;
        $conf{MOCAT_mapping_mode} = "bwa$tmp";
      }
      if ( $conf{MOCAT_mapping_mode} eq "bowtie2" )
      {
        unless ( $conf{bowtie2_options} ) { $conf{bowtie2_options} = "" }
        my $tmp = $conf{bowtie2_options};
        $tmp =~ s/ //g;
        $tmp =~ s/\t//g;
        $conf{MOCAT_mapping_mode} = "bowtie2$tmp";
      }

      # Define output files
      $output_gene_folder = "coverage.$databases.$conf{MOCAT_data_type}";

      #$output_folder        = "$profiling_mode.profiles.$databases.$conf{MOCAT_data_type}";
      # original v1.3 # $output_gene_file     = "$sample.filtered.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
      my ( $input, $input_old );
      if ( $do_profiling_bwa[0] )
      {
        $output_gene_file     = "$sample.gene.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$BWA.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
        $output_file          = "$sample.$profiling_mode.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$BWA.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.z$conf{bwa_min_perc_query_aligned}";
        $output_file_rownames = "$sample_file_basename.$profiling_mode.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$BWA.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.z$conf{bwa_min_perc_query_aligned}";
        $input                = "$sample.filtered.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$BWA.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.z$conf{bwa_min_perc_query_aligned}";
        $input_old            = "$sample.filtered.$read_type.$reads.on.$databases_old.$conf{MOCAT_data_type}.$BWA.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.z$conf{bwa_min_perc_query_aligned}";
        $input_folder         = "$cwd/$sample/reads.filtered.$databases.$conf{MOCAT_data_type}";
        $input_file           = "$input_folder/$input.bam";
        $input_file_old       = "$input_folder/$input.bam";
      } else
      {
        $output_gene_file     = "$sample.gene.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
        $output_file          = "$sample.$profiling_mode.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
        $output_file_rownames = "$sample_file_basename.$profiling_mode.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
        $input                = "$sample.filtered.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
        $input_old            = "$sample.filtered.$read_type.$reads.on.$databases_old.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
        $input_folder         = "$cwd/$sample/reads.filtered.$databases.$conf{MOCAT_data_type}";
        $input_folder_old     = "$cwd/$sample/reads.filtered.$databases_old.$conf{MOCAT_data_type}";
        if ( $conf{MOCAT_mapping_mode} eq "allbest" )
        {
          $input_file     = "$input_folder/$input.bam";
          $input_file_old = "$input_folder_old/$input_old.bam";
        } elsif (    $conf{MOCAT_mapping_mode} eq "unique"
                  || $conf{MOCAT_mapping_mode} eq "random" )
        {
          $input_file = "$input_folder/$input.soap.gz";
          $input_file = "$input_folder_old/$input_old.soap.gz";
        } else
        {
          die "ERROR & EXIT: Unrecognized MOCAT mapping mode '$conf{MOCAT_mapping_mode}'";
        }

      }
      $columnName = "";    # set to nothing, because only applies to scaftigs, contigs, ...
    }
    #############################################################################################################################################

    #############################################################################################################################################
    # CONTINUE CREATING JOBS
    #############################################################################################################################################

    # Define other files
    if ( $do_profiling_bwa[0] )
    {
      $out_stats_file = "$cwd/$sample/stats/$sample.coverage.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$BWA.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.z$conf{bwa_min_perc_query_aligned}.stats";
      if ( $reads eq 'reads.processed' )
      {
        $out_PE_stats_file = "$cwd/$sample/stats/$sample.extracted.$databases.after.PE.filter.and.within.padded.region.$BWA.$conf{MOCAT_data_type}.stats";
      } else
      {
        $out_PE_stats_file = "$cwd/$sample/stats/$sample.extracted.screened.$reads.on.$databases.after.PE.filter.and.within.padded.region.$BWA.$conf{MOCAT_data_type}.stats";
      }
    } else
    {
     $out_stats_file = "$cwd/$sample/stats/$sample.coverage.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.stats";
      if ( $reads eq 'reads.processed' )
      {
        $out_PE_stats_file = "$cwd/$sample/stats/$sample.extracted.$databases.after.PE.filter.and.within.padded.region.$conf{MOCAT_mapping_mode}.$conf{MOCAT_data_type}.stats";
      } else
      {
        $out_PE_stats_file = "$cwd/$sample/stats/$sample.extracted.screened.$reads.on.$databases.after.PE.filter.and.within.padded.region.$conf{MOCAT_mapping_mode}.$conf{MOCAT_data_type}.stats";
      }
    }

    # Check input file
    unless ( -e $input_file )
    {
      if ( -e $input_file_old )
      {
        print localtime() . ": WARNING! Could not find new format input file '$input_file', but using old format input file '$input_file_old'\n";
        $input_file = $input_file_old;
      } else
      {
        die "\nERROR & EXIT: Missing filtered mapping results file: $input_file";
      }
    }

    #		# CReate output folders
    #		MOCATCore::mkdir_or_die("$cwd/$sample/base.$output_gene_folder");
    #		MOCATCore::mkdir_or_die("$cwd/$sample/insert.$output_gene_folder");
    #		unless ( $profiling_mode eq 'gene' ) {
    #			MOCATCore::mkdir_or_die("$cwd/$sample/$output_folder");
    #		}
    #
    # Set temp output folder for staoring temp files
    if ($SHM)
    {
      MOCATCore::mkdir_or_die("/dev/shm/$username/MOCAT_temp/$sample/");
      $temp_file   = "/dev/shm/$username/MOCAT_temp/$sample/$sample.$job.$date";
      $temp_folder = "/dev/shm/$username/MOCAT_temp/$sample";
    } else
    {
      $temp_file   = "$temp_dir/TMP.$sample.$job.$date";
      $temp_folder = "$temp_dir/MOCAT.tmp";
      MOCATCore::mkdir_or_die("$temp_folder");
    }

    #my $covfile = "$cwd/$sample/stats/$sample.coverage.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.stats";
    #$covfile = "$cwd/$sample/stats/$sample.coverage.$reads.on.$end.$assembly_type.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.stats";

    # Define input stats file

    if ( $do_profiling_bwa[0] )
    {
      if ($calculateTaxonomy_previous_calc_coverage_stats_file)
      {
        if ( $reads eq 'reads.processed' )
        {
          $stats_file = "$cwd/$sample/stats/$sample.readtrimfilter.after.PE.filter.and.within.padded.region.$BWA.$conf{MOCAT_data_type}.stats";
        } else
        {
          $stats_file = "$cwd/$sample/stats/$sample.$read_type.$reads.after.PE.filter.and.within.padded.region.$BWA.$conf{MOCAT_data_type}.stats";
        }
      } elsif ($calculateTaxonomy_manual_stats_file)
      {
        $stats_file = "$calculateTaxonomy_manual_stats_file.$sample";
      } else
      {
        if ( $reads eq 'reads.processed' )
        {
          $stats_file = "$cwd/$sample/stats/$sample.readtrimfilter.$conf{MOCAT_data_type}.stats";
        } else
        {
          $stats_file = "$cwd/$sample/stats/$sample.$read_type.$reads.$conf{MOCAT_data_type}.stats";
        }
      }
    } else
    {

      if ($calculateTaxonomy_previous_calc_coverage_stats_file)
      {
        if ( $reads eq 'reads.processed' )
        {
          $stats_file = "$cwd/$sample/stats/$sample.readtrimfilter.after.PE.filter.and.within.padded.region.$conf{MOCAT_mapping_mode}.$conf{MOCAT_data_type}.stats";
        } else
        {
          $stats_file = "$cwd/$sample/stats/$sample.$read_type.$reads.after.PE.filter.and.within.padded.region.$conf{MOCAT_mapping_mode}.$conf{MOCAT_data_type}.stats";
        }
      } elsif ($calculateTaxonomy_manual_stats_file)
      {
        $stats_file = "$calculateTaxonomy_manual_stats_file.$sample";
      } else
      {
        if ( $reads eq 'reads.processed' )
        {
          $stats_file = "$cwd/$sample/stats/$sample.readtrimfilter.$conf{MOCAT_data_type}.stats";
        } else
        {
          $stats_file = "$cwd/$sample/stats/$sample.$read_type.$reads.$conf{MOCAT_data_type}.stats";
        }
      }
    }

    ( -e $stats_file ) or die "ERROR & EXIT: Missing input stats file $stats_file";

    if ($profiling_SAM)
    {
      $profiling_SAM = "-sam";
    } else
    {
      $profiling_SAM = "";
    }

    # get length files, let's hope no ones names them stupidly
    if ($VERBOSE)
    {
      $VERBOSE = " -verbose ";
    } else
    {
      $VERBOSE = "";
    }

    my @levels;
    my $print_levels = "";
    if ( $conf{profiling_PRINT_FUNCTIONAL_LEVELS} )
    {
      my $line = $conf{profiling_PRINT_FUNCTIONAL_LEVELS};
      $line =~ s/ //g;
      if ( $line ne 'undef' && $line ne '' )
      {
        @levels = split ",", $line;
      }
    }
    if ( scalar @levels > 0 )
    {
      $print_levels = " -levels " . join( " ", @levels ) . " ";
    }
    if ($CALCULATE_HORIZONTAL_COVERAGE)
    {
      $CALCULATE_HORIZONTAL_COVERAGE = " -horizon ";
    } else
    {
      $CALCULATE_HORIZONTAL_COVERAGE = "";
    }

    # Print job
    print JOB "exec 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log && ";
    print JOB " $scr_dir/MOCATProfiling.pl " . " -output_file $cwd/$sample/$output_file " . " -input_file $input_file " . " -stats_file $stats_file " . " -out_stats_file $out_stats_file " . " -out_PE_stats_file $out_PE_stats_file " . " -mode $profiling_mode " . " -map_file '$profiling_map' " . " -PE_filter $conf{profiling_paired_end_filtering} " . " -samtools_executable $bin_dir/samtools " . " -sample_name $sample " . " -temp_file $temp_file " . " -zip 'gzip -5'" . " -date $date " . " -threads $processors " . " -blocksize $blocksize " . " -coord_files $coord_files " . " -rownames $rownames " . " $print_rownames " . " $profiling_SAM "    # adding this will force format to SAM
      . " -temp_folder $temp_folder " . " -bin '$bin_dir' " . " -firstColumnName '$columnName' " . " -version 'MOCAT$INTERNAL_MOCAT_VERSION:$MOCAT_ID' " . " $print_levels " . " $CALCULATE_HORIZONTAL_COVERAGE " . " $VERBOSE $LOG \n";
    #############################################################################################################################################
  }    # End per sample
  close JOB;
}
1;
