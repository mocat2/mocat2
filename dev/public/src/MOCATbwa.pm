package MOCATbwa;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;
use File::Basename;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

sub create_job
{

  # Set variables
  my $job        = $_[0];
  my $processors = $_[1];
  $ZIP =~ s/pigz.*/pigz -p $processors/;
  open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?";

  print localtime() . ": Creating $job jobs...\n";
  my ( $max, $avg, $kmer );
  my $screen_source;
  my $screen_save;
  my $basename;
  my $temp_dir      = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
  my $read_type     = 'screened';
  my $assembly_type = 'assembly';
  my $end;
  my $processors2 = $processors - 2;

  if ( $processors2 < 1 )
  {
    $processors2 = 1;
  }

  my $OSX = "";
  if ( $systemType =~ m/Darwin/ )
  {
    $OSX = "_OSX";
  }
  if ($use_extracted_reads)
  {
    $read_type = 'extracted';
  }

  unless ( $conf{bwa_postprocess_mapfile} )
  {
    $conf{bwa_postprocess_mapfile} = 'NOMAP';
  }
  unless ( $conf{bwa_options} ) { $conf{bwa_options} = "" }
  my $tmp = $conf{bwa_options};
  $tmp =~ s/ //g;
  $tmp =~ s/\t//g;
  $conf{MOCAT_mapping_mode} = "bwa$tmp";

  #  if ( $conf{MOCAT_mapping_mode} eq 'unique' )
  #  {
  #    die "ERROR & EXIT: bwa mapping mode only allows 'allbest' mapping mode";
  #  } elsif ( $conf{MOCAT_mapping_mode} eq 'random' )
  #  {
  #    die "ERROR & EXIT: bwa mapping mode only allows 'allbest' mapping mode";
  #  } elsif ( $conf{MOCAT_mapping_mode} eq 'allbest' )
  #  {
  #    #$mapping_mode = "-a";
  #  } else
  #  {
  #    die "ERROR & EXIT: Unknown MOCAT_mapping_mode";
  #  }
  if ( $reads eq 'reads.processed' )
  {
    $screen_source = "reads.processed";
    $screen_save   = "reads.processed";
  } else
  {
    $screen_source = "reads.$read_type.$reads";
    $screen_save   = "$read_type.$reads";
  }

  # Check DB, only if not equal to s,c,f,r AND NOT as fasta file
  if (    !$SCREEN_FASTA_FILE
       && !( $screen[0] eq 's' || $screen[0] eq 'c' || $screen[0] eq 'f' || $screen[0] eq 'r' ) )
  {
    my $counter = -1;

    my $new_db_format = MOCATCore::checkAndReturnDB( \@screen );    # this is the new DB format in v1.4+, it is only first used in Filter, but we make sure before screening that it is a valid format
    foreach my $screen (@screen)
    {

      if ( $screen =~ m/\// )
      {
        unless ( $screen =~ m/^\// )
        {
          die "ERROR & EXIT: You specified a path to a DB, (not only a single name indicating it exists in the $data_dir folder).\nBUT, the path you specified was not an ABSOLUTE path. Please Specify the path from the root like: /usr/home/data/database";
        }
      }

      $counter++;
      if ( -e "$data_dir/$screen.bwt" )
      {
        print localtime() . ": DATABASE status: Found $screen.bwt in $data_dir = ALL OK!\n";
        print localtime() . ": Continuing creating $job jobs...\n";
        unless ( -e "$data_dir/$screen.len" )
        {
          print localtime() . ": Creating length file $data_dir/$screen.len...";
          system "$scr_dir/MOCATFilter_falen.pl -infile $data_dir/$screen -outfile $data_dir/$screen.len";
          print " OK!\n";
          print localtime() . ": Continuing creating $job jobs...";
        }
      } else
      {
        if ( -e "$data_dir/$screen" )
        {
          print localtime() . ": DATABASE status: Found $screen, but not $screen.bwa in $data_dir\n";
          die "\nERROR & EXIT: $data_dir/$screen.bwt does not exist, please index the database by running: $bin_dir/bwa index $data_dir/$screen $data_dir/$screen";
        } elsif ( -e "$screen.bwa" )
        {
          print localtime() . ": DATABASE status: Found external database $screen, but it wasn't imported into $data_dir.\n";
          print localtime() . ": DATABASE status: Importing $screen into $data_dir by symlinks...";
          chomp( my $base = `basename $screen` );
          my $dirn = dirname($screen);
          system "ln -fs $dirn/$base.bwa $data_dir/ ; ln -fs $dirn/$base $data_dir/";
          if ( -e "$data_dir/$base.bwa" )
          {
            print " OK!\n";
            print "\n\nNEXT TIME: Specify this database only as '$base'\n\n\n";
            print localtime() . ": Continuing creating $job jobs...\n";
            $screen = $base;
            $screen[$counter] = $screen;
          } else
          {
            $screen = $base;
            $screen[$counter] = $screen;
            die "\nERROR & EXIT: Could not import $screen into $data_dir ($data_dir/$screen.bwa does not exists, after attempt to create it)";
          }
        } elsif ( -e "$screen" )
        {
          print localtime() . ": DATABASE status: Found external database $screen, but it wasn't indexed.\n";
          chomp( my $base = `basename $screen` );
          my $dirn = dirname($screen);
          die "\nERROR & EXIT: $data_dir/$screen.bwa does not exist, please index and link the database by running: $bin_dir/bwa index $screen $screen.bwa && ln -fs $dirn/$base.bwa $data_dir/ && ln -fs $dirn/$base $data_dir/";
        } else
        {
          die "\nERROR & EXIT: $screen is not a valid database path";
        }
      }
    }
  }    # End DB check

  # Loop over samples, and create jobs
  foreach my $sample (@samples)
  {
    my $LOG       = "2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log";
    my $LOG2      = "2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log";
    my $inputfile = "$temp_dir/$sample/temp/SDB_INPUT.$date.gz";

    my $stats;
    if ( $reads eq 'reads.processed' )
    {
      $stats = "$cwd/$sample/stats/$sample.readtrimfilter.$conf{MOCAT_data_type}.stats";
    } else
    {
      $stats = "$cwd/$sample/stats/$sample.$read_type.$reads.$conf{MOCAT_data_type}.stats";
    }

    print JOB "exec 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log && ";
    print JOB "mkdir -p $temp_dir/$sample/temp; ";

    # Make folder, but should exist
    print JOB "mkdir -p $cwd/$sample/stats";
    foreach my $screen (@screen)
    {

      # Define variables
      my $db_on_db;
      if ( $reads eq 'reads.processed' )
      {
        if ($SCREEN_FASTA_FILE)
        {
          $db_on_db = $basename;
        } else
        {
          $db_on_db = $screen;
        }
      } else
      {
        if ($SCREEN_FASTA_FILE)
        {
          $db_on_db = "$read_type.$reads.on.$basename";
        } else
        {
          $db_on_db = "$read_type.$reads.on.$screen";
        }
      }
      my $sf            = "$cwd/$sample/reads.screened.$db_on_db.$conf{MOCAT_data_type}";
      my $ef            = "$cwd/$sample/reads.extracted.$db_on_db.$conf{MOCAT_data_type}";
      my $mf            = "$cwd/$sample/reads.mapped.$db_on_db.$conf{MOCAT_data_type}";
      my $file_output   = "$mf/$sample.mapped.$screen_save.on.$screen.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}";
      my $database      = "$data_dir/$screen";
      my $output_folder = "$cwd/$sample/reads.filtered.$screen.$conf{MOCAT_data_type}";
      my $output_file   = "$output_folder/$sample.filtered.$read_type.$reads.on.$screen.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{bwa_length_cutoff}.p$conf{bwa_percent_cutoff}.z$conf{bwa_min_perc_query_aligned}";

      my $ids_file;
      if ($SCREEN_FASTA_FILE)
      {
        $ids_file = "$temp_dir/$sample/temp/$sample.aligned.$reads.on.$basename.ids";
      }
      if ( !$SCREEN_FASTA_FILE )
      {
        $ids_file = "$temp_dir/$sample/temp/$sample.aligned.$reads.on.$screen.ids";
      }
      my $stats_file  = "$cwd/$sample/stats/$sample.screened.$db_on_db.$conf{MOCAT_data_type}.stats";
      my $estats_file = "$cwd/$sample/stats/$sample.extracted.$db_on_db.$conf{MOCAT_data_type}.stats";

      my $len_file = "$temp_dir/$sample/temp/lengths.$date";

      # Redefine variables if screening against a scaftig, contig, scaffold or revised scaftig
      if (    $screen eq 's'
           || $screen eq 'c'
           || $screen eq 'f'
           || $screen eq 'r' )
      {
        ( $max, $avg, $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-r" );
        if ( $screen eq 's' )
        {
          $end = 'scaftig';
        }
        if ( $screen eq 'c' )
        {
          $end = 'contig';
        }
        if ( $screen eq 'f' )
        {
          $end = 'scafSeq';
        }
        if ( $screen eq 'r' )
        {
          $assembly_type = 'assembly.revised';
          $end           = 'scaftig';
        }
        my $addon = "";
        if ( $reads eq 'reads.processed' )
        {
          $addon = "";
        }
        $sf          = "$cwd/$sample/reads.screened.$end.$assembly_type.K$kmer.$addon$conf{MOCAT_data_type}";
        $ef          = "$cwd/$sample/reads.extracted.$end.$assembly_type.K$kmer.$addon$conf{MOCAT_data_type}";
        $mf          = "$cwd/$sample/reads.mapped.$end.$assembly_type.K$kmer.$addon$conf{MOCAT_data_type}";
        $database    = "$cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end";
        $file_output = "$mf/$sample.mapped.$reads.on.$end.$assembly_type.K$kmer.$addon$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}";
        $ids_file    = "$temp_dir/$sample/temp/$sample.aligned.$reads.on.$assembly_type.$end.ids";
        $stats_file  = "$cwd/$sample/stats/$sample.screened.$end.$assembly_type.K$kmer$addon.stats";
        $estats_file = "$cwd/$sample/stats/$sample.extracted.$end.$assembly_type.K$kmer$addon.stats";

        $output_folder = "$cwd/$sample/reads.filtered.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}";
        $output_file   = "$output_folder/$sample.filtered.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{bwa_length_cutoff}.p$conf{bwa_percent_cutoff}.z$conf{bwa_min_perc_query_aligned}";

        # Check if files exist
        unless ( -e "$cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end" || -e "$cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end.gz" )
        {
          die "ERROR & EXIT: Missing $cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end(.gz)\nDid you enter the correct -r option?";
        }

        # If unzipped files doesn't exist
        unless ($only_regenerate_reads)
        {
          unless ( -e "$cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end" )
          {
            print JOB " && $ZIP -dc $cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end.gz > $cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end ";
          }
        }

        # If file not indexed, index it
        unless ($only_regenerate_reads)
        {
          unless ( -e "$database.bwt" )
          {
            print JOB " && $bin_dir/bwa index $cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end $cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end.bwa $LOG ";
          }
        }

        if ( $conf{bwa_make_unique_sorted} eq 'yes' )
        {

          # Get input file, and check if it exists, also create .len file
          my $assembly_file = "$cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end";
          $len_file = "$assembly_file.len";
          ( -e "$assembly_file.gz" ) or die "\nERROR & EXIT: Missing $end file: $assembly_file.gz";
          unless ( -e $len_file )
          {
            print JOB "$scr_dir/MOCATFilter_falen.pl -infile $assembly_file.gz -outfile $len_file -zip && ";
          }
        }
      }

      # Get lane IDs
      my @lanes = `ls -1 $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/*pair.1.fq.gz`;
      foreach my $i ( 0 .. ( scalar @lanes - 1 ) )
      {
        chomp( $lanes[$i] );
        $lanes[$i] =~ s/$cwd\/$sample\/$screen_source.$conf{MOCAT_data_type}//;
        $lanes[$i] =~ s/.pair.1.fq.gz//;
      }

      # Make folders
      if ($screened_files)
      {
        print JOB " && mkdir -p $sf ";
      }
      if ($extracted_files)
      {
        print JOB " && mkdir -p $ef ";
      }

      print JOB " mkdir -p $output_folder && cat $data_dir/" . join( ".len $data_dir/", @screen ) . ".len | sort -u > $len_file ";

      # bwa mapping
      chomp( my @ONES = `ls $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/*pair.1.fq.gz` );
      my $ONES = join ",", @ONES;
      chomp( my @TWOS = `ls $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/*pair.2.fq.gz` );
      my $TWOS = join ",", @TWOS;
      chomp( my @UNIQUES = `ls $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/*single.fq.gz` );
      my $UNIQUES = join ",", @UNIQUES;

      print JOB " && $bin_dir/bwa.fixed mem $conf{bwa_options} -t $processors2 $database <($src_dir/MOCATbwa_load.pl 1 $ONES $UNIQUES) <($src_dir/MOCATbwa_load.pl 2 $TWOS $UNIQUES) $LOG2 ";
      print JOB " | grep -v '\^\@' | awk '{if(\$6 != \"*\"){print}}' ";

      # This was initially added for Bowtie2, but because we hacked the source code and added bwa_load it's not needed here anymore
      #      if ( $conf{bwa_paired_end_filtering} eq 'yes' && $conf{MOCAT_paired_end} eq "yes")
      #      {
      #        print JOB "| $scr_dir/MOCATBowtie2_process.pl ";
      #      } else
      #      {
      #        print JOB "| perl -F\"\\t\" -lane 'if (\$F[2] ne \"*\"){ unless (\$F[0] =~ m/\\/1\$/){\$F[0]=\"\$F[0]/1\"}; print join \"\\t\", \@F}' ";
      #      }

      if ( -e "$database.coord" )
      {
        print JOB " | perl $scr_dir/MOCATFilter_remove_in_padded.pl -db $database -sam $LOG2";
      }

      #print JOB " | perl $src_dir/MOCATbwa_filter.pl $conf{screen_percent_cutoff} $conf{screen_length_cutoff} $LOG2 ";
      print JOB " | $bin_dir/msamtools$OSX -S --besthit -m filter -l $conf{bwa_length_cutoff} -p $conf{bwa_percent_cutoff} -z $conf{bwa_min_perc_query_aligned} -t $len_file - $LOG2 ";

      if ( $conf{bwa_paired_end_filtering} eq 'yes' && $conf{MOCAT_paired_end} eq "yes" )
      {
        print JOB "| $scr_dir/MOCATFilter_filterPE.pl $LOG2 ";
      }

      #print JOB " | $scr_dir/MOCATFilter_besthit.pl $LOG2 ";
      print JOB " | $scr_dir/MOCATFilter_stats_bwa.pl -length $conf{filter_length_cutoff} -identity $conf{filter_percent_cutoff} -db $screen -format SAM -stats $stats $LOG2 ";
      print JOB " | $bin_dir/msamtools$OSX -Sb -m merge -t $len_file - $LOG2 > $output_file.bam.tmp";
      print JOB " && sync  && test -e $output_file.bam.tmp && mv $output_file.bam.tmp $output_file.bam && test -e $output_file.bam ";
      
      # shouldn't we make the $ids_file? NO this is generated inside the MOCATScreen_filter.pl file 

      if ( $conf{bwa_postprocess} eq 'yes' )
      {
        print JOB " && $src_dir/MOCATbwa_postprocess.sh $output_file.bam $len_file $conf{bwa_postprocess_mapfile} $bin_dir $LOG2 ";
      }

      # Create the screened and extracted files
      print JOB " && $scr_dir/MOCATScreen_filter.pl -zip \'$ZIP\' -ziplevel $ziplevel ";
      print JOB " -in ";
      foreach my $lane (@lanes)
      {
        print JOB " $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/$lane ";
      }
      if ($screened_files)
      {
        print JOB " -out ";
        foreach my $lane (@lanes)
        {
          print JOB " $sf/$lane.screened.$screen ";
        }
      }
      if ($extracted_files)
      {
        print JOB " -ex ";
        foreach my $lane (@lanes)
        {
          print JOB " $ef/$lane.extracted.$screen ";
        }
      }
      print JOB " -bindir $bin_dir -bwa $output_file.bam -toremove not_used_in_BWA_mode -stats $stats_file -estats $estats_file -identity $conf{bwa_percent_cutoff} -length $conf{bwa_length_cutoff} -soapmaxmm NA $LOG ";

      # Remove ids file, it can get really big
      # NOT GENERATED print JOB " && rm -f $ids_file";

    }    # End, for each DB

    # If not scfr, remove temp file
    if (    !( $screen[0] eq 's' || $screen[0] eq 'c' || $screen[0] eq 'f' || $screen[0] eq 'r' )
         && !$SCREEN_FASTA_FILE )
    {
      print JOB " && rm -f $inputfile";
    }

    # If contig, scaftig... remove the unzipped file
    if (    $screen[0] eq 's'
         || $screen[0] eq 'c'
         || $screen[0] eq 'f'
         || $screen[0] eq 'r' )
    {
      print JOB " && rm -f $cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end";
    }
    print JOB "\n";
  }    # End, each sample
  print localtime() . ": Created jobs!\n";
  print localtime() . ": Temp directory is $temp_dir/SAMPLE/temp\n";
  close JOB;
}

sub pre_check_files
{

  # Define variables
  print localtime() . ": Checking files...";
  my $read_type = 'screened';
  if ($use_extracted_reads)
  {
    $read_type = 'extracted';
  }
  foreach my $sample (@samples)
  {

    # Define variables
    my $screen_source;
    if ( $reads eq 'reads.processed' )
    {
      $screen_source = "reads.processed.$conf{MOCAT_data_type}";
    } else
    {
      $screen_source = "reads.$read_type.$reads.$conf{MOCAT_data_type}";
    }

    # Check files
    my @F          = `ls -1 $cwd/$sample/$screen_source/*pair*fq.gz $cwd/$sample/$screen_source/*single*fq.gz 2>/dev/null`;
    my @statsfiles = `ls -1 $cwd/$sample 2>/dev/null`;
    my $line       = "";
    foreach my $sf (@statsfiles)
    {
      if ( $sf =~ m/^reads.screened/ || $sf =~ m/^reads.processed/ )
      {
        chomp($sf);
        $sf =~ s/^reads.screened.//;
        $sf =~ s/^reads.extracted.//;
        $sf =~ s/\.solexaqa$//;
        $sf =~ s/\.fastx$//;
        $line = $line . "- $sf\n";
      }

    }
    unless ( ( scalar @F % 3 ) == 0 && scalar @F >= 3 )
    {
      die "\nERROR & EXIT: Missing $screen_source files for $sample.\nPossible reason: Read Trim Filter has not been performed. Or you specified the wrong database.\nSpecify -r to be either 'FASTA FILE', 'DATABASE' or 'reads.processed'\nPerhaps you meant either of these as -r:\n$line";
    }
  }
  print " OK!\n";
}

sub post_check_files
{

  # Define variables
  print localtime() . ": Checking files... ";
  my $p         = 1;
  my $read_type = 'screened';
  my $screen_before;
  my $db_on_db;
  my $basename;
  if ($use_extracted_reads)
  {
    $read_type = 'extracted';
  }
  if ( $reads eq 'reads.processed' )
  {
    $screen_before = "reads.processed";
  } else
  {
    $screen_before = "$read_type.$reads";
  }
  if ($SCREEN_FASTA_FILE)
  {
    chomp( $basename = `basename $screen[0]` );
  }
  foreach my $sample (@samples)
  {
    foreach my $screen (@screen)
    {

      # Define variables
      my $end;
      my $assembly_type = "assembly";
      my $db_on_db;
      if ( $screen eq 's' )
      {
        $end = 'scaftig';
      } elsif ( $screen eq 'c' )
      {
        $end = 'contig';
      } elsif ( $screen eq 'f' )
      {
        $end = 'scafSeq';
      } elsif ( $screen eq 'r' )
      {
        $assembly_type = 'assembly.revised';
        $end           = 'scaftig';
      } else
      {
        if ( $reads eq 'reads.processed' )
        {

          $db_on_db = $screen;

        } else
        {
          $db_on_db = "$read_type.$reads.on.$screen";
        }
      }
      my ( $file_output, $sf, $mf, $ef );
      if (    $screen eq 's'
           || $screen eq 'c'
           || $screen eq 'f'
           || $screen eq 'r' )
      {
        ( my $max, my $avg, my $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-r" );
        my $addon = "";
        if ( $reads eq 'reads.processed' )
        {
          $addon = "";
        }
        $sf = "$cwd/$sample/reads.screened.$end.$assembly_type.K$kmer.$addon$conf{MOCAT_data_type}";
        $ef = "$cwd/$sample/reads.extracted.$end.$assembly_type.K$kmer.$addon$conf{MOCAT_data_type}";

        #$mf          = "$cwd/$sample/reads.mapped.$end.$assembly_type.K$kmer.$addon$conf{MOCAT_data_type}";
        #$file_output = "$mf/$sample.mapped.$reads.on.$end.$assembly_type.K$kmer.$addon$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.soap.gz";
        $mf          = "$cwd/$sample/reads.filtered.$screen.$conf{MOCAT_data_type}";
        $file_output = "$mf/$sample.filtered.$read_type.$reads.on.$screen.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{bwa_length_cutoff}.p$conf{bwa_percent_cutoff}.z$conf{bwa_min_perc_query_aligned}.bam";
      } else
      {
        $sf = "$cwd/$sample/reads.screened.$db_on_db.$conf{MOCAT_data_type}";
        $ef = "$cwd/$sample/reads.extracted.$db_on_db.$conf{MOCAT_data_type}";

        #$mf          = "$cwd/$sample/reads.mapped.$db_on_db.$conf{MOCAT_data_type}";
        #$file_output = "$mf/$sample.mapped.$screen_before.on.$screen.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.soap.gz";
        $mf          = "$cwd/$sample/reads.filtered.$screen.$conf{MOCAT_data_type}";
        $file_output = "$mf/$sample.filtered.$read_type.$reads.on.$screen.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{bwa_length_cutoff}.p$conf{bwa_percent_cutoff}.z$conf{bwa_min_perc_query_aligned}.bam";

      }

      # Check map file
      if ( !$SCREEN_FASTA_FILE )
      {
        unless ( -e "$file_output" )
        {
          print "\nERROR: Missing $file_output";
          $p = 0;
        }
      }

      # If produced, check screened and extracted files
      if ($screened_files)
      {
        my @F = `ls -1 $sf/*pair*.fq.gz $sf/*single*.fq.gz 2>/dev/null`;
        unless ( ( scalar @F % 3 ) == 0
                 && scalar @F >= 3 )
        {
          print "\nERROR: Missing screened pair and/or single files for $sample";
          $p = 0;
        }
      }
      if ($extracted_files)
      {
        my @F = `ls -1 $ef/*pair*.fq.gz $ef/*single*.fq.gz 2>/dev/null`;
        unless ( ( scalar @F % 3 ) == 0
                 && scalar @F >= 3 )
        {
          print "\nERROR: Missing extracted pair and/or single files for $sample";
          $p = 0;
        }
      }

    }
  }
  if ( $p == 1 )
  {
    print " OK!\n";
  } else
  {
    die "\nERROR & EXIT: Missing one or more files!";
  }
}
1;
