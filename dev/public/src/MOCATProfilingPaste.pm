package MOCATProfilingPaste;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2014
# This code is released under GNU GPL v3.

sub system_
{
  my $cmd = shift;
  system($cmd) == 0 or die("\nERROR & EXIT: system($cmd) failed: $!\n");
}

sub system_2
{
  my $cmd = shift;

  #print "\n\n$cmd\n\n";
  system($cmd) == 0 or die("\nERROR & EXIT: system($cmd) failed: $!\n\nPlease check that the count type, norm type and levels in the config file are valid,\nbecause it seems like the file that was supposed to be extracted did not exist in the zip file.\n");
}

sub paste2
{

### BWA SUPPORT ###
  my $bwaAddon = "";
  if ( $do_profiling_bwa[0] )
  {
    unless ( $conf{bwa_options} ) { $conf{bwa_options} = "" }
    my $tmp = $conf{bwa_options};
    $tmp =~ s/ //g;
    $tmp =~ s/\t//g;
    $conf{MOCAT_mapping_mode} = "bwa$tmp";
    $bwaAddon = ".z$conf{bwa_min_perc_query_aligned}";
  }
### BWA SUPPORT ###

  ### DEFINE VARIABLES AND OPEN OUTPUT FILE ###
  my $job        = $_[0];
  my $processors = $_[1];
  $ZIP =~ s/pigz.*/pigz -p $processors/;
  my $read_type = 'screened';
  if ($use_extracted_reads)
  {
    $read_type = 'extracted';
  }
  my $uniq;
  my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
  if ($SHM)
  {
    $temp_dir = "/dev/shm/$username/MOCAT_temp/ProfilingPaste.$date";
  } else
  {
    $temp_dir = "$temp_dir/ProfilingPaste.$date";
  }

  print localtime() . ": PASTING PROFILES...\n";
  print localtime() . ": temp dir = $temp_dir\n";

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

  my $databases_old = join( "_AND_", @databases );
  my $databases     = MOCATCore::checkAndReturnDB( \@databases );
  my $assembly_type = 'assembly';
  my $end;
  my $assembly = 0;
  my $profiling_mode2;
  my $profiling_mode3;    # this is just set to '' for gene and functional
  chomp( my $bn = `basename $sample_file` );
  my $profiling_mode_folder;
  my $OUTFOLDER;
  my %LEVELS;
  my %VERSION;

  # Here it gets a bit tricky, this main loop should be run ONCE if not any assmelby type, but has to be run one by one for each sample if assembly...
  my $LOOP             = 0;
  my @original_samples = @samples;
  if (    $databases[0] eq 's'
       || $databases[0] eq 'c'
       || $databases[0] eq 'f'
       || $databases[0] eq 'r' )
  {
    $LOOP = scalar @samples - 1;
  }

  for my $LOOPER ( 0 .. $LOOP )
  {
    MOCATCore::mkdir_or_die("$temp_dir/in");
    MOCATCore::mkdir_or_die("$temp_dir/out");
    if (    $databases[0] eq 's'
         || $databases[0] eq 'c'
         || $databases[0] eq 'f'
         || $databases[0] eq 'r' )
    {
      @samples = $original_samples[$LOOPER];
      if ( $databases[0] eq 's' ) { $end = 'scaftig' }
      if ( $databases[0] eq 'c' ) { $end = 'contig' }
      if ( $databases[0] eq 'f' ) { $end = 'scafSeq' }
      if ( $databases[0] eq 'r' )
      {
        $assembly_type = 'assembly.revised';
        $end           = 'scaftig';
      }
      $assembly       = 1;
      $databases      = "$assembly_type.$end";
      $profiling_mode = 'gene';
    }

    my @BOI;
    my @NORM;
    my @LOOP;
    my @LEVELS;
    my @levels;

    unless ($profiling_mode)
    {
      $profiling_mode = '';
    }

    if ( $profiling_mode eq 'RefMG' )
    {
      $profiling_mode = 'NCBI';
    }
    if ( ( $profiling_mode eq 'mOTU' || $profiling_mode eq 'gene' ) )
    {
      if ( $profiling_mode eq 'mOTU' )
      {
        @LEVELS                = ('mOTU');
        $profiling_mode2       = "motu";
        $profiling_mode3       = "mOTU";
        $profiling_mode_folder = "taxonomic";
        MOCATCore::mkdir_or_die("$temp_dir/out/COGs");
      } elsif ( $profiling_mode eq 'gene' )
      {
        @LEVELS                = ('gene');
        $profiling_mode2       = "gene";
        $profiling_mode3       = "";
        $profiling_mode_folder = "gene";
      }
      my $ct = $conf{profiling_PRINT_COUNT_TYPE};
      unless ($ct)
      {
        $ct = "";
      }
      $ct =~ s/ //g;
      my $nt = $conf{profiling_PRINT_NORM_TYPES};
      unless ($nt)
      {
        $nt = "";
      }
      $nt =~ s/ //g;
      if ( $ct ne '' && $nt ne '' && $ct ne 'undef' && $nt ne 'undef' )
      {
        @BOI  = split ",", $ct;
        @NORM = split ",", $nt;
        foreach my $boi (@BOI)
        {
          foreach my $norm (@NORM)
          {
            foreach my $i (@LEVELS)
            {
              push @LOOP, "$boi:$norm:$i";
              if ($CALCULATE_HORIZONTAL_COVERAGE)
              {
                push @LOOP, "horizontal:x:$i";
              }

            }
          }
        }
      } else
      {
        @BOI = ( 'base', 'insert' );
        @NORM = ( 'raw', 'norm', 'only.unique.raw', 'only.unique.norm', 'scaled', 'only.unique.scaled', 'mm.dist.among.unique.raw', 'mm.dist.among.unique.norm', 'mm.dist.among.unique.scaled' );
        @levels = @LEVELS;
        foreach my $boi (@BOI)
        {
          foreach my $norm (@NORM)
          {
            foreach my $i (@LEVELS)
            {
              push @LOOP, "$boi:$norm:$i";
              if ($CALCULATE_HORIZONTAL_COVERAGE)
              {
                push @LOOP, "horizontal:x:$i";
              }

            }
          }
        }
      }
    } elsif ( $profiling_mode eq 'functional' || $profiling_mode eq 'NCBI' )
    {
      my $NAME1;
      my $NAME2;
      if ( $profiling_mode eq 'functional' )
      {
        $NAME1 = "profiling_PRINT_ONLY_THESE_FUNCTIONAL_LEVELS";
        $NAME2 = "profiling_PRINT_FUNCTIONAL_LEVELS";
        my $databases                = MOCATCore::checkAndReturnDB( \@do_profiling );
        my $profiling_functional_map = "$data_dir/$databases.functional.map";
        my $MAP;
        if ( -e "$profiling_functional_map.gz" )
        {
          open $MAP, "gunzip -c $profiling_functional_map.gz | ";
        } else
        {
          open $MAP, "<", "$profiling_functional_map" or die "ERROR & EXIT: Cannot load functional map $profiling_functional_map";
        }
        chomp( my $levels = <$MAP> );
        $levels =~ s/^#//;
        @LEVELS = split /\s+/, $levels;
        close $MAP;
        $profiling_mode2       = "functional";
        $profiling_mode3       = "";
        $profiling_mode_folder = "functional";

        if ( $conf{profiling_PRINT_FUNCTIONAL_LEVELS} )
        {
          $conf{profiling_PRINT_FUNCTIONAL_LEVELS} =~ s/\s+//g;
          my @lc = split ",", $conf{profiling_PRINT_FUNCTIONAL_LEVELS};
          print localtime() . ": These levels defined in config:   " . join( ", ", @lc ) . "\n";
          print localtime() . ": These levels defined in map file: " . join( ", ", @LEVELS ) . "\n";
          my @intersection =
            grep { defined } @{ { map { lc, => $_ } @LEVELS } }{ map { lc } @lc };
          @LEVELS = @intersection;
          print localtime() . ": Summarizing these levels:         " . join( ", ", @LEVELS ) . "\n";
        } else
        {
          print localtime() . ": Summarizing these levels: " . join( ", ", @LEVELS ) . "\n";
        }
      } else
      {
        @LEVELS                = ( 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'specI_cluster', 'taxaid' );
        $profiling_mode2       = "ncbi";
        $profiling_mode3       = "NCBI";
        $profiling_mode_folder = "taxonomic";
        $NAME1                 = "profiling_PRINT_ONLY_THESE_NCBI_LEVELS";
        $NAME2                 = "profiling_PRINT_NCBI_LEVELS";
      }
      foreach my $k (@LEVELS)
      {
        $LEVELS{$k} = 1;
      }

      my $line = $conf{$NAME1};
      $line =~ s/ //g;
      if ( $line ne 'undef' && $line ne '' )
      {
        my @line = split ",", $line;
        foreach my $unit (@line)
        {
          my @unit = split ":", $unit;
          if ( $unit[0] && $unit[1] && $unit[2] )
          {
            push @LOOP, "$unit[0]:$unit[1]:$unit[2]";
            if ($CALCULATE_HORIZONTAL_COVERAGE)
            {
              push @LOOP, "horizontal:x:$unit[2]";
            }
          }
        }
      } else
      {
        my $ct = $conf{profiling_PRINT_COUNT_TYPE};
        $ct =~ s/ //g;
        my $nt = $conf{profiling_PRINT_NORM_TYPES};
        $nt =~ s/ //g;
        my $l = $conf{$NAME2};
        unless ($l)
        {
          $l = "";
        }
        $l =~ s/ //g;
        if ( $ct ne '' && $nt ne '' && $l ne '' && $ct ne 'undef' && $nt ne 'undef' && $l ne 'undef' )
        {
          @BOI    = split ",", $ct;
          @NORM   = split ",", $nt;
          @levels = split ",", $l;
          foreach my $boi (@BOI)
          {
            foreach my $norm (@NORM)
            {
              foreach my $i (@levels)
              {
                push @LOOP, "$boi:$norm:$i";
                if ($CALCULATE_HORIZONTAL_COVERAGE)
                {
                  push @LOOP, "horizontal:x:$i";
                }
              }
            }
          }
        } else
        {
          @BOI = ( 'base', 'insert' );
          @NORM = ( 'raw', 'norm', 'only.unique.raw', 'only.unique.norm', 'scaled', 'only.unique.scaled', 'mm.dist.among.unique.raw', 'mm.dist.among.unique.norm', 'mm.dist.among.unique.scaled' );
          foreach my $boi (@BOI)
          {
            foreach my $norm (@NORM)
            {
              foreach my $i (@LEVELS)
              {
                push @LOOP, "$boi:$norm:$i";
                if ($CALCULATE_HORIZONTAL_COVERAGE)
                {
                  push @LOOP, "horizontal:x:$i";
                }
              }
            }
          }
        }
      }
    } else
    {
      die "ERROR & EXIT: Unknown profiling mode. Please specify a correct -mode (NCBI, mOTU, gene or functional)";
    }

    # set mode back from gene
    if ( $assembly == 1 )
    {
      if ( $databases[0] eq 's' )
      {
        $profiling_mode  = 'gene';
        $profiling_mode2 = 'scaftigs';
        $profiling_mode3 = $samples[0];
      }
      if ( $databases[0] eq 'c' )
      {
        $profiling_mode  = 'gene';
        $profiling_mode2 = 'contigs';
        $profiling_mode3 = $samples[0];
      }
      if ( $databases[0] eq 'f' )
      {
        $profiling_mode  = 'gene';
        $profiling_mode2 = 'scafSeqs';
        $profiling_mode3 = $samples[0];
      }
      if ( $databases[0] eq 'r' )
      {
        $profiling_mode  = 'gene';
        $profiling_mode2 = 'revised.scaftigs';
        $profiling_mode3 = $samples[0];
      }
      $profiling_mode_folder = $profiling_mode2;
      $sample_file_basename  = $samples[0];
    }

    if ( scalar @LOOP == 0 )
    {
      @BOI = ( 'base', 'insert' );
      @NORM = ( 'raw', 'norm', 'only.unique.raw', 'only.unique.norm', 'scaled', 'only.unique.scaled', 'mm.dist.among.unique.raw', 'mm.dist.among.unique.norm', 'mm.dist.among.unique.scaled' );
      @levels = @LEVELS;
      foreach my $boi (@BOI)
      {
        foreach my $norm (@NORM)
        {
          foreach my $i (@LEVELS)
          {
            push @LOOP, "$boi:$norm:$i";
            if ($CALCULATE_HORIZONTAL_COVERAGE)
            {
              push @LOOP, "horizontal:x:$i";
            }
          }
        }
      }
    }

    if ( $assembly == 1 )
    {
      $OUTFOLDER = "$cwd/PROFILES/$profiling_mode_folder.profiles";
    } else
    {
      $OUTFOLDER = "$cwd/PROFILES/$profiling_mode_folder.profiles/$sample_file_basename/$profiling_mode3";
    }
    MOCATCore::mkdir_or_die("$OUTFOLDER");
    my $printError = 1;
    my $VERSION;

    my %tmp;
    foreach my $l (@LOOP)
    {
      $tmp{$l} = 1;
    }
    @LOOP = ();
    foreach my $l ( sort keys %tmp )
    {
      push @LOOP, $l;
    }

    foreach my $loop (@LOOP)
    {
      my @loop = split ":", $loop;
      my $boi  = $loop[0];
      my $norm = $loop[1];
      my $i    = $loop[2];
      my $full_end;
      if ( $boi eq 'horizontal' )
      {
        $full_end = "$boi.$i";
      } else
      {
        $full_end = "$boi.$norm.$i";
      }
      print localtime() . ": Preparing & processing $full_end\n";
      my $to_paste = '';
      my $rownames;
      my $rownamesOut;
      my $rownames_old;
      my $folder;
      my $folder_old;
      my $file;
      my $fileOut;
      my $file_old;
      my %hash;
      $VERSION = "";
      my $counter1 = 0;
      my $counter2 = 0;

      foreach my $sample (@samples)
      {

        if ( $pcf_uniq && $i eq 'gene' )
        {
          $uniq = ".uniq";
        } else
        {
          $uniq = "";
        }

        $counter1++;
        
        
  if ( !$assembly )
        {
          $folder_old   = "$cwd/$sample/$sample.$profiling_mode.profile.$read_type.$reads.on.$databases_old.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon";
          $file         = "$sample.$profiling_mode.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$full_end";
          $file_old     = "$sample.$profiling_mode.profile.$read_type.$reads.on.$databases_old.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$full_end";
          $folder       = "$cwd/$sample/$sample.$profiling_mode.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon";
          $fileOut      = "$sample.$profiling_mode.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$full_end";
          $rownamesOut  = "$sample.$profiling_mode.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$i.rownames$uniq";
          $rownames_old = "$sample.$profiling_mode.profile.$read_type.$reads.on.$databases_old.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$i.rownames$uniq";
          $rownames     = "$sample.$profiling_mode.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$i.rownames$uniq";

        } elsif ($assembly)
        {
          ( my $max, my $avg, my $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-r" );
          $folder      = "$cwd/$sample/$sample.$profiling_mode.profile.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon";
          $file        = "$sample.$profiling_mode.profile.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$full_end";
          $rownamesOut = "$sample.$profiling_mode.profile.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$i.rownames$uniq";
          $fileOut     = "$sample.$profiling_mode.profile.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$full_end";
          $rownames    = "$sample.$profiling_mode.profile.$reads.on.$end.$assembly_type.K$kmer.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$i.rownames$uniq";
          $folder_old  = "";
          $file_old    = "";
        }
   


        # In MOCAT v1.4+ the data files are in zip format. Two options: 1. extract all files for all samples,
        # but this would use A LOT of disk space. Instead we extract the files one by one. This is slower, but hopefully still feasible.

        unless ( -e "$folder.zip" )
        {
          if ( -e "$folder_old.zip" && $assembly == 0 )
          {
            if ($printError)
            {
              print localtime() . ": WARNING! Could not find new format input file '$folder.zip', but using old format input file '$folder_old.zip'\n";
            }
            $folder     = $folder_old;
            $rownames   = $rownames_old;
            $file       = $file_old;
            $printError = undef;
          } else
          {
            die "\nERROR & EXIT: Missing profiling file: $folder.zip";
          }
        }

        unless ( -e "$temp_dir/in/$rownamesOut" )
        {
          system_2("$bin_dir/unzip -p $folder.zip $rownames > $temp_dir/in/$rownamesOut");
          unless ( -e "$temp_dir/in/$rownamesOut" )
          {
            die "\nERROR & EXIT: Missing $temp_dir/in/$rownamesOut";
          }
        }
        system_("$bin_dir/unzip -p $folder.zip $file > $temp_dir/in/$fileOut");
        unless ( -e "$temp_dir/in/$fileOut" )
        {
          die "\nERROR & EXIT: Missing $temp_dir/in/$fileOut";
        }
        unless ( $hash{$counter2} )
        {
          $hash{$counter2} = "";
        }
        $hash{$counter2} = $hash{$counter2} . " $temp_dir/in/$fileOut ";

        $VERSION{$version} = 1;
        my $version = "MOCAT version unknown";
        if ( -e "$folder.zip.version" )
        {
          chomp( $version = `head -1 $folder.zip.version` );
        }
        $VERSION .= "\t$version";

        if ( $counter1 == 100 )
        {
          $counter1 = 0;
          $counter2++;
        }
      }
      my $norm2 = $norm;

      #my $name          = "$temp_dir/out/$profiling_mode2.profiles/$databases/$read_type.$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}/$sample_file_basename.$profiling_mode2.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.$boi.$norm.$i";
      my $name = "$temp_dir/out/$sample_file_basename.$profiling_mode2.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$full_end";
      if ( $assembly == 1 )
      {
        $name = "$temp_dir/out/$sample_file_basename.$profiling_mode2.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$boi.$norm.$profiling_mode_folder";
      }
      my $original_name = "$OUTFOLDER/$sample_file_basename.$profiling_mode2.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$full_end";
      if ( $assembly == 1 )
      {
        $original_name = "$OUTFOLDER/$sample_file_basename.$profiling_mode2.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$boi.$norm.$profiling_mode2";
      }
      if ( $profiling_mode eq 'functional' )
      {
        $original_name = "$OUTFOLDER*/$sample_file_basename.$profiling_mode2.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$full_end";
      }

      if ( $i eq 'gene' && $profiling_mode eq 'functional' )
      {
        $name          = "$temp_dir/out/$sample_file_basename.gene.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$full_end";
        $original_name = "$OUTFOLDER*/$sample_file_basename.gene.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$full_end";
      }

      my $COG_name          = "$temp_dir/out/COGs/$sample_file_basename.$profiling_mode2.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$full_end";
      my $original_COG_name = "$OUTFOLDER/$sample_file_basename.$profiling_mode2.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.$full_end";
      system_("cp `cat $temp_dir/in/$rownamesOut` $name.tmp1");

      foreach my $f ( sort { $a <=> $b } keys %hash )
      {
        system_("paste $name.tmp1 $hash{$f} > $name.tmp2 && mv $name.tmp2 $name.tmp1 && rm $hash{$f}");
      }
      chomp( my $currentDate = `date` );
      my $nice = "$boi.$norm";
      if ( $nice eq 'horizontal.x' )
      {
        $nice = "horizontal";
      }
      open C, ">$name.comments" or die "\nERROR & EXIT: Cannot write to $!";
      print C "# MOCAT $version ($MOCAT_ID) $profiling_mode2 profile ($nice) summarized at $i level\n";
      print C "# Created $currentDate by $username @ $hostname with identifier $date\n";
      print C "# $cwd/MOCAT.pl " . join( " ", @args ) . "\n";
      print C "# Original filename: $original_name\n";
      unless ( scalar keys %VERSION == 1 )
      {
        print C "# Versions$VERSION\n";
      }
      close C;
      system_("cat $name.comments $name.tmp1 > $name && rm -f $name.comments $name.tmp1");

      # if mOTU mode
      if ( $profiling_mode eq 'mOTU' )
      {
        my @mg = ( 'COG0012', 'COG0049', 'COG0052', 'COG0048', 'COG0016', 'COG0018', 'COG0080', 'COG0088', 'COG0081', 'COG0087', 'COG0090', 'COG0085', 'COG0091', 'COG0092', 'COG0093', 'COG0094', 'COG0096', 'COG0097', 'COG0098', 'COG0099', 'COG0100', 'COG0102', 'COG0103', 'COG0124', 'COG0172', 'COG0184', 'COG0185', 'COG0186', 'COG0197', 'COG0200', 'COG0201', 'COG0202', 'COG0215', 'COG0256', 'COG0522', 'COG0495', 'COG0533', 'COG0525', 'COG0552', 'COG0541' );
        print localtime() . ": Grouping $i by COGs...";
        my $nice = "$boi.$norm";
        if ( $nice eq 'horizontal.x' )
        {
          $nice = "horizontal";
        }
        foreach my $i (@mg)
        {
          open IN,  "$name"         or die "\nERROR & EXIT: Missing $name.tmp";
          open OUT, ">$COG_name.$i" or die "\nERROR & EXIT: Cannot write to $COG_name.$i";
          print OUT "# MOCAT $version ($MOCAT_ID) $profiling_mode2 profile ($nice) summarized for $i\n";
          print OUT "# Created $currentDate by $username @ $hostname with identifier $date\n";
          print OUT "# $cwd/MOCAT.pl " . join( " ", @args ) . "\n";
          print OUT "# Original filename: $original_COG_name.$i\n";
          unless ( scalar keys %VERSION == 1 )
          {
            print OUT "# Versions$VERSION\n";
          }

          while (<IN>)
          {
            if ( $. == 5 || $_ =~ /^$i\./ )
            {    ### WE HAVE 5 HEADER ROWS ###
              print OUT $_;
            }
          }
          close IN;
          close OUT;

          #system "$ZIP -$ziplevel -f $COG_name.$i.tmp && mv $COG_name.$i.tmp.gz $COG_name.$i.gz";
        }
        print " OK!\n";
      }

      #system_("$scr_dir/MOCATFraction.pl -in $name -out $name.fraction");

      #system "$ZIP -f -$ziplevel $name.tmp $name.fraction.tmp && mv $name.tmp.gz $name.gz && mv $name.fraction.tmp.gz $name.fraction.gz";
    }

    # remove temp files
    print localtime() . ": removing $temp_dir/in\n";
    system_("rm -r $temp_dir/in");

    # create final zip file
    my $name = "$OUTFOLDER/$sample_file_basename.$profiling_mode2.profiles.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon";
    print localtime() . ": zipping output files\n";

    # split gene files and functional files
    if ( $profiling_mode eq 'functional' )
    {

      my $files;
      my $ran = 0;

      # KEGG
      my %LEVELS2 = %LEVELS;

      # Previously we split it into COG and EGG profiles, but now we will place all of them in the same file, will make everyone's lives easier.
      #			if ( $LEVELS{ko} || $LEVELS{module} || $LEVELS{pathway} ) {
      #				delete $LEVELS2{ko};
      #				delete $LEVELS2{module};
      #				delete $LEVELS2{pathway};
      #				$ran = 1;
      #				$OUTFOLDER = "$cwd/PROFILES/$profiling_mode_folder.profiles/$sample_file_basename/KEGG";
      #				system_("mkdir -p $OUTFOLDER");
      #				my $name = "$OUTFOLDER/$sample_file_basename.KEGG.profiles.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
      #				chomp(my @files = `cd $temp_dir/out/; ls *.ko *.module *.pathway 2>/dev/null`);
      #				#system_("rm -f $name.zip && cd $temp_dir/out && $bin_dir/zip -q -y -D -4 $name.zip $files && rm -f $files");
      #				#system_("rm -f $name.zip && cd $temp_dir/out && echo \"$files\" | sed 's/ /\\n/g' > ../filelist && cat ../filelist | $bin_dir/zip -@ -q -y -D -4 $name.zip && for f in `cat ../filelist`; do rm -f \$f; done && rm -fr ../filelist");
      #				system_("rm -f $name.zip && cd $temp_dir/out && echo \"". join("\n", @files) ."\" > ../filelist && cat ../filelist | $bin_dir/zip -@ -q -y -D -4 $name.zip && for f in `cat ../filelist`; do rm -f \$f; done && rm -fr ../filelist");
      #			}
      #
      #			# COG
      #			if ( $LEVELS{cog} ) {
      #				delete $LEVELS2{cog};
      #				$ran = 1;
      #				$OUTFOLDER = "$cwd/PROFILES/$profiling_mode_folder.profiles/$sample_file_basename/COG";
      #				system_("mkdir -p $OUTFOLDER");
      #				$name  = "$OUTFOLDER/$sample_file_basename.COG.profiles.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}";
      #				chomp(my @files = `cd $temp_dir/out/; ls *.cog 2>/dev/null`);
      #				#system_("rm -f $name.zip && cd $temp_dir/out && $bin_dir/zip -q -y -D -4 $name.zip $files && rm -f $files");
      #				#system_("rm -f $name.zip && cd $temp_dir/out && echo \"$files\" | sed 's/ /\\n/g' > ../filelist && cat ../filelist | $bin_dir/zip -@ -q -y -D -4 $name.zip && for f in `cat ../filelist`; do rm -f \$f; done && rm -fr ../filelist");
      #				system_("rm -f $name.zip && cd $temp_dir/out && echo \"". join("\n", @files) ."\" > ../filelist && cat ../filelist | $bin_dir/zip -@ -q -y -D -4 $name.zip && for f in `cat ../filelist`; do rm -f \$f; done && rm -fr ../filelist");
      #			}

      # gene
      if ( $LEVELS{gene} )
      {
        delete $LEVELS2{gene};
        $ran       = 1;
        $OUTFOLDER = "$cwd/PROFILES/gene.profiles/$sample_file_basename/";
        system_("mkdir -p $OUTFOLDER");
        $name = "$OUTFOLDER/$sample_file_basename.gene.profiles.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon";
        system_("rm -f $name.zip && cd $temp_dir/out && ls *.gene > ../filelist && cat ../filelist | $bin_dir/zip -@ -q -y -D -4 $name.zip && for f in `cat ../filelist`; do rm -f \$f; done && rm -fr ../filelist");
      }

      #other
      if ( scalar keys %LEVELS2 > 0 )
      {
        $OUTFOLDER = "$cwd/PROFILES/$profiling_mode_folder.profiles/$sample_file_basename/";
        system_("mkdir -p $OUTFOLDER");
        my $name = "$OUTFOLDER/$sample_file_basename.functional.profiles.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon";
        system_("rm -f $name.zip && cd $temp_dir/out && ls > ../filelist && cat ../filelist | $bin_dir/zip -@ -q -y -D -4 $name.zip");

      }

      # reset for nice print
      $OUTFOLDER = "$cwd/PROFILES/$profiling_mode_folder.profiles/$sample_file_basename AND $cwd/PROFILES/gene.profiles/$sample_file_basename"

    } else
    {
      system_("rm -f $name.zip && cd $temp_dir/out && ls > ../filelist && cat ../filelist | $bin_dir/zip -@ -q -y -D -4 $name.zip && rm -fr ../filelist");
    }

    print localtime() . ": output files zipped\n";

    # Make nice output files
    if ( $OUTPUT_FOLDER =~ m/^\// )
    {
    } else
    {
      $OUTPUT_FOLDER = "$cwd/$OUTPUT_FOLDER";
    }
    MOCATCore::mkdir_or_die($OUTPUT_FOLDER);
    my @input_files1 = ();
    my @input_files2 = ();
    my @output_files = ();

    # Run perl script
    if ( $profiling_mode eq 'mOTU' )
    {
      my $LOG = "$cwd/logs/motu_profiles/jobs/MOCATJob_motu_species_profiles\_$date.log";
      MOCATCore::mkdir_or_die("$cwd/logs/motu_profiles/jobs/");
      print localtime() . ": mode is mOTO, generating mOTU profiles...\n";
      my $name = "$temp_dir/out/$sample_file_basename.$profiling_mode2.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.insert.mm.dist.among.unique.scaled.mOTU";
      my $pre  = "";
      my $file = "$sample_file_basename.$profiling_mode2.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon";
      my $db   = "$temp_dir/out";
      if ( -e "$data_dir/$databases.motu.linkage.map" )
      {
        system_("echo -e \"perl $scr_dir/MOCATPasteTaxonomyCoverageFiles_generate_mOTU_tables.pl -full-output -wd '$cwd' -map '$data_dir/$databases.motu.linkage.map' -table '$name' -prefix '$db/$pre/$file.insert.mm.dist.among.unique.scaled' >>$LOG 2>>$LOG \n\n\" >>$LOG");
        print localtime() . ": Generating mOTU tables...\n";
        system_("perl $scr_dir/MOCATPasteTaxonomyCoverageFiles_generate_mOTU_tables.pl -full-output -wd '$cwd' -map '$data_dir/$databases.motu.linkage.map' -table '$name' -prefix '$db/$pre/$file.insert.mm.dist.among.unique.scaled' >>$LOG 2>>$LOG");
        print localtime() . ": Post processing mOTU tables...\n";
        my @tables = <$db/$pre/$file*tab>;

        foreach my $table (@tables)
        {
          my $tableNoTab = $table;
          $tableNoTab =~ s/.tab$//;
          my $annot = "";
          if ( $table =~ m/annotated.mOTU.clusters/ )
          {
            $annot = " annotated";
          }
          my $frac = "";
          if ( $tableNoTab =~ m/.fraction$/ )
          {
            $frac = " (fractions)";
          }
          my $tableNoTab2 = $tableNoTab;    # tricky here, this #2 is used only for the name of the file inside the file
          $tableNoTab2 =~ s/$temp_dir\/out/$OUTFOLDER/;
          chomp( my $currentDate = `date` );
          open IN,  "$table"     or die "ERROR & EXIT: Missing $table";
          open OUT, ">$table.p2" or die "ERROR & EXIT: Cannot write to $table.p2\n";
          print OUT "# MOCAT $version ($MOCAT_ID)$annot mOTU clusters abundance table (equivalent of species level)\n";
          print OUT "# Created $currentDate by $username @ $hostname with identifier $date\n";
          print OUT "# $cwd/MOCAT.pl " . join( " ", @args ) . "\n";
          print OUT "# Original filename: $tableNoTab2\n";

          unless ( scalar keys %VERSION == 1 )
          {
            print OUT "# Versions$VERSION\n";
          }

          while (<IN>)
          {
            print OUT $_;
          }
          close IN;
          close OUT;
          system_("rm -f $table && mv $table.p2 $tableNoTab");
        }
      } else
      {
        system_("echo \"ERROR: Missing $data_dir/$databases.motu.linkage.map - could not create mOTU linkage groups.\" >>$LOG");
      }

      # Here we need to extract only these specific files...
      #		print localtime() . ": Zipping mOTU tables...\n";
      #		for my $i ( 0 .. scalar @tables - 1 ) {
      #			$tables[$i] =~ s/.tab$//;
      #		}
      #		system "$ZIP -f -$ziplevel " . join( " ", @tables );

      # Make easy output
      my @input_files0 = ( "$temp_dir/out/$file.insert.mm.dist.among.unique.scaled.annotated.mOTU.clusters", "$temp_dir/out/$file.insert.mm.dist.among.unique.scaled.mOTU.clusters" );
      foreach my $file (@input_files0)
      {
        if ( -e $file )
        {
          system_("$scr_dir/MOCATFraction.pl -in $file -out $file.fraction");
          push @input_files1, "$file.fraction";
          push @input_files1, $file;
        } else
        {
          push @input_files1, "";
          push @input_files1, "";
        }
      }
      @input_files2 = ( "$OUTFOLDER/$file.insert.mm.dist.among.unique.scaled.annotated.mOTU.clusters.fraction", "$OUTFOLDER/$file.insert.mm.dist.among.unique.scaled.annotated.mOTU.clusters", "$OUTFOLDER/$file.insert.mm.dist.among.unique.scaled.mOTU.clusters.fraction", "$OUTFOLDER/$file.insert.mm.dist.among.unique.scaled.mOTU.clusters" );
      @output_files = ( 'annotated.mOTU.abundances', 'mOTU.abundances', 'annotated.mOTU.counts', 'mOTU.counts' );

      my @p1 = ( "$data_dir/mOTU.v1.map.txt", "$data_dir/mOTU-LG.v1.annotations.txt" );
      my @p2 = ( 'mOTU.v1.map.txt',           'mOTU-LG.v1.annotations.txt' );
      for my $i ( 0 .. scalar @p1 - 1 )
      {
        system_("ln -sf $p1[$i] $OUTPUT_FOLDER/$p2[$i]");
      }
      print localtime() . ": zipping COG output files\n";
      $name = "$OUTFOLDER/$sample_file_basename.$profiling_mode2.profiles.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon";
      system_("rm -f $name.COGs.zip && cd $temp_dir/out/COGs && $bin_dir/zip -q -y -D -4 $name.COGs.zip *");
      print localtime() . ": COG output files zipped\n";
    }

    if ( $profiling_mode eq 'NCBI' )
    {
      print localtime() . ": creating fraction files\n";
      my $file = "$temp_dir/out/$sample_file_basename.$profiling_mode2.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.insert.mm.dist.among.unique.scaled.species";
      if ( -e $file )
      {
        system_("$scr_dir/MOCATFraction.pl -in $file -out $file.fraction");
        @input_files1 = ( "$file.fraction",                                                                                                                                                                                                                                            $file );
        @input_files2 = ( "$OUTFOLDER/$sample_file_basename.$profiling_mode2.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}$bwaAddon.insert.mm.dist.among.unique.scaled.species.fraction", "$OUTFOLDER/$sample_file_basename.$profiling_mode2.profile.$read_type.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$conf{filter_length_cutoff}.p$conf{filter_percent_cutoff}.insert.mm.dist.among.unique.scaled.species" );
        @output_files = ( 'NCBI.species.abundances',                                                                                                                                                                                                                                   'NCBI.species.counts' );
      }
    }

    if ( scalar @input_files1 > 0 )
    {
      for my $i ( 0 .. scalar @input_files1 - 1 )
      {
        unless ( $input_files1[$i] eq '' )
        {
          system_("cp $input_files1[$i] $input_files2[$i] && ln -sf $input_files2[$i] $OUTPUT_FOLDER/$output_files[$i]");
        }
      }
    }
    print localtime() . ": removing $temp_dir\n";

    system_("rm -r $temp_dir");
  }

  print localtime() . ": done\n";
  if ( $assembly == 1 )
  {
    $OUTFOLDER = "$cwd/PROFILES/$profiling_mode_folder.profiles/$sample_file_basename";
  }
  print localtime() . ": RESULTS SAVED IN $OUTFOLDER\n";
}

1;
