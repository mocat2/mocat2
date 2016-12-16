package MOCATScreen;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;
use File::Basename;

# Reverted to old screen implmentation

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
  my $mapping_mode;
  my $screen_source;
  my $screen_save;
  my $basename;
  my $temp_dir      = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
  my $read_type     = 'screened';
  my $assembly_type = 'assembly';
  my $end;

  if ($use_extracted_reads)
  {
    $read_type = 'extracted';
  }
  if ( $conf{MOCAT_mapping_mode} eq 'unique' )
  {
    $mapping_mode = "-r 0";
  } elsif ( $conf{MOCAT_mapping_mode} eq 'random' )
  {
    $mapping_mode = "-r 1";
  } elsif ( $conf{MOCAT_mapping_mode} eq 'allbest' )
  {
    $mapping_mode = "-r 2";
  } else
  {
    die "ERROR & EXIT: Unknown MOCAT_mapping_mode";
  }
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
      if ( -e "$data_dir/$screen.index.sai" )
      {
        print localtime() . ": DATABASE status: Found $screen.index* in $data_dir = ALL OK!\n";
        print localtime() . ": Continuing creating $job jobs...\n";
      } else
      {
        if ( -e "$data_dir/$screen" )
        {
          print localtime() . ": DATABASE status: Found $screen, but not $screen.index* in $data_dir\n";
          die "\nERROR & EXIT: $data_dir/$screen.index* does not exist, please index the database by running: $bin_dir/2bwt-builder $data_dir/$screen";

          #					print localtime() . ": DATABASE status: Indexing database now...";
          #					system "$bin_dir/2bwt-builder $data_dir/$screen >/dev/null";
          #					if ( -e "$data_dir/$screen.index.sai" ) {
          #						print " OK!\n";
          #						print localtime() . ": Continuing creating $job jobs...\n";
          #					}
          #					else {
          #						die "\nERROR & EXIT: $data_dir/$screen.index* does not exist, and could not be created.";
          #					}
        } elsif ( -e "$screen.index.sai" )
        {
          print localtime() . ": DATABASE status: Found external database $screen, but it wasn't imported into $data_dir.\n";
          print localtime() . ": DATABASE status: Importing $screen into $data_dir by symlinks...";
          chomp( my $base = `basename $screen` );
          my $dirn = dirname($screen);
          system "ln -fs $dirn/$base.index.* $data_dir/ ; ln -fs $dirn/$base $data_dir/";
          if ( -e "$data_dir/$base.index.sai" )
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
            die "\nERROR & EXIT: Could not import $screen into $data_dir ($data_dir/$screen.index.sai does not exists, after attempt to create it)";
          }
        } elsif ( -e "$screen" )
        {
          print localtime() . ": DATABASE status: Found external database $screen, but it wasn't indexed.\n";

          #print localtime() . ": DATABASE status: Indexing and importing $screen into $data_dir by symlinks...";
          chomp( my $base = `basename $screen` );
          my $dirn = dirname($screen);
          die "\nERROR & EXIT: $data_dir/$screen.index* does not exist, please index and link the database by running: $bin_dir/2bwt-builder $screen && ln -fs $dirn/$base.index.* $data_dir/ && ln -fs $dirn/$base $data_dir/";

          #					system "$bin_dir/2bwt-builder $screen >/dev/null && ln -fs $dirn/$base.index.* $data_dir/ ; ln -fs $dirn/$base $data_dir/";
          #					if ( -e "$data_dir/$base.index.sai" ) {
          #						print " OK!\n";
          #						print "\n\nNEXT TIME: Specify the this database only as '$base'\n\n\n";
          #
          #						# I got tired of pressing enter each time...
          #						#print "Press <enter> to continue.";
          #						#chomp( my $key = <STDIN> );
          #						print "\n" . localtime() . ": Continuing creating $job jobs...\n";
          #						$screen = $base;
          #						$screen[$counter] = $screen;
          #					}
          #					else {
          #						die "\nERROR & EXIT: Could not index and import $screen into $data_dir";
          #					}
        } else
        {
          die "\nERROR & EXIT: $screen is not a valid database path";
        }
      }
    }
  }    # End DB check

  # Fasta file check
  if ($SCREEN_FASTA_FILE)
  {
    if ( $screen[0] =~ m/\// )
    {
      unless ( $screen[0] =~ m/^\// )
      {
        die "ERROR & EXIT: The path you specify to the fasta file have to be either of:\n1. full path, eg. /usr/home/files/fastafile/fna\n2. single name, then MOCAT will look into your\n   data folder ($data_dir) and\n   current path ($cwd)\n   for the file.\nNOTE! FIRST THE DATA FOLDER WILL BE SEARCHED for THE FILE AND THE THE CURRENT FOLDER!";
      }
    }
    chomp( $basename = `basename $screen[0]` );
    if ( $screen[0] =~ m/^\// )
    {
      unless ( -e $screen[0] )
      {
        die "ERROR & EXIT: Specified fasta file does not exist in data folder, current folder or in the specified path.";
      }
      print localtime() . ": File used: $screen[0]\n";
      unless ( "$screen[0]" eq "$data_dir/$basename" )
      {
        print localtime() . ": THIS FASTA FILE IS IMPORTED INTO YOU DATA FOLDER BY A LINK!\n";
        system "ln -sf $screen[0] $data_dir";
      }
      print "\n\nNEXT TIME: Specify the this fasta file only as '$basename'\n\n\n";

      #print "Press <enter> to continue.";
      #chomp( my $key = <STDIN> );
      print "\n" . localtime() . ": Continuing creating $job jobs...\n";

    } else
    {
      if ( -e "$data_dir/$screen[0]" && -e "$cwd/$screen[0]" )
      {
        print localtime() . ": WARNING! FOUND FASTA FILE IN DATA FOLDER AND CURRENT FOLDER!\n";
        print localtime() . ": Fasta file used: $data_dir/$screen[0]\n";
      } elsif ( -e "$data_dir/$screen[0]" )
      {
        print localtime() . ": FOUND FASTA FILE IN DATA FOLDER:\n";
        print localtime() . ": File used: $data_dir/$screen[0]\n";
      } elsif ( -e "$cwd/$screen[0]" )
      {
        print localtime() . ": FOUND FASTA FILE IN CURRENT FOLDER\n";
        print localtime() . ": THIS FASTA FILE IS IMPORTED INTO YOU DATA FOLDER BY A LINK!\n";
        print localtime() . ": File used: $data_dir/$screen[0]\n";
        system "ln -sf $cwd/$screen[0] $data_dir";
        unless ( -e "$data_dir/$screen[0]" )
        {
          die "ERROR & EXIT: Import FAILED! Please manually copy the $screen[0] file into $data_dir";
        }
        print "\n\nNEXT TIME: Specify the this fasta file only as '$basename'\n\n\n";

        # I got tired of pressing enter each time...
        #print "Press <enter> to continue.";
        #chomp( my $key = <STDIN> );
        print "\n" . localtime() . ": Continuing creating $job jobs...\n";
      } else
      {
        die "ERROR & EXIT: Specified fasta file does not exist in data folder, current folder or in the specified path.";
      }
      $screen[0] = "$data_dir/$screen[0]";
    }

  }

  # Loop over samples, and create jobs
  foreach my $sample (@samples)
  {
    my $LOG       = "2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log";
    my $LOG2      = "2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log";
    my $inputfile = "$temp_dir/$sample/temp/SDB_INPUT.$date.gz";

    print JOB "exec 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log && ";
    print JOB "mkdir -p $temp_dir/$sample/temp; ";

    # If not Fasta file, save input files into a temp file
    if ( !$SCREEN_FASTA_FILE )
    {

      unless ($only_regenerate_reads)
      {
        chomp( our $systemType = `uname -s` );
        if ( $systemType =~ m/Darwin/ )
        {
          print JOB "$ZCAT $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/*pair*gz $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/*single*gz > $inputfile && ";
        } else
        {
          print JOB "cat $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/*pair*gz $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/*single*gz > $inputfile && ";
        }
      }
    }

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
      my $sf          = "$cwd/$sample/reads.screened.$db_on_db.$conf{MOCAT_data_type}";
      my $ef          = "$cwd/$sample/reads.extracted.$db_on_db.$conf{MOCAT_data_type}";
      my $mf          = "$cwd/$sample/reads.mapped.$db_on_db.$conf{MOCAT_data_type}";
      my $file_output = "$mf/$sample.mapped.$screen_save.on.$screen.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}";
      my $database    = "$data_dir/$screen";
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
          unless ( -e "$database.index.sai" )
          {
            print JOB " && $bin_dir/2bwt-builder $cwd/$sample/$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly_type.$reads.$conf{MOCAT_data_type}.K$kmer.$end $LOG ";
          }
        }
      }

      # Get lane IDs
      my @lanes = `ls -1 $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/*pair.1.fq.gz`;
      foreach my $i ( 0 .. ( scalar @lanes - 1 ) )
      {
        chomp( $lanes[$i] );
        my @tmp = split /\//, $lanes[$i];
        $lanes[$i] = $tmp[-1];
        $lanes[$i] =~ s/.pair.1.fq.gz//;
      }

      # Get lanes in the sample folder, for fasta file screen
      %use3files = ();
      my $end;
      

      # make directory
      print JOB " && mkdir -p $mf ";

      # Make folders
      if ($screened_files)
      {
        print JOB " && mkdir -p $sf ";
      }
      if ($extracted_files)
      {
        print JOB " && mkdir -p $ef ";
      }

      # Usearch, if fasta file
      my @fqs;
      if ($SCREEN_FASTA_FILE)
      {
        @fqs = @lanes;
        print JOB " && echo 'FASTA FILE USED IN SCREEN: $screen[0]' >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log ";
        if ( scalar @fqs == 0 )
        {
          die "ERROR & EXIT: There is no RAW DATA in the sample folder $sample!";
        }
        print JOB " && rm -f $ids_file";    # Because read ids are appended to this file below, we need to first delete the file
        foreach my $lane (@fqs)
        {
          my $fa      = "$temp_dir/$sample/temp/$lane.tmp";
          my $blast   = "$mf/$lane.usearch";
          my $logfile = "$mf/$lane.usearch.log";
          my $stats   = "$cwd/$sample/stats/$lane.usearch.on.$basename.$conf{MOCAT_data_type}.total.stats";
          my $stats2  = "$cwd/$sample/stats/$lane.usearch.on.$basename.$conf{MOCAT_data_type}.match.stats";

          unless ($only_regenerate_reads)
          {
            print JOB " && gunzip -c $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/$lane*fq.gz | $scr_dir/MOCATScreenFastaFile_aux.pl fq2fa > $fa";

            #<(gunzip -c $cwd/$sample/$screen_source.$conf{MOCAT_data_type}/$lane*fq.gz | $scr_dir/MOCATScreenFastaFile_aux.pl fq2fa)

            if ( $conf{screen_fasta_file_usearch_version} eq '5' )
            {
              print JOB " && $conf{MOCAT_dir}/ext/usearch/$conf{screen_fasta_file_usearch_version_5_exe} $conf{screen_fasta_file_additional_usearch_cmd} --query $fa --minlen $conf{screen_fasta_file_blast_read_min_length} --db $screen --blast6out $blast --evalue $conf{screen_fasta_file_blast_e_value} --log $logfile $LOG";
            } elsif ( $conf{screen_fasta_file_usearch_version} eq '6' )
            {
              print JOB " && $conf{MOCAT_dir}/ext/usearch/$conf{screen_fasta_file_usearch_version_6_exe} $conf{screen_fasta_file_additional_usearch_cmd} --threads $processors --strand both --maxaccepts 0 --ublast $fa --db $screen --blast6out $blast --evalue $conf{screen_fasta_file_blast_e_value} --log $logfile $LOG";
            } else
            {
              die "ERROR & EXIT: Invalid screen_fasta_file_usearch_version in config. It should be 5 or 6.";
            }
            print JOB " && faEntries=\`grep -c \'>\' $fa\` && totHits=\`wc -l $blast | cut -f 1 -d\" \"\` && uniqHits=\`cut -f 1 $blast | sort -u | wc -l\`";
            print JOB " && echo -e \"total_reads\\ttotal_hits\\tunique_hits\\n\$faEntries\\t\$totHits\\t\$uniqHits\" > $stats";
          }

          # in MOCAT v2.2+; we're loading not from the raw FastQ files but rather from the processed ones, and dont need the conversion
          print JOB " && rm $fa $LOG && cut -f 1 $blast  | sed 's|/[12]\$||' >> $ids_file ";
          print JOB " && perl -lane \'BEGIN{sub fun {\$h{\$b}<=>\$h{\$a}}}; \$h{\$F[1]}++; END{foreach \$k (sort fun(keys(\%h))){print \"\$k\\t\$h{\$k}\"}}\' $blast > $stats2";
        }    # End, for each lane
      }    # End, if fasta file

      # Mapping, if not fasta file
      else
      {
        unless ($only_regenerate_reads)
        {
          print JOB " && $bin_dir/soap2.21 -a $inputfile -D $database.index $mapping_mode $conf{screen_soap_cmd} -l " . $conf{screen_soap_seed_length} . " -v " . $conf{screen_soap_max_mm} . " -p $processors $LOG2 -o /dev/stdout ";
        }
      }

      # Do filtering if not fasta file, otherwise not
      if ( !$SCREEN_FASTA_FILE )
      {
        unless ($only_regenerate_reads)
        {

          #print JOB "&& cat $temp_soap_output";
        }
        if ($only_regenerate_reads)
        {
          print JOB "&& zcat $file_output.soap.gz ";
        }
        if ( -e "$database.coord" )
        {
          print JOB " | perl $scr_dir/MOCATFilter_remove_in_padded.pl -db $database 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log";
        }

        # OLD # print JOB " | perl -F\"\\t\" -lane '\$len=\$F[5];\$mm = \$F[-1] =~ tr/[a-zA-Z]/[a-zA-Z]/;\$as = 100-(\$mm/\$len)*100; if (\$as >= $conf{screen_percent_cutoff} && \$len >= $conf{screen_length_cutoff}){\$F[0] =~ /(.+)\\/[12]/;\$ins = \$1;\$inserts{\$ins}++; print STDERR \"\$_\"; }END{foreach my \$i (keys %inserts){print \"\$i\/1\\n\$i\/2\"}}' >> $ids_file 2> $file_output.soap.tmp2 && mv $file_output.soap.tmp2 $file_output.soap && rm $temp_soap_output";
        print JOB " | perl -F\"\\t\" -lane '\$len=\$F[5];\$mm = \$F[-1] =~ tr/[a-zA-Z]/[a-zA-Z]/;\$as = 100-(\$mm/\$len)*100; if (\$as >= $conf{screen_percent_cutoff} && \$len >= $conf{screen_length_cutoff}){\$F[0] =~ /(.+)\\/[12]/;\$ins = \$1;\$inserts{\$ins}++; print STDERR \"\$_\"; print STDOUT \"\$ins\"}' >> $ids_file  2> $file_output.soap.tmp2 && mv $file_output.soap.tmp2 $file_output.soap ";
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
          if ($SCREEN_FASTA_FILE)
          {
            print JOB " $sf/$lane.screened.$basename ";
          } else
          {
            print JOB " $sf/$lane.screened.$screen ";
          }
        }
      }
      if ($extracted_files)
      {
        print JOB " -ex ";
        foreach my $lane (@lanes)
        {
          if ($SCREEN_FASTA_FILE)
          {
            print JOB " $ef/$lane.extracted.$basename ";
          } else
          {
            print JOB " $ef/$lane.extracted.$screen ";
          }
        }
      }

      print JOB " -toremove $ids_file -stats $stats_file -estats $estats_file -identity $conf{screen_percent_cutoff} -length $conf{screen_length_cutoff} -soapmaxmm $conf{screen_soap_max_mm} $LOG ";

      # If saving in SAM format
      unless ($only_regenerate_reads)
      {
        if ( $conf{screen_save_format} eq 'sam' && !$SCREEN_FASTA_FILE )
        {

          #print JOB " && $scr_dir/MOCATExternal_soap2sam.pl '$ZIP -$ziplevel -c -f' '$file_output.sam' < $file_output.soap $LOG ";
          print JOB " && $scr_dir/MOCATExternal_soap2sam.pl < $file_output.soap | $ZIP -$ziplevel -c -f > $file_output.sam.gz $LOG2 ";
        }
      }

      # Zip soap file
      unless ($only_regenerate_reads)
      {
        if ( !$SCREEN_FASTA_FILE )
        {
          print JOB " && $ZIP -$ziplevel -f $file_output.soap $LOG";
        }
      }

      # Remove ids file, it can get really big
      print JOB " && rm -f $ids_file";

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

  # Fasta file specific
  if ($SCREEN_FASTA_FILE)
  {
    unless ( -s "$conf{MOCAT_dir}/ext/usearch/usearch" )
    {
      die "ERROR & EXIT: Missing external program 'usearch' in $conf{MOCAT_dir}/ext/usearch/usearch.";
      print "1. Please Download Usearch from http://www.drive5.com/usearch/nonprofit_form.html\n";
      print "2. Rename the downloaded file to 'usearch'\n";
      print "3. Copy it to $conf{MOCAT_dir}/ext/usearch\n";
      print "4. Run chmod u+x $conf{MOCAT_dir}/ext/usearch/usearch to make it executable!\n";
    }
  }

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
        $sf          = "$cwd/$sample/reads.screened.$end.$assembly_type.K$kmer.$addon$conf{MOCAT_data_type}";
        $ef          = "$cwd/$sample/reads.extracted.$end.$assembly_type.K$kmer.$addon$conf{MOCAT_data_type}";
        $mf          = "$cwd/$sample/reads.mapped.$end.$assembly_type.K$kmer.$addon$conf{MOCAT_data_type}";
        $file_output = "$mf/$sample.mapped.$reads.on.$end.$assembly_type.K$kmer.$addon$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.soap.gz";
      } else
      {
        $sf          = "$cwd/$sample/reads.screened.$db_on_db.$conf{MOCAT_data_type}";
        $ef          = "$cwd/$sample/reads.extracted.$db_on_db.$conf{MOCAT_data_type}";
        $mf          = "$cwd/$sample/reads.mapped.$db_on_db.$conf{MOCAT_data_type}";
        $file_output = "$mf/$sample.mapped.$screen_before.on.$screen.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.soap.gz";
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
