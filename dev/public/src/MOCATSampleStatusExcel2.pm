package MOCATSampleStatusExcel2;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;
use Excel::Writer::XLSX;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

sub run
{

############## EXCEL SUMMARY CODE ###############

  my $original_sample_file = $sample_file;
  if ( $original_sample_file =~ m/\// )
  {
    my @a = split '/', $original_sample_file;
    $original_sample_file = $a[-1];
  }
  my $sample_file        = "$cwd/SUMMARIES/$original_sample_file/$original_sample_file";
  my $sample_file_folder = "$cwd/SUMMARIES/$original_sample_file/";
  system "rm -f $sample_file.summary.xlsx";

  my $WB         = Excel::Writer::XLSX->new("$sample_file.summary.xlsx");
  my $WS_support = $WB->add_worksheet('Support');
  $WS_support->add_write_handler( qr[\w], \&store_string_widths );

  $WS_support->write( 0,  0, "Color codes" );
  $WS_support->write( 1,  0, "Read QC Statistics" );
  $WS_support->write( 1,  1, "Blue indicates raw seq depth. Dark blue is deeper." );
  $WS_support->write( 2,  1, "Red, yellow, green indicates whether above cutoff 1 (green), b/w cutoff 1 and 2 (yellow) or below cutoff 2." );
  $WS_support->write( 3,  1, "If column name is 'trimmed': cutoff 1: 70%, cutoff 2: 50% left from previous step (raw)" );
  $WS_support->write( 4,  1, "If column name is 'adapter screened' OR 'Adp_seq.fasta': cutoff 1: 90%, cutoff 2: 80% left from previous step (trimmed normally)" );
  $WS_support->write( 5,  1, "If column name is 'hg19 screened'   : cutoff 1: 97%, cutoff 2: 90% left from previous step (adapter normally)" );
  $WS_support->write( 6,  1, "If column name is other    : cutoff 1: 80%, cutoff 2: 50% left from previous step (unknown)" );
  $WS_support->write( 8,  0, "Revised Assembly Statistics & Gene Prediction Statistics" );
  $WS_support->write( 8,  1, "Blue to white gradient. Blue is higher, white lower." );
  $WS_support->write( 10, 0, "DB X" );
  $WS_support->write( 10, 1, "Green to red gradient. Green is higher, red lower." );

  my @name = ('Read QC Statistics');
  my @map = ( "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "AA", "AB", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AJ" );
  my @assembly_files;
  my $assembly_counter = 0;
  $assembly_files[0] = "";
  chomp( my @revised_assembly_file    = `ls -1 $sample_file.assembly.revised.*.scaftig.stats.summary 2>/dev/null` );
  chomp( my @nonrevised_assembly_file = `ls -1 $sample_file.assembly.*.scaftig.stats.summary 2>/dev/null | grep -v '$sample_file.assembly.revised.'` );
  unless ( $nonrevised_assembly_file[0] )
  {
    $nonrevised_assembly_file[0] = '';
  }
  unless ( $revised_assembly_file[0] )
  {
    $revised_assembly_file[0] = '';
  }

  if ( $conf{sample_status_force_cfg_hq_reads} eq 'no' )
  {
    if ( scalar @nonrevised_assembly_file > 1 || scalar @revised_assembly_file > 1 )
    {
      print STDOUT "WARNING: Found multiple assemblies (either revised or non-revised), defaulting to HQ reads defined as sample_status_default_hq_reads=$conf{sample_status_default_hq_reads} in the config file.\n";
      if ( $conf{sample_status_default_hq_reads} )
      {
        $assembly_files[0] = "$sample_file.assembly.$conf{sample_status_default_hq_reads}.$conf{MOCAT_data_type}.scaftig.stats.summary";
      }
    } else
    {
      unless ( $nonrevised_assembly_file[0] eq '' )
      {
        $assembly_files[$assembly_counter] = $nonrevised_assembly_file[0];
        $assembly_counter++;
        push( @name, 'Non-revised Assembly Statistics' );
      }
      unless ( $revised_assembly_file[0] eq '' )
      {
        $assembly_files[$assembly_counter] = $revised_assembly_file[0];
        push( @name, 'Revised Assembly Statistics' );
      }

      if ( $assembly_files[0] eq '' )
      {
        $assembly_files[0] = "$sample_file.assembly.unassembled.$conf{MOCAT_data_type}.scaftig.stats.summary";

        if ( $conf{sample_status_default_hq_reads} )
        {
          $assembly_files[0] = "$sample_file.assembly.$conf{sample_status_default_hq_reads}.$conf{MOCAT_data_type}.scaftig.stats.summary";
        }

      }
    }
  } else
  {
 if ( $conf{sample_status_default_hq_reads} )
      {
        $assembly_files[0] = "$sample_file.assembly.$conf{sample_status_default_hq_reads}.$conf{MOCAT_data_type}.scaftig.stats.summary";
      }
      else {
       die "ERROR & EXIT: Config 'sample_status_force_cfg_hq_reads' set not set to no, but sample_status_default_hq_reads not defined"
      }
  }
  my $assembly_file = $assembly_files[0];
  my $CC            = 1;
  print "$sep\n";
  print localtime() . ": Generating second Excel summary";
  my $HQ_name = "";
  my @Trimmed_array;
  my $chartSeries = 2;
  my $assembly_name;
  my $short_assembly_name;

  if ( $assembly_file ne '' )
  {
    $assembly_file =~ m/$sample_file.assembly.(.*).scaftig.stats.summary/;
    $assembly_name       = $1;
    $short_assembly_name = $1;
  } else
  {
    $assembly_name       = "N/A";
    $short_assembly_name = "N/A";
  }
  $short_assembly_name =~ s/revised.//;

  $assembly_name =~ m/(.*).(solexaqa|fastx)/;
  if ($1)
  {
    $HQ_name = ucfirst("$1");
  } else
  {
    $HQ_name = "N/A";
  }
  if ( $short_assembly_name =~ m/(screened|extracted).*.on.(.*).(solexaqa|fastx)/ )
  {
    $HQ_name = ucfirst("$2 $1");
  }

  chomp( my @files = `ls -1 $sample_file.*.summary | grep -v '.xls'` );
  my %R;
  my %C;

  my %checked;

  foreach my $assembly_file (@assembly_files)
  {

    $assembly_file =~ m/$sample_file.assembly.(.*).scaftig.stats.summary/;
    
    my $assembly_name       = $1;
    my $short_assembly_name = $1;
    $short_assembly_name =~ s/revised.//;
    my $revised = "";
    if ( $assembly_name =~ m/revised\./ )
    {
      $revised = 'Revised ';
    }

    my $RTF = 1;
    foreach my $file (@files)
    {
      if ( $file =~ m/$sample_file.screened.$short_assembly_name.summary/ )
      {
        $RTF = 0;
      }
    }

    #die join ("\n", @files);
    my %opened;
    foreach my $file (@files)
    {
      chomp $file;

      # DB coverage file
      unless ( $checked{$file} )
      {
        if ( $file =~ m/$sample_file.coverage.summary/ )
        {
          open IN, $file;

          #print "6 opening $file : \n";
          chomp( my $header = <IN> );
          while (<IN>)
          {
            chomp;
            my @l = split "\t", $_;
            $C{ $l[1] }{ $l[0] }{'Total inserts'}    = $l[2];
            $C{ $l[1] }{ $l[0] }{'Mapped inserts'}   = $l[3];
            $C{ $l[1] }{ $l[0] }{'Unmapped inserts'} = $l[2] - $l[3];
            $C{ $l[1] }{ $l[0] }{'% mapped inserts'} = $l[4];
            $C{ $l[1] }{ $l[0] }{'Total bases'}      = $l[5];
            $C{ $l[1] }{ $l[0] }{'Mapped bases'}     = $l[6];
            $C{ $l[1] }{ $l[0] }{'Unmapped bases'}   = $l[5] - $l[6];
            $C{ $l[1] }{ $l[0] }{'% mapped bases'}   = $l[7];
          }
          close IN;
        }
      }

      # Trimmed data, this could be any number of files... Painful.
      # It's inside this loop we process the LQ reads

      unless ( $opened{$file} )
      {
        unless ( $checked{$file} )
        {
          if ( $file =~ m/$sample_file.(.*).(solexaqa|fastx).summary/ || ( $file =~ m/$sample_file.readtrimfilter.(solexaqa|fastx).summary/ ) )
          {
            my $file_part = $1;

            if ( ( $short_assembly_name =~ m/$file_part/ && $short_assembly_name ne $file_part ) || $file =~ m/$sample_file.readtrimfilter.(solexaqa|fastx).summary/ )
            {

              #print "file_part = $file_part | short_assembly_name == $short_assembly_name :: processing $file\n";

              $chartSeries++;
              if ( $file =~ m/$sample_file.readtrimfilter.(solexaqa|fastx).summary/ )
              {
                $file_part = "Trimmed";
              } else
              {
                $file_part =~ s/screened.(.*)/$1 screened/;
                $file_part =~ s/extracted.(.*)/$1 extracted/;
                $file_part = ucfirst($file_part);
              }

              open IN, $file;

              #print "5 opening $file : \n";
              chomp( my $header = <IN> );
              while (<IN>)
              {
                chomp;
                my @l = split "\t", $_;
                $R{ $l[0] }{"$file_part reads"}   = $l[1];
                $R{ $l[0] }{"$file_part bases"}   = $l[2];
                $R{ $l[0] }{"$file_part inserts"} = $l[6];
              }
              if ( $file =~ m/$sample_file.readtrimfilter.(solexaqa|fastx).summary/ )
              {
                @Trimmed_array = ( "$file_part bases", "$file_part reads", "$file_part inserts", @Trimmed_array );
              } else
              {
                push @Trimmed_array, ( "$file_part bases", "$file_part reads", "$file_part inserts" );
              }
              close IN;
              $opened{$file} = 1;
            }
          }
        }
      }

      # HQ data
      unless ( $checked{$file} )
      {
        if ( $file =~ m/$sample_file.screened.$short_assembly_name.summary/ )
        {

          #print "\n1 opening $file : \n";
          open IN, $file;
          chomp( my $header = <IN> );
          while (<IN>)
          {
            chomp;
            my @l = split "\t", $_;
            $R{ $l[0] }{"$HQ_name reads"}  = $l[1];
            $R{ $l[0] }{"$HQ_name bases"}  = $l[2];
            $R{ $l[0] }{'Max read length'} = $l[3];
            $R{ $l[0] }{'Avg read length'} = $l[4];
            $R{ $l[0] }{'K-mer'}           = $l[5];
            $R{ $l[0] }{'Inserts'}         = $l[6];
          }
          close IN;
        }
      }

      # HQ data
      unless ( $opened{$file} )
      {
        unless ( $checked{$file} )
        {
          if ( $file =~ m/$sample_file.readtrimfilter.(solexaqa|fastx).summary/ )
          {
            if ( $RTF == 1 )
            {
              #print "2 opening $file : \n";
              $opened{$file} = 1;
              open IN, $file;
              chomp( my $header = <IN> );
              while (<IN>)
              {
                chomp;
                my @l = split "\t", $_;
                $R{ $l[0] }{"$HQ_name reads"}  = $l[1];
                $R{ $l[0] }{"$HQ_name bases"}  = $l[2];
                $R{ $l[0] }{'Max read length'} = $l[3];
                $R{ $l[0] }{'Avg read length'} = $l[4];
                $R{ $l[0] }{'K-mer'}           = $l[5];
                $R{ $l[0] }{'Inserts'}         = $l[6];
              }
              close IN;
            }
          }
        }
      }

      # Raw data
      unless ( $checked{$file} )
      {
        if ( $file =~ m/$sample_file.raw.reads.lanes.summary/ )
        {
          open IN, $file;

          #print "3 opening $file : \n";
          chomp( my $header = <IN> );
          while (<IN>)
          {
            chomp;
            my @l = split "\t", $_;
            unless ( $R{ $l[0] }{'Raw reads'} )
            {
              $R{ $l[0] }{'Raw reads'} = 0;
            }
            unless ( $R{ $l[0] }{'Raw bases'} )
            {
              $R{ $l[0] }{'Raw bases'} = 0;
            }
            $R{ $l[0] }{'Raw reads'} = $R{ $l[0] }{'Raw reads'} + $l[2];
            $R{ $l[0] }{'Raw bases'} = $R{ $l[0] }{'Raw bases'} + $l[3];
          }
          close IN;
        }
      }

      # Insert sizes
      if ( $file =~ m/$sample_file.assembly.$short_assembly_name.inserts.summary/ )
      {
        open IN, $file;
        chomp( my $header = <IN> );
        while (<IN>)
        {
          chomp;
          my @l = split "\t", $_;
          unless ( $R{ $l[0] }{"${revised}Avg insert size"} )
          {
            $R{ $l[0] }{"${revised}Avg insert size"} = 0;
          }
          unless ( $R{ $l[0] }{"${revised}number_insert"} )
          {
            $R{ $l[0] }{"${revised}number_insert"} = 0;
          }
          $R{ $l[0] }{"${revised}Avg insert size"} = $R{ $l[0] }{"${revised}Avg insert size"} + $l[2];
          $R{ $l[0] }{"${revised}number_insert"}++;
        }
        close IN;
      }

      # Assembly data
      if ( $file =~ m/$sample_file.assembly.$assembly_name.scaftig.stats.summary/ )
      {
        open IN, $file;
        chomp( my $header = <IN> );
        while (<IN>)
        {
          if (m/\t500$/)
          {
            chomp;
            my @l = split "\t", $_;
            $R{ $l[0] }{"${revised}Scaftigs"}     = $l[1];
            $R{ $l[0] }{"${revised}Total length"} = $l[2];
            $R{ $l[0] }{"${revised}N50"}          = $l[3];
            $R{ $l[0] }{"${revised}N90"}          = $l[4];
            $R{ $l[0] }{"${revised}Max length"}   = $l[5];
            $R{ $l[0] }{"${revised}Min length"}   = $l[6];
          }
        }
        close IN;
      }

      # gene prediction data
      if ( $file =~ m/$sample_file.gene.prediction.assembly.$assembly_name(.*).summary/ )
      {
        my $gp = $1;
        $gp =~ s/^\.//;
        $gp =~ s/MetaGeneMark/MGM/;
        $gp =~ s/Prodigal/P/;
        push @name, "${revised}GP $gp";
        open IN, $file;
        chomp( my $header = <IN> );
        while (<IN>)
        {
          chomp;
          my @l = split "\t", $_;
          $R{"$l[0]"}{"$revised$gp Genes"}                   = $l[1];
          $R{"$l[0]"}{"$revised$gp Complete genes"}          = $l[2];
          $R{"$l[0]"}{"$revised$gp Incomplete genes"}        = $l[3];
          $R{"$l[0]"}{"$revised$gp Total gene length"}       = $l[6];
          $R{"$l[0]"}{"$revised$gp Length complete genes"}   = $l[7];
          $R{"$l[0]"}{"$revised$gp Length incomplete genes"} = $l[8];
        }
        close IN;
      }
      $checked{$file} = 1;
    }

  }    # End foreach assembly file

  # Header
  my $heading = $WB->add_format( bold => 1, align => 'center' );
  $heading->set_bottom(2);
  $heading->set_left(2);
  $heading->set_right(2);

  # Empty cell 0,0
  my $empty = $WB->add_format();
  $empty->set_bottom(2);
  $empty->set_right(2);
  $empty->set_left(2);

  # Sample
  my $sample = $WB->add_format( bold => 1 );
  $sample->set_bottom(2);
  $sample->set_right(2);
  $sample->set_left(2);

  # Var
  my $Var = $WB->add_format( bold => 1, align => 'center' );
  $Var->set_bottom(2);

  # Sum
  my $Sum = $WB->add_format();
  $Sum->set_right(2);
  $Sum->set_left(2);

  # Left
  my $left = $WB->add_format();
  $left->set_right(2);

  # Top
  my $top = $WB->add_format();
  $top->set_top(2);

  #Numbers
  my $numbers = $WB->add_format( align => 'right' );
  $numbers->set_num_format('#,##');

  #NumbersC
  my $numbersC = $WB->add_format( align => 'center' );
  $numbersC->set_num_format('#,##');
  $numbersC->set_bottom(2);

  #NumbersC2
  my $numbersC2 = $WB->add_format( align => 'center' );
  $numbersC2->set_num_format('#,##');

  #Numbers with line
  my $numbersB = $WB->add_format( align => 'right' );
  $numbersB->set_num_format('#,##');
  $numbersB->set_bottom(2);

  for ( my $k = 0 ; $k <= scalar @name - 1 ; $k++ )
  {
    $name[$k] = substr $name[$k], 0, 31;
    my $WS = $WB->add_worksheet( $name[$k] );
    $WS->add_write_handler( qr[\w], \&store_string_widths );
    $WS->set_zoom(200);
    $WS->freeze_panes('B3');
    my $i = 3;
    my @VAR;

    if ( $name[$k] eq 'Read QC Statistics' )
    {
      $WS->activate();
      @VAR = ( "Raw bases", "Raw reads" );
      if ( scalar @Trimmed_array > 0 )
      {
        push @VAR, @Trimmed_array;
      }
      push @VAR, ( "$HQ_name bases", "$HQ_name reads", "Inserts", "Max read length", "Avg read length" );
    }

    if ( $name[$k] eq 'Non-revised Assembly Statistics' )
    {
      @VAR = ( "K-mer", "Avg insert size", "Scaftigs", "Total length", "Max length", "Min length", "N50", "N90" );
    }
    if ( $name[$k] eq 'Revised Assembly Statistics' )
    {
      @VAR = ( "Revised Avg insert size", "Revised Scaftigs", "Revised Total length", "Revised Max length", "Revised Min length", "Revised N50", "Revised N90" );
    }

    if ( $name[$k] =~ /^Revised GP/ )
    {
      $name[$k] =~ m/.*GP (.*)/;
      my $name = $1;
      @VAR = ( "Revised $name Genes", "Revised $name Complete genes", "Revised $name Incomplete genes", "Revised $name Total gene length", "Revised $name Length complete genes", "Revised $name Length incomplete genes" );
    }
    if ( $name[$k] =~ /^GP/ )
    {
      $name[$k] =~ m/.*GP (.*)/;
      my $name = $1;
      @VAR = ( "$name Genes", "$name Complete genes", "$name Incomplete genes", "$name Total gene length", "$name Length complete genes", "$name Length incomplete genes" );
    }

    $WS->print_area( 0, 0, scalar @samples + 3, scalar @VAR );

    $WS->merge_range( 0, 1, 0, scalar @VAR, $name[$k], $heading );
    $WS->write( 0, 0, "",        $empty );
    $WS->write( 1, 0, "Sample",  $sample );
    $WS->write( 2, 0, "Sum",     $Sum );
    $WS->write( 3, 0, "Average", $empty );

    $i = 3;
    my %samples;

    #foreach my $sample (@samples) {
    #	$samples{$sample} = 1;
    #}
    foreach my $sample ( sort keys %R )
    {
      $i++;
      $WS->write( $i, 0, $sample, $Sum );
    }
    for ( my $t = 0 ; $t <= scalar @VAR ; $t++ )
    {
      $WS->write( scalar @samples + 4, $t, "", $top );
    }

    my $j = 0;
    foreach my $var (@VAR)
    {
      print ".";
      $j++;
      $i = 3;

      if ( $j == scalar @VAR )
      {

        # Var
        $Var = $WB->add_format( bold => 1, align => 'center' );
        $Var->set_bottom(2);
        $Var->set_right(2);

        #Numbers
        $numbers = $WB->add_format( align => 'right' );
        $numbers->set_num_format('#,##');
        $numbers->set_right(2);

        #NumbersC
        $numbersC = $WB->add_format( align => 'center' );
        $numbersC->set_num_format('#,##');
        $numbersC->set_bottom(2);
        $numbersC->set_right(2);

        #NumbersC2
        $numbersC2 = $WB->add_format( align => 'center' );
        $numbersC2->set_num_format('#,##');
        $numbersC2->set_right(2);

        #Numbers with line
        $numbersB = $WB->add_format( align => 'right' );
        $numbersB->set_num_format('#,##');
        $numbersB->set_bottom(2);
        $numbersB->set_right(2);
      } else
      {

        # Var
        $Var = $WB->add_format( bold => 1, align => 'center' );
        $Var->set_bottom(2);

        #Numbers
        $numbers = $WB->add_format( align => 'right' );
        $numbers->set_num_format('#,##');

        #NumbersC
        $numbersC = $WB->add_format( align => 'center' );
        $numbersC->set_num_format('#,##');
        $numbersC->set_bottom(2);

        #NumbersC2
        $numbersC2 = $WB->add_format( align => 'center' );
        $numbersC2->set_num_format('#,##');

        #Numbers with line
        $numbersB = $WB->add_format( align => 'right' );
        $numbersB->set_num_format('#,##');
        $numbersB->set_bottom(2);

      }

      my $var_fixed = $var;
      $var_fixed =~ s/(MetaGeneMark|Prodigal)\.\d+\.//;
      $WS->write( 1, $j, $var_fixed, $Var );

      if ( $var eq 'Max read length' || $var eq 'Avg read length' || $var eq 'Max length' || $var eq 'Min length' || $var eq 'N50' || $var eq 'N90' || $var eq 'K-mer' || $var eq 'Avg insert size' || $var eq 'Revised Max read length' || $var eq 'Revised Avg read length' || $var eq 'Revised Max length' || $var eq 'Revised Min length' || $var eq 'Revised N50' || $var eq 'Revised N90' || $var eq 'Revised K-mer' || $var eq 'Revised Avg insert size' )
      {
        $WS->write( 2, $j, "-", $numbersC2 );
      } else
      {
        $WS->write( 2, $j, "=SUM($map[($j)]5:$map[($j)]" . scalar @samples + 4 . ")", $numbers );
      }
      if ( $var eq 'Max read length' || $var eq 'Avg read length' || $var eq 'K-mer' || $var eq 'Avg insert size' || $var eq 'Min length' )
      {
        $WS->write( 3, $j, "=AVERAGE($map[($j)]5:$map[($j)]" . scalar @samples + 4 . ")", $numbersC );
      } else
      {
        $WS->write( 3, $j, "=AVERAGE($map[($j)]5:$map[($j)]" . scalar @samples + 4 . ")", $numbersB );
      }

      # THIS SECTION IS REMOVED BECAUSE THE GENERATED EXPRESSION IS TOO BIG IF MORE THAN 161 SAMPLES. IT IS ALSO NOT USED.
      #				# Support for conditional formatting below, written to special WS
      #				if ( $var =~ m/^Trimmed bases$/ || $var =~ m/^Trimmed reads$/ || $var =~ m/ (screened|extracted) reads$/ || $var =~ m/ (screened|extracted) bases$/ ) {
      #					my $ii = $i + 1;
      #					my $jj;
      #					my $min;
      #					my $max;
      #					if ( $var =~ m/^Trimmed bases$/ || $var =~ m/^Trimmed reads$/ ) {
      #						$jj = $j - 2;
      #					}
      #					else {
      #						$jj = $j - 3;
      #					}
      #
      #					my $expr = "";
      #					for ( my $o = 5 ; $o < scalar @samples + 4 ; $o++ ) {
      #						$expr = $expr . "'Read QC Statistics'!$map[$j]$o/'Read QC Statistics'!$map[$jj]$o,";
      #					}
      #					$expr = $expr . "'Read QC Statistics'!$map[$j]" . ( scalar @samples + 4 ) . "/'Read QC Statistics'!$map[$jj]" . ( scalar @samples + 4 );
      #					$WS_support->write( 1, $j, "STDEV of ($map[$j]/$map[$jj]) from Read QC Statistics" );
      #					$WS_support->write( 2, $j, "=STDEV($expr)" );
      #
      #					$WS_support->write( 4, $j, "AVERAGE of ($map[$j]/$map[$jj]) from Read QC Statistics" );
      #					$WS_support->write( 5, $j, "=AVERAGE($expr)" );
      #				}

      foreach my $kk ( sort keys %R )
      {
        $i++;
        if ( $R{$kk}{$var} )
        {
          if ( $var eq 'Avg insert size' )
          {
            $WS->write_number( $i, $j, $R{$kk}{$var} / $R{$kk}{'number_insert'}, $numbersC2 );
          } elsif ( $var eq 'Max read length' || $var eq 'Avg read length' || $var eq 'K-mer' || $var eq 'Min length' )
          {
            $WS->write_number( $i, $j, $R{$kk}{$var}, $numbersC2 );
          } elsif ( $var eq 'Revised Avg insert size' )
          {
            $WS->write_number( $i, $j, $R{$kk}{$var} / $R{$kk}{'Revised number_insert'}, $numbersC2 );
          } elsif ( $var eq 'Revised Max read length' || $var eq 'Revised Avg read length' || $var eq 'Revised Min length' )
          {
            $WS->write_number( $i, $j, $R{$kk}{$var}, $numbersC2 );
          } else
          {
            $WS->write_number( $i, $j, $R{$kk}{$var}, $numbers );
          }

          if ( $var =~ m/^Trimmed bases$/ || $var =~ m/^Trimmed reads$/ || $var =~ m/ (screened|extracted) reads$/ || $var =~ m/ (screened|extracted) bases$/ || $var =~ m/^\S+ bases$/ || $var =~ m/^\S+ reads/ )
          {

            # Conditional formatting
            # Light red fill with dark red text.
            my $red = $WB->add_format( bg_color => '#ff798a', );

            # Green fill with dark green text.
            my $green = $WB->add_format( bg_color => '#8aff79', );

            # Yellow fill with dark yellow text.
            my $yellow = $WB->add_format(
              bg_color => '#ffee79',

            );

            my $ii = $i + 1;
            my $jj;
            my $min;
            my $max;
            if ( $var =~ m/^Trimmed bases$/ || $var =~ m/^Trimmed reads$/ )
            {
              $jj = $j - 2;
            } else
            {
              $jj = $j - 3;
            }
            if ( $var =~ m/^Trimmed bases$/ || $var =~ m/^Trimmed reads$/ )
            {
              $min = 0.5;
              $max = 0.7;
            } elsif ( $var =~ m/^Adapter screened bases$/ || $var =~ m/^Adapter screened reads$/ )
            {
              $min = 0.8;
              $max = 0.9;
            } elsif ( $var =~ m/^Adp_seq.fasta bases$/ || $var =~ m/^Adp_seq.fasta reads$/ )
            {
              $min = 0.8;
              $max = 0.9;
            } elsif ( $var =~ m/^Hg19 screened bases$/ || $var =~ m/^Hg19 screened reads$/ )
            {
              $min = 0.9;
              $max = 0.97;
            } else
            {
              $min = 0.5;
              $max = 0.8;
            }

            $WS->conditional_formatting(
              "$map[$j]$ii",
              {
                type     => 'formula',
                criteria => "=\$$map[$j]\$$ii/\$$map[$jj]\$$ii >= $max",

                #criteria => "=\$$map[$j]\$$ii/\$$map[$jj]\$$ii > 'Support'!$map[$j]6  - 1*'Support'!$map[$j]3",
                format => $green,
              }
            );
            $WS->conditional_formatting(
              "$map[$j]$ii",
              {
                type     => 'formula',
                criteria => "=\$$map[$j]\$$ii/\$$map[$jj]\$$ii >= $min ",

                #criteria => "=\$$map[$j]\$$ii/\$$map[$jj]\$$ii > 'Support'!$map[$j]6  - 2*'Support'!$map[$j]3",
                format => $yellow,
              }
            );
            $WS->conditional_formatting(
              "$map[$j]$ii",
              {
                type     => 'formula',
                criteria => "=\$$map[$j]\$$ii/\$$map[$jj]\$$ii < $min ",

                #criteria => "=\$$map[$j]\$$ii/\$$map[$jj]\$$ii > 'Support'!$map[$j]6  - 2*'Support'!$map[$j]3",
                format => $red,
              }
            );

          }

        }
      }

      if ( $k == 0 )
      {
        my @A = ( 'B', 'C' );
        foreach my $v (@A)
        {
          $WS->conditional_formatting(
                                       "$v" . "5:$v" . ( scalar @samples + 4 ),
                                       {
                                         type      => '3_color_scale',
                                         min_color => "#d5e4f5",
                                         mid_color => "#2b6bb4",
                                         max_color => "#183c67",
                                       }
          );
        }

      }

      if ( $k == 1 || $k >= 2 )
      {
        my @A;
        if ( $k == 1 )
        {
          @A = ( 'B', 'C', 'D', 'E', 'F', 'H', 'I' );
        } else
        {
          @A = ( 'B', 'C', 'D', 'E', 'F', 'G' );
        }
        foreach my $v (@A)
        {
          $WS->conditional_formatting(
                                       "$v" . "5:$v" . ( scalar @samples + 4 ),
                                       {
                                         type      => '3_color_scale',
                                         min_color => "#d5e4f5",
                                         mid_color => "#2b6bb4",
                                         max_color => "#183c67",
                                       }
          );
        }
      }

      autofit_columns($WS);
      $WS->fit_to_pages( 1, 0 );
    }
  }

  # Process the coverage file
  my $DBc = 0;
  my @VAR = ( 'Total inserts', 'Mapped inserts', 'Unmapped inserts', '% mapped inserts', 'Total bases', 'Mapped bases', 'Unmapped bases', '% mapped bases' );
  foreach my $DB ( sort keys %C )
  {
    $DBc++;
    my $WS = $WB->add_worksheet("DB $DBc");
    $WS->freeze_panes('B3');
    $WS->add_write_handler( qr[\w], \&store_string_widths );
    $WS->set_zoom(200);
    $WS->write( 0, 0, "", $empty );
    $WS->merge_range( 0, 1, 0, scalar @VAR, $DB, $heading );
    $WS->write( 1, 0, "Sample",  $sample );
    $WS->write( 2, 0, "Sum",     $Sum );
    $WS->write( 3, 0, "Average", $empty );

    for ( my $t = 0 ; $t <= scalar @VAR ; $t++ )
    {
      $WS->write( scalar @samples + 4, $t, "", $top );
    }

    my $j = 0;
    foreach my $var (@VAR)
    {
      $j++;
      if ( $j == scalar @VAR )
      {

        # Var
        $Var = $WB->add_format( bold => 1, align => 'center' );
        $Var->set_bottom(2);
        $Var->set_right(2);

        #Numbers with line
        $numbersB = $WB->add_format( align => 'right' );
        $numbersB->set_num_format('#,##');
        $numbersB->set_right(2);

        #NumbersC
        $numbersC = $WB->add_format( align => 'right' );
        $numbersC->set_num_format('#,##');
        $numbersC->set_bottom(2);
        $numbersC->set_right(2);
      } else
      {

        # Var
        $Var = $WB->add_format( bold => 1, align => 'center' );
        $Var->set_bottom(2);

        #Numbers with line
        $numbersB = $WB->add_format( align => 'right' );
        $numbersB->set_num_format('#,##');

        #NumbersC
        $numbersC = $WB->add_format( align => 'right' );
        $numbersC->set_num_format('#,##');
        $numbersC->set_bottom(2);

      }
      $WS->write( 1, $j, $var, $Var );
      if ( $var eq '% mapped inserts' || $var eq '% mapped bases' )
      {
        $numbersB->set_num_format(9);
        $numbersB->set_align('center');
        $numbersC->set_num_format(9);
        $numbersC->set_align('center');
      }
      if ( $var ne '% mapped inserts' && $var ne '% mapped bases' )
      {
        $WS->write( 2, $j, "=SUM($map[($j)]5:$map[($j)]" . scalar @samples + 4 . ")", $numbersB );
      } else
      {
        $WS->write( 2, $j, "-", $numbersB );
      }

      $WS->write( 3, $j, "=AVERAGE($map[($j)]5:$map[($j)]" . scalar @samples + 4 . ")", $numbersC );
    }

    my $row = 3;
    foreach my $sample ( sort keys %{ $C{$DB} } )
    {
      print ".";
      $row++;

      my $col = 0;
      $WS->write( $row, 0, $sample, $Sum );
      foreach my $var (@VAR)
      {
        $col++;
        if ( $col == scalar @VAR )
        {

          # Var
          $Var = $WB->add_format( bold => 1, align => 'center' );
          $Var->set_bottom(2);
          $Var->set_right(2);

          #Numbers
          $numbers = $WB->add_format( align => 'right' );
          $numbers->set_num_format('#,##');
          $numbers->set_right(2);

          #NumbersC
          $numbersC = $WB->add_format( align => 'center' );
          $numbersC->set_num_format('#,##');
          $numbersC->set_bottom(2);
          $numbersC->set_right(2);

          #NumbersC2
          $numbersC2 = $WB->add_format( align => 'center' );
          $numbersC2->set_num_format('#,##');
          $numbersC2->set_right(2);

          #Numbers with line
          $numbersB = $WB->add_format( align => 'right' );
          $numbersB->set_num_format('#,##');
          $numbersB->set_bottom(2);
          $numbersB->set_right(2);
        } else
        {

          # Var
          $Var = $WB->add_format( bold => 1, align => 'center' );
          $Var->set_bottom(2);

          #Numbers
          $numbers = $WB->add_format( align => 'right' );
          $numbers->set_num_format('#,##');

          #NumbersC
          $numbersC = $WB->add_format( align => 'center' );
          $numbersC->set_num_format('#,##');
          $numbersC->set_bottom(2);

          #NumbersC2
          $numbersC2 = $WB->add_format( align => 'center' );
          $numbersC2->set_num_format('#,##');

          #Numbers with line
          $numbersB = $WB->add_format( align => 'right' );
          $numbersB->set_num_format('#,##');
          $numbersB->set_bottom(2);
        }
        if ( $var eq '% mapped inserts' || $var eq '% mapped bases' )
        {
          $numbers->set_num_format(9);
          $numbers->set_align('center');
        }

        $WS->write_number( $row, $col, $C{$DB}{$sample}{$var}, $numbers );
      }

    }

    # Conditional formatting for a DB sheet
    my @A = ( 'B', 'E', 'F', 'I' );
    foreach my $v (@A)
    {
      $WS->conditional_formatting(
                                   "$v" . "5:$v" . ( scalar @samples + 4 ),
                                   {
                                     type      => '3_color_scale',
                                     min_color => "#ff798a",
                                     mid_color => "#ffee79",
                                     max_color => "#8aff79",
                                   }
      );
    }

    #autofit_columns($WS);

    # Coverage charts

    print ".";
    my $chart = $WB->add_chart( type => 'area', embedded => 1, subtype => 'stacked' );
    my $num   = ( scalar @samples + 4 );
    my $name  = "DB$DBc";
    $chart->add_series(
                        name       => "='DB $DBc'!\$C\$2",
                        categories => "='DB $DBc'!\$A\$5:\$A\$$num",
                        values     => "='DB $DBc'!\$C\$5:\$C\$$num",
    );
    $chart->add_series(
                        name       => "='DB $DBc'!\$D\$2",
                        categories => "='DB $DBc'!\$A\$5:\$A\$$num",
                        values     => "='DB $DBc'!\$D\$5:\$D\$$num",
    );
    $chart->set_title( name => 'Mapped & Unmapped Inserts' );
    $chart->set_y_axis( name => 'Number of inserts' );
    $chart->set_legend( position => 'bottom' );
    $CC++;
    $chart->set_style(10);

    print ".";
    $WS->insert_chart( 'J3', $chart, 20, 5 );

    $chart = $WB->add_chart( type => 'area', embedded => 1, subtype => 'stacked' );
    $chart->add_series(
      name       => "='DB $DBc'!\$G\$2",
      categories => "='DB $DBc'!\$A\$5:\$A\$$num",
      values     => "='DB $DBc'!\$G\$5:\$G\$$num",

    );
    $chart->add_series(
                        name       => "='DB $DBc'!\$H\$2",
                        categories => "='DB $DBc'!\$A\$5:\$A\$$num",
                        values     => "='DB $DBc'!\$H\$5:\$H\$$num",
    );
    $chart->set_title( name => 'Mapped & Unmapped Bases' );
    $chart->set_y_axis( name => 'Number of bases' );
    $chart->set_legend( position => 'bottom' );
    $CC++;
    $chart->set_style(10);
    $WS->insert_chart( 'J19', $chart, 20, 5 );

    # Done with coverage file
    $WS->fit_to_pages( 1, 0 );

  }

  # Plots
  my $num = scalar @samples + 4;

  #Create a new chart object. In this case an embedded chart.
  my $WS2    = $WB->add_worksheet('Graphics');
  my $chart  = $WB->add_chart( type => 'column', embedded => 1, name => 'Base Seq Depth' );
  my $chart2 = $WB->add_chart( type => 'column', embedded => 1, name => 'Read Seq Depth' );

  #Configure series
  for ( my $i = 1 ; $i <= $chartSeries ; $i++ )
  {
    print ".";
    my $col;
    my $col2;
    if ( $i == 1 )
    {
      $col  = 1;
      $col2 = 2;
    } else
    {
      $col  = $i * 3 - 3;
      $col2 = $i * 3 - 2;
    }
    $chart->add_series(
      name       => "='Read QC Statistics'!\$$map[($col)]\$2",
      categories => "='Read QC Statistics'!\$A\$5:\$A\$$num",
      values     => "='Read QC Statistics'!\$$map[($col)]\$5\$:$map[($col)]\$$num",

    );
    $chart2->add_series(
                         name       => "='Read QC Statistics'!\$$map[($col2)]\$2",
                         categories => "='Read QC Statistics'!\$A\$5:\$A\$$num",
                         values     => "='Read QC Statistics'!\$$map[($col2)]\$5\$:$map[($col2)]\$$num",
    );
  }

  for ( my $k = 1 ; $k <= scalar @name - 1 ; $k++ )
  {
    if ( $name[$k] eq 'Non-revised Assembly Statistics' || $name[$k] eq 'Revised Assembly Statistics' )
    {
      my $Letter;
      my $Revised;
      if ( $name[$k] eq 'Non-revised Assembly Statistics' )
      {
        $Letter  = "I";
        $Revised = "";
      }
      if ( $name[$k] eq 'Revised Assembly Statistics' )
      {
        $Letter  = "Q";
        $Revised = "Revised ";
      }
      my $j = -14;
      my $chart3 = $WB->add_chart( type => 'column', embedded => 1 );
      $chart3->add_series(
                           name       => "='$name[$k]'!\$C\$2",
                           categories => "='$name[$k]'!\$A\$5:\$A\$$num",
                           values     => "='$name[$k]'!\$C\$5:\$C\$$num",
      );
      $chart3->set_title( name => "${Revised}Average Insert Size" );
      $chart3->set_y_axis( name => "${Revised}Average insert size" );

      my $chart4 = $WB->add_chart( type => 'column', embedded => 1 );
      $chart4->add_series(
                           name       => "='$name[$k]'!\$D\$2",
                           categories => "='$name[$k]'!\$A\$5:\$A\$$num",
                           values     => "='$name[$k]'!\$D\$5:\$D\$$num",
      );
      $chart4->set_title( name => "${Revised}Number of Scaftigs" );
      $chart4->set_y_axis( name => "${Revised}Number of scaftigs" );

      my $chart5 = $WB->add_chart( type => 'column', embedded => 1 );
      $chart5->add_series(
                           name       => "='$name[$k]'!\$E\$2",
                           categories => "='$name[$k]'!\$A\$5:\$A\$$num",
                           values     => "='$name[$k]'!\$E\$5:\$E\$$num",
      );
      $chart5->set_title( name => "${Revised}Total Scaftig Length" );
      $chart5->set_y_axis( name => "${Revised}Total scaftig length" );

      my $chart6 = $WB->add_chart( type => 'column', embedded => 1 );
      $chart6->add_series(
                           name       => "='$name[$k]'!\$F\$2",
                           categories => "='$name[$k]'!\$A\$5:\$A\$$num",
                           values     => "='$name[$k]'!\$F\$5:\$F\$$num",
      );
      $chart6->set_title( name => "${Revised}Max Scaftig Length" );
      $chart6->set_y_axis( name => "${Revised}Max scaftig length" );

      my $chart7 = $WB->add_chart( type => 'column', embedded => 1 );
      $chart7->add_series(
                           name       => "='$name[$k]'!\$H\$2",
                           categories => "='$name[$k]'!\$A\$5:\$A\$$num",
                           values     => "='$name[$k]'!\$H\$5:\$H\$$num",
      );
      $chart7->set_title( name => "${Revised}N50" );
      $chart7->set_y_axis( name => "${Revised}N50" );

      my $chart8 = $WB->add_chart( type => 'column', embedded => 1 );
      $chart8->add_series(
                           name       => "='$name[$k]'!\$I\$2",
                           categories => "='$name[$k]'!\$A\$5:\$A\$$num",
                           values     => "='$name[$k]'!\$I\$5:\$I\$$num",
      );
      $chart8->set_title( name => "${Revised}N90" );
      $chart8->set_y_axis( name => "${Revised}N90" );
      $chart3->set_style(10);
      $chart4->set_style(10);
      $chart5->set_style(10);
      $chart6->set_style(10);
      $chart7->set_style(10);
      $chart8->set_style(10);
      $chart3->set_legend( position => 'bottom' );
      $chart4->set_legend( position => 'bottom' );
      $chart5->set_legend( position => 'bottom' );
      $chart6->set_legend( position => 'bottom' );
      $chart7->set_legend( position => 'bottom' );
      $chart8->set_legend( position => 'bottom' );
      $WS2->insert_chart( "${Letter}1",  $chart3, 20, 5 );
      $WS2->insert_chart( "${Letter}16", $chart4, 20, 5 );
      $WS2->insert_chart( "${Letter}31", $chart5, 20, 5 );
      $WS2->insert_chart( "${Letter}46", $chart6, 20, 5 );
      $WS2->insert_chart( "${Letter}61", $chart7, 20, 5 );
      $WS2->insert_chart( "${Letter}76", $chart8, 20, 5 );
    }

    if ( $name[$k] =~ m/^GP/ || $name[$k] =~ m/Revised GP/ )
    {
      my $Letter;
      if ( $name[$k] =~ m/^GP/ )
      {
        $Letter = "I";
      }
      if ( $name[$k] =~ m/^Revised GP/ )
      {
        $Letter = "Q";
      }
      my $j = 76;
      $name[$k] =~ m/^(.*)GP(.*)/;
      my $N = "$1$2";
      $N =~ s/\s*$//;
      my $chart9 = $WB->add_chart( type => 'area', embedded => 1, subtype => 'stacked' );
      $chart9->add_series(
                           name       => "='$name[$k]'!\$C\$2",
                           categories => "='$name[$k]'!\$A\$5:\$A\$$num",
                           values     => "='$name[$k]'!\$C\$5:\$C\$$num",
      );
      $chart9->add_series(
                           name       => "='$name[$k]'!\$D\$2",
                           categories => "='$name[$k]'!\$A\$5:\$A\$$num",
                           values     => "='$name[$k]'!\$D\$5:\$D\$$num",
      );
      $chart9->set_title( name => "Number of Predicted Genes ($N)" );
      $chart9->set_y_axis( name => 'Predicted genes' );

      my $chart10 = $WB->add_chart( type => 'column', embedded => 1 );
      $chart10->add_series(
                            name       => "='$name[$k]'!\$F\$2",
                            categories => "='$name[$k]'!\$A\$5:\$A\$$num",
                            values     => "='$name[$k]'!\$F\$5:\$F\$$num",
      );
      $chart10->add_series(
                            name       => "='$name[$k]'!\$G\$2",
                            categories => "='$name[$k]'!\$A\$5:\$A\$$num",
                            values     => "='$name[$k]'!\$G\$5:\$G\$$num",
      );
      $chart10->set_title( name => "Total Length of Predicted Genes ($N)" );
      $chart10->set_y_axis( name => 'Total length' );

      $chart9->set_style(10);
      $chart10->set_style(10);
      $chart9->set_legend( position => 'bottom' );
      $chart10->set_legend( position => 'bottom' );
      $j += 15;
      $WS2->insert_chart( "$Letter$j", $chart9, 20, 5 );
      $j += 15;
      $WS2->insert_chart( "$Letter$j", $chart10, 20, 5 );

    }
  }

  # Add a chart title and some axis labels.
  $chart->set_title( name => 'Base Sequencing Depth' );
  $chart->set_y_axis( name => 'Number of bases' );
  $chart2->set_title( name => 'Read Sequencing Depth' );
  $chart2->set_y_axis( name => 'Number of reads' );

  # Insert the chart into the worksheet (with an offset).
  $chart->set_style(10);
  $chart2->set_style(10);

  $chart->set_legend( position => 'bottom' );
  $chart2->set_legend( position => 'bottom' );

  $WS2->insert_chart( 'A1',  $chart,  20, 5 );
  $WS2->insert_chart( 'A16', $chart2, 20, 5 );

  $WS2->fit_to_pages( 1, 0 );
  $WS2->set_zoom(170);

  autofit_columns($WS_support);
  $WS_support->set_zoom(200);
  $WS_support->fit_to_pages( 1, 0 );

  $WB->close();
  print " OK!\n";
}

############### EXCEL SUMMARY CODE ###############

###############################################################################
#
# Functions used for Autofit.
#
###############################################################################

###############################################################################
#
# Adjust the column widths to fit the longest string in the column.
#
sub autofit_columns
{

  my $worksheet = shift;
  my $col       = 0;

  for my $width ( @{ $worksheet->{__col_widths} } )
  {

    $worksheet->set_column( $col, $col, $width ) if $width;
    $col++;
  }
}

###############################################################################
#
# The following function is a callback that was added via add_write_handler()
# above. It modifies the write() function so that it stores the maximum
# unwrapped width of a string in a column.
#
sub store_string_widths
{

  my $worksheet = shift;
  my $col       = $_[1];
  my $token     = $_[2];

  # Ignore some tokens that we aren't interested in.
  return if not defined $token;       # Ignore undefs.
  return if $token eq '';             # Ignore blank cells.
  return if ref $token eq 'ARRAY';    # Ignore array refs.
  return if $token =~ /^=/;           # Ignore formula

  # Ignore numbers
  return if $token =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/;

  # Ignore various internal and external hyperlinks. In a real scenario
  # you may wish to track the length of the optional strings used with
  # urls.
  return if $token =~ m{^[fh]tt?ps?://};
  return if $token =~ m{^mailto:};
  return if $token =~ m{^(?:in|ex)ternal:};

  # We store the string width as data in the Worksheet object. We use
  # a double underscore key name to avoid conflicts with future names.
  #
  my $old_width    = $worksheet->{__col_widths}->[$col];
  my $string_width = string_width($token);

  if ( not defined $old_width or $string_width > $old_width )
  {

    # You may wish to set a minimum column width as follows.
    #return undef if $string_width < 10;

    $worksheet->{__col_widths}->[$col] = $string_width;
  }

  # Return control to write();
  return undef;
}

###############################################################################
#
# Very simple conversion between string length and string width for Arial 10.
# See below for a more sophisticated method.
#
sub string_width
{

  return 0.9 * length $_[0];
}

###############################################################################
#
# This function uses an external module to get a more accurate width for a
# string. Note that in a real program you could "use" the module instead of
# "require"-ing it and you could make the Font object global to avoid repeated
# initialisation.
#
# Note also that the $pixel_width to $cell_width is specific to Arial. For
# other fonts you should calculate appropriate relationships. A future version
# of S::WE will provide a way of specifying column widths in pixels instead of
# cell units in order to simplify this conversion.
#

1;
