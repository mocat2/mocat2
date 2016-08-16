package MOCATCore;
use warnings;
use strict;
use MOCATVariables;
use MOCATUsage2;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.
### DETERMINE TEMP DIRECOTRY ###

sub determine_temp
{
  ### NOTE ###
  # When calling the temp function, it's currently run: my $temp_dir = MOCATCore::determine_temp((200*1024*1024));
  # But the required disk space is not parsed. The user have to make sure there's enough space.
  # This is because it's hard to know exactly how much space is needed in each step.
  my $required = $_[0];
  my $temp_dir = $cwd;
  my $hostname = `hostname`;
  my @temp_info;
  if ($no_temp)
  {
  } else
  {
    chomp($hostname);
    my $number = 0;
    for my $key ( sort keys %conf )
    {
      if ( $key =~ m/^MOCAT_temp_dir\_/ )
      {
        $number++;
      }
    }
    my $set_temp = 1;
    for my $i ( 1 .. $number )
    {
      @temp_info = split( /\|/, $conf{"MOCAT_temp_dir_$i"} );
      if ( ( ( $hostname eq $temp_info[0] ) || $temp_info[0] eq 'any' )
           && $set_temp )
      {

        #	  my @space = `df -P $temp_info[2]`;
        #	  chomp($space[-1]);
        #	  $_ = $space[-1];
        #	  my @free_space = split;
        #	  if ($free_space[3] * 0.9 >= $required && $temp_info[1] >= $required) {
        $temp_dir = $temp_info[1];

        #	  }
        $set_temp = 0;
      }
    }
  }
  if ( $RUN_MODULE ne "1" )
  {
    if ( $temp_dir ne $cwd )
    {
      $temp_dir = "$temp_dir/$username/MOCAT_temp";
      foreach my $sample (@samples)
      {
        system "mkdir -p $temp_dir/$sample/temp";
      }

    } else
    {
      foreach my $sample (@samples)
      {
        system "mkdir -p $temp_dir/$sample/temp";
      }
    }
  }
  return ($temp_dir);
}

sub mkdir_or_die
{
  my $dir = shift;
  ( -d $dir ) or system "mkdir -p $dir";
  ( -d $dir ) or die "Could not create $dir.\nPlease check your permissions.";
}

### COUNTS NUMBER OF FQ FILES IN A SAMPLE FOLDER ###
sub count_fq
{
  my $counter = 0;
  foreach my $sample (@samples)
  {
    my @fqs = <$cwd/$sample/*.fq $cwd/$sample/*.fq.gz>;
    my @fqs2;
    for my $i ( 0 .. ( scalar @fqs - 1 ) )
    {
      if ( $fqs[$i] =~ m/pair\..*\.fq.gz/ || $fqs[$i] =~ m/single\.fq.gz/ )
      {
      } else
      {
        push @fqs2, $fqs[$i];
      }
    }
    @fqs = @fqs2;
    foreach my $fq (@fqs)
    {
      $counter++;
    }
  }
  return ($counter);
}

sub read_config
{
  my $config = shift;
  my $print  = shift;

  #  my %conf = ();
  if ( $print eq "print_config" )
  {
    print "GLOBAL SETTINGS\n";
  }
  if ( -e $config )
  {
    open( IN, "<$config" );
    while (<IN>)
    {
      unless ( $_ =~ /#/ || $_ !~ /:/ )
      {
        chomp;
        my @line = split( /\s+:\s+/, $_ );
        chomp @line;
        if ( $line[0] ne "screen_fasta_file_additional_usearch_cmd" )
        {
          if ( $line[0] =~ /qsub_add_param/ && !defined $line[1] )
          {
            $line[1] = "";
          }
          unless ( defined $line[1] && $line[0] )
          {
            die "\nERROR & EXIT: The field for $line[0] is not set in the config file.";
          }
        }
        my @line1;
        unless ( $line[0] =~ m/qsub_add_param/ )
        {
          @line1 = split( /\[|\(/, $line[1] );
        } else
        {
          $line1[0] = $line[1];
        }
        $line[0] =~ s/\s+//;
        $line1[0] =~ s/\s+$//;
        if ( $print eq "print_config" )
        {
          printf( "%-42s : %0s\n", $line[0], $line1[0] );
        }
        $conf{ $line[0] } = $line1[0];
      }
    }
    close IN;

    # SET ADDITIONAL CONF
    unless ( $conf{MOCAT_dir} )
    {
      die "ERROR & EXIT: Missing MOCAT_dir in config file. Have you specified a correct config file?";
    }
    unless ( $conf{MOCAT_data_dir} )
    {
      $conf{MOCAT_data_dir} = "$conf{MOCAT_dir}/data";
    }
    unless ( $conf{MOCAT_rownames_dir} )
    {
      $conf{MOCAT_rownames_dir} = "$conf{MOCAT_dir}/data";
    }
    return (%conf);
  } else
  {
    die "ERROR: Expected config file in $config.\nPlease copy the config file to this folder by executing 'cp MOCAT_DIRECTORY/MOCAT.cfg .'\n";
    exit;
  }
}

sub read_samples
{
  my $sample = shift;
  if ( -s $sample )
  {
    open( SAMPLE, '<', $sample );
    while (<SAMPLE>)
    {
      chomp $_;
      if ( $_ =~ m/^\S+$/ )
      {
        push( @samples, $_ );
      } else
      {
        print "WARNING: Misformatted line in sample file ($sample): $_\n";
      }
    }
    close SAMPLE;
    return @samples;
  } else
  {
    die "ERROR & EXIT: Sample file missing or incorrect. Specify sample file using --sample_file";
  }
}

sub check_files
{
  my $folder     = $_[0];
  my $file_post  = $_[1];
  my $file_post2 = "";
  print localtime() . ": Checking files...";
  for my $sample (@samples)
  {
    unless ( -s "$cwd/$sample/$folder/$sample.$file_post" )
    {
      die "\nERROR & EXIT: Missing file: $cwd/$sample/$folder/$sample.$file_post";
    }
  }
  print " OK!\n";
}

sub check_databases
{
  my @dbs = @{ $_[0] };
  print localtime() . ": Checking databases...";
  foreach my $db (@dbs)
  {
    unless ( -s "$data_dir/$db" || -s "$data_dir/$db.index.bwt" )
    {
      die "\nERROR & EXIT: Missing database '$db' in $data_dir";
    }
  }
  print " OK!\n";
}

sub execute_job
{
  my $job        = shift;
  my $processors = shift;
  $ZIP =~ s/pigz.*/pigz -p $processors/;

  my $numberJobs = shift;
  my $memory     = shift;
  my $ID         = "MC_$job";
  my $jobfile    = "$jobdir/MOCATJob_$job\_$date";
  my $now        = localtime();

  print "$now: Job Identifier is $date\n";
  print "$now: Executing $job jobs...\n";
  system "mkdir -p $cwd/logs $cwd/logs/other $cwd/logs/$job $cwd/logs/status";
  system "mkdir -p $cwd/logs/$job/startstop $cwd/logs/$job/samples $cwd/logs/$job/jobs $cwd/logs/$job/commands";
  my $temp_dir      = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
  my $range         = 9999999999;
  my $random_number = rand($range);
  $ID = "MC" . substr( $random_number, 0, 8 );
  chomp( my $hm = `hostname` );
  chomp( my $me = `whoami` );
  my $cmdargs = join( ' ', @args );
  my $on_host = '';

  if ( $run_on_host ne "" )
  {
    $on_host = "EXECUTE ON HOST:      $run_on_host\n";
  }
  open OUT, ">", "$cwd/logs/$job/commands/MOCATJob_$job.$date.command.log" or die "ERROR & EXIT: Cannot write $cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
  my $xa = join( ' ', @args );
  print OUT <<"EOF";

$sep
                   MOCAT - Metagenomics Analysis Toolkit                $version
 by Jens Roat Kultima, Luis Pedro Coelho, Shinichi Sunagawa @ Bork Group, EMBL
$sep
MOCAT ID:             $MOCAT_ID
COMMAND EXECUTED:     MOCAT.pl $xa
EXECUTION START TIME: $now
JOB IDENTIFIER:       $date
USERNAME:             $me
PROCESSORS:           $processors
HOSTNAME:             $hm
$on_host
USER:                 $me
QUEUE ID:             $ID
CURRENT DIRECTTORY:   $cwd
TEMP DIRECTTORY:      $temp_dir/<SAMPLE>/temp


LOG FILE:             $cwd/logs/$job/samples/MOCATJob_$job.<SAMPLE.log
START STOP LOG FILE:  $cwd/logs/$job/startstop/MOCATJob_$job.startstop.log
EOF

  chomp( our $systemType = `uname -s` );

  unless ( $systemType =~ m/Darwin/ )
  {
    print OUT "FREE MEMORY (in GB):\n";
    system "free -g >> $cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
  }
  close OUT;

  open OUT, ">>", "$cwd/logs/$job/commands/MOCATJob_$job.$date.command.log" or die "ERROR & EXIT: Cannot write $cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
  print OUT "\n\nFREE DISK:\n";
  close OUT;
  system "df -h $cwd/ >> $cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
  open OUT, ">>", "$cwd/logs/$job/commands/MOCATJob_$job.$date.command.log" or die "ERROR & EXIT: Cannot write $cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
  print OUT "\n\nFREE TEMP DISK:\n";
  close OUT;
  system "df -h $temp_dir/ >> $cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
  open OUT, ">>", "$cwd/logs/$job/commands/MOCATJob_$job.$date.command.log" or die "ERROR & EXIT: Cannot write $cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
  print OUT "\n\nSAMPLES:\n";
  my $counter = 0;

  foreach my $sample (@samples)
  {
    $counter++;
    print OUT "$counter -> $sample\n";
  }
  print OUT "\n\n";
  print OUT "MOCAT SETTINGS:\n";
  print MOCATCore::print_settings_OUT("MOCAT");
  print OUT "\n\n";
  print OUT "JOB SETTINGS:\n";
  print MOCATCore::print_settings_OUT("$job");
  print OUT "\n\n";
  print OUT "COMMAND TO EXECUTE:\n";

  # Edit pre execute command
  if ( $systemType =~ m/Darwin/ )
  {
    if ( $conf{MOCAT_pre_execute} ne '' )
    {
      system "sed -i '' 's%^%$conf{MOCAT_pre_execute}; %' $jobfile";
    }
  } else
  {
    if ( $conf{MOCAT_pre_execute} ne '' )
    {
      system "sed -i 's%^%$conf{MOCAT_pre_execute}; %' $jobfile";
    }
  }

  if ( $qsub eq "SGE" )
  {
    print localtime() . ": Queue Identifier is $ID\n";
    my $SGEEXTRAPARAMS = $conf{MOCAT_SGE_qsub_add_param};
    unless ( defined $SGEEXTRAPARAMS )
    {
      $SGEEXTRAPARAMS = "";
    }
    print OUT "qsub -sync y $SGEEXTRAPARAMS -V -t 1-" . $numberJobs . " $cwd/MOCATJobArraySGE.$job.$date.sh\n";
  } elsif ( $qsub eq "PBS" )
  {
    print localtime() . ": Queue Identifier is $ID\n";
    my $SGEEXTRAPARAMS = $conf{MOCAT_PBS_qsub_add_param};
    unless ( defined $SGEEXTRAPARAMS )
    {
      $SGEEXTRAPARAMS = "";
    }
    print OUT "qsub $SGEEXTRAPARAMS -W block=true -V $cwd/MOCATJobArrayPBS.$job.$date.sh\n";
  } elsif ( $qsub eq "LSF" )
  {
    my $JOB = substr( "$ID", 0, 15 );
    my $SGEEXTRAPARAMS = $conf{MOCAT_LSF_qsub_add_param};
    unless ( defined $SGEEXTRAPARAMS )
    {
      $SGEEXTRAPARAMS = "";
    }
    if ( $conf{MOCAT_LSF_memory_limit} ne "" )
    {
      $SGEEXTRAPARAMS .= " -M $conf{MOCAT_LSF_memory_limit} ";
      $memory = $conf{MOCAT_LSF_memory_limit};
    } else
    {
      $memory =~ s/gb/1000/;
      $memory =~ s/G/1000/;
    }
    if ( $conf{MOCAT_LSF_queue} ne "" )
    {
      $SGEEXTRAPARAMS .= " -q $conf{MOCAT_LSF_queue} ";
    }
    print localtime() . ": Queue Identifier is $ID\n";
    unless ( defined $SGEEXTRAPARAMS )
    {
      $SGEEXTRAPARAMS = "";
    }
    if ( $conf{MOCAT_LSF_max_concurrent} ne '' )
    {
      $conf{MOCAT_LSF_max_concurrent} = "\%$conf{MOCAT_LSF_max_concurrent}";
    } else
    {
      $conf{MOCAT_LSF_max_concurrent} = "";
    }

    print OUT "bsub $SGEEXTRAPARAMS " . "-J \"$JOB\[1-$numberJobs\]$conf{MOCAT_LSF_max_concurrent}\" " . "-cwd \"$cwd\" " . " -M $memory " .    # the captial K flag waits for jobs to finish, or so they say...
      "-R \"select[(mem >= $memory)] rusage[mem=$memory] span[hosts=1]\" " . "-n $processors " . "-o $cwd/logs/other/MOCATJob_$job\_$date.\%I.log " . "-e $cwd/logs/other/MOCATJob_$job\_$date.\%I.err " . "$cwd/MOCATJobArrayLSF.$job.$date.sh\n";
  } elsif ( $qsub eq "none" )
  {
    print OUT "bash $jobfile\n";
  }
  my $error = 0;
  $counter = 0;
  unless ( $mod_single_job eq 'TRUE' )
  {
    # Currently removed from MOCAT because we now use one log file per sample per run instead of one log file per sample
    #	foreach my $sample (@samples) {
    #		$counter++;
    #		open my $X, '>>', "$cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log" or die "ERROR & EXIT: Cannot write $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log";
    #		print $X "\n\n$sep\n";
    #		print $X "                   MOCAT - Metagenomics Analysis Toolkit                $version\n";
    #		print $X " by Jens Roat Kultima, Luis Pedro Coelho, Shinichi Sunagawa @ Bork Group, EMBL\n$sep\n";
    #		print $X "EXECUTION START TIME: " . localtime() . "\n";
    #		print $X "USERNAME:             $me\n";
    #		print $X "JOB IDENTIFIER:       $date\n";
    #		print $X "PROCESSORS:           $processors\n";
    #		print $X "HOSTNAME:             $hm\n";
    #		if ( $run_on_host ne "" ) {
    #			print $X "EXECUTE ON HOST:      $run_on_host\n";
    #		}
    #		print $X "QUEUE ID:             $ID.$counter\n";
    #		print $X "SAMPLE FILE:          $sample_file\n";
    #		print $X "CURRENT DIRECTTORY:   $cwd\n";
    #		print $X "COMMAND EXECUTED:     MOCAT.pl " . join( ' ', @args ) . "\n";
    #		print $X "MOCAT ID:             $MOCAT_ID\n";
    #		print $X "TEMP DIRECTTORY:      $temp_dir/$sample/temp\nOUTPUT:\n";
    #		close $X;
    #	}
  }
  open X, '>>', "$cwd/logs/$job/startstop/MOCATJob_$job.startstop.log" or die "ERROR & EXIT: Cannot write $cwd/logs/$job/startstop/MOCATJob_$job.startstop.log";
  print X << "EOF";

$sep
                   MOCAT - Metagenomics Analysis Toolkit                $version
 by Jens Roat Kultima, Luis Pedro Coelho, Shinichi Sunagawa @ Bork Group, EMBL
$sep
MOCAT ID:             $MOCAT_ID
EXECUTION START TIME: $now
JOB IDENTIFIER:       $date
USERNAME:             $me
PROCESSORS:           $processors
HOSTNAME:             $hm
$on_host
QUEUE ID:             $ID
CURRENT DIRECTTORY:   $cwd
COMMAND EXECUTED:     MOCAT.pl $cmdargs
TEMP DIRECTTORY:      $temp_dir/SAMPLE/temp
SAMPLES:
EOF
  $counter = 0;

  foreach my $sample (@samples)
  {
    $counter++;
    print X "$counter -> $sample\n";
  }
  print X "\nSTATUS:\n";
  close X;
  if ( $qsub eq "SGE" )
  {
    my $SGEEXTRAPARAMS = $conf{MOCAT_SGE_qsub_add_param};
    unless ( defined $SGEEXTRAPARAMS )
    {
      $SGEEXTRAPARAMS = "";
    }
    if ( $SGEEXTRAPARAMS =~ m/(.*-l h_vmem=)(\d+)(G.*)/ )
    {
      my $two = $2 / $processors;
      $SGEEXTRAPARAMS = "$1$two$3";
    }

    my $CWD = $cwd;
    $CWD =~ s/\//\\\//g;
    if ( $systemType =~ m/Darwin/ )
    {
      system <<"EOF";
                cp $scr_dir/MOCATExecuteJobSGE.pl $cwd/MOCATExecuteJob.$date.pl &&
                sed -i '' 's/RJSS/MOCATJob_$job.$date.\$i.sh/' $cwd/MOCATExecuteJob.$date.pl &&
                $cwd/MOCATExecuteJob.$date.pl $jobfile &&
                cp $scr_dir/MOCATJobArraySGE.sh $cwd/MOCATJobArraySGE.$job.$date.sh &&
                sed -i '' 's/pe smp 8/pe $conf{MOCAT_SGE_parallell_env} $processors/' $cwd/MOCATJobArraySGE.$job.$date.sh &&
                sed -i '' 's/-N MOCATJobArray/-N $ID/' $cwd/MOCATJobArraySGE.$job.$date.sh &&
                sed -i '' 's/RJB/MOCATJob_$job.$date/g' $cwd/MOCATJobArraySGE.$job.$date.sh &&
                sed -i '' 's/SGEJOBID/$date/' $cwd/MOCATJobArraySGE.$job.$date.sh &&
                sed -i '' 's/SGEARRAY/$job/g' $cwd/MOCATJobArraySGE.$job.$date.sh &&
                sed -i '' 's/MCsaMPLE/sample/g' $cwd/MOCATJobArraySGE.$job.$date.sh &&
                sed -i '' 's/CWD/$CWD/g' $cwd/MOCATJobArraySGE.$job.$date.sh &&
                rm $cwd/MOCATExecuteJob.$date.pl
EOF
    } else
    {
      system <<"EOF";
                cp $scr_dir/MOCATExecuteJobSGE.pl $cwd/MOCATExecuteJob.$date.pl &&
                sed -i 's/RJSS/MOCATJob_$job.$date.\$i.sh/' $cwd/MOCATExecuteJob.$date.pl &&
                $cwd/MOCATExecuteJob.$date.pl $jobfile &&
                cp $scr_dir/MOCATJobArraySGE.sh $cwd/MOCATJobArraySGE.$job.$date.sh &&
                sed -i 's/pe smp 8/pe $conf{MOCAT_SGE_parallell_env} $processors/' $cwd/MOCATJobArraySGE.$job.$date.sh &&
                sed -i 's/-N MOCATJobArray/-N $ID/' $cwd/MOCATJobArraySGE.$job.$date.sh &&
                sed -i 's/RJB/MOCATJob_$job.$date/g' $cwd/MOCATJobArraySGE.$job.$date.sh &&
                sed -i 's/SGEJOBID/$date/' $cwd/MOCATJobArraySGE.$job.$date.sh &&
                sed -i 's/SGEARRAY/$job/g' $cwd/MOCATJobArraySGE.$job.$date.sh &&
                sed -i 's/MCsaMPLE/sample/g' $cwd/MOCATJobArraySGE.$job.$date.sh &&
                sed -i 's/CWD/$CWD/g' $cwd/MOCATJobArraySGE.$job.$date.sh &&
                rm $cwd/MOCATExecuteJob.$date.pl
EOF
    }
    my $cmd = "qsub -sync y $SGEEXTRAPARAMS -V -t 1-" . $numberJobs . " $cwd/MOCATJobArraySGE.$job.$date.sh | tee -a $cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
    if ( $run_on_host ne "" )
    {
      $cmd = "ssh $run_on_host \"cd $cwd && $cmd \"";
    }

    if ($only_print)
    {
      push @all_commands, "$cmd";
      close OUT;
      print "TO RUN YOUR JOBS, EXECUTE: $cmd\n";
    } else
    {
      print OUT "\n\nQSUB STATUS:\n";
      close OUT;
      my @ERR;
      if ( $realtime_status || $conf{realtime_status_use} eq 'yes' )
      {
        if ($realtime_status)
        {
        } else
        {
          $realtime_status = $conf{realtime_status_timer};
        }
        $error = MOCATUsage2::usage( $cmd, $ID, $job, $date );
        if ( $error == 0 )
        {
          if ( $conf{realtime_status_print} eq 'yes' )
          {
            print "\n";
          }
          print "$sep\n";
        }
      } else
      {
        @ERR = `$cmd`;
      }
      open OUT, ">>", "$cwd/logs/$job/commands/MOCATJob_$job.$date.command.log" or die "ERROR & EXIT: Cannot write $cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
      open STATUS, ">>", "$cwd/logs/status/$job.failed";
      print STATUS "$sep\n";
      open STATUS2, ">>", "$cwd/logs/status/$job.successful";
      print STATUS2 "$sep\n";
      print OUT "\n";
      foreach my $err (@ERR)
      {
        $err =~ m/Job (.*) exited.*exit code (\d+)./;
        my $nr1 = $1;
        my $two = "undefed";
        if ( defined $2 )
        {
          $two = $2;
        } else
        {
          if ( $err =~ m/Job (.*) exited because of signal SIGKILL/ )
          {
            $two = "999";
          } elsif ( $err =~ m/Unable to run job (.*)/ )
          {
            $two = "999";
          }
          $nr1 = $1;
        }
        unless ( $two eq "undefed" )
        {
          if ( $two != 0 )
          {
            if ( !$realtime_status )
            {
              my $sample = $samples[0];
              if ( scalar @samples > 1 )
              {
                my @id = split /\./, $nr1;
                $sample = $samples[ $id[1] - 1 ];
              }

              my @E;
              if ( $sample eq $date )
              {    # has to do with module support
                @E = `tail -3 $cwd/logs/$job/samples/MOCATJob_$job.$date.log`;
                print "\n\nSAMPLE $sample FAILED. Here are the last 3 lines of the log file\n($cwd/logs/$job/samples/MOCATJob_$job.$date.log):\n$sep\n";
              } else
              {
                @E = `tail -3 $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log`;
                print "\n\nSAMPLE $sample FAILED. Here are the last 3 lines of the log file\n($cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log):\n$sep\n";
              }
              print OUT "FAILED SAMPLE: $sample\n";
              print STATUS "$sample\t$date\tMOCAT.pl " . join( " ", @args ) . "\n";
              foreach my $e (@E)
              {
                print $e;
              }
              print "$sep\n\n";
            }
            $error = 1;
          } elsif ( $two == 0 )
          {
            if ( !$realtime_status )
            {
              my $sample = $samples[0];
              if ( scalar @samples > 1 )
              {
                my @id = split /\./, $nr1;
                $sample = $samples[ $id[1] - 1 ];
              }
              print STATUS2 "$sample\t$date\tMOCAT.pl " . join( " ", @args ) . "\n";
            }
          }
        }
      }
      close STATUS;
      close STATUS2;
      print OUT "\n";
      close OUT;
    }
  }

  # PBS
  elsif ( $qsub eq "PBS" )
  {
    my $SGEEXTRAPARAMS = $conf{MOCAT_PBS_qsub_add_param};
    unless ( defined $SGEEXTRAPARAMS )
    {
      $SGEEXTRAPARAMS = "";
    }
    my $JOB = substr( "$ID", 0, 15 );
    if ( $systemType =~ m/Darwin/ )
    {
      system <<"EOF";
			    cp $scr_dir/MOCATExecuteJobPBS.pl $cwd/MOCATExecuteJob.$date.pl &&
                sed -i '' 's/RJSS/MOCATJob_$job.$date.\$i.sh/' $cwd/MOCATExecuteJob.$date.pl &&
                $cwd/MOCATExecuteJob.$date.pl $cwd/MOCATJob_$job\_$date &&
                cp $scr_dir/MOCATJobArrayPBS.sh $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i '' 's/ncpus=CPUS/ncpus=$processors/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i '' 's/-N MOCATJobArray/-N $JOB/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i '' 's/RJB/MOCATJob_$job.$date/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i '' 's/PBSJOBID/$date/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i '' 's/PBSARRAY/$job/g' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i '' 's/MAXJOBS/$numberJobs/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i '' 's/SAMPLE/$date/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i '' 's/PBSDATE/$date/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i '' 's/#PBS -J 1-1\$/##PBS -J 1-1/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i '' 's/MCsaMPLE/sample/g' $cwd/MOCATJobArrayPBS.$job.$date.sh && rm $cwd/MOCATExecuteJob.$date.pl
EOF
      if ( $numberJobs == 1 )
      {
        system "sed -i '' 's/\$PBS_ARRAY_INDEX/1/' $cwd/MOCATJobArrayPBS.$job.$date.sh";
      }
    } else
    {
      system <<"EOF";
			    cp $scr_dir/MOCATExecuteJobPBS.pl $cwd/MOCATExecuteJob.$date.pl &&
                sed -i 's/RJSS/MOCATJob_$job.$date.\$i.sh/' $cwd/MOCATExecuteJob.$date.pl &&
                $cwd/MOCATExecuteJob.$date.pl $cwd/MOCATJob_$job\_$date &&
                cp $scr_dir/MOCATJobArrayPBS.sh $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i 's/ncpus=CPUS/ncpus=$processors/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i 's/-N MOCATJobArray/-N $JOB/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i 's/RJB/MOCATJob_$job.$date/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i 's/PBSJOBID/$date/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i 's/PBSARRAY/$job/g' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i 's/MAXJOBS/$numberJobs/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i 's/PBSDATE/$date/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i 's/SAMPLE/$date/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i 's/#PBS -J 1-1\$/##PBS -J 1-1/' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                sed -i 's/MCsaMPLE/sample/g' $cwd/MOCATJobArrayPBS.$job.$date.sh &&
                rm $cwd/MOCATExecuteJob.$date.pl
EOF
      if ( $numberJobs == 1 )
      {
        system "sed -i 's/\$PBS_ARRAY_INDEX/1/' $cwd/MOCATJobArrayPBS.$job.$date.sh";
      }
    }
    my $cmd = "qsub $SGEEXTRAPARAMS -W block=true -V $cwd/MOCATJobArrayPBS.$job.$date.sh | tee -a $cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
    if ( $run_on_host ne "" )
    {
      $cmd = "ssh $run_on_host \"cd $cwd && $cmd \"";
    }

    if ($only_print)
    {
      push @all_commands, "$cmd";
      print "TO RUN YOUR JOBS, EXECUTE: $cmd\n";
      close OUT;
    } else
    {
      print OUT "\n\nQSUB STATUS:\n";
      close OUT;
      my @ERR = `$cmd`;
      open IN, '<', "$cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
      while (<IN>)
      {
        if ( $_ =~ m/.*PBSERRORBYMOCAT!.*/ )
        {
          $error = 1;
        }
      }
      close IN;
    }
  }

  # LSF
  elsif ( $qsub eq "LSF" )
  {
    system "mkdir -p $cwd/logs/LSF_specific";
    system "touch $cwd/logs/LSF_specific/$job.$date.status";
    my $SGEEXTRAPARAMS = $conf{MOCAT_LSF_qsub_add_param};
    unless ( defined $SGEEXTRAPARAMS )
    {
      $SGEEXTRAPARAMS = "";
    }
    if ( $conf{MOCAT_LSF_memory_limit} ne "" )
    {
      $SGEEXTRAPARAMS .= " -M $conf{MOCAT_LSF_memory_limit} ";
      $memory = $conf{MOCAT_LSF_memory_limit};
    } else
    {
      $memory =~ s/gb/1000/;
      $memory =~ s/G/1000/;
    }
    if ( $conf{MOCAT_LSF_queue} ne "" )
    {
      $SGEEXTRAPARAMS .= " -q $conf{MOCAT_LSF_queue} ";
    }

    my $JOB = substr( "$ID", 0, 15 );
    if ( $systemType =~ m/Darwin/ )
    {
      system <<"EOF";
			    cp $scr_dir/MOCATExecuteJobLSF.pl $cwd/MOCATExecuteJob.$date.pl &&
                sed -i '' 's/RJSS/MOCATJob_$job.$date.\$i.sh/' $cwd/MOCATExecuteJob.$date.pl &&
                $cwd/MOCATExecuteJob.$date.pl $cwd/MOCATJob_$job\_$date &&
                cp $scr_dir/MOCATJobArrayLSF.sh $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i '' 's/BSUB -n CPUS/BSUB -n $processors/' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i '' 's/-J MOCATJobArray/-J $JOB/' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i '' 's/RJB/MOCATJob_$job.$date/g' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i '' 's/LSFJOBID/$date/' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i '' 's/LSFARRAY/$job/g' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i '' 's/MAXJOBS/$numberJobs/' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i '' 's/SAMPLE/$date/' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i '' 's/LSFDATE/$date/' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i '' s|CWD|$cwd|g' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                rm -f $cwd/MOCATExecuteJob.$date.pl
EOF
    } else
    {
      system <<"EOF";
			    cp $scr_dir/MOCATExecuteJobLSF.pl $cwd/MOCATExecuteJob.$date.pl &&
                sed -i 's/RJSS/MOCATJob_$job.$date.\$i.sh/' $cwd/MOCATExecuteJob.$date.pl &&
                $cwd/MOCATExecuteJob.$date.pl $cwd/MOCATJob_$job\_$date &&
                cp $scr_dir/MOCATJobArrayLSF.sh $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i 's/BSUB -n CPUS/BSUB -n $processors/' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i 's/-J \"MOCATJobArray/-J \"$JOB/' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i 's/RJB/MOCATJob_$job.$date/g' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i 's/LSFJOBID/$date/' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i 's/LSFARRAY/$job/g' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i 's/MAXJOBS/$numberJobs/' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i 's/LSFDATE/$date/' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i 's/SAMPLE/$date/' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                sed -i 's|CWD|$cwd|g' $cwd/MOCATJobArrayLSF.$job.$date.sh &&
                rm -f $cwd/MOCATExecuteJob.$date.pl
EOF
    }

    # $conf{MOCAT_LSF_max_concurrent} a;ready extended with a % in the beginning above
    my $cmd = "bsub $SGEEXTRAPARAMS -J \"$JOB\[1-$numberJobs\]$conf{MOCAT_LSF_max_concurrent}\" -cwd \"$cwd\" -R \"select[(mem > $memory)] rusage[mem=$memory] span[hosts=1]\" -n $processors -o $cwd/logs/other/MOCATJob_$job\_$date.\%I.log -e $cwd/logs/other/MOCATJob_$job\_$date.\%I.err $cwd/MOCATJobArrayLSF.$job.$date.sh | tee -a $cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
    if ( $run_on_host ne "" )
    {
      $cmd = "ssh $run_on_host \"cd $cwd && $cmd \"";
    }

    if ($only_print)
    {
      push @all_commands, "$cmd";
      print "TO RUN YOUR JOBS, EXECUTE: $cmd\n";
      close OUT;
    } else
    {
      print OUT "QSUB STATUS:\n";
      print localtime() . ": Submitting LSF jobs...\n";
      print "SUBMIT: $cmd\n";
      my $return = `$cmd`;
      $return =~ m/Job <(.*)> is.*/;
      my $id = $1;
      print localtime() . ": $return";
      my $continue = 1;

      while ($continue)
      {

        unless ( -e "$cwd/logs/LSF_specific/$job.$date.status" )
        {
          $continue = 0;
          $error    = 1;
        }

        my @lines = `qstat -u $username | grep $username | grep -E ' $id\[[0-9]*\] '`;
        if ( scalar @lines == 0 )
        {
          $continue = 0;
        }
        chomp( my @finished = `cut -f 1 -d'.' $cwd/logs/LSF_specific/$job.$date.status 2>/dev/null` );
        my %finished;
        my @missing;
        foreach my $finished (@finished)
        {
          $finished{$finished} = 1;
        }
        my $total = 0;
        for my $i ( 1 .. $numberJobs )
        {
          if ( $finished{$i} )
          {
            $total++;
          } else
          {
            push @missing, $i;
          }
        }
        if ( $total == $numberJobs )
        {
          $continue = 0;
        }

        #print localtime() . ": waiting for $JOB | qID=$id | queued=" . scalar @lines . " | finished=$total | total=$numberJobs | unfinished=";
        #foreach my $i (@missing) {
        #	print "$i:" . $samples[ ( $i - 1 ) ] . " ";
        #}
        #print "\n";
        last if $continue == 0;
        sleep 180;
      }
      open STATUS, ">>", "$cwd/logs/status/$job.failed";
      print STATUS "$sep\n";
      chomp( my @finished = `grep 'failed' $cwd/logs/LSF_specific/$job.$date.status 2>/dev/null | cut -f 1 -d'.'` );
      foreach my $finished (@finished)
      {
        print OUT "FAILED SAMPLE: " . $samples[ ( $finished - 1 ) ] . "\n";
        my $sample = $samples[ ( $finished - 1 ) ];

        my @E;
        if ( $date eq $sample )
        {    # this if has to do with the module support
          @E = `tail -3 $cwd/logs/$job/samples/MOCATJob_$job.$date.log`;
          print "\n\nSAMPLE $sample FAILED. Here are the last 3 lines of the log file\n($cwd/logs/$job/samples/MOCATJob_$job.$date.log):\n$sep\n";
        } else
        {
          @E = `tail -3 $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log`;
          print "\n\nSAMPLE $sample FAILED. Here are the last 3 lines of the log file\n($cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log):\n$sep\n";
        }

        print STATUS "$sample\t$date\tMOCAT.pl " . join( " ", @args ) . "\n";
        foreach my $e (@E)
        {
          print $e;
        }
        print "$sep\n\n";

      }
      if ( scalar @finished > 0 )
      {
        $error = 1;
      }
      close OUT;
      close STATUS;
      open STATUS2, ">>", "$cwd/logs/status/$job.successful";
      print STATUS2 "$sep\n";
      @finished = ();
      chomp( @finished = `grep 'completed' $cwd/logs/LSF_specific/$job.$date.status 2>/dev/null | cut -f 1 -d'.'` );
      foreach my $finished (@finished)
      {
        my $sample = $samples[ ( $finished - 1 ) ];
        print STATUS2 "$sample\t$date\tMOCAT.pl " . join( " ", @args ) . "\n";
      }
      close STATUS2;
      open IN, '<', "$cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
    }
  }

  # SLURM
  elsif ( $qsub eq "SLURM" )
  {
    system "mkdir -p $cwd/logs/SLURM_specific/";
    system "touch $cwd/logs/SLURM_specific/$job.$date.status";

    if ( $conf{MOCAT_SLURM_max_concurrent} ne '' )
    {
      $conf{MOCAT_SLURM_max_concurrent} = "\%$conf{MOCAT_SLURM_max_concurrent}";
    } else
    {
      $conf{MOCAT_SLURM_max_concurrent} = "";
    }

    my $SGEEXTRAPARAMS = $conf{MOCAT_SLURM_qsub_add_param};
    unless ( defined $SGEEXTRAPARAMS )
    {
      $SGEEXTRAPARAMS = "";
    }
    my $JOB = substr( "$ID", 0, 15 );
    if ( $systemType =~ m/Darwin/ )
    {
      die "Unsupported. Sorry.";
    } else
    {
      system <<"EOF";
			    cp $scr_dir/MOCATExecuteJobSLURM.pl $cwd/MOCATExecuteJob.$date.pl &&
                sed -i 's/RJSS/MOCATJob_$job.$date.\$i.sh/' $cwd/MOCATExecuteJob.$date.pl &&
                $cwd/MOCATExecuteJob.$date.pl $cwd/MOCATJob_$job\_$date &&
                cp $scr_dir/MOCATJobArraySLURM.sh $cwd/MOCATJobArraySLURM.$job.$date.sh &&
                sed -i 's/--mincpus=1/--mincpus=$processors/' $cwd/MOCATJobArraySLURM.$job.$date.sh &&
                sed -i 's/--cpus-per-task=1/--cpus-per-task=$processors/' $cwd/MOCATJobArraySLURM.$job.$date.sh &&
                sed -i 's/--job-name=NAME/--job-name=$JOB/' $cwd/MOCATJobArraySLURM.$job.$date.sh &&
                sed -i 's/RJB/MOCATJob_$job.$date/' $cwd/MOCATJobArraySLURM.$job.$date.sh &&
                sed -i 's/SLURMJOBID/$date/' $cwd/MOCATJobArraySLURM.$job.$date.sh &&
                sed -i 's/SLURMARRAY/$job/g' $cwd/MOCATJobArraySLURM.$job.$date.sh &&
                sed -i 's/MAXJOBS/$numberJobs$conf{MOCAT_SLURM_max_concurrent}/' $cwd/MOCATJobArraySLURM.$job.$date.sh &&
                sed -i 's/SLURMDATE/$date/' $cwd/MOCATJobArraySLURM.$job.$date.sh &&
                sed -i 's|CWD|$cwd|g' $cwd/MOCATJobArraySLURM.$job.$date.sh &&
                rm $cwd/MOCATExecuteJob.$date.pl
EOF
      my $cmd = "sbatch $SGEEXTRAPARAMS --output=$cwd/logs/other/$job.$date.log --error=$cwd/logs/other/$job.$date.log $cwd/MOCATJobArraySLURM.$job.$date.sh | tee -a $cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
      if ( $run_on_host ne "" )
      {
        $cmd = "ssh $run_on_host \"cd $cwd && $cmd \"";
      }

      if ($only_print)
      {
        push @all_commands, "$cmd";
        print "TO RUN YOUR JOBS, EXECUTE: $cmd\n";
      } else
      {
        print OUT "$cmd\n";
        print OUT "\n\nQSUB STATUS:\n";
        chomp( my $return = `$cmd` );
        $return =~ m/.* (\d+)/;
        my $id = $1;
        unless ($id) {
         die "ERROR & EXIT: No SLURM ID available. Did the queuing of the job fail?";
        }
        print localtime() . ": SLURM Job ID is $id - $ID\n";
        system("touch $cwd/logs/SLURM_specific/$job.$date.status");
        sleep 10;
        print localtime() . ": Waiting for jobs to finish...\n";
        my $continue = 1;
        $error = 0;

        while ($continue)
        {
          unless ( -e "$cwd/logs/SLURM_specific/$job.$date.status" )
          {
            $continue = 0;
            $error    = 1;
          }
          chomp( my $FAILED    = `sacct -j $id | grep $JOB | grep -c 'FAILED'` );
          chomp( my $COMPLETED = `sacct -j $id | grep $JOB | grep -c 'COMPLETED'` );
          if ( ($FAILED) + ($COMPLETED) == $numberJobs )
          {
            $continue = 0;
          }

          chomp( my @finished = `cut -f 1 -d'.' $cwd/logs/SLURM_specific/$job.$date.status 2>/dev/null` );
          my %finished;
          my @missing;
          foreach my $finished (@finished)
          {
            $finished{$finished} = 1;
          }
          my $total = 0;
          for my $i ( 1 .. $numberJobs )
          {
            if ( $finished{$i} )
            {
              $total++;
            } else
            {
              push @missing, $i;
            }
          }
          if ( $total == $numberJobs )
          {
            $continue = 0;
          }
          last if $continue == 0;
          sleep 180;
        }
        open STATUS, ">>", "$cwd/logs/status/$job.failed";
        print STATUS "$sep\n";
        chomp( my @finished = `grep 'failed' $cwd/logs/SLURM_specific/$job.$date.status 2>/dev/null | cut -f 1 -d'.'` );
        foreach my $finished (@finished)
        {
          print OUT "FAILED SAMPLE: " . $samples[ ( $finished - 1 ) ] . "\n";
          my $sample = $samples[ ( $finished - 1 ) ];

          my @E;
          if ( $date eq $sample )
          {    # this if has to do with the module support
            @E = `tail -3 $cwd/logs/$job/samples/MOCATJob_$job.$date.log`;
            print "\n\nSAMPLE $sample FAILED. Here are the last 3 lines of the log file\n($cwd/logs/$job/samples/MOCATJob_$job.$date.log):\n$sep\n";
          } else
          {
            @E = `tail -3 $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log`;
            print "\n\nSAMPLE $sample FAILED. Here are the last 3 lines of the log file\n($cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log):\n$sep\n";
          }

          print STATUS "$sample\t$date\tMOCAT.pl " . join( " ", @args ) . "\n";
          foreach my $e (@E)
          {
            print $e;
          }
          print "$sep\n\n";

        }
        if ( scalar @finished > 0 )
        {
          $error = 1;
        }
        close OUT;
        close STATUS;
        open STATUS2, ">>", "$cwd/logs/status/$job.successful";
        print STATUS2 "$sep\n";
        @finished = ();
        chomp( @finished = `grep 'completed' $cwd/logs/SLURM_specific/$job.$date.status 2>/dev/null | cut -f 1 -d'.'` );
        foreach my $finished (@finished)
        {
          my $sample = $samples[ ( $finished - 1 ) ];
          print STATUS2 "$sample\t$date\tMOCAT.pl " . join( " ", @args ) . "\n";
        }
        close STATUS2;
        open IN, '<', "$cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
      }
    }
  }

  # NONE
  elsif ( $qsub eq "none" )
  {
    close OUT;
    system "mkdir -p $cwd/logs";
    if ($only_print)
    {
      push @all_commands, "$cwd/MOCATJob_$job\_$date";
      system "cp $cwd/MOCATJob_$job\_$date $cwd/logs/$job/jobs/MOCATJob_$job.$date.sh";
      my $cmd = "bash $cwd/MOCATJob_$job\_$date";
      if ( $run_on_host ne "" )
      {
        $cmd = "ssh $run_on_host \"cd $cwd && $cmd \"";
      }
      print "TO RUN YOUR JOBS, EXECUTE: $cmd\n";
    } else
    {
      system "cp $cwd/MOCATJob_$job\_$date $cwd/logs/$job/jobs/MOCATJob_$job.$date.sh";
      my $cmd = "bash $cwd/MOCATJob_$job\_$date";
      if ( $run_on_host ne "" )
      {
        $cmd = "ssh $run_on_host \"cd $cwd && $cmd \"";
      }
      system "$cmd";
    }
  } else
  {
    print OUT "\n\nEXECUTION END TIME:\n" . localtime() . "\n\n\n";
    print OUT "EXIT STATUS:\nEXECUTION FAIL (Not correct queuing system in config file)\n\n";
    close OUT;
    die "ERROR & EXIT: Expected 'qsub_system' in config file to be either 'SGE', 'PBS', 'LSF', 'SLURM' or 'none'.";
  }

  if ($only_print)
  {
    open OUT, ">>", "$cwd/logs/$job/commands/MOCATJob_$job.$date.command.log" or die "ERROR & EXIT: Cannot write $cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
    print OUT "\n\nEXECUTION END TIME:\n" . localtime() . "\n\n\n";
    print OUT "EXIT STATUS:\nEXECUTION NOT PERFORMED\n\n";
    close OUT;
    print localtime() . ": EXECUTION NOT PERFORMED. ONLY WROTE JOB FILES.\n";
  } else
  {

    #Print status of filesystem to log files
    #foreach my $sample (@samples) {
    #	system "cd $cwd/$sample/ && ls -lhR >> $cwd/logs/$job/filesystem/MOCATJob_$job.$date.filesystem.log";
    #}

    if ( $? eq 0 && $error == 0 )
    {
      open OUT, ">>", "$cwd/logs/$job/commands/MOCATJob_$job.$date.command.log" or die "ERROR & EXIT: Cannot write $cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
      system "rm -f $cwd/MOCATJob_$job\_$date 2>/dev/null";
      system "rm -f $cwd/MOCATJobArray*.$job.$date.sh 2>/dev/null";
      system "mv $cwd/MOCATJob_$job.$date* $cwd/logs/$job/jobs/ 2>/dev/null";
      system "mv $cwd/$ID.po* $cwd/$ID.pe* $cwd/$ID.o* $cwd/$ID.e* $cwd/logs/other/ 2>/dev/null";
      chomp( my $cd = `date +'\%F \%T'` );
      chomp( my $h  = `hostname` );
      system "echo 'MOCAT.pl " . join( " ", @args ) . " @ $cd ($h) ID: $job.$date' >> MOCAT.successful";
      print OUT "\n\nEXECUTION END TIME:\n" . localtime() . "\n\n\n";
      print OUT "EXIT STATUS:\nEXECUTION SUCCESS\n\n";
      close OUT;
      open STATUS, ">>", "$cwd/logs/status/$job.successful";
      print STATUS "$sep\n";

      foreach my $sample (@samples)
      {
        print STATUS "$sample\t$date\tMOCAT.pl " . join( " ", @args ) . "\n";
      }
      close STATUS;
      print localtime() . ": Execution completed successfully.\n";
    } else
    {
      open OUT, ">>", "$cwd/logs/$job/commands/MOCATJob_$job.$date.command.log" or die "ERROR & EXIT: Cannot write $cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
      system "rm -f $cwd/MOCATJob_$job\_$date 2>/dev/null";
      system "rm -f $cwd/MOCATJobArray*.$job.$date.sh";
      system "mv $cwd/MOCATJob_$job.$date* $cwd/logs/$job/jobs/ 2>/dev/null";
      system "mv $cwd/$ID.po* $cwd/$ID.pe* $cwd/$ID.o* $cwd/$ID.e* $cwd/logs/other/ 2>/dev/null";
      chomp( my $cd = `date +'\%F \%T'` );
      chomp( my $h  = `hostname` );
      system "echo 'MOCAT.pl " . join( " ", @args ) . " @ $cd ($h) ID: $job.$date' >> MOCAT.failed";
      print OUT "\n\nEXECUTION END TIME:\n" . localtime() . "\n\n\n";
      print OUT "EXIT STATUS:\nEXECUTION FAIL\n\n";
      close OUT;
      die "ERROR & EXIT: Execution failed.\nLog files from this run are stored in the folder $cwd/logs/$job";
    }
  }
}

sub pre_check
{
  my $assembly_type = shift;
  my $reads         = shift;
  my $command1      = shift;
       ( $assembly_type eq "" )
    or ( $assembly_type eq 'assembly' )
    or ( $assembly_type eq 'assembly.revised' )
    or die "\nERROR & EXIT: Option -$command1 'ASSEMBLY ORIGIN' must be either of 'assembly' or 'assembly.revised'";

  if ( $reads eq '' )
  {
    die "\nERROR & EXIT: Option -r 'READS ORIGIN' must be either 'reads.processed', 'FASTA FILE' or 'DATABASE NAME'";
  }
}

sub get_kmer
{
  my $sample  = shift;
  my $read    = shift;
  my $command = shift;
  my $statsFile;
  if ( $read eq "reads.processed" )
  {
    $statsFile = "$cwd/$sample/stats/$sample.readtrimfilter.$conf{MOCAT_data_type}.stats";
  } else
  {
    $statsFile = "$cwd/$sample/stats/$sample.screen.$read.$conf{MOCAT_data_type}.stats";
    if ( !( -e $statsFile ) )
    {
      $statsFile = "$cwd/$sample/stats/$sample.screened.$read.$conf{MOCAT_data_type}.stats";
    }
  }
  unless ( -e $statsFile )
  {
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
    die <<"EOF";

ERROR & EXIT: Missing statistics file $statsFile.
Have you specified the correct database using $command 'DATABASE', 'reads.processed' or 'FASTA FILE'?

Perhaps you meant either of these:
$line
EOF

  }
  open IL, '<', $statsFile;
  my $lastline;
  while (<IL>)
  {
    chomp;
    $lastline = $_;
  }
  close IL;
  my @tokens = split( /\t/, $lastline );
  return ( $tokens[2], $tokens[3], $tokens[4] );
}

sub fake_data_init
{
  my $move     = $_[0];
  my $movename = $move;
  $movename =~ s/.index//;
  print localtime() . ": Copying $movename into memory...";
  chomp( my $user = `whoami` );
  system "mkdir -p /dev/shm/$user/MOCAT/data; rsync -avL $data_dir/$move* /dev/shm/$user/MOCAT/data/ 2>/dev/null >/dev/null";
  $data_dir_bkup = $data_dir;
  $data_dir      = "/dev/shm/$user/MOCAT/data";
  print " OK!\n";
  print localtime() . ": Data dir is now $data_dir\n";
}

sub fake_data_done
{
  my $move     = $_[0];
  my $movename = $move;
  $movename =~ s/.index//;
  print localtime() . ": Removing $movename from memory...";
  chomp( my $user = `whoami` );
  system "rm -r /dev/shm/$user/MOCAT/data/$move*";
  $data_dir = $data_dir_bkup;
  print " OK!\n";
  print localtime() . ": Data dir is now $data_dir\n";
}

sub print_settings
{
  my $text      = $_[0];
  my $maxLength = 0;
  my $printed   = 0;
  for my $key ( sort keys %conf )
  {
    if ( $key =~ m/^$text\_/ )
    {
      if ( $maxLength < length($key) )
      {
        $maxLength = length($key);
      }
    }
  }
  for my $key ( sort keys %conf )
  {
    if ( $key =~ m/^$text\_/ )
    {
      $printed = 1;
      print "$key";
      print " " x ( $maxLength - length($key) );
      print " : $conf{$key}\n";
    }
  }
  if ( $conf{MOCAT_prompt_before_run} eq 'yes' )
  {
    print "Please double check settings, and press <enter>.";
    chomp( my $key = <STDIN> );
  }
  if ($printed)
  {
    print "$sep\n";
  }
}

sub print_settings_OUT
{
  my $text      = $_[0];
  my $maxLength = 0;
  for my $key ( sort keys %conf )
  {
    if ( $key =~ m/^$text\_/ )
    {
      if ( $maxLength < length($key) )
      {
        $maxLength = length($key);
      }
    }
  }
  for my $key ( sort keys %conf )
  {
    if ( $key =~ m/^$text\_/ )
    {
      print OUT "$key";
      print OUT " " x ( $maxLength - length($key) );
      print OUT " : $conf{$key}\n";
    }
  }
}

sub record_script_version
{
  my $cmdargs = join( ' ', @args );
  my $now = localtime();
  print "Storing md5sum information about .ncbi.map .motu.map.gz .pl .pm .R files\n";
  print "$sep\n";

  #chomp(my @stat = `stat $data_dir/*.ncbi.map* $data_dir/*.motu.map* $data_dir/*.coord* $data_dir/*.functional.map* $scr_dir/*.pl $scr_dir/*.pm $scr_dir/*.sh $scr_dir/MOCATProfiling*.pm 2>/dev/null | grep -P 'File:|Size:|Modify:' | sed ':a;N;$!ba;s/\n/ /g' | sed 's/File:/\nFile:/g' | sort`);

  my @md5sums = `md5sum $data_dir/*.ncbi.map $data_dir/*.ncbi.map.gz md5sum $data_dir/*.motu.map $data_dir/*.motu.map.gz $scr_dir/*.pl $scr_dir/*.pm $conf{MOCAT_DEV_DIR}/*.pm $conf{MOCAT_DEV_DIR}/*.pl $conf{MOCAT_DEV_DIR}/*.R $conf{MOCAT_DEV_DIR}/analyze/*.R $conf{MOCAT_DEV_DIR}/ResistanceScreen/*.pl $conf{MOCAT_DEV_DIR}/VirulenceScreen/*.pl 2>/dev/null`;
  system "mkdir -p $cwd/logs/file_versions/scripts $cwd/logs/file_versions/logs $cwd/logs/file_versions/MOCAT_IDs";

  open X, '>>', "$cwd/logs/file_versions/logs/MOCAT_session.$date.log" or die "ERROR & EXIT: Cannot write $cwd/logs/file_versions/logs/MOCAT_session.$date.log";
  open Y, '>>', "$cwd/logs/file_versions/MOCAT_IDs/tmp.$date"          or die "ERROR & EXIT: Cannot write $cwd/logs/file_versions/MOCAT_IDs/tmp.$date";
  print X << "EOF";

$sep
                   MOCAT - Metagenomics Analysis Toolkit                 $version
 by Jens Roat Kultima, Luis Pedro Coelho, Shinichi Sunagawa @ Bork Group, EMBL
$sep
COMMAND EXECUTED:     MOCAT.pl $cmdargs
JOB IDENTIFIER:       $date
CURRENT DIRECTTORY:   $cwd
EXECUTION START TIME: $now
MOCAT ID:             $MOCAT_ID
SCRIPT VERSIONS:
EOF

  foreach my $md5sum (@md5sums)
  {
    chomp $md5sum;
    $md5sum =~ m/(\S+)\s+(\S+)/;
    print X "$md5sum\n";
    print Y "$md5sum\n";
    unless ( -e "$cwd/logs/file_versions/scripts/$1" )
    {
      system "cp $2 $cwd/logs/file_versions/scripts/$1";
    }
  }
  close X;
  close Y;
  my $mc_id = `md5sum $cwd/logs/file_versions/MOCAT_IDs/tmp.$date`;
  chomp $mc_id;
  $mc_id =~ m/(\S+)\s+(\S+)/;
  system "mv $cwd/logs/file_versions/MOCAT_IDs/tmp.$date $cwd/logs/file_versions/MOCAT_IDs/$1";
  $MOCAT_ID = "$1";
  print "MOCAT ID:          $MOCAT_ID\n";
  print "MOCAT version:     $INTERNAL_MOCAT_VERSION\n";
  print "MOCAT description: $MOCAT_DESCRIPTION\n";
  print "MOCAT desc file:   $MOCAT_DESCRIPTION_FILE\n";

  print "$sep\n";
  system "echo -e \"$MOCAT_ID\t$INTERNAL_MOCAT_VERSION\t$INTERNAL_TAXONOMIC_VERSION\t$INTERNAL_FUNCTIONAL_VERSION\t$INTERNAL_MOTU_VERSION\t$INTERNAL_GENE_VERSION\t$INTERNAL_RESISTANCE_VERSION\t$INTERNAL_VIRULENCE_VERSION\t$INTERNAL_X_VERSION\t$INTERNAL_Y_VERSION\t$INTERNAL_Z_VERSION\t$INTERNAL_A_VERSION\t$INTERNAL_B_VERSION\t$INTERNAL_C_VERSION\t$INTERNAL_GENERAL_VERSION\t$INTERNAL_NCBI_VERSION\t$MOCAT_DESCRIPTION\" >> $MOCAT_DESCRIPTION_FILE 2>/dev/null";
}

sub checkAndReturnDB
{
  my $db = shift;
  my @db = @{$db};
  if ( scalar @db > 1 )
  {
    ( $db[0] =~ m/(.*)\.1/ ) or die "ERROR & EXIT: database name '$db[0]' must be on the format NAME.1";
    my $base    = $1;
    my $counter = 0;
    foreach my $db (@db)
    {
      $counter++;
      ( $db eq "$base.$counter" ) or die "ERROR & EXIT: Expected '$db' to be '$base.$counter',\nplease rename database(s) on the format NAME.1 NAME.2 ... NAME.X";
    }
    return "$base.1-$counter";
  }
  return $db[0];
}

1;
