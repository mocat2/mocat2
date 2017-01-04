package MOCATUsage2;
use warnings;
use strict;
use POSIX;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

sub nearest {
	my $targ = abs(shift);
	my @res  = map {
		if ( $_ >= 0 )
		{
			$targ * int( ( $_ + 0.50000000000008 * $targ ) / $targ );
		}
		else {
			$targ * POSIX::ceil( ( $_ - 0.50000000000008 * $targ ) / $targ );
		}
	} @_;

	return (wantarray) ? @res : $res[0];
}

sub usage {
	my $temp_dir = MOCATCore::determine_temp( ( 200 * 1024 * 1024 ) );
	chomp( my $dateS = `date +'\%F \%T'` );
	my $cmd  = $_[0];
	my $ID   = $_[1];
	my $job  = $_[2];
	my $date = $_[3];
	my @ERR;
	my $error = 0;
	$| = 1;
	my $pid = fork();

	if ( not defined $pid ) {
		print "resources not avilable.\n";
	}
	elsif ( $pid == 0 ) {
		@ERR = `$cmd`;
		my @errors;
		foreach my $e (@ERR) {
			my $one;
			if ( $e =~ m/.*exit code (\d+)./ ) {
				$one = $1;
			}
			if ( $e =~ m/Job.*SIGKILL/ ) {
				$one = "999";
			}
			if ( $e =~ m/Unable to run job/ ) {
				$one = "999";
			}
			if ( defined $one ) {
				if ( $one != 0 ) {
					if ( $e =~ m/Job \d+ / ) {
						push @errors, 1;
					}
					elsif ( $e =~ m/Job \d+\.(\d+) / ) {
						push @errors, $1;
					}
					$error = 1;
				}
			}
		}
		open WRITER0, ">$cwd/logs/other/$job\_$date.errcode";
		print WRITER0 "$error\n" . join( "\n", @errors );
		close WRITER0;
		exit(0);

	}
	else {
	}
	my %ever_max_mem;
	my %ever_max_pmem;
	my %ever_max_cpu;
	my %qhost;
	my %oldhost;
	my %main_id;
	my %kill_id;
	my %temp_size;
	my %ever_frac_mem_used;

	foreach my $sample (@samples) {
		$ever_max_mem{$sample}       = 0;
		$ever_max_pmem{$sample}      = 0;
		$ever_max_cpu{$sample}       = 0;
		$temp_size{$sample}          = 0;
		$ever_frac_mem_used{$sample} = 0;
		$qhost{$sample}              = "NOT STARTED";
	}
	chomp( my $current_host = `hostname` );
	chomp( my $current_user = `whoami` );

	my $number = 0;
	my %fixes;
	for my $key ( sort keys %conf ) {
		if ( $key =~ m/^realtime_status_fix_/ ) {
			$number++;
		}
	}
	for my $i ( 1 .. $number ) {
		my @temp_info = split( /\|/, $conf{"realtime_status_fix_$i"} );
		$fixes{ $temp_info[0] } = $temp_info[0];
	}

	while ( -e "/proc/$pid/status" && !( `grep 'State:' /proc/$pid/status` =~ m/Z/ ) ) {
		my $str = "$sep\n";
		$str = $str . "                   MOCAT - Metagenomics Analysis Toolkit                $version\n";
		$str = $str . "  by Jens Roat Kultima, Luis Pedro Coelho, Shinichi Sunagawa @ Bork Group, EMBL\n$sep\n";
		chomp( my $dateX = `date +'\%F \%T'` );
		$str = $str . "Executed MOCAT.pl " . join( ' ', @args ) . " @ $cwd\n";
		$str = $str . "Currently running $job job $date on $current_host:$ID, screen updated at $dateX\n$sep\n";
		my $max = scalar @samples;
		for my $nr ( 1 .. scalar @samples ) {
			my $percent = nearest( 1, 100 * $nr / $max );
			if ( $conf{realtime_status_print} eq 'yes' ) {
				system "echo -ne \"\r[Updating process data: $percent% complete]\"";
			}
			my $sample = $samples[ $nr - 1 ];

			my $qhost = $qhost{$sample};
			if ( $qhost{$sample} eq "NOT STARTED" ) {
				my $cmd = "qstat -t";
				if ( $run_on_host ne "" ) {
					$cmd = "ssh $run_on_host \"cd $cwd && $cmd\"";
				}
				my @qstat = `$cmd`;
				foreach my $qstat (@qstat) {
					$qstat =~ s/^\s+//;
					my @line = split( /\s+/, $qstat );
					if ( $line[2] && $line[9] ) {
						if ( $line[2] eq $ID && $line[9] eq $nr ) {
							$qhost = $line[7];
							$qhost =~ s/.*@//;
							$qhost{$sample} = $qhost;
							if ( $fixes{$current_host} ) {
								if ( $conf{realtime_status_log} eq 'yes' ) {
									print "Fixing, host $current_host exists in log\ncurrent host: $current_host\nqhost: $qhost\nfixed qhost: $fixes{$current_host}\n";
								}
								$qhost{$sample} = $fixes{$current_host};
								$qhost = $fixes{$current_host};
							}
						}
					}
				}
			}

			my @possible_matches;
			my @lines;
			if ( $qhost ne $current_host && $qhost ne "FINISHED" && $qhost ne "NOT STARTED" && !$main_id{$sample} ) {
				my $cmd1 = "ssh $qhost \"pstree -ap \" 2>/dev/null ";
				my $cmd2 = "| grep 'MOCATJob_$job.$date.$nr.sh'";
				if ( $run_on_host ne "" ) {
					$cmd = "ssh $run_on_host \"cd $cwd && $cmd1\" $cmd2 ";
				}
				else {
					$cmd = "$cmd1 $cmd2";
				}
				@possible_matches = `$cmd`;
				if ( $conf{realtime_status_log} eq 'yes' ) {
					print "$cmd\n";
					print "MATCHES:\n";
					foreach my $r (@possible_matches) {
						print ": $r";
					}
				}
				my $child;
				foreach my $m (@possible_matches) {
					if ( $m =~ m/.*-MOCATJob_.*,(\d+)\s+.*/ ) {
						$child = $1;
						my $cmd = "ssh $qhost \"ps -o \"pgrp,pid\" $child\" 2>/dev/null";
						if ( $run_on_host ne "" ) {
							$cmd = "ssh $run_on_host \"cd $cwd && $cmd\"";
						}
						if ( $conf{realtime_status_log} eq 'yes' ) {
							print "line = $cmd\n";
						}
						my @line = `$cmd`;
						if ( scalar @line >= 2 ) {
							unless ( $line[-2] ) {
								$line[-2] = '';
							}
						}
						if ( $line[-2] =~ m/PGRP/ ) {
							$line[-1] =~ s/^\s+//;
							my @lines = split /\s+/, $line[-1];
							$main_id{$sample} = $lines[0];
							$kill_id{$sample} = $child;
							if ( $conf{realtime_status_log} eq 'yes' ) {
								print "Child $child : Parent $main_id{$sample}\n";
							}
						}
						last;
					}
				}
			}
			elsif ( $qhost eq $current_host && $qhost ne "FINISHED" && $qhost ne "NOT STARTED" && !$main_id{$sample} ) {
				my $cmd1 = "pstree -ap 2>/dev/null";
				my $cmd2 = "| grep 'MOCATJob_$job.$date.$nr.sh'";
				if ( $run_on_host ne "" ) {
					$cmd = "ssh $run_on_host \"cd $cwd && $cmd1\" $cmd2";
				}
				else {
					$cmd = "$cmd1 $cmd2";
				}
				@possible_matches = `$cmd`;
				if ( $conf{realtime_status_log} eq 'yes' ) {
					print "$cmd\n";
					print "MATCHES:\n";
					foreach my $r (@possible_matches) {
						print ": $r";
					}
				}
				foreach my $m (@possible_matches) {
					if ( $m =~ m/.*-MOCATJob_.*,(\d+)\s+.*/ ) {
						my $child = $1;
						my $cmd   = "ps -o \"pgrp,pid\" $child";
						if ( $run_on_host ne "" ) {
							$cmd = "ssh $run_on_host \"cd $cwd && $cmd\"";
						}
						if ( $conf{realtime_status_log} eq 'yes' ) {
							print "line = $cmd\n";
						}
						my @line = `$cmd`;
						if ( $line[-2] =~ m/PGRP/ ) {
							$line[-1] =~ s/^\s+//;
							my @lines = split /\s+/, $line[-1];
							$main_id{$sample} = $lines[0];
							$kill_id{$sample} = $child;
							if ( $conf{realtime_status_log} eq 'yes' ) {
								print "Child $child : Parent $main_id{$sample}\n";
							}
						}
						last;
					}
				}
			}

			my @memory_usage;
			if ( $qhost ne $current_host && $main_id{$sample} ) {
				my $cmd1 = "ssh $qhost \"ps -eo \'pgrp,pid,pcpu,pmem,vsize,comm,stat\' \" 2>/dev/null";
				my $cmd2 = "| grep -P '^(\\s+|)$main_id{$sample}'";
				if ( $run_on_host ne "" ) {
					$cmd = "ssh $run_on_host \"cd $cwd && $cmd1\" $cmd2";
				}
				else {
					$cmd = "$cmd1 $cmd2";
				}
				@lines = `$cmd`;
				if ( $conf{realtime_status_log} eq 'yes' ) {
					print "$cmd\n";
				}
				$cmd = "ssh $qhost \"free\" 2>/dev/null";
				if ( $run_on_host ne "" ) {
					$cmd = "ssh $run_on_host \"cd $cwd && $cmd\"";
				}
				@memory_usage = `$cmd`;
				if ( $conf{realtime_status_log} eq 'yes' ) {
					print "$cmd\n";
				}
			}
			elsif ( $qhost eq $current_host && $main_id{$sample} ) {
				my $cmd1 = "ps -eo \"pgrp,pid,pcpu,pmem,vsize,comm,stat\" 2>/dev/null";
				my $cmd2 = "| grep -P '^(\\s+|)$main_id{$sample}'";
				if ( $run_on_host ne "" ) {
					$cmd = "ssh $run_on_host \"cd $cwd && $cmd1\" $cmd2";
				}
				else {
					$cmd = "$cmd1 $cmd2";
				}
				if ( $conf{realtime_status_log} eq 'yes' ) {
					print "$cmd\n";
				}
				@lines = `$cmd`;
				$cmd   = "free 2>/dev/null";
				if ( $run_on_host ne "" ) {
					$cmd = "ssh $run_on_host \"cd $cwd && $cmd\"";
				}
				@memory_usage = `$cmd`;
				if ( $conf{realtime_status_log} eq 'yes' ) {
					print "$cmd\n";
				}
			}

			my $total;
			my $used;
			my $cached;
			my $frac_mem_used = 0;
			if ( $main_id{$sample} ) {
				foreach my $mem (@memory_usage) {
					if ( $mem =~ m/^Mem:/ ) {
						my @mem = split /\s+/, $mem;
						$total         = $mem[1];
						$used          = $mem[2];
						$cached        = $mem[6];
						$frac_mem_used = ( $used - $cached ) / ($total);
						$frac_mem_used = $frac_mem_used * 100;
						if ( $frac_mem_used > 99 ) {
							system "echo 'MOCAT KILLS $sample:PID=$kill_id{$sample} because free memory on $qhost is less than 1%' >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log";
							my $cmd = "kill -15 $kill_id{$sample}";
							if ( $run_on_host ne "" ) {
								$cmd = "ssh $run_on_host \"cd $cwd && $cmd\"";
							}
							system "$cmd";
						}
						$frac_mem_used = nearest( .1, $frac_mem_used );
					}
				}
			}

			if ( scalar @lines == 0 && $qhost ne "NOT STARTED" && $ever_max_cpu{$sample} > 0 && $qhost ne "FINISHED" ) {
				$oldhost{$sample} = $qhost;
				$qhost            = "FINISHED";
				$qhost{$sample}   = $qhost;
			}

			$str = $str . $sample . " \@ $qhost ";
			my $mem_tot  = 0;
			my $pmem_tot = 0;
			my $cpu_tot  = 0;
			my $str2     = "";
			for my $line (@lines) {
				$line =~ s/^\s+//;
				my @line  = split /\s+/, $line;
				my $cpu   = $line[2];
				my $pmem  = $line[3];
				my $mem   = $line[4];
				my @names = split '/', $line[5];
				my $name  = substr( $names[-1], 0, 6 );
				$pmem_tot = $pmem_tot + $pmem;
				$cpu_tot  = $cpu_tot + $cpu;
				$mem_tot  = $mem_tot + $mem;
				my $state = $line[6];
				my $fix   = "k";

				if ( $mem > 1024 ) {
					$mem = $mem / 1024;
					$fix = "M";
					if ( $mem > 1024 ) {
						$mem = $mem / 1024;
						$fix = "G";
					}
				}
				if ( $fix eq "k" ) {
					$mem = nearest( 1, $mem );
				}
				if ( $fix eq "M" ) {
					$mem = nearest( 1, $mem );
				}
				if ( $fix eq "G" ) {
					$mem = nearest( .1, $mem );
				}
				$cpu = nearest( 1, $cpu );
				if ( $name =~ m/^bash$/ || $name =~ m/^sh$/ ) {
				}
				else {
					$str2 = $str2 . "$name($state $mem$fix $cpu%) ";
				}
			}

			if ( $mem_tot > $ever_max_mem{$sample} ) {
				$ever_max_mem{$sample} = $mem_tot;
			}
			if ( $cpu_tot > $ever_max_cpu{$sample} ) {
				$ever_max_cpu{$sample} = $cpu_tot;
			}
			if ( $pmem_tot > $ever_max_pmem{$sample} ) {
				$ever_max_pmem{$sample} = $pmem_tot;
			}
			if ( $frac_mem_used > $ever_frac_mem_used{$sample} ) {
				$ever_frac_mem_used{$sample} = $frac_mem_used;
			}
			my $max_mem     = $ever_max_mem{$sample};
			my $max_mem_fix = "k";
			if ( $max_mem > 1024 ) {
				$max_mem     = $max_mem / 1024;
				$max_mem_fix = "M";
				if ( $max_mem > 1024 ) {
					$max_mem     = $max_mem / 1024;
					$max_mem_fix = "G";
				}
			}
			my $mem = $mem_tot;
			my $fix = "k";
			if ( $mem > 1024 ) {
				$mem = $mem / 1024;
				$fix = "M";
				if ( $mem > 1024 ) {
					$mem = $mem / 1024;
					$fix = "G";
				}
			}
			if ( $fix eq "k" ) {
				$mem = nearest( 1, $mem );
			}
			if ( $fix eq "M" ) {
				$mem = nearest( 1, $mem );
			}
			if ( $fix eq "G" ) {
				$mem = nearest( .1, $mem );
			}
			$mem_tot = $mem;

			if ( $max_mem_fix eq "k" ) {
				$max_mem = nearest( 1, $max_mem );
			}
			if ( $max_mem_fix eq "M" ) {
				$max_mem = nearest( 1, $max_mem );
			}
			if ( $max_mem_fix eq "G" ) {
				$max_mem = nearest( .1, $max_mem );
			}

			my $cmd = "du -s $temp_dir/$sample/temp/ 2>/dev/null";

			# We assume that all hosts have access to the same filesystem
			#if ( $run_on_host ne "" ) {
			#	$cmd = "ssh $run_on_host \"cd $cwd && $cmd\"";
			#}
			my $temp_sizes = `$cmd`;
			my @temp_sizes = split /\s+/, $temp_sizes;
			my $temp_size  = $temp_sizes[0];
			unless ($temp_size) {
			   $temp_size = 0;
			}
			if ( $temp_size > $temp_size{$sample} ) {
				$temp_size{$sample} = $temp_size;
			}
			my $temp_fix = "k";
			if ( $temp_size > 1024 ) {
				$temp_size = $temp_size / 1024;
				$temp_fix  = "M";
			}
			if ( $temp_size > 1024 ) {
				$temp_size = $temp_size / 1024;
				$temp_fix  = "G";
			}
			if ( $temp_fix eq "k" ) {
				$temp_size = nearest( 1, $temp_size );
			}
			if ( $temp_fix eq "M" ) {
				$temp_size = nearest( 1, $temp_size );
			}
			if ( $temp_fix eq "G" ) {
				$temp_size = nearest( .1, $temp_size );
			}
			my $temp_size_now = $temp_size;
			my $temp_fix_now  = $temp_fix;
			$temp_size = $temp_size{$sample};
			$temp_fix  = "k";
			if ( $temp_size > 1024 ) {
				$temp_size = $temp_size / 1024;
				$temp_fix  = "M";
			}
			if ( $temp_size > 1024 ) {
				$temp_size = $temp_size / 1024;
				$temp_fix  = "G";
			}
			if ( $temp_fix eq "k" ) {
				$temp_size = nearest( 1, $temp_size );
			}
			if ( $temp_fix eq "M" ) {
				$temp_size = nearest( 1, $temp_size );
			}
			if ( $temp_fix eq "G" ) {
				$temp_size = nearest( .1, $temp_size );
			}

			$str2 =~ s/ $//;
			$cpu_tot = nearest( 1, $cpu_tot );
			$ever_max_cpu{$sample} = nearest( 1, $ever_max_cpu{$sample} );
			my ( $l1, $l2 );
			if ( $temp_fix_now eq $temp_fix ) {
				$l1 = "$temp_size_now($temp_size)$temp_fix";
			}
			else {
				$l1 = "$temp_size_now$temp_fix_now($temp_size$temp_fix)";
			}
			if ( $fix eq $max_mem_fix ) {
				$l2 = "$mem_tot($max_mem)$max_mem_fix";
			}
			else {
				$l2 = "$mem_tot$fix($max_mem$max_mem_fix)";
			}
			$str2 = "[$l1] [$l2 $pmem_tot($ever_max_pmem{$sample})% $frac_mem_used($ever_frac_mem_used{$sample})% $cpu_tot($ever_max_cpu{$sample})%] [$str2]\n";
			$str  = $str . $str2;
		}
		unless ( $conf{realtime_status_log} eq 'yes' ) {
			if ( $conf{realtime_status_print} eq 'yes' ) {
				system "clear";
			}
		}
		if ( $conf{realtime_status_print} eq 'yes' ) {
			print $str;
		}
		if ( $conf{realtime_status_print} eq 'yes' ) {
			system "echo -ne \"\r[Sleeping $realtime_status seconds]\"";
		}
		sleep($realtime_status);
	}
	system "mkdir -p $cwd/logs/resources";
	open OUT,  ">>", "$cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
	open OUT2, ">>", "$cwd/logs/resources/used_resources_for_jobs.txt";
	print OUT "\n\nUSED RESOURCES:\n";
	my $counter = 0;
	chomp( my $dateF = `date +'\%F \%T'` );
	foreach my $sample (@samples) {
		my $tot    = $ever_max_mem{$sample};
		my $totfix = "k";
		if ( $tot > 1024 ) {
			$tot    = $tot / 1024;
			$totfix = "M";
			if ( $tot > 1024 ) {
				$tot    = $tot / 1024;
				$totfix = "G";
			}
		}
		my $temp_size = $temp_size{$sample};
		my $temp_fix  = "k";
		if ( $temp_size > 1024 ) {
			$temp_size = $temp_size / 1024;
			$temp_fix  = "M";
		}
		if ( $temp_size > 1024 ) {
			$temp_size = $temp_size / 1024;
			$temp_fix  = "G";
		}
		$tot       = nearest( .1, $tot );
		$temp_size = nearest( .1, $temp_size );

		unless ( $oldhost{$sample} ) {
			$oldhost{$sample} = $qhost{$sample};
		}

		print OUT "Resources for $job and sample $sample @ $oldhost{$sample} mem: $tot $totfix mem: $ever_max_pmem{$sample}% max_frac_mem_used: $ever_frac_mem_used{$sample}% cpu: $ever_max_cpu{$sample}% temp_size: $temp_size $temp_fix\n";
		print OUT2 "$job\t$date\t$sample\t$current_host\t$oldhost{$sample}\t$tot $totfix\t$ever_max_pmem{$sample} %\t$ever_frac_mem_used{$sample} %\t$ever_max_cpu{$sample} %\t$temp_size $temp_fix\t$dateS\t$dateF\t" . join( " ", @args ) . "\n";
	}
	close OUT;
	close OUT2;
	open WRITER0, "<$cwd/logs/other/$job\_$date.errcode";
	$error = <WRITER0>;
	open OUTx,   ">>", "$cwd/logs/$job/commands/MOCATJob_$job.$date.command.log";
	open STATUS, ">>", "$cwd/logs/status/$job.failed";
	print STATUS "$sep\n";
	print OUTx "\n";
	my %failed;

	while (<WRITER0>) {
		chomp;
		my $sample = $samples[ $_ - 1 ];
		$failed{$sample} = 1;

		my @E;
		if ( $date eq $sample ) {
			@E = `tail -3 $cwd/logs/$job/samples/MOCATJob_$job.$date.log`;
			print "\n\nSAMPLE $sample FAILED. Here are the last 3 lines of the log file\n($cwd/logs/$job/samples/MOCATJob_$job.$date.log):\n$sep\n";
		}
		else {
			@E = `tail -3 $cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log`;
			print "\n\nSAMPLE $sample FAILED. Here are the last 3 lines of the log file\n($cwd/logs/$job/samples/MOCATJob_$job.$sample.$date.log):\n$sep\n";
		}
		print OUTx "FAILED SAMPLE: $sample\n";
		print STATUS "$sample\t$date\tMOCAT.pl " . join( " ", @args ) . "\n";
		foreach my $e (@E) {
			print $e;
		}
		print "$sep\n\n";
	}
	print OUTx "\n";
	close OUTx;
	close STATUS;
	close WRITER0;
	open STATUS, ">>", "$cwd/logs/status/$job.successful";
	print STATUS "$sep\n";
	foreach my $sample (@samples) {

		unless ( $failed{$sample} ) {
			print STATUS "$sample\t$date\tMOCAT.pl " . join( " ", @args ) . "\n";
		}
	}
	close STATUS;

	system "rm -r $cwd/logs/other/$job\_$date.errcode";
	return $error;
}

1;
