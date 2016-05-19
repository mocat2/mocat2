package MOCATProfiling::MainQueuer;
use strict;
use warnings;
use MOCATProfiling::Variables;
use Storable;

# this is the routine started for each thred
sub mainQueuer {
	#############################################################################################################################################
	# DEFINE VARIABLES
	#############################################################################################################################################
	my $line;
	$prev_insert = "0";
	my $blocks      = 1;
	my $fileCounter = 1;
	my $counter;
	my %sendHash;    # this is used for sending a hash to the threads
	my $fileOut;
	my $j;
	my %sendHashLevels;
	#############################################################################################################################################
	#############################################################################################################################################
	# OPEN FILE
	#############################################################################################################################################
	#| $zip_execute -c
	open $fileOut, " | gzip -c -1 > $temp_file.job.$fileCounter.data" or die "ERROR & EXIT: Cannot write to OUTPUT PIPE | gzip -c -1 > $temp_file.job.$fileCounter.data";
	if ($VERBOSE) {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [QUEUE] line 1: write to $temp_file.job.$fileCounter.data\n";
	}
	#############################################################################################################################################

	#############################################################################################################################################
	# MAIN LOOP FOR SUBMITTING JOBS
	#############################################################################################################################################
	while ( $line = <$READ> ) {
		my @currentLine = split "\t", $line;
		if ( $input_file_format eq 'BAM' || $input_file_format eq 'SAM' ) {
			$ref_id = $currentLine[2];
		}
		elsif ( $input_file_format eq 'SOAP' ) {
			$ref_id = $currentLine[7];
		}
		if ( $mode eq 'NCBI' ) {
			$ref_id =~ m/^(\d+)\..*/;
			chomp $line;
			$line .= "\t" . $MAP{$1} . "\n";
		}
		$counter = 0;
		while ( exists $HASH{"$ref_id.$counter"} ) {
			$sendHash{"$ref_id.$counter"} = $HASH{"$ref_id.$counter"};
			$counter++;
		}

		$currentLine[0] =~ s/\/[12]$//;
		if ( $currentLine[0] eq $prev_insert || $. == 1 ) {    # $. is the line number
			print {$fileOut} $line;
		}
		else {
			if ( $blocks == $BLOCKSIZE ) {
				close $fileOut or die "ERROR & EXIT: Could not close $fileOut";
				store \%sendHash, "$temp_file.job.$fileCounter.data.hash" or die "ERROR & EXIT: Error saving sendHash to $temp_file.job.$fileCounter.data.hash";
				$request_q->enqueue("$temp_file.job.$fileCounter.data");
				if ($VERBOSE) {
					print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [QUEUE] line $.: QUEUED $temp_file.job.$fileCounter.data\n";
				}
				%sendHash       = ();
				%sendHashLevels = ();
				$fileCounter++;
				my $waiting  = $request_q->pending();
				my $waiting2 = $return_q->pending();
				if ($VERBOSE) {
					print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [QUEUE] line $.: start writing to $temp_file.job.$fileCounter.data [$waiting waiting for submission | $waiting2 waiting for merging]\n";
				}
				open $fileOut, " | gzip -c -5 > $temp_file.job.$fileCounter.data" or die "ERROR & EXIT: Cannot write to OUTPUT PIPE | gzip -c -5 > $temp_file.job.$fileCounter.data.gz";
				$blocks = 1;
				print {$fileOut} $line;
				$counter = 0;
				while ( exists $HASH{"$ref_id.$counter"} ) {
					$sendHash{"$ref_id.$counter"} = $HASH{"$ref_id.$counter"};
					$counter++;
				}

				#############################################################################################################################################
				# MERGE RESULTS IF AVAILABLE
				#############################################################################################################################################
				if ( $waiting >= $THREADS ) {
					while ( $request_q->pending() > 1 ) {
						my $file = $return_q->dequeue();
						MOCATProfiling::Threading::merge($file) if $file;
					}
				}
				else {
					if ($VERBOSE) {
						print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [QUEUE] $waiting queued. Require $THREADS queued before merging HASH(_LEVELS)\n";
					}
				}
				#############################################################################################################################################

				#############################################################################################################################################
				# WAIT TO NOT HAVE TOO MANY JOBS
				#############################################################################################################################################
				my $pending = $request_q->pending();
				if ( $pending > ( $THREADS * 2 ) ) {    # we don't want to fill up the queue with requests too much, uses a lot of memory
					my $death_counter = 0;
					my $pending       = $request_q->pending();
					while ( $pending > ($THREADS) ) {
						if ($VERBOSE) {
							print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [QUEUE] waiting... ($death_counter)\n";
						}
						my $file = $return_q->dequeue();
						if ($file) {
							MOCATProfiling::Threading::merge($file);
						}
						else {
							sleep 60;
							$death_counter++;
						}
						if ( $death_counter >= 60 ) {
							die "ERROR & EXIT: Waited 60 minutes for any thread to finish, and then gave up. Please change -blocksize to a lower number or just restart the program.";
						}
						$pending = $request_q->pending();
					}
				}
				#############################################################################################################################################
			}
			else {
				$blocks++;
				print {$fileOut} $line;
			}
		}
		$prev_insert = $currentLine[0];
	}
	#############################################################################################################################################

	#############################################################################################################################################
	# PROCESS THE LAST LINE TO ENQUEUE
	#############################################################################################################################################
	my $waiting  = $request_q->pending();
	my $waiting2 = $return_q->pending();
	if ($VERBOSE) {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [QUEUE] line $. (last): start writing to $temp_file.job.$fileCounter.data [$waiting waiting for submission | $waiting2 waiting for merging]\n";
	}
	close $fileOut or die "ERROR & EXIT: Could not close $fileOut";
	store \%sendHash, "$temp_file.job.$fileCounter.data.hash" or die "ERROR & EXIT: Error saving sendHash to $temp_file.job.$fileCounter.data.hash";
	if ( $mode eq 'mOTU' || $mode eq 'NCBI' ) {
		store \%sendHashLevels, "$temp_file.job.$fileCounter.data.hashlevels" or die "ERROR & EXIT: Error saving sendHash to $temp_file.job.$fileCounter.data.hashlevels";
	}
	$request_q->enqueue("$temp_file.job.$fileCounter.data");    # do the last line
	if ($VERBOSE) {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [QUEUE] line $. (last): QUEUED $temp_file.job.$fileCounter.data\n";
	}

	# Signal to threads that there is no more work.
	$request_q->enqueue(undef) for 1 .. ( $THREADS * 2 );
	if ($VERBOSE) {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [QUEUE]: QUEUED " . ( $THREADS * 2 ) . " undef jobs to close threads\n";
	}
	#############################################################################################################################################

	#############################################################################################################################################
	# END
	#############################################################################################################################################
	# clear hashes from memory;
	%sendHash = ();
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [QUEUE] FINISHED ENQUEUING ALL LINES\n";
	return (1);
	#############################################################################################################################################
}

1;
