package MOCATProfiling::Load;
use strict;
use warnings;
use MOCATProfiling::Variables;
use threads;
use threads::shared;
use Storable;

#############################################################################################################################################
# LOAD MAP FILES
#############################################################################################################################################
sub loadMapFiles {
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Loading $map_file\n";
	my $MAP;
	if ( -e "$map_file.gz" ) {
		open $MAP, "gunzip -c $map_file.gz | " or die;
	}
	else {
		open $MAP, "<", "$map_file" or die;
	}
	my $count = 0;
	while (<$MAP>) {
		next if m/^#/;
		chomp;
		my $line = $_;
		my @line = split "\t";
		if ( $mode eq 'mOTU' ) {
			if ( scalar @line != 4 ) {
				die "ERROR & EXIT: map file ($map_file) does not have correct number of fields. Perhaps different mode or map file should be specified?";
			}
			@{ $HASH{"$line[0].0"} }[$INIT_NUM] = "$line[2].$line[3]";

			# for mOTU we store the length in field 1, we'll see how this works out...
			@{ $HASH_LEVELS{'mOTU'}{"$line[2].$line[3]"} }[0] = 1;
			@{ $HASH_LEVELS{'mOTU'}{"$line[2].$line[3]"} }[1] += $line[1];    # length
		}
		if ( $mode eq 'NCBI' ) {
			if ( scalar @line != 10 ) {
				die "ERROR & EXIT: map file ($map_file) does not have correct number of fields. Perhaps different mode or map file should be specified?";
			}
			@{ $HASH_LEVELS{'kingdom'}{ $line[1] } }[ 0 .. 1 ]       = ( 1, 1 );
			@{ $HASH_LEVELS{'phylum'}{ $line[2] } }[ 0 .. 1 ]        = ( 1, 1 );
			@{ $HASH_LEVELS{'class'}{ $line[3] } }[ 0 .. 1 ]         = ( 1, 1 );
			@{ $HASH_LEVELS{'order'}{ $line[4] } }[ 0 .. 1 ]         = ( 1, 1 );
			@{ $HASH_LEVELS{'family'}{ $line[5] } }[ 0 .. 1 ]        = ( 1, 1 );
			@{ $HASH_LEVELS{'genus'}{ $line[6] } }[ 0 .. 1 ]         = ( 1, 1 );
			@{ $HASH_LEVELS{'species'}{ $line[7] } }[ 0 .. 1 ]       = ( 1, 1 );
			@{ $HASH_LEVELS{'specI_cluster'}{ $line[8] } }[ 0 .. 1 ] = ( 1, 1 );
			@{ $HASH_LEVELS{'taxaid'}{ $line[0] } }[ 0 .. 1 ]        = ( 1, 1 );
			$NCBI_TAXAID_LENGTH{ $line[0] } = $line[9];
			$MAP{ $line[0] }                = $line;
			$ALL_TAXA_MEAN_LENGTH += $line[9];
			$count++;
		}
	}
	if ( $mode eq 'NCBI' ) {
		$ALL_TAXA_MEAN_LENGTH = $ALL_TAXA_MEAN_LENGTH / $count;
	}
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Finished loading $map_file\n";
	close $MAP or die "ERROR & EXIT: Closing pipe failed, this means something is wrong. Corrput files?";
}
#############################################################################################################################################

#############################################################################################################################################
# LOAD COORD FILES
#############################################################################################################################################
sub loadCoordFile {

	my $thread = shift;
	my $i      = $thread - 1;

	# Load coordinates
	# this is a bit fancy, if we have multiple files, MOCAT filtering also creates a reference.gz file with the references that are in that map
	# file, this way we can load only the coords that we need for this specific thread
	my $text = " [THREAD $thread]";
	if ( $thread == 0 || $thread == -1 ) {
		$text = "";
	}

	# [0] because right now we only support one input file!
	my %IDs;
	my $use_ref = 0;
	if ( -e "$input_files[0].references.gz" && $thread != 0 ) {    # if the thread is zero, it means we want all entries
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") :$text Using reference file $input_files[0].references.gz to decrease coord hash size\n";
		open my $IN, "gunzip -c $input_files[0].references.gz |";
		while (<$IN>) {
			chomp;
			$IDs{$_} = 1;
		}
		close $IN or die "ERROR & EXIT: Closing pipe failed, this means something is wrong. Corrput files?";
		$use_ref = 1;
	}

	# used for checking if we have previously created HASH for these coord files
	my $name = MOCATProfiling::Misc::checkAndReturnDB( \@coord_files );
	$name = "$name.MOCAT_DATA";

	# we have stored the data form a previous run
	if ( -e "$name.HASH" && -e "$name.avg" ) {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Loading coordinates from MOCAT DATA file $name.HASH\n";
		%HASH                    = %{ retrieve("$name.HASH") };
		$db_average_entry_length = ${ retrieve("$name.avg") };
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Finished loading $name.HASH\n";

	}

	# we have not stored the data before
	else {
		foreach my $coord_file (@coord_files) {
			my $IN;
			if ( -e "$coord_file.gz" ) {
				open $IN, "gunzip -c $coord_file.gz | " or die;
				print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") :$text Loading coordinates from $coord_file.gz\n";

			}
			else {
				open $IN, "<", "$coord_file" or die;
				print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") :$text Loading coordinates from $coord_file\n";

			}
			my $prev_id = "0";
			my $last_counter;    # this is used for keeping track of the numbering when threads=0
			my $tot_counter = 0;
			while (<$IN>) {
				next if m/^#/;
				(m/^(\S+)\s+(\d+)\s+(\d+)$/) or die "ERROR & EXIT: $coord_file has line:\n$_\nExpected the line to be of format: NONSPACECHARS\\tDIGIT\\tDIGIT";
				my $counter = 0;
				if ( ( $use_ref == 1 && $IDs{$1} ) || $use_ref == 0 ) {    # we only load the entry into hash if it exists in IDs, if we have an ID hash
					if ( $thread != 0 ) {
						while ( exists $HASH{"$1.$counter"} ) { ++$counter; }
						my @array;                             # : shared;
						@array = ( $2, $3 );
						$HASH{"$1.$counter"} = \@array;
					}
					else {                                           #  thread=0 means we are repopulating the hash, but we only want to add it if it doesn't exist already
						if ( $prev_id eq $1 ) {
							$last_counter++;
							@{ $HASH{"$1.$last_counter"} }[ 0 .. 1 ] = ( $2, $3 );
							$counter++;

						}
						else {
							@{ $HASH{"$1.$counter"} }[ 0 .. 1 ] = ( $2, $3 );
						}
					}
					$db_average_entry_length += $3 - $2 + 1;
					$tot_counter++;
				}
				$prev_id      = $1;
				$last_counter = $counter;
			}
			$db_average_entry_length = $db_average_entry_length / $tot_counter;
			print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") :$text Finished loading $. coordinates from $coord_file\n";
			close $IN or die "ERROR & EXIT: Closing pipe failed, this means something is wrong. Corrput files?";
		}

		# store data for fast access
		my $name = MOCATProfiling::Misc::checkAndReturnDB( \@coord_files );
		$name = "$name.MOCAT_DATA";
		if ( -e "$name.HASH.tmp" && "$name.avg.tmp" ) {
			print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") :$text $name.HASH.tmp already exist, skipped creating $name.HASH\n";
		}
		else {
			print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") :$text saving coordinates to MOCAT DATA file $name.HASH\n";
			store \%HASH,                    "$name.HASH.tmp" or die "ERROR & EXIT: Cannot save HASH for future use in $name.HASH.tmp";
			store \$db_average_entry_length, "$name.avg.tmp"  or die "ERROR & EXIT: Cannot save avg for future use in $name.HASH.tmp";
			system_("mv $name.HASH.tmp $name.HASH && mv $name.avg.tmp $name.avg");
			print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") :$text saved coordinates\n";
		}
	}

	%IDs = ();    # this is probably unnecessary, perl will delete it anyway, but why not...
}
#############################################################################################################################################

#############################################################################################################################################
# Load stats file
#############################################################################################################################################
sub loadStatsFile {
	if ( -e "$stats_file" ) {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Loading stats file $stats_file\n";
		open my $IN, "$stats_file" or die "ERROR & EXIT: Missing stats file $stats_file";
		my $line;
		$line = <$IN>;
		$line = <$IN>;
		close $IN or die "ERROR & EXIT: Closing pipe failed, this means something is wrong. Corrput files?";
		chomp $line;
		my @line = split /\t/, $line;
		@STATS_FILE_DATA                      = @line;
		$TOTAL{'base'}                        = $line[1];
		$TOTAL{'insert'}                      = $line[5];
		$TOTAL{'base.only.unique'}            = $line[1];
		$TOTAL{'base.mm.dist.among.unique'}   = $line[1];
		$TOTAL{'insert.only.unique'}          = $line[5];
		$TOTAL{'insert.mm.dist.among.unique'} = $line[5];
	}
	else {

		# I decided not to have any support if the stats file is missing, it causes too much head ache with divisions by 0 and sclaed files won't work
		# but for now I'll keep it
		print STDERR "====== ERROR ====== : Missing stats file $stats_file - continuing anyway, some stats will be incorrect\n";
		@STATS_FILE_DATA = ( 0, 0, 0, 0, 0, 0, 0, 0, 0 );
		$TOTAL{'base'}   = 0;
		$TOTAL{'insert'} = 0;
		$TOTAL{'base.only.unique'}            = 0;
		$TOTAL{'base.mm.dist.among.unique'}   = 0;
		$TOTAL{'insert.only.unique'}          = 0;
		$TOTAL{'insert.mm.dist.among.unique'} = 0;
	}
}
#############################################################################################################################################

#############################################################################################################################################
# SYSTEM COMMAND
#############################################################################################################################################
sub system_ {
	my $cmd = shift;
	system($cmd) == 0 or die("ERROR & EXIT: system($cmd) failed: $!\n");
}
#############################################################################################################################################

#############################################################################################################################################
# LOAD FUNCITONAL LEVELS
#############################################################################################################################################
sub loadFunctLevels {

	my $levels;
	@LEVELS = ();
	my %levels;
	my %levels2;
	my @levels;
	my $MAP;
	if ( -e "$map_file.gz" ) {
		open $MAP, "gunzip -c $map_file.gz | " or die;
	}
	else {
		open $MAP, "<", "$map_file" or die;
	}
	$levels = <$MAP>;
	chomp $levels;
	$levels =~ s/#gene\s+//;
	@levels = split /\s+/, $levels;
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Loading functional levels from map file header: " . ( join " ", @levels ) . "\n";
	foreach my $k (@levels) {
		$levels{$k} = 1;
	}
	close $MAP;
	if (@CONFIG_LEVELS) {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Loading functional levels from config file:";
		foreach my $l (@CONFIG_LEVELS) {
			unless ( $l eq 'gene' ) {
				$levels2{$l} = 1;
				print STDERR " $l";
			}
		}
		print "\n";
		foreach my $k ( keys %levels2 ) {
			if ( $levels{$k} && $k ne 'gene' ) {

				push @LEVELS, $k;
			}
		}
	}
	else {
		foreach my $k ( keys %levels ) {
			push @LEVELS, $k;
		}
	}
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : Set functional levels (except gene) to: " . join( " ", @LEVELS ) . "\n";
	if ( scalar @LEVELS < 1 ) {
		die "ERROR & EXIT: There are no levels that are defined in both the map file and the config";
	}
}
#############################################################################################################################################

1;
