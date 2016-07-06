package MOCATProfiling::Threading;
use strict;
use warnings;
use MOCATProfiling::Variables;
use Storable;
#use File::Sync qw(fsync sync);

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2014
# This code is released under GNU GPL v3.

sub merge {
	my $file = shift;
	if ($VERBOSE) {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [QUEUE] merging $file.HASH(_LEVELS)\n";
	}

	# we use this array because this is the code we had from the Threading.pm module
	my @returned;
	sync();
	open my $IN, "<", "$file.HASH";
	fsync($IN);
	close $IN;
	open $IN, "<", "$file.HASH_LEVELS";
	fsync($IN);
	close $IN;

	$returned[6] = retrieve("$file.HASH");
	$returned[7] = retrieve("$file.HASH_LEVELS");

	# HASH
	my $LAST = 1;
	if ($CALCULATE_HORIZONTAL_COVERAGE) {
		$LAST = 2;    # If we calculate horizontal coverage, the last entry shouldn't be summarized
		foreach my $key1 ( keys %{ $returned[6] } ) {
			my $string1 = @{ ${ $returned[6] }{$key1} }[ $INIT_NUM + 1 ];
			my $string2 = @{ $HASH{$key1} }[ $INIT_NUM + 1 ];
			unless ( length $string2 ) {
				$string2 = "0" x 100;
			}
			unless ( length $string1 ) {
				$string1 = "0" x 100;
			}
			for my $i ( 0 .. 99 ) {
				if ( substr( $string1, $i, 1 ) eq '1') {
					substr($string2, $i, 1) = "1"; 
				}
			}
			@{ $HASH{$key1} }[ $INIT_NUM + 1 ] = $string2;
		}
	}
	foreach my $key1 ( keys %{ $returned[6] } ) {
		my @array = @{ ${ $returned[6] }{$key1} };
		for my $i ( 2 .. scalar @array - $LAST ) {
			if ( $array[$i] ) {
				@{ $HASH{$key1} }[$i] += $array[$i];
			}
		}
	}

	# HASH_LEVELS
	foreach my $key1 ( keys %{ $returned[7] } ) {
		foreach my $key2 ( keys %{ $returned[7]{$key1} } ) {
			my @array = @{ ${ $returned[7]{$key1} }{$key2} };
			for my $i ( 2 .. scalar @array - 1 ) {
				if ( $array[$i] ) {
					@{ $HASH_LEVELS{$key1}{$key2} }[$i] += $array[$i];
				}
			}
		}
	}

	unlink "$file.HASH"        or warn "ERROR: Could not unlink $file.HASH: $!\n";
	unlink "$file.HASH_LEVELS" or warn "ERROR: Could not unlink $file.HASH_LEVELS: $!\n";
	if ($VERBOSE) {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [QUEUE] finished merging $file.HASH(_LEVELS)\n";
	}
}

sub parseThreadResult {
	my $returned_ref = shift;
	my $i            = shift;
	my @returned     = @$returned_ref;

	if ($VERBOSE) {
		print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : [THREAD " . ( $i + 1 ) . "] parsing final thread results\n";
	}

	# return values
	$TOTAL_LENGTH                                   += $returned[0];
	$TOTAL_LENGTH_MAPPED                            += $returned[1];
	$reference_for_read_do_not_exists_in_coord_file += $returned[4];
	$reference_for_read_exists_in_coord_file        += $returned[5];

	# TOTAL_MAPPED
	foreach my $key1 ( keys %{ $returned[2] } ) {
		foreach my $key2 ( keys %{ ${ $returned[2] }{$key1} } ) {
			$TOTAL_MAPPED{$key1}{$key2} += ${ ${ $returned[2] }{$key1} }{$key2};
		}
	}

	# TOTAL_MAPPED_LENGTH_NORM
	foreach my $key1 ( keys %{ $returned[3] } ) {
		foreach my $key2 ( keys %{ ${ $returned[3] }{$key1} } ) {
			$TOTAL_MAPPED_LENGTH_NORM{$key1}{$key2} += ${ ${ $returned[3] }{$key1} }{$key2};
		}
	}
}

sub signal_handler {
	print(">>> Terminating <<<\n");
	$TERMINATE = 1;
}

sub initThreads {
	my @initThreads;
	for ( my $i = 1 ; $i <= $THREADS ; $i++ ) {
		push( @initThreads, $i );
	}
	return @initThreads;
}

1;
