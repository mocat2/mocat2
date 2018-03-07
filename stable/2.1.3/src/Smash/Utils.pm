package Smash::Utils;
use strict;
use warnings;

# clip/trim coordinate conversion

############################################
# Forge format: ($start, $end) as clip coordinates means: 
# (in 0-based coordinates)
# [0,$start] + [$start+1,$end-1] + [$end,]
# [0,$start]         is vector,   length = $start+1
# [$start+1, $end-1] is sequence, length = $end-$start-1
# [$end,]            is vector,   length = $seq_length - $end
# maui format: ($start, $length) means:
# [0,$start-1] + [$start,$start+$length-1] + [$start+$length,]
#
# case Pass:
# 	depending on the mode, print the right coordinates
# case Fail:
#	if mode=quality or both, DONT print a line
#	else leave them in tact, so print a line
############################################

sub forge2maui {
	my $this   = shift;
	my $input  = shift;
	my $output = shift;
	my $mode   = shift || "both";

	open(IN, "<$input") || die "Cannot open $input: $!";
	open(OUT, ">$output") || die "Cannot open $output: $!";
	while (<IN>) {
		chomp();
		my @words = split(/\s/);
		if ($words[0] eq "Pass") {
			my ($def, $start, $end);
			if ($mode eq "quality") {
				($def, $start, $end) = @words[1,2,3];
			} elsif ($mode eq "both") {
				($def, $start, $end) = @words[1,6,7];
			} else { # vector or lucy
				($def, $start, $end) = @words[1,4,5];
			}
			printf OUT "%s\t%d\t%d\t%d\n", $def, 1, $start+1, $end-$start-1;
		} elsif ($words[0] eq "Fail") {
			if ($mode ne "quality" && $mode ne "both") { # The whole sequence remains
				my $def = $words[1];
				printf OUT "%s\t%d\t%d\t%d\n", $def, 1, 0, -1;
			}
		}
	}
	close(OUT);
	close(IN);
}

############################################
# Execute a given command. 
# Returns the return status of the process
# Ignores output and error streams
############################################

sub execute {
	my $this = shift;
	my ($command) = @_;
	my $status;

	#print "# $command\n";
	$status = system($command);
	if ($status != 0) {
		warn "Error executing command:\n\t$command\nError Status:\n\t$status\n";
	}
	return $status;
}

############################################
# Execute a given command. 
# Returns the output as a string
# Ignores error stream
############################################

sub execute_capture {
	my $this      = shift;
	my ($command) = @_;
	my $output;

	#warn "# $command\n";
	$output = `$command`;
	return $output;
}

1;
