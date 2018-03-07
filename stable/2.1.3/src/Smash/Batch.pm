package Smash::Batch;

use strict;
use warnings;
use Fcntl;
use File::Temp;
use Smash::Global qw(:all);

sub new {
	my $class   = shift;
	my %options = @_;
	my $this    = bless {%options}, $class;

	# Calculate the prefix

	my $prefix  = $this->label;

	# output

	my $outfile = "$prefix.out";
	my $OUTFH;
	open($OUTFH, ">$outfile") || die "Cannot open $outfile: $!";
	$this->{OUTFILE}  = $outfile;
	$this->{OUTFILE_FH} = $OUTFH;

	my $rerunfile = "$prefix.rerun";
	my $RERUNFH;
	open($RERUNFH, ">$rerunfile") || die "Cannot open $rerunfile: $!";
	$this->{RERUNFILE}= "$prefix.rerun";
	$this->{RERUNFILE_FH} = $RERUNFH;

	$this->{PROGRESS} = "$prefix.progress";

	$this->{INDENT_STR} = " " unless $this->{INDENT_STR};
	
	printf STDERR "Smash::Batch instance: Output of this run will be reported in $outfile\n";
	printf STDOUT "Smash::Batch instance: Output of this run will be reported in $outfile\n";
	return $this;
}

sub label      {shift->{LABEL}}
sub outfile    {shift->{OUTFILE}}
sub outfile_fh {shift->{OUTFILE_FH}}
sub rerunfile  {shift->{RERUNFILE}}
sub rerunfile_fh  {shift->{RERUNFILE_FH}}
sub progress   {shift->{PROGRESS}}
sub indent_str {shift->{INDENT_STR}}

sub process_batch {
	my $this  = shift;
	my $flat  = shift;
	my @batch = split(/\n+/, $flat);
	my $outfile_fh = $this->outfile_fh;
	print $outfile_fh "<SmashBatch>\n";
	my $quit;
	$this->{VARIABLES}->{SMASHBINDIR} = $SMASH_SCRIPT_LOCATION;
	STEP: foreach my $step (@batch) {

		# Transform this step with current variable substitutions

		next STEP if $step =~ /^\s*#/ || $step =~ /^\s*$/;
		$step = $this->transform($step);

		# If there was error before, write the remaining commands to be processed

		if ($quit) {
			$this->write_rerun_command($step);
			next STEP;
		}

		my ($mode, @rest) = split(/\s+/, $step);
		my $command;
		if ($mode eq "SET") {
			my $target  = shift @rest;
			$this->{VARIABLES}->{$target} = join(" ", @rest);
			$this->write_rerun_command($step);
		} elsif ($mode eq "CAPTURE") {
			my $target  = shift @rest;
			my $script  = shift @rest;
			   $script  =~ s|^[\./]||;
			   $script  = "$SMASH_SCRIPT_LOCATION/$script";
			   $command = join(" ", $script, @rest);
			my $capture = $this->get_output($command);
			if ($capture) {
				$this->{VARIABLES}->{$target} = $capture;
				$this->write_rerun_command("SET $target $capture");
			} else {
				$quit = "$command failed to generate output\n";
				$this->write_rerun_command($step);
			}
		} elsif ($mode eq "SUCCEED") { # Successfully run a Smash command
			my $script  = shift @rest;
			   $script  =~ s|^[\./]||;
			   $script  = "$SMASH_SCRIPT_LOCATION/$script";
			   $command = join(" ", $script, @rest);
			my $status  = $this->check_success($command);
			if ($status == 1) {
				$this->write_rerun_command("#$command");
			} else {
				$quit = "$command failed with status $status\n";
				$this->write_rerun_command($step);
			}
		} elsif ($mode eq "SYSTEM") { # Successfully run a system command
			   $command = join(" ", @rest);
			my $status  = $this->run_system_command($command);
			if ($status == 0) {
				$this->write_rerun_command("#$command");
			} else {
				$quit = "$command failed with status $status\n";
				$this->write_rerun_command($step);
			}
		} elsif ($mode eq "EXIT") {
			print $outfile_fh "\tExiting\n";
			last STEP;
		}
	}
	if ($quit) {
		die $quit;
	}
	print $outfile_fh "</SmashBatch>\n";
}

sub write_rerun_command {
	my $this = shift;
	my $string = shift;
	my $fh = $this->rerunfile_fh;
	print $fh "$string\n";
}

=item C<transform($string)>

transforms a given string by replacing the current values of all variables.

=cut

sub transform {
	my $this   = shift;
	my $string = shift;
	while (my ($key, $value) = each %{$this->{VARIABLES}}) {
		$string =~ s/\@${key}@/$value/g;
	}
	return $string;
}

=item C<check_success($command)>

runs the command and checks if it runs successfully

=item C<get_output($command)>

runs the command and captures and returns the output

=item C<run_system_command($command)>

runs the given command without capturing anything - STDOUT and STDERR
go to the STDOUT and STDERR of the calling script

=cut

sub check_success {
	my $this    = shift;
	my $command = shift;
	my $output  = $this->run_smash_step($command);
	return 1 if $output and lc($output) eq "success";
	return 0;
}

sub get_output {
	my $this    = shift;
	my $command = shift;
	return $this->run_smash_step($command);
}

# No output/error capturing!

sub run_system_command {
	my $this    = shift;
	my $command = shift;
	print "$command\n";
	return system($command);
}

=item C<run_smash_step($command)>

runs the given command and prints progress to the output file. It creates the
following block in the batch output file:

<SmashStep command='$command'>
<CapturedStdOut>
...
</CapturedStdOut>
<CapturedStdErr>
...
</CapturedStdErr>
<ReturnValue>
...
</ReturnValue>
</SmashStep>

=cut

sub run_smash_step {
	use File::Temp;

	my $this    = shift;
	my $command = shift;

	my $output;

	my $STEP_ERR     = new File::Temp();

	$this->output_indented("<SmashStep>", 1);
	$this->output_indented("<Command>", 2);
	$this->output_indented($command, 3);
	$this->output_indented("</Command>", 2);
	$this->output_indented("<StartTime>", 2);
	$this->output_indented(scalar(localtime), 3);
	$this->output_indented("</StartTime>", 2);
	$this->output_indented("<CapturedStdOut>", 2);

	# Run the command and start capturing

	print "$command\n";

	open(STEPOUT, "$command 2>$STEP_ERR |") || die "Cannot open pipe for step: $!";

	# Process STDOUT

	LINE: while (<STEPOUT>) {
		$this->output_indented("$_", 3);
		if (m|<output>(.*)</output>|) {
			$output = $1;
		}
	}
	$this->output_indented("</CapturedStdOut>", 2);

	# Process STDERR

	seek $STEP_ERR, SEEK_SET, 0;
	$this->output_indented("<CapturedStdErr>", 2);
	while (<$STEP_ERR>) {
		$this->output_indented("$_", 3);
	}
	$this->output_indented("</CapturedStdErr>", 2);

	$this->output_indented("<ReturnValue>", 2);
	$this->output_indented($output || "NULL", 3);
	$this->output_indented("</ReturnValue>", 2);
	$this->output_indented("<FinishTime>", 2);
	$this->output_indented(scalar(localtime), 3);
	$this->output_indented("</FinishTime>", 2);
	$this->output_indented("</SmashStep>", 1);
	return $output;
}

sub die_gracefully {
	my $this    = shift;
}

sub output_indented {
	my $this = shift;
	$this->print_indented($this->outfile_fh, @_);
}

sub print_indented {
	my $this = shift;
	my ($FH, $message, $indent) = @_;

	$indent = 0 unless $indent;
	chomp($message);

	printf $FH "%s%s\n", $this->indent_str x $indent, $message;
}

sub output_indented_multiline {
	my $this = shift;
	$this->print_indented_multiline($this->outfile_fh, @_);
}

sub print_indented_multiline {
	my $this = shift;
	my ($FH, $message, $indent) = @_;

	$indent = 0 unless $indent;
	chomp($message);
	$message =~ s/\t/  /g;
	
	# Calculate how much can be written per line

	my $line_max = 80 - $indent*length($this->indent_str) - 1; # subtract 1 for the prefix character

	SUBSTR:while (1) {
		my $substr = substr($message, 0, $line_max);
		printf $FH "%s|%s\n", $this->indent_str x $indent, $substr;
		if (length($message) > $line_max) {
			$message = substr($message, $line_max);
		} else {
			last SUBSTR;
		}
	}
}

1;
