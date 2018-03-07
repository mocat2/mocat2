package Thread::Deadlock;

# Make sure we have version info for this module
# Make sure we do everything by the book from now on

our $VERSION : unique = '0.01';
use strict;

# Make sure we have threads
# Make sure we can lock
# Make sure signals will have END executed

use threads;
use threads::shared ();
use sigtrap qw(die normal-signals);

# Make sure we can cluck
# Report flag (so only (n)one thread will report)
# Report from each thread

use Carp ();
my $report : shared = 'STDERR';
my %report : shared;

# Save current coderefs

my $cond_wait = \&threads::shared::cond_wait;
my $cond_signal = \&threads::shared::cond_signal;
my $cond_broadcast = \&threads::shared::cond_broadcast;

# Install hi-jacked coderefs, we can't do lock() yet ;-(

*threads::shared::cond_wait =
 sub (\[$@%]) { _remember( 'cond_wait()' ); goto &$cond_wait };
*threads::shared::cond_signal =
 sub (\[$@%]) { _remember( 'cond_signal()' ); goto &$cond_signal };
*threads::shared::cond_broadcast =
 sub (\[$@%]) { _remember( 'cond_broadcast()' ); goto &$cond_broadcast };

# Satisfy -require-

1;

#---------------------------------------------------------------------------

# class methods

#---------------------------------------------------------------------------
#  IN: 1 class (ignored)
#      2 name of file to write to (no change)
# OUT: 1 current setting

sub report { $report = $_[1] if @_ >1; $report } #report

#---------------------------------------------------------------------------

# internal routines

#---------------------------------------------------------------------------
#  IN: 1 what we're remembering

sub _remember {

# Obtain the thread we're in
# Obtain the stacktrace
# Remove this call
# Add what we're remembering
# Indicate start of thread if appropriate

    my $tid = threads->tid;
    my @cluck = split( m#(?<=$/)#,Carp::longmess() );
    shift( @cluck );
    $cluck[0] =~ s#.*?called#"Thread $tid: ".shift#e;
    $cluck[-1] =~ s#eval \{\.\.\.\} called#thread started#;

# Initialize local report
# For all of the lines,
#  Remove the thread information (if any)
#  Add to the report for this thread
# Save the report in the shared hash

    my $report;
    foreach (@cluck) {
        s/, thread #(\d+)$//;
        $report .= $_;
    }
    $report{$tid} = $report;
} #_remember

#---------------------------------------------------------------------------

# routines for standard perl features

#---------------------------------------------------------------------------
#  IN: 1 class (ignored)
#      2 file to which to save (optional: default)

sub import { goto &report } #import

#---------------------------------------------------------------------------

END {

# Attempt to lock the report flag
# Return now if we don't need to report

    lock( $report );
    return unless $report;

# Initialize handle to write to
# If we have the default value
#  Set to write to standard error

    my $handle;
    if ($report eq 'STDERR') {
        $handle = *STDERR;
    } elsif ($report eq 'STDOUT') {
        $handle = *STDOUT;

# Elseif successful in opening it as a file (no action)
# Else (not successful in opening file)
#  Set to use standard error
#  And let the world know why

    } elsif (open( $handle,'>',$report )) {
    } else {
        $handle = *STDERR;
	print $handle <<EOD;
Could not report to $report ($!)
Writing to STDERR instead
EOD
    }

# Tell the world what it is
# For all of the keys in the report
#  Tell the world about it
# Indicate that no-one else needs to report

    print $handle '*** '.__PACKAGE__." report ***\n";
    foreach (sort {$a <=> $b} keys %report) {
        print $handle "$report{$_}\n";
    }
    $report = '';
} #END

#---------------------------------------------------------------------------

__END__

=head1 NAME

Thread::Deadlock - report deadlocks with stacktrace

=head1 SYNOPSIS

    perl -MThread::Deadlock program            # report to STDERR
    perl -MThread::Deadlock=filename program   # report to file

    use Thread::Deadlock;                      # report to STDERR
    use Thread::Deadlock 'filename';           # report to file

    Thread::Deadlock->report( 'filename' );    # report to file
    Thread::Deadlock->report( '' );            # do not report

=head1 DESCRIPTION

                  *** A note of CAUTION ***

 This module only functions on Perl versions 5.8.0 and later.
 And then only when threads are enabled with -Dusethreads.  It
 is of no use with any version of Perl before 5.8.0 or without
 threads enabled.

                  *************************

The Thread::Deadlock module allows you to find out B<where> your threaded
application may be deadlocked.  It does B<not> prevent any deadlocks, nor is
it capable of resolving any deadlocks.

If you use the Thread::Deadlock module, all occurences of cond_wait(),
cond_signal() and cond_broadcast() in the source are checkpointed to
remember where it was exactly in your source and where it was in the execution
stack.  When your program finishes (either as intended or after you killed the
program, e.g. by pressing Control-C), then a report will be generated for each
thread, indicating where each thread had its last checkpoint.  By default, this
report is written to STDERR, but can be redirected to a file of your choice.

=head1 CLASS METHODS

There is one class method.

=head2 report

 Thread::Deadlock->report( 'filename' );  # write to specific file

 Thread::Deadlock->report( '' );          # disable reporting

 $report = Thread::Deadlock->report;      # obtain current setting

The "report" class method returns the current setting for the thread
checkpoint report.  It can also be used to set the name of the file to which
the report will be written, or to disable report writing altogether.

=head1 CAVEATS

This module was originally conceived as hi-jacking the core lock() function.
However, this proved to be impossible, at least with Perl 5.8.0.  It was
therefore re-written to hi-jack the cond_wait(), cond_signal() and
cond_broadcast() routines from threads::shared.pm.  This is not exactly the
same, but since most deadlocking problems are caused by mixups of cond_wait()
and cond_signal()/cond_broadcast(), this seems to be as good a solution.

=head1 AUTHOR

Elizabeth Mattijsen, <liz@dijkmat.nl>.

Please report bugs to <perlbugs@dijkmat.nl>.

=head1 COPYRIGHT

Copyright (c) 2002 Elizabeth Mattijsen <liz@dijkmat.nl>. All rights
reserved.  This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=head1 SEE ALSO

L<threads>, L<threads::shared>.

=cut
