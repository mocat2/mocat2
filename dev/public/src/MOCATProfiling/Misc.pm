package MOCATProfiling::Misc;
use strict;
use warnings;
use MOCATProfiling::Variables;
use File::Copy;
use File::Basename;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2014
# This code is released under GNU GPL v3.

sub getLoggingTime {
	my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime(time);
	my $nice_timestamp = sprintf( "%04d%02d%02d %02d:%02d:%02d", $year + 1900, $mon + 1, $mday, $hour, $min, $sec );
	return $nice_timestamp;
}

sub rename_or_die {
	my $orig = shift;
	my $dest = shift;
	move( $orig, $dest ) or die("ERROR & EXIT: Move failed ['$orig' -> '$dest']: $!");
}

sub zipFiles {
	my $output = shift;

	# zip, we used to move, but now we zip (note that using gz and tar isn't good because zip has a uilt in index for fetching files faster than tar, which has to process entire file)
	system_ ("rm -f $output && cd $temp_file.output/ && $bin_dir/zip -q -y -D -4 $output * ");
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : saved final output file: $output\n";
	system_("echo '$version' > $output.version");
	print STDERR "(" . MOCATProfiling::Misc::getLoggingTime() . ") : saved final output version file: $output.version\n";
}

sub mkdir_or_die {
	my $dir = shift;
	( -d $dir ) or system "mkdir -p $dir";
	( -d $dir ) or die "Could not create $dir.\nPlease check your permissions.";
}

sub system_ {
	my $cmd = shift;
	system($cmd) == 0 or die("ERROR & EXIT: system($cmd) failed: $!\n");
}

sub checkAndReturnDB {
	my $db = shift;
	my @db = @{$db};
	if (scalar @db > 1) {
		($db[0] =~ m/(.*)\.1.coord/) or die "ERROR & EXIT: database name '$db[0]' must be on the format NAME.1.coord";
		my $base = $1;
		my $basename = basename ($db[0]);
		$basename =~ m/(.*)\.1\.coord/;
		my $basebasename = $1;
		my $counter = 0;
		foreach my $db (@db) {
			$counter++;
			($db eq "$base.$counter.coord") or die "ERROR & EXIT: Expected '$db' to be '$base.$counter.coord',\nplease rename database coord files to the format NAME.1.coord NAME.2.coord ... NAME.X.coord";
		}
		return ("$base.1-$counter.coord");
	}
	return $db[0];
}


1;
