package Smash::Databases::MetaGenomeDB::Loader;

use strict;
use warnings;

our @ISA = qw(Smash::Core Exporter);

use File::Path;
use Smash::Databases::MetaGenomeDB::MetaGenomeLoader;
use Smash::Databases::MetaGenomeDB::AssemblyLoader;
use Smash::Databases::MetaGenomeDB::GenePredictionLoader;

our @EXPORT_OK = qw(load_smash);

############################################
# Constructor
############################################

sub new {
	my $class  = shift;
	my %params = @_;

	bless {
		%params
	}, $class;
}

sub _overload {shift->{LOADERTYPE}}

sub get_correct_type {
	my $type  = shift;
	my @types = qw(MetaGenome Assembly GenePrediction);
	my $correct;
	TYPE: foreach (@types) {
		if (lc($_) eq lc($type)) {
			$correct = $_;
			last TYPE;
		}
	}
	if (!$correct) {
		die "Smash does not know how to load $type!";
	}
	$correct  = "${correct}Loader";
	return $correct;
}

sub load_smash {
	my $type    = shift;
	my %options = @_;

	# Get the right loader

	my $class = get_correct_type($type);
	$options{loadertype} = $class;

	# Hash keys used in making Loader objects are upper-case equivalents of the lower case options
	# obtained from command line. So the hash keys are uppercased here.

	my $loader = "Smash::Databases::MetaGenomeDB::$class"->new(map {uc($_) => $options{$_}} keys %options);

	# run the (un)loader
	$loader->init();
	if (!$options{unload}) {
		$loader->run();
	} else {
		$loader->unload_db();
		if ($options{wipeout}) {
			$loader->wipe_out();
		}
	}
	$loader->finish();
}

1;
