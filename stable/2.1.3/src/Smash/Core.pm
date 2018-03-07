############################################
#  Core methods for Smash
############################################

package Smash::Core;
use strict;
use Carp;
use Config;
use DBI;
use FAlite;
use File::Path;
use File::Basename;
use Smash::Global qw(:all);
use overload '""' => '_overload';

use base 'Exporter';
our @EXPORT_OK = qw(suicide safe_open_file get_median max2);
our %EXPORT_TAGS = ('all' => [@EXPORT_OK]); 

############################################
############################################
##    Member variables through methods    ##
############################################
############################################

=head1 NAME

Smash::Core - Basic interface to the Smash library

=head1 SYNOPSIS

	use Smash::Core;
	my $smash = new Smash::Core(GENEPRED => "MC1.MG1.AS1.GP1");
	$smash->init();

	my $collection = $smash->collection; # "MC1"
	my $metagenome = $smash->metagenome; # "MC1.MG1"
	my $assembly   = $smash->assembly;   # "MC1.MG1.AS1"
	my $genepred   = $smash->genepred;   # "MC1.MG1.AS1.GP1"

	my ($dir, $dbh);

	$dir = $smash->read_dir($metagenome);
	$dir = $smash->assembly_dir($assembly);
	$dir = $smash->genepred_dir($genepred);
	
	my $contig_id   = $smash->get_id_by_name("contig", "MC2.MG1.AS1.C34");
	my $lib_id      = $smash->get_id_by_name("library", "MC2.BAAU");

	my $fasta = new FAlite(\*STDIN);
	if (my $entry = $fasta->nextEntry) {
		my $gc = $smash->get_gc_percent($fasta->seq, 0, 100);
		print $smash->pretty_fasta($fasta->seq);
	}

	$dbh = $smash->get_db_handle();        # collection specific DB
	$dbh = $smash->get_smashdb_handle();   # general SmashDB
	$smash->finish(); # closes these two handles

=head1 DESCRIPTION

Smash::Core is the core module for the Smash perl codebase. It provides all the necessary
functionalities for most of the modules in L<Smash::Analyses>, L<Smash::Utils> and L<Smash::Databases>. 
Almost all the functionalities of Smash::Core are object-oriented: you cannot just call a function of
Smash::Core like

	Smash::Core::function($arg1, $arg2); # illegal

since it expects to be called as:

	Smash::Core->function($arg1, $arg2); # legal, yet incorrect

Moreover, calling it in an object-oriented fashion just doesn't cut it. The object must have been
created and initialized properly, like so:

	my $smash = new Smash::Core(COLLECTION => "MC1");
	$smash->init();
	$smash->function($arg1, $arg2);

When an instance of Smash::Core is created and its C<init()> method is called, it does the
following:

=over 4

1. parses the configuration file and registers all the configuration details

2. makes two database connections (SmashDB and collection-specific)

=back

When you are done with the Smash::Core object, be sure to destroy it using the mandatory method C<finish()>, like so:

	$smash->finish();

The C<finish()> method commits any pending changes to the database and closes all open database connections.

Some important roles of Smash::Core are:

=over 4

=item *

parsing the configuration file and registering all the configuration details

=item *

providing data location for all levels of data (collection, metagenome, reads, assembly, gene prediction)

=item *

providing software package locations

=item *

making database connections depending on the database engine used

=item *

querying database for internal ids using external names of sample, library, assembly, gene prediction, contig, gene, etc.

=back

There are groups of functions in Smash::Core that perform the above roles.

=head1 OBJECT CREATION AND DESTRUCTION

=over 4

=item B<new>

Returns a new Smash::Core object. Must be called as follows:

	my $smash = Smash::Core->new(GENEPRED => "MC1.MG1.AS1.GP1");

Valid keys include (but are not limited to) COLLECTION, METAGENOME, ASSEMBLY and GENEPRED.
However, one of these four must be passed to the constructor to use collection specific
databases. Calling C<new> without any of these keys will initialize connections to the
general SmashDB database and a subsequent call to C<get_db_handle> results in an error.
See L<"DBI interface with database engine">.

=item B<init>

Parses config file and initializes database connections.

=item B<finish>

Closes all open database connections.

=back

=cut

############################################
# new
############################################

sub new {
	my $class  = shift;
	my %params = @_;
	bless {
		%params
	}, $class;
}

############################################
# init
############################################

sub init {
	my $this = shift;

	# read in conf file
	$this->parse_config();
	$this->set_globals();

	# parse id's

	my $id = (shift || $this->{GENEPRED} || $this->{ASSEMBLY} || $this->{METAGENOME} || $this->{COLLECTION});
	my ($collection, $metagenome, $assembly, $genepred) = $this->parse_concat_id($id);
	$this->{GENEPRED}   = $genepred;
	$this->{ASSEMBLY}   = $assembly;
	$this->{METAGENOME} = $metagenome;
	$this->{COLLECTION} = $collection;

	####
	# File permissions for everything created by Smash
	# User and group get everything. Others get nothing
	####

	umask 0007; 
	$this->{FILE_PERM} = 0770;

	# Host
	my $host = `hostname -s`;
	chomp($host);
	$host =~ s/\..*//;
	$this->{HOSTNAME} = $host;

	# Architecture

	my $software_dir = $this->get_conf_value("Software", "software_dir");
	my $arch = `$software_dir/config.guess`;
	chomp($arch);
	$this->{ARCHITECTURE} = $arch;
	return $this;
}

sub set_globals {
	my $this = shift;
	$TAXONOMY_LOCAL_LOCATION  = $this->get_conf_value("Taxonomy", "local_repository");
	$NCBI_TAXONOMY_REMOTE_LOCATION = $this->get_conf_value("Taxonomy", "remote_repository");
	($SMASH_SCRIPT_NAME, $SMASH_SCRIPT_LOCATION) = fileparse($0);
	$SMASH_PERL = $Config{perlpath};
}

sub finish {
	my $this = shift;
	if ($this->_smashdb_handle) {
		$this->_smashdb_handle->commit();
		$this->close_smashdb_handle();
	}
	#$this->_refgenomedb_handle->commit();
	#$this->_refgenomedb_handle->disconnect();
	if ($this->{COLLECTION}) {
		if ($this->_db_handle) {
			$this->_db_handle->commit();
			$this->close_db_handle();
		}
	}
}

sub _overload {"Smash"};

############################################
# Members
############################################

=head1 MEMBER VARIABLES

These are member variables of the Smash object, once it is initialized.

=over 4

=item B<collection>

collection id

=item B<metagenome>

metagenome id

=item B<assembly>

assembly id

=item B<genepred>

genepred id

=item B<config>

hash containing the configuration

=item B<host>

short name of the execution host, set using C<`hostname -s`>

=item B<cluster>

name to identify a group of hosts with the same architecture, so a program built
on one of them can run on all of them

=back

=cut

sub common_name  {shift->{COMMON_NAME}}
sub collection   {shift->{COLLECTION}}
sub metagenome   {shift->{METAGENOME}}
sub assembly     {shift->{ASSEMBLY}}
sub genepred     {shift->{GENEPRED}}
sub config       {shift->{CONFIG}}
sub config_file  {shift->{CONFIG_FILE}}
sub file_perm    {shift->{FILE_PERM}}
sub host         {shift->{HOSTNAME}}
sub arch         {shift->{ARCHITECTURE}}
sub cluster      {shift->{CLUSTER}}

############################################
# Members
############################################

############################################
# Parser methods
############################################

=head1 FUNCTIONS

=head2 Smash-specific configuration file parsing

=over 4

=item B<parse_config>

Parses the config file and sets the relevant variables. The following locations
are searched for F<smash.conf> in the given order: current working directory 
and F<$HOME/.smash>.

=item B<get_smash_conf_value($key)>

Get a value from the config file for any 
key C<$key> under C<[Smash]> section. For example,

	$smash->get_smash_conf_value("data_dir");

returns the value for C<data_dir> under C<[Smash]> section.

=item B<get_conf_value($section, $key)>

Get a value from the config file for any 
key C<$key> under C<[$section]> section. For example,

	$smash->get_conf_value("Taxonomy", "data_dir");

returns the value for C<data_dir> under C<[Taxonomy]> section.

=back

=cut

############################################
# Config Parser
############################################

sub parse_config {
	use Env "HOME";
	use Smash::Config::ConfigParser;
	my $this = shift;

	if (defined($this->{CONFIG})) {
		return;
	}

	# Parse the config file
	my @config_files = ("smash.conf", "$HOME/.smash/smash.conf");
	my $config = {};
	my $parsed = 0;
	foreach my $file (@config_files) {
		if (-f $file) {
			my $parser = new Smash::Config::ConfigParser($file);
			$config    = $parser->parse();
			$parsed    = 1;
			$this->{CONFIG_FILE} = $file;
			last;
		}
	}

	# No config file in the places I expect

	if (!$parsed) {
		warn "FATAL ERROR:\n";
		warn "Missing config file: smash.conf\n";
		warn "I checked: ".join(", ", @config_files)."\n";
		confess "Did you forget to place it in one of these places?\n";
	}

	# [Smash] section missing in config file

	if (!$config->{"Smash"}) {
		confess "[Smash] section missing in config file\n";
	}

	# [SmashDB] section missing in config file

	if (!$config->{"SmashDB"}) {
		confess "[SmashDB] section missing in config file\n";
	}

	$this->{CONFIG} = $config;

	return $config;
}

############################################
# Get a value from the config file for any 
# key under [Smash] section
############################################

sub get_smash_conf_value {
	my $this   = shift;
	my $key    = shift;
	return $this->get_conf_value("Smash", $key);
}

sub get_conf_value {
	my $this = shift;
	my ($section, $key) = @_;
	if (!$this->config) {
		confess "Config not parsed!";
	}
	my %config = %{$this->config};
	return $config{$section}{$key};
	
}

=head2 Smash-specific string parsing

=over 4

=item B<parse_metagenome_id>

Parses a metagenome id and returns the metagenome collection id. For example,

	$smash->parse_metagenome_id("MC2.MG1");

returns C<"MC2">;


=item B<parse_assembly_id>

Parses an assembly id and returns the metagenome collection id and metagenome id. For example,

	$smash->parse_assembly_id("MC2.MG1.AS1");

returns C<("MC2", "MC2.MG1")>;

=item B<parse_genepred_id>

Parses a genepred id and returns the metagenome collection id, metagenome id and assembly id. For example,

	$smash->parse_assembly_id("MC2.MG1.AS1.GP1");

returns C<("MC2", "MC2.MG1", "MC2.MG1.AS1")>;

=back

=cut

############################################
# ID parser
############################################

sub parse_metagenome_id {
	my $this       = shift;
	my $metagenome = shift;
	my ($collection) = $this->parse_concat_id($metagenome);
	if (!defined($collection)) {
		confess "Invalid metagenome id: $metagenome";
	}
	return $collection;
}

sub parse_assembly_id {
	my $this       = shift;
	my $assembly   = shift;
	my ($collection, $metagenome) = $this->parse_concat_id($assembly);
	if (!defined($metagenome)) {
		confess "Invalid assembly id: $assembly";
	}
	return ($collection, $metagenome);
}

sub parse_genepred_id {
	my $this       = shift;
	my $genepred   = shift;
	my ($collection, $metagenome, $assembly) = $this->parse_concat_id($genepred);
	if (!defined($assembly)) {
		confess "Invalid genepred id: $genepred";
	}
	return ($collection, $metagenome, $assembly);
}

sub parse_concat_id {
	my $this       = shift;
	my $concat_id  = shift;
	my $collection = undef;
	my $metagenome = undef;
	my $assembly   = undef;
	my $genepred   = undef;
	my $collection_prefix = $this->get_smash_conf_value("collection_prefix");
	my $metagenome_prefix = $this->get_smash_conf_value("metagenome_prefix");
	my $assembly_prefix   = $this->get_smash_conf_value("assembly_prefix");
	my $genepred_prefix   = $this->get_smash_conf_value("genepred_prefix");

	my $buffer     = $concat_id;

	if ($buffer =~ /^(${collection_prefix}[0-9]+)/) {
		$collection   = $1;
	} else {
		return ($collection, $metagenome, $assembly, $genepred);
	}

	if ($buffer =~ /^$collection\.(${metagenome_prefix}[0-9]+)/) {
		$metagenome   = "$collection.$1";
	} else {
		return ($collection, $metagenome, $assembly, $genepred);
	}

	if ($buffer =~ /^$metagenome\.(${assembly_prefix}[0-9]+)/) {
		$assembly   = "$metagenome.$1";
	} else {
		return ($collection, $metagenome, $assembly, $genepred);
	}

	if ($buffer =~ /^$assembly\.(${genepred_prefix}[0-9]+)/) {
		$genepred   = "$assembly.$1";
	}
	return ($collection, $metagenome, $assembly, $genepred);
}

############################################
# Record locator
############################################

=head2 Data location

=over 4

=item B<data_dir>

Returns the location of Smash data repository.

=item B<read_dir>

Returns the raw data location for the given metagenome id.

=item B<assembly_dir>

Returns the assembled data location for the given assembly id.

=item B<genepred_dir>

Returns the predicted protein/gene data location for the given gene prediction id.

=item B<analyses_dir>

Returns the analyses data location for the given metagenome id.

=item B<get_blastdb>

Returns the full path of a given blast database in the repository. For example,

	$smash->get_blastdb("STRING7");

will return the full path like F</home/smash/data_repos/databases/STRING7>.

=item B<maui_dir>

Returns the path of the maui program installation for this host or cluster.

=back

=cut

sub data_dir {
	shift->get_smash_conf_value("data_dir");
}

sub analyses_dir {
	my $this         = shift;
	my $metagenome   = shift;
	my $data_dir     = $this->data_dir;
	my ($collection) = $this->parse_metagenome_id($metagenome);
	return "$data_dir/metagenome_collections/$collection/metagenomes/$metagenome/analyses";
}

sub bootstrap_dir {
	my $this         = shift;
	my $metagenome   = shift;
	return $this->analyses_dir($metagenome)."/bootstrap";
}

sub collection_dir {
	my $this         = shift;
	my $collection   = shift;
	my $data_dir     = $this->data_dir;
	return "$data_dir/metagenome_collections/$collection";
}

sub read_dir {
	my $this         = shift;
	my $metagenome   = shift;
	my $data_dir     = $this->data_dir;
	my ($collection) = $this->parse_metagenome_id($metagenome);
	return "$data_dir/metagenome_collections/$collection/metagenomes/$metagenome/reads";
}

sub assembly_dir {
	my $this         = shift;
	my $assembly     = shift;
	my $data_dir     = $this->data_dir;
	my ($collection, $metagenome) = $this->parse_assembly_id($assembly);
	return "$data_dir/metagenome_collections/$collection/metagenomes/$metagenome/assemblies/$assembly";
}

sub genepred_dir {
	my $this         = shift;
	my $genepred     = shift;
	my $data_dir     = $this->data_dir;
	my ($collection, $metagenome) = $this->parse_genepred_id($genepred);
	return "$data_dir/metagenome_collections/$collection/metagenomes/$metagenome/gene_predictions/$genepred";
}

sub get_blastdb {
	use File::Glob;
	my $this         = shift;
	my $db           = shift;
	my $data_dir     = $this->data_dir;
	my @files        = <$data_dir/databases/$db.*>;
	if (scalar(@files) == 0) {
		confess "Cannot find $db in $data_dir";
	}
	return "$data_dir/databases/$db";
}

sub maui_dir {
	my $this         = shift;
	my $machine      = shift;
	my $software_dir = $this->get_conf_value("Software", "software_dir");
	my $multi_arch   = $this->get_conf_value("Software", "multi_arch");
	my $maui_dir;

	if ($machine) {
		$maui_dir = "$software_dir/$machine/maui";
	} elsif ($multi_arch eq "no") {
		$maui_dir = "$software_dir/maui";
	} else {
		MACHINE:foreach my $arch ($this->cluster, $this->host, $this->arch, "generic") {
			my $dir = "$software_dir/$arch/maui";
			if (-d $dir) {
				$maui_dir = $dir;
				last MACHINE;
			}
		}
	}
	if (!$maui_dir || ! -d $maui_dir) {
		confess sprintf "Cannot find an installation of 'maui' for host '%s'", $this->host;
	}
	return $maui_dir;
}

sub software_dir {
	my $this         = shift;
	my $software     = shift;
	my $version      = shift;
	my @names        = ($software);
	my $multi_arch   = $this->get_conf_value("Software", "multi_arch");
	my $software_dir = $this->get_conf_value("Software", "software_dir");
	my $alias        = $this->get_conf_value("Software", $software);
	my $pkg_dir;
	my @locations;

	if (!$alias) {
		warn "FATAL ERROR:\n";
		confess sprintf "\tSoftware '%s' not found in your config file '%s'\n", $software, $this->config_file;
	}

	push(@names, $alias);

	if ($version eq "current") {
		my $maybe = $this->get_conf_value("Current Version", $software);
		$version = $maybe if $maybe;
	}

	if ($multi_arch eq "no") {
		 LOCATION:foreach my $location (@names) {
			$pkg_dir = "$software_dir/$location/$version";
			push(@locations, $pkg_dir);
			if (-d $pkg_dir || -l $pkg_dir) {
				last LOCATION;
			}
		}
	} else {
		ARCH:foreach my $arch ($this->cluster, $this->host, $this->arch, "generic") {
			next ARCH unless defined($arch);
			LOCATION:foreach my $location (@names) {
				$pkg_dir = "$software_dir/$arch/$location/$version";
				push(@locations, $pkg_dir);
				if (-d $pkg_dir || -l $pkg_dir) {
					last ARCH;
				}
			}
		}
	}

	# Now figure out the version as well.
	# May be if it is a link, like current pointing to some location.

	if ($version eq "current" && -l $pkg_dir) { 
		$version = readlink($pkg_dir);
	}

	if (!$pkg_dir || ! -d $pkg_dir || !$version) {
		warn "FATAL ERROR:\n";
		warn sprintf "\tUsing your config file: %s,\n", $this->config_file;
		warn sprintf "\tCannot find version '%s' of %s for host '%s'\n", $version, $software, $this->host;
		warn "\tI tried the following locations:\n";
		confess join("\n", map {"\t\t$_"} @locations)."\n";
	}
	return ($pkg_dir, $version);
}

############################################
# DBI interface
############################################

=head2 DBI interface with database engine

=over 4

=item B<get_db_handle>

Returns a database handle to the given metagenome collection.

=item B<get_smashdb_handle>

Closes the database handle to the SmashDB meta-database.

=item B<last_db_insert_id>

Returns the autoincrement value of an autoincrement primary key in collection specific DB.
Encapsulates the database specific functions this way.

=back

=cut

sub get_smashdb_sqlite_file {
	my $this    = shift;
	return sprintf("%s/%s.sqlite", $this->data_dir, $this->get_conf_value("SmashDB", "database_name"));
}

sub get_collection_sqlite_file {
	my $this       = shift;
	my $collection = shift;
	return sprintf("%s/metagenome_collections/%s/%s.sqlite", $this->data_dir, $collection, $collection);
}

sub get_refproteindb_sqlite_file {
	my $this    = shift;
	return sprintf("%s/%s.sqlite", $this->get_conf_value("RefProteinDB", "data_dir"), $this->get_conf_value("RefProteinDB", "database_name"));
}

sub get_reforganismdb_sqlite_file {
	my $this    = shift;
	return sprintf("%s/%s.sqlite", $this->get_conf_value("RefOrganismDB", "data_dir"), $this->get_conf_value("RefOrganismDB", "database_name"));
}

####
# These will be created when make_db_handle is called
# and closed when finish is called
####

sub _db_engine {shift->{DB_ENGINE}}
sub _db_handle {shift->{DB_HANDLE}}

sub make_db_handle {
	my $this = shift;
	my $collection = $this->collection;
	my %config     = %{$this->config};
	my $engine     = $config{SmashDB}{database_engine};
	my $host       = $config{SmashDB}{host};
	my $port       = $config{SmashDB}{port};
	my $user       = $config{SmashDB}{user} || "";
	my $passwd     = $config{SmashDB}{pass} || "";
	my $dbh;

	if (!$collection) {
		confess "COLLECTION not defined";
	}

	if ($engine eq "mysql") {
		$dbh = DBI->connect("DBI:mysql:database=$collection;host=$host;port=$port", $user, $passwd, {AutoCommit => 0, RaiseError => 1}) || confess "Couldn't connect to database: " .DBI->errstr;
	} elsif ($engine eq "sqlite3") {
		my $db_file = $this->get_collection_sqlite_file($collection);
		$dbh    = DBI->connect("DBI:SQLite:dbname=$db_file", $user, $passwd, {AutoCommit => 0, RaiseError => 1}) || confess "Couldn't connect to database: " .DBI->errstr;
		$dbh->do("PRAGMA foreign_keys = ON");
		$dbh->func(600000, 'busy_timeout');
	} else {
		confess "Database engine $engine is not supported by Smash";
	}

	$this->{DB_HANDLE} = $dbh;
	$this->{DB_ENGINE} = $engine;
	return $dbh;
}

sub get_db_handle {
	my $this = shift;
	return ($this->_db_handle || $this->make_db_handle);
}

sub close_db_handle {
	my $this = shift;
	my $dbh  = $this->_db_handle;
	if ($this->_db_engine eq "sqlite3") {
		$dbh->commit();
		undef $dbh;
	} else {
		$dbh->disconnect();
	}
	delete $this->{DB_HANDLE};
}

sub _smashdb_engine {shift->{SMASHDB_ENGINE}}
sub _smashdb_handle {shift->{SMASHDB_HANDLE}}

sub make_smashdb_handle {
	my $this   = shift;

	my %config     = %{$this->config};
	my $engine     = $config{SmashDB}{database_engine};
	my $database   = $config{SmashDB}{database_name};
	my $host       = $config{SmashDB}{host};
	my $port       = $config{SmashDB}{port};
	my $user       = $config{SmashDB}{user} || "";
	my $passwd     = $config{SmashDB}{pass} || "";
	my $dbh;

	if ($engine eq "mysql") {
		$dbh    = DBI->connect("DBI:mysql:database=$database;host=$host;port=$port", $user, $passwd, {AutoCommit => 0, RaiseError => 1}) || confess "Couldn't connect to database: " .DBI->errstr;
	} elsif ($engine eq "sqlite3") {
		my $db_file = $this->get_smashdb_sqlite_file();
		$dbh    = DBI->connect("DBI:SQLite:dbname=$db_file", $user, $passwd, {AutoCommit => 0, RaiseError => 1}) || confess "Couldn't connect to database: " .DBI->errstr;
		$dbh->do("PRAGMA foreign_keys = ON");
		$dbh->func(600000, 'busy_timeout');
	} else {
		confess "Database engine $engine is not supported by Smash";
	}
	$this->{SMASHDB_HANDLE} = $dbh;
	$this->{SMASHDB_ENGINE} = $engine;
	return $dbh;
}

sub get_smashdb_handle {
	my $this = shift;
	return ($this->_smashdb_handle || $this->make_smashdb_handle);
}

sub close_smashdb_handle {
	my $this = shift;
	my $dbh  = $this->_smashdb_handle;
	if ($this->_smashdb_engine eq "sqlite3") {
		$dbh->commit();
		undef $dbh;
	} else {
		$dbh->disconnect();
	}
	delete $this->{SMASHDB_HANDLE};
}

sub _refgenomedb_engine {shift->{REFGENOMEDB_ENGINE}}
sub _refgenomedb_handle {shift->{REFGENOMEDB_HANDLE}}

sub init_refgenomedb_database {
	my $this   = shift;
	my $dbh    = shift;
	my %config = %{$this->config};
	my $engine = $config{RefOrganismDB}{database_engine};
	my @tables = (
		"CREATE TABLE IF NOT EXISTS \
			taxonomy(\
				taxonomy_id BIGINT PRIMARY KEY, \
				organism VARCHAR(255) NOT NULL \
			)",
		"CREATE TABLE IF NOT EXISTS \
			project(\
				project_id BIGINT PRIMARY KEY, \
				taxonomy_id BIGINT NOT NULL REFERENCES taxonomy(taxonomy_id), \
				source VARCHAR(255) NOT NULL, \
				accession VARCHAR(255) DEFAULT NULL \
			)",
		"CREATE INDEX IF NOT EXISTS \
			p_tid\
				ON project(taxonomy_id)",
		"CREATE TABLE IF NOT EXISTS \
			sequence(\
				sequence_id VARCHAR(255) PRIMARY KEY, \
				project_id BIGINT NOT NULL REFERENCES project(project_id), \
				seq_type VARCHAR(10) NOT NULL, \
				mol_type VARCHAR(10) NOT NULL, \
				length INTEGER NOT NULL, \
				definition VARCHAR(255) NOT NULL
			)",
		"CREATE INDEX IF NOT EXISTS \
			s_pid\
				ON sequence(project_id)",
		"CREATE TABLE IF NOT EXISTS \
			feature(\
				sequence_id VARCHAR(255) REFERENCES sequence(sequence_id), \
				start INTEGER NOT NULL, \
				end   INTEGER NOT NULL, \
				type  INTEGER NOT NULL, \
				info  VARCHAR(512), \
				PRIMARY KEY(sequence_id, start, end, type) \
			)",
		"CREATE TABLE IF NOT EXISTS \
			gene(\
				gene_id INTEGER PRIMARY KEY AUTO_INCREMENT, \
				external_id VARCHAR(255), \
				sequence_id VARCHAR(255) REFERENCES sequence(sequence_id), \
				type VARCHAR(32), \
				gene_info VARCHAR(512), \
				length INTEGER NOT NULL, \
				start INTEGER NOT NULL, \
				end INTEGER NOT NULL, \
				strand CHAR(1) NOT NULL, \
				start_codon INTEGER NOT NULL, \
				stop_codon INTEGER NOT NULL, \
				gc FLOAT(3,1) NOT NULL, \
				UNIQUE (external_id) \
			)",
		"CREATE INDEX IF NOT EXISTS \
			g_sid\
				ON gene(sequence_id)",
		"CREATE TABLE IF NOT EXISTS \
			mature_rna(\
				gene_id INTEGER REFERENCES gene(gene_id), \
				start INTEGER NOT NULL, \
				end INTEGER NOT NULL, \
				PRIMARY KEY(gene_id, start) \
			)",
		"CREATE TABLE IF NOT EXISTS \
			gene2og(\
				gene VARCHAR(255) REFERENCES gene(external_id), \
				og VARCHAR(32), \
				PRIMARY KEY (gene, og) \
			)",
		"CREATE INDEX IF NOT EXISTS \
			g2o_og\
				ON gene2og(og)",
		"CREATE TABLE IF NOT EXISTS \
			marker_genes(\
				gene VARCHAR(255) REFERENCES gene(external_id), \
				og VARCHAR(32), \
				taxonomy_id BIGINT REFERENCES taxonomy(taxonomy_id), \
				project_id BIGINT NOT NULL REFERENCES project(project_id), \
				sequence_id VARCHAR(255) REFERENCES sequence(sequence_id), \
				start INTEGER NOT NULL, \
				end INTEGER NOT NULL, \
				strand CHAR(1) NOT NULL, \
				PRIMARY KEY (gene, og) \
			)",
		"CREATE INDEX IF NOT EXISTS \
			mg_og\
				ON marker_genes(og)",
		"CREATE INDEX IF NOT EXISTS \
			mg_tax\
				ON marker_genes(taxonomy_id)",
		"CREATE INDEX IF NOT EXISTS \
			mg_proj\
				ON marker_genes(project_id)",
	);
	$this->execute_statements($dbh, $engine, @tables);
	$dbh->commit();
}

sub make_refgenomedb_handle {
	my $this   = shift;

	my %config     = %{$this->config};
	my $data_dir   = $config{RefOrganismDB}{data_dir};
	my $engine     = $config{RefOrganismDB}{database_engine};
	my $database   = $config{RefOrganismDB}{database_name};
	my $host       = $config{RefOrganismDB}{host};
	my $port       = $config{RefOrganismDB}{port};
	my $user       = $config{RefOrganismDB}{user} || "";
	my $passwd     = $config{RefOrganismDB}{pass} || "";
	my $dbh;

	if ($engine eq "mysql") {
		$dbh    = DBI->connect("DBI:mysql:database=$database;host=$host;port=$port", $user, $passwd, {AutoCommit => 0, RaiseError => 1}) || confess "Couldn't connect to database: " .DBI->errstr;
	} elsif ($engine eq "sqlite3") {
		my $db_file = $this->get_reforganismdb_sqlite_file();
		if (! -f $db_file) {
			$dbh = DBI->connect("DBI:SQLite:dbname=$db_file", $user, $passwd, {AutoCommit => 0, RaiseError => 1}) || confess "Couldn't connect to database: " .DBI->errstr;
			$this->init_refgenomedb_database($dbh);
		} else {
			$dbh = DBI->connect("DBI:SQLite:dbname=$db_file", $user, $passwd, {AutoCommit => 0, RaiseError => 1}) || confess "Couldn't connect to database: " .DBI->errstr;
		}
		$dbh->do("PRAGMA foreign_keys = ON");
		$dbh->func(600000, 'busy_timeout');
	} else {
		confess "Database engine $engine is not supported by Smash";
	}
	$this->{REFGENOMEDB_HANDLE} = $dbh;
	$this->{REFGENOMEDB_ENGINE} = $engine;
	return $dbh;
}

sub get_refgenomedb_handle {
	my $this = shift;
	return ($this->_refgenomedb_handle || $this->make_refgenomedb_handle);
}

sub close_refgenomedb_handle {
	my $this = shift;
	my $dbh  = $this->_refgenomedb_handle;
	if ($this->_refgenomedb_engine eq "sqlite3") {
		$dbh->commit();
		undef $dbh;
	} else {
		$dbh->disconnect();
	}
	delete $this->{REFGENOMEDB_HANDLE};
}

sub _refproteindb_engine {shift->{REFPROTEINDB_ENGINE}}
sub _refproteindb_handle {shift->{REFPROTEINDB_HANDLE}}

sub init_refproteindb_database {
	my $this   = shift;
	my $dbh    = shift;
	my %config = %{$this->config};
	my $engine = $config{RefProteinDB}{database_engine};
	my @tables = (
		"CREATE TABLE IF NOT EXISTS \
			kegg_details(\
				protein VARCHAR(255) NOT NULL, \
				length INTEGER NOT NULL, \
				PRIMARY KEY (protein) \
			)",
		"CREATE TABLE IF NOT EXISTS \
			protein2ko(\
				protein VARCHAR(255) NOT NULL, \
				ko VARCHAR(31) NOT NULL, \
				PRIMARY KEY (protein, ko) \
			)",
		"CREATE INDEX IF NOT EXISTS \
			p2k_ko\
				ON protein2ko(ko)",
		"CREATE TABLE IF NOT EXISTS \
			ko2name(\
				ko VARCHAR(31) NOT NULL, \
				name VARCHAR(511) NOT NULL, \
				PRIMARY KEY (ko) \
			)",
		"CREATE TABLE IF NOT EXISTS \
			ko2go(\
				ko VARCHAR(31) NOT NULL, \
				go VARCHAR(127) NOT NULL, \
				PRIMARY KEY (ko, go) \
			)",
		"CREATE INDEX IF NOT EXISTS \
			k2g_go\
				ON ko2go(go)",
		"CREATE TABLE IF NOT EXISTS \
			ko2module(\
				ko VARCHAR(31) NOT NULL, \
				module VARCHAR(31) NOT NULL, \
				PRIMARY KEY (ko, module) \
			)",
		"CREATE INDEX IF NOT EXISTS \
			k2m_mo\
				ON ko2module(module)",
		"CREATE TABLE IF NOT EXISTS \
			ko2pathway(\
				ko VARCHAR(31) NOT NULL, \
				pathway VARCHAR(31) NOT NULL, \
				PRIMARY KEY (ko, pathway) \
			)",
		"CREATE INDEX IF NOT EXISTS \
			k2p_pw\
				ON ko2pathway(pathway)",
		"CREATE TABLE IF NOT EXISTS \
			module2name(\
				module VARCHAR(31) NOT NULL, \
				name VARCHAR(511) NOT NULL, \
				PRIMARY KEY (module) \
			)",
		"CREATE TABLE IF NOT EXISTS \
			module2def(\
				module VARCHAR(31) NOT NULL, \
				def VARCHAR(2047) NOT NULL, \
				PRIMARY KEY (module) \
			)",
		"CREATE TABLE IF NOT EXISTS \
			module2pathway(\
				module VARCHAR(31) NOT NULL, \
				pathway VARCHAR(31) NOT NULL, \
				PRIMARY KEY (module, pathway) \
			)",
		"CREATE INDEX IF NOT EXISTS \
			m2p_pw\
				ON module2pathway(pathway)",
		"CREATE TABLE IF NOT EXISTS \
			pathway2name(\
				pathway VARCHAR(31) NOT NULL, \
				name VARCHAR(511) NOT NULL, \
				PRIMARY KEY (pathway) \
			)",
		"CREATE TABLE IF NOT EXISTS \
			pathway2go(\
				pathway VARCHAR(31) NOT NULL, \
				go VARCHAR(127) NOT NULL, \
				PRIMARY KEY (pathway, go) \
			)",
		"CREATE INDEX IF NOT EXISTS \
			p2g_go\
				ON pathway2go(go)",
	);
	$this->execute_statements($dbh, $engine, @tables);
	$dbh->commit();
}

sub make_refproteindb_handle {
	my $this   = shift;

	my $config     = $this->config->{RefProteinDB};
	my $data_dir   = $config->{data_dir};
	my $engine     = $config->{database_engine};
	my $database   = $config->{database_name};
	my $host       = $config->{host};
	my $port       = $config->{port};
	my $user       = $config->{user} || "";
	my $passwd     = $config->{pass} || "";
	my $dbh;

	if ($engine eq "mysql") {
		$dbh    = DBI->connect("DBI:mysql:database=$database;host=$host;port=$port", $user, $passwd, {AutoCommit => 0, RaiseError => 1}) || confess "Couldn't connect to database: " .DBI->errstr;
	} elsif ($engine eq "sqlite3") {
		my $db_file = $this->get_refproteindb_sqlite_file();
		if (! -f $db_file) {
			$dbh = DBI->connect("DBI:SQLite:dbname=$db_file", $user, $passwd, {AutoCommit => 0, RaiseError => 1}) || confess "Couldn't connect to database: " .DBI->errstr;
			$this->init_refproteindb_database($dbh);
		} else {
			$dbh = DBI->connect("DBI:SQLite:dbname=$db_file", $user, $passwd, {AutoCommit => 0, RaiseError => 1}) || confess "Couldn't connect to database: " .DBI->errstr;
		}
		$dbh->do("PRAGMA foreign_keys = ON");
		$dbh->func(600000, 'busy_timeout');
	} else {
		confess "Database engine $engine is not supported by Smash";
	}
	$this->{REFPROTEINDB_HANDLE} = $dbh;
	$this->{REFPROTEINDB_ENGINE} = $engine;
	return $dbh;
}

sub get_refproteindb_handle {
	my $this = shift;
	return ($this->_refproteindb_handle || $this->make_refproteindb_handle);
}

sub close_refproteindb_handle {
	my $this = shift;
	my $dbh  = $this->_refproteindb_handle;
	if ($this->_refproteindb_engine eq "sqlite3") {
		$dbh->commit();
		undef $dbh;
	} else {
		$dbh->disconnect();
	}
	delete $this->{REFPROTEINDB_HANDLE};
}

sub last_db_insert_id {
	my $this    = shift;
	my $dbh     = $this->get_db_handle;
	my $engine  = $this->_db_engine;
	if ($engine eq "mysql") {
		return $dbh->{mysql_insertid};
	} elsif ($engine eq "sqlite3") {
		return $dbh->func("last_insert_rowid");
	} else {
		confess "last_db_insert_id() not implemented for $engine";
	}
}

sub last_refgenomedb_insert_id {
	my $this    = shift;
	my $dbh     = $this->get_refgenomedb_handle;
	my $engine  = $this->_refgenomedb_engine;
	if ($engine eq "mysql") {
		return $dbh->{mysql_insertid};
	} elsif ($engine eq "sqlite3") {
		return $dbh->func("last_insert_rowid");
	} else {
		confess "last_refgenomedb_insert_id() not implemented for $engine";
	}
}

sub last_refproteindb_insert_id {
	my $this    = shift;
	my $dbh     = $this->get_refproteindb_handle;
	my $engine  = $this->_refproteindb_engine;
	if ($engine eq "mysql") {
		return $dbh->{mysql_insertid};
	} elsif ($engine eq "sqlite3") {
		return $dbh->func("last_insert_rowid");
	} else {
		confess "last_refproteindb_insert_id() not implemented for $engine";
	}
}

sub last_smashdb_insert_id {
	my $this    = shift;
	my $dbh     = $this->get_smashdb_handle;
	my $engine  = $this->_smashdb_engine;
	if ($engine eq "mysql") {
		return $dbh->{mysql_insertid};
	} elsif ($engine eq "sqlite3") {
		return $dbh->func("last_insert_rowid");
	} else {
		confess "last_smashdb_insert_id() not implemented for $engine";
	}
}

=head2 Methods accessing database

=over 4

=item B<get_metagenomes_for_collection>

Returns a list of metagenomes present in the given collection.

=item B<get_metagenome_label>

Returns the external label for an internal metagenome id.

=item B<get_metagenome_description>

Returns the description (specified when creating the metagenome) for an internal metagenome id.

=item B<get_refsequence_details($seq_id)>

returns (taxonomy_id, definition, length, display_tax_id) as a hash given a sequence
identifier.

=back

=cut

sub get_metagenomes_for_collection {
	my $this        = shift;
	my $collection  = shift;
	my $dbh         = $this->get_smashdb_handle();
	my $sth         = $dbh->prepare('SELECT metagenome_id FROM metagenome WHERE collection_id=? ORDER BY title');
	my @metagenomes = ();
	$sth->execute($collection);
	while (my ($id) = $sth->fetchrow_array()) {
		push(@metagenomes, $id);
	}
	$sth->finish();
	$this->close_smashdb_handle();
	return @metagenomes;
}

sub get_metagenome_label {
        my $this        = shift;
        my $metagenome  = shift;
        my $dbh         = $this->get_smashdb_handle();
	my $label;
	{
		my $sth         = $dbh->prepare("SELECT title FROM metagenome WHERE metagenome_id=?");
		$sth->execute($metagenome);
		($label)     = $sth->fetchrow_array();
		$sth->finish();
	}
	$this->close_smashdb_handle();
        return $label;
}

sub get_metagenome_description {
        my $this        = shift;
        my $metagenome  = shift;
        my $dbh         = $this->get_smashdb_handle();
        my $sth         = $dbh->prepare("SELECT description FROM metagenome WHERE metagenome_id=?");
        $sth->execute($metagenome);
        my ($label)     = $sth->fetchrow_array();
        $sth->finish();
	$this->close_smashdb_handle();
        return $label;
}

sub get_refsequence_details {
        my $this        = shift;
        my $seq_id      = shift;
        my $dbh         = $this->get_refgenomedb_handle();
        my $sth         = $dbh->prepare("SELECT taxonomy_id, definition, length, display_tax_id FROM sequence INNER JOIN taxonomy USING (taxonomy_id) WHERE sequence_id=?");
        $sth->execute($seq_id);
        my ($tax_id, $def, $length, $display_tax_id) = $sth->fetchrow_array();
        $sth->finish();
	$this->close_refgenomedb_handle();
        return {taxonomy_id => $tax_id, definition => $def, length => $length, display_tax_id => $display_tax_id};
}

=head2 TO DO

clean up here

=over 4

=item B<get_metagenome_files($metagenome, $extension)>

Returns an array of all files in the read directory with the given extension. This function is called by
C<fasta_files($metagenome)>, C<qual_files($metagenome)> and C<xml_files($metagenome)> as 
C<get_metagenome_files($metagenome, "fasta")>,
C<get_metagenome_files($metagenome, "qual")> and
C<get_metagenome_files($metagenome, "xml")> respectively. Can be extended by the user to retrieve any kind of
files.

=item B<fasta_files($metagenome)>

Returns a list of read fasta files corresponding to this metagenome.

=item B<qual_files($metagenome)>

Returns a list of read quality files corresponding to this metagenome.

=item B<xml_files($metagenome)>

Returns a list of read xml files corresponding to this metagenome.

=item B<get_read_lengths($metagenome)>

Returns a hash containing read length information grouped into templates. For example,
if paired end reads for template C<GFTRSDA> are C<GFTRSDA.b> from the forward primer and
C<GFTRSDA.z> from the reverse primer, the values in the hash will be:

	{"GFTRSDA" => {"GFTRSDA.b" => 656, "GFTRSDA.z" => 643}}

=back

=cut

sub get_metagenome_files {
	my $this        = shift;
	my $metagenome  = shift;
	my $extension   = shift;
	my $read_dir    = $this->read_dir($metagenome);
	my @files       = ();
	foreach my $file (<$read_dir/*.$extension>) {
		push(@files, $file);
	}
	return @files;
}

sub fasta_files {
	shift->get_metagenome_files(shift, "fasta");
}

sub qual_files {
	shift->get_metagenome_files(shift, "qual");
}

sub xml_files {
	shift->get_metagenome_files(shift, "xml");
}

sub get_read_lengths {
	my $this       = shift;
	my $metagenome = shift;
	my $dbh        = $this->get_db_handle();

	my $ReadLength = {};

	# Get read lengths

	{
		my $sth = $dbh->prepare_cached("SELECT read_id, template_id, length FROM readinfo r INNER JOIN library l USING (library_id) INNER JOIN sample s USING (sample_id) WHERE metagenome_id=?");
		$sth->execute($metagenome);
		my ($name, $template, $read_length);
		$sth->bind_columns( \( $name, $template, $read_length ));
		while ($sth->fetchrow_arrayref()) {
			$ReadLength->{$template}->{$name} = $read_length;
		}
		$sth->finish();
	}

	return $ReadLength;
}

=head2 Generic functions to map between internal and external ids

These functions get the relevant internal(external) ids given external(internal) ids. Works for 
a few elements where the table structure is simple. The list of element types where these functions
work is listed below.

=over 4

=item B<get_id_by_name>

Returns the internal id for a given element from its corresponding table. Supported element
types are: C<assembly>, C<gene_prediction>, C<contig>, C<scaffold>, C<gene>, C<library>. For
example,

	$smash->get_id_by_name("contig", "MC1.MG1.AS1.C2");

will return the integer primary key in the C<contig> table.

=cut

sub get_id_by_name_generic {
	my $this  = shift;
	my $dbh   = shift;
	my $table = shift;
	my $value = shift;

	my $id;
	{
		my $sth = $dbh->prepare_cached("SELECT ${table}_id FROM $table WHERE external_id=?");
		$sth->execute($value);
		($id) = $sth->fetchrow_array();
		while ($sth->fetchrow_array()){}; # perl DBI shush
		$sth->finish();
	}
	return $id;
}

sub get_id_by_name_smashdb {
	my $this  = shift;
	my $dbh   = shift;
	my $table = shift;
	my $value = shift;

	my $id;
	{
		my $sth = $dbh->prepare("SELECT ${table}_id FROM $table WHERE external_id=?");
		$sth->execute($value);
		($id) = $sth->fetchrow_array();
		#$sth->fetchrow_array(); # perl DBI shush
		$sth = undef;
	}
	return $id;
}

sub get_id_by_name {
	my $this  = shift;
	my $table = shift;
	my $value = shift;
	my $id;
	my @smashdb_tables = qw(assembly gene_prediction);
	my @mc_tables      = qw(contig scaffold gene library);
	if ((grep {$_ eq $table} @smashdb_tables) > 0) {
		my $dbh = $this->get_smashdb_handle();
		$id = $this->get_id_by_name_smashdb($dbh, $table, $value);
		$this->close_smashdb_handle();
		return $id;
	} elsif ((grep {$_ eq $table} @mc_tables) > 0) {
		my $dbh = $this->get_db_handle();
		$id = $this->get_id_by_name_generic($dbh, $table, $value);
		$this->close_db_handle();
		return $id;
	} else {
		confess "Do not know how to get id for $value!";
	}
}

=item B<get_name_by_id>

Returns the external id for a given element from its corresponding table. Supported element
types are: C<assembly>, C<gene_prediction>, C<contig>, C<scaffold>, C<gene>, C<library>. For
example,

	$smash->get_name_by_id("contig", 1532);

where C<1532> is the integer primary key in the C<contig> table, returns the external string
id of that contig, such as C<"MC1.MG1.AS1.C2">.

=back

=cut

sub get_name_by_id_generic {
	my $this  = shift;
	my $dbh   = shift;
	my $table = shift;
	my $value = shift;

	my $name;
	{
		my $sth   = $dbh->prepare_cached("SELECT external_id FROM $table WHERE ${table}_id=?");
		$sth->execute($value);
		($name) = $sth->fetchrow_array();
		while ($sth->fetchrow_array()){}; # perl DBI shush
		$sth->finish();
	}
	return $name;
}

sub get_name_by_id {
	my $this  = shift;
	my $table = shift;
	my $value = shift;
	my $name;
	my @smashdb_tables = qw(assembly gene_prediction);
	my @mc_tables      = qw(contig scaffold gene library);
	if ((grep {$_ eq $table} @smashdb_tables) > 0) {
		my $dbh = $this->get_smashdb_handle();
		$name = $this->get_name_by_id_generic($dbh, $table, $value);
		$this->close_smashdb_handle();
		return $name;
	} elsif ((grep {$_ eq $table} @mc_tables) > 0) {
		my $dbh = $this->get_db_handle();
		$name = $this->get_name_by_id_generic($dbh, $table, $value);
		$this->close_db_handle();
		return $name;
	} else {
		confess "Do not know how to get id for $value!";
	}
}

=head2 Specific functions to get internal ids for external information

These functions get the relevant internal ids given external information. They create an entry if it 
does not exist, and return the newly created internal id.

=over 4

=item B<get_program_id>

Returns an internal id for a software given its name, version and parameters.
Creates an entry in the database if necessary. For example,

	$smash->get_program_id("GeneMark", "0.96", "heu_11_gc");

checks if there is a record of C<GeneMark> version C<0.96> using parameter key C<heu_11_gc>, and makes
one if necessary.

=item B<get_sample_id>

Returns an internal id for a sample given the metagenome id and metadata of this
specific sample. Creates an entry in the database if necessay. For example,

	$smash->get_sample_id("MC2.MG1", "Depth:50m, Temp:14F");

checks if there is a record of sample C<"MC2.MG1"> with associated metadata C<"Depth:50m, Temp:14F">, and
makes one if necessary.

=item B<get_library_id>

Returns an internal id for a library given the sample id, type of reads, and
if applicable insert length and standard deviation. Creates an entry in the 
database if necessary. For example,

	$smash->get_library_id("MC2.BAAU", 23, "sanger", 10000, 3000);

checks if there is a record of library C<"MC2.BAAU"> associated with sample C<23>, containing C<sanger> reads
of insert size C<10000> and insert standard deviation C<3000>.

=back

=cut

sub get_program_details {
	my $this       = shift;
	my $id         = shift;
	my $dbh        = $this->get_smashdb_handle();

	my (undef, undef, $asm, $gp) = $this->parse_concat_id($id);
	my @vals;
	if ($gp) {
		my $sth = $dbh->prepare('SELECT p.name, p.version, p.parameters FROM gene_prediction gp INNER JOIN program p ON gp.predictor = p.program_id AND gp.external_id=?');
		$sth->execute($gp);
		if (@vals = $sth->fetchrow_array()) {
		} else {
			die "Cannot find entry $id in database!";
		}
	} elsif ($asm) {
		my $sth = $dbh->prepare('SELECT p.name, p.version, p.parameters FROM assembly a INNER JOIN program p ON a.assembler = p.program_id AND a.external_id=?');
		$sth->execute($asm);
		if (@vals = $sth->fetchrow_array()) {
		} else {
			die "Cannot find entry $id in database!";
		}
	} else {
		die "Need assembly or gene_prediction id to retrieve software details!";
	}
	$this->close_smashdb_handle();
	return @vals;
}

sub get_program_id {
	my $this       = shift;
	my $name       = shift;
	my $version    = shift;
	my $parameters = shift;
	my $dbh        = $this->get_smashdb_handle();
	my $program;
	{
		my $sth        = $dbh->prepare('SELECT program_id FROM program WHERE name=? and version=? and parameters=?');

		      $name =~ s/^\s*//;       $name =~ s/\s*$//;
		   $version =~ s/^\s*//;    $version =~ s/\s*$//;
		$parameters =~ s/^\s*//; $parameters =~ s/\s*$//;

		$sth->execute($name, $version, $parameters);
		if (my ($id) = $sth->fetchrow_array()) {
			$program = $id;
		} else {
			my $sth2 = $dbh->prepare('INSERT INTO program(name, version, parameters) VALUES(?, ?, ?)');
			$sth2->execute($name, $version, $parameters);
			$program = $this->last_smashdb_insert_id($dbh);
			$sth2->finish();
			$dbh->commit();
		}
		$sth->finish();
	}
	$this->close_smashdb_handle();
	return $program;
}

sub get_sample_id {
	my $this       = shift;
	my $metagenome = shift;
	my $metadata   = shift;

	my $sample;
	my $dbh        = $this->get_db_handle();
	my $sth        = $dbh->prepare('SELECT sample_id FROM sample WHERE metagenome_id=? and metadata=?');

	$sth->execute($metagenome, $metadata);
	if (my ($id) = $sth->fetchrow_array()) {
		$sample = $id;
	} else {
		my $sth2 = $dbh->prepare('INSERT INTO sample(metagenome_id, metadata) VALUES(?, ?)');
		$sth2->execute($metagenome, $metadata);
		$sample = $this->last_db_insert_id($dbh);
		$sth2->finish();
	}
	$sth->finish();
	$this->close_db_handle();
	return $sample;
}

sub get_library_id {
	my $this       = shift;
	my ($name, $sample_id, $type, $insert_length, $insert_stdev) = @_;

	my $library;
	my $dbh        = $this->get_db_handle();
	my $sth;
	if (defined($insert_length) && defined($insert_stdev)) {
		$sth        = $dbh->prepare('SELECT library_id FROM library WHERE external_id=? AND sample_id=? AND type=? AND insert_length=? AND insert_stdev=?');
		$sth->execute($name, $sample_id, $type, $insert_length, $insert_stdev);
	} elsif (defined($insert_length)) {
		$sth        = $dbh->prepare('SELECT library_id FROM library WHERE external_id=? AND sample_id=? AND type=? AND insert_length=? AND insert_stdev IS NULL');
		$sth->execute($name, $sample_id, $type, $insert_length);
	} else {
		$sth        = $dbh->prepare('SELECT library_id FROM library WHERE external_id=? AND sample_id=? AND type=? AND insert_length IS NULL AND insert_stdev IS NULL');
		$sth->execute($name, $sample_id, $type);
	}

	if (my ($id) = $sth->fetchrow_array()) {
		$library = $id;
	} else {
		my $sth2 = $dbh->prepare('INSERT INTO library(external_id, sample_id, type, insert_length, insert_stdev) VALUES(?, ?, ?, ?, ?);');
		$sth2->execute($name, $sample_id, $type, $insert_length, $insert_stdev);
		$library = $this->last_db_insert_id($dbh);
		$sth2->finish();
	}
	$sth->finish();
	$this->close_db_handle();
	return $library;
}

=head2 Creating new entries and removing them

=over 4

=item B<make_new_assembly>

Makes a new assembly. Any instance of subclasses of L<Smash::Analyses::Assembler> should call this function to get a new assembly
id from Smash. It creates an entry in the C<assembly> table with the C<program_id> of this instance.

=item B<remove_assembly_by_name>

Removes the assembly from the C<assembly> table given the external id.

=item B<safe_make_new_entry>

Safely adds a new entry into the SmashDB database. First selects the entries
in the table using C<$select_st> using

	$select_sth->execute(@$st1_args);

finds the maximum value after removing C<$prefix>, increments it by one,
makes the new entry using 

	$entry = "$prefix$max";

and then inserts into the table using

	$insert_sth->execute($entry, @$st2_args);

For example, C<make_new_assembly()> calls this function as follows:

	sub make_new_assembly {
		my $this       = shift;
		my $metagenome = shift;
		my $assembler  = shift;
		my $as         = $this->get_smash_conf_value("assembly_prefix");
		my $prefix     = "$metagenome.$as";

		my $st1        = 'SELECT external_id FROM assembly WHERE metagenome_id=?';
		my $st2        = 'INSERT INTO assembly(external_id, metagenome_id, assembler) VALUES(?, ?, ?)';
		my $assembly   = $this->safe_make_new_entry($prefix, $st1, $st2, [$metagenome], [$metagenome, $assembler]);

		return $assembly;
	}

=cut

sub safe_make_new_entry {
	my $this = shift;
	my ($prefix, $select_st, $insert_st, $st1_args, $st2_args) = @_;

	my $entry;
	my $dbh = $this->get_smashdb_handle();
	{
		my $attempt    = 0;
		TRY:while (1) { # Keep trying until you succeed
			my $max = 0;
			{
				my $select_sth = $dbh->prepare($select_st);
				   $select_sth->execute(@$st1_args);
				while (my ($id) = $select_sth->fetchrow_array()) {
					$id =~ s/^${prefix}//;
					if ($id > $max) {
						$max = $id;
					}
				}
			}
			$max++;
			$entry = "$prefix$max";

			$attempt++;
			my $insert_sth = $dbh->prepare($insert_st);
			eval {
				$insert_sth->execute($entry, @$st2_args);
			};

			if ($@) {
				my $sleep_time = 30+int(rand(30));
				warn "NOTE: safe_make_new_entry attempt $attempt failed. Will retry in $sleep_time seconds!\n";
				$this->close_smashdb_handle();
				sleep($sleep_time);
				$dbh = $this->get_smashdb_handle();
			} else {
				warn "NOTE: safe_make_new_entry attempt $attempt succeeded!\n" if ($attempt > 1);
				$dbh->commit();
				last TRY;
			}
		}
	}
	$this->close_smashdb_handle();
	return $entry;
}

sub execute_statements {
	my $this       = shift;
	my $dbh        = shift;
	my $engine     = shift;
	my @statements = @_;
	STATEMENT:foreach my $statement (@statements) {
		if ($engine eq "sqlite3") {
			$statement =~ s/AUTO_INCREMENT/AUTOINCREMENT/g;
		}
		print "Executing: '$statement'\n";
		my $sth    = $dbh->prepare($statement);
		my $status = $sth->execute();
		if (!$status) {
			warn "WARNING: Execution failed!\n";
		}
	}
}

sub init_smash_environment {
	my $this = shift;
	mkpath($this->get_smash_conf_value("data_dir"), {mode => $this->file_perm});
	$this->init_global_database();
}

sub init_global_database {
	my $this       = shift;

	# for sqlite, just making a connection to the MC database creates it if it doesnt exist.
	# So we can directly call get_smashdb_handle().
	# For MySQL though, we need to make this database!
	# Since that turned out to be complicated, making the SmashDB database is part of the 
	# installation process now, since there is a manual step of adding a user anyway.
	# Now, we can just call get_smashdb_handle() with MySQL as well.

	my $dbh     = $this->get_smashdb_handle();
	my %config  = %{$this->config};
	my $engine  = $this->get_conf_value("SmashDB", "database_engine");
	my @queries = (
			"CREATE TABLE IF NOT EXISTS \
				metagenome_collection (\
					collection_id VARCHAR(5) PRIMARY KEY, \
					title VARCHAR(40) NOT NULL, \
					description TEXT \
				)",
			"CREATE TABLE IF NOT EXISTS \
				metagenome (\
					metagenome_id VARCHAR(10) PRIMARY KEY, \
					collection_id VARCHAR(5) NOT NULL, \
					title VARCHAR(40) NOT NULL, \
					default_assembly_id INTEGER DEFAULT NULL, \
					description TEXT \
				)",
			"CREATE TABLE IF NOT EXISTS \
				program (\
					program_id INTEGER PRIMARY KEY AUTO_INCREMENT, \
					name VARCHAR(64) NOT NULL, \
					version VARCHAR(64) NOT NULL, \
					parameters VARCHAR(1024) NOT NULL\
				)",
			"CREATE TABLE IF NOT EXISTS \
				assembly (\
					assembly_id INTEGER PRIMARY KEY AUTO_INCREMENT, \
					external_id VARCHAR(32) NOT NULL UNIQUE, \
					metagenome_id VARCHAR(32) NOT NULL, \
					assembler INTEGER NOT NULL REFERENCES program(program_id)\
				)",
			"CREATE TABLE IF NOT EXISTS \
				gene_prediction (\
					gene_prediction_id INTEGER PRIMARY KEY AUTO_INCREMENT, \
					external_id VARCHAR(32) NOT NULL UNIQUE, \
					assembly_id INTEGER NOT NULL REFERENCES assembly(assembly_id), \
					predictor INTEGER NOT NULL REFERENCES program(program_id)\
				)"

	);
	$this->execute_statements($dbh, $engine, @queries);
	$this->close_smashdb_handle();
}

sub init_collection_database {
	my $this       = shift;
	my $collection = shift;

	# for sqlite, just making a connection to the MC database creates it if it doesnt exist.
	# So we can directly init() the object with COLLECTION attribute.
	# For MySQL though, we need to make this database!

	my $engine     = $this->get_conf_value("SmashDB", "database_engine");
	if ($engine eq "mysql") {
		my $dbh  = $this->get_smashdb_handle();
		my $user = $this->get_conf_value("SmashDB", "user");
		$dbh->func("createdb", $collection, 'admin');
		$this->close_smashdb_handle();
		$dbh = $this->get_smashdb_handle();
		$dbh->do("GRANT ALL PRIVILEGES ON $collection.* TO `$user`@`localhost` WITH GRANT OPTION;");
		$dbh->do("GRANT ALL PRIVILEGES ON $collection.* TO `$user`@`%` WITH GRANT OPTION;");
		$dbh->commit();
		$dbh->func("reload", 'admin');
		$this->close_smashdb_handle();
	}

	# Now reinit everything

	$this->init($collection);

	# Make the tables in the new database

	my $dbh = $this->get_db_handle;
	my @queries = (
		"CREATE TABLE IF NOT EXISTS \
				sample(\
					sample_id INTEGER PRIMARY KEY AUTO_INCREMENT, \
					metagenome_id VARCHAR(10) NOT NULL, \
					metadata VARCHAR(255)\
				)",
		"CREATE TABLE IF NOT EXISTS \
				library(\
					library_id INTEGER PRIMARY KEY AUTO_INCREMENT, \
					external_id VARCHAR(255) NOT NULL, \
					sample_id INTEGER NOT NULL REFERENCES sample(sample_id), \
					type VARCHAR(10) NOT NULL, \
					insert_length INTEGER, \
					insert_stdev INTEGER\
				)",
		"CREATE TABLE IF NOT EXISTS \
				readinfo(\
					read_id VARCHAR(255) PRIMARY KEY, \
					defline VARCHAR(255) NOT NULL, \
					template_id VARCHAR(255), \
					direction VARCHAR(1), \
					library_id INTEGER REFERENCES library(library_id) ON UPDATE CASCADE, \
					clip_left INTEGER, \
					clip_right INTEGER, \
					length INTEGER NOT NULL\
				)",
		"CREATE TABLE IF NOT EXISTS \
				scaffold(\
					scaffold_id INTEGER PRIMARY KEY AUTO_INCREMENT, \
					external_id VARCHAR(255) NOT NULL UNIQUE, \
					contig_count INTEGER NOT NULL, \
					length INTEGER NOT NULL, \
					assembly_id INTEGER NOT NULL\
				)",
		"CREATE TABLE IF NOT EXISTS \
				scaffold2contig(\
					scaffold_id INTEGER NOT NULL REFERENCES scaffold(scaffold_id) ON DELETE CASCADE, \
					contig_id INTEGER NOT NULL REFERENCES contig(contig_id), \
					start INTEGER NOT NULL, \
					end INTEGER NOT NULL, \
					strand CHAR(1) NOT NULL, \
					PRIMARY KEY(scaffold_id, contig_id)\
				)",
		"CREATE TABLE IF NOT EXISTS \
				contig(\
					contig_id INTEGER PRIMARY KEY AUTO_INCREMENT, \
					external_id VARCHAR(255) NOT NULL UNIQUE, \
					read_count INTEGER NOT NULL, \
					length INTEGER NOT NULL, \
					assembly_id INTEGER NOT NULL\
				)",
		"CREATE TABLE IF NOT EXISTS \
				contig2read(\
					contig_id INTEGER NOT NULL REFERENCES contig(contig_id) ON DELETE CASCADE, \
					read_id VARCHAR(255) NOT NULL REFERENCES readinfo(read_id), \
					start INTEGER NOT NULL, \
					end INTEGER NOT NULL, \
					strand CHAR(1) NOT NULL, \
					PRIMARY KEY(contig_id, read_id, start)\
				)",
		"CREATE TABLE IF NOT EXISTS \
				gene(\
					gene_id INTEGER PRIMARY KEY AUTO_INCREMENT, \
					external_id VARCHAR(255) NOT NULL, \
					gene_prediction_id INTEGER NOT NULL, \
					contig_id INTEGER NOT NULL REFERENCES contig(contig_id), \
					length INTEGER NOT NULL, \
					start INTEGER NOT NULL, \
					end INTEGER NOT NULL, \
					strand CHAR(1) NOT NULL, \
					start_codon INTEGER NOT NULL, \
					stop_codon INTEGER NOT NULL, \
					gc FLOAT(3,1) NOT NULL \
				)",
		"CREATE TABLE IF NOT EXISTS \
				read_refgenome_map (\
					read_id VARCHAR(255) NOT NULL REFERENCES readinfo(read_id), \
					ref_taxonomy_id BIGINT NOT NULL, \
					confidence FLOAT(6,4) DEFAULT NULL, \
					PRIMARY KEY(read_id, ref_taxonomy_id) \
				)",
		"CREATE INDEX \
				s_metagenome_id \
				ON sample(metagenome_id)",
		"CREATE INDEX \
				l_sample_id \
				ON library(sample_id)",
		"CREATE INDEX \
				r_library_id \
				ON readinfo(library_id)",
		"CREATE INDEX \
				c_external_id \
				ON contig(external_id)",
		"CREATE INDEX \
				c_assembly_id \
				ON contig(assembly_id)",
		"CREATE INDEX \
				c2r_read_id \
				ON contig2read(read_id)",
		"CREATE INDEX \
				s_external_id \
				ON scaffold(external_id)",
		"CREATE INDEX \
				s_assembly_id \
				ON scaffold(assembly_id)",
		"CREATE INDEX \
				s2c_contig_id \
				ON scaffold2contig(contig_id)",
		"CREATE INDEX \
				g_contig_id \
				ON gene(contig_id)",
		"CREATE INDEX \
				g_gp_id \
				ON gene(gene_prediction_id)",
		"CREATE INDEX \
				g_external_id \
				ON gene(external_id)",
		"CREATE TABLE IF NOT EXISTS \
				gene2og(\
					gene_name VARCHAR(255), \
					string_protein VARCHAR(255), \
					string_version FLOAT(5,2) NOT NULL, \
					og VARCHAR(63) NOT NULL, \
					placement_start INTEGER NOT NULL, \
					placement_end INTEGER NOT NULL, \
					bitscore FLOAT(16,2) NOT NULL, \
					PRIMARY KEY(gene_name, string_protein, string_version, og, placement_start, placement_end)\
				)"
	);
	$this->execute_statements($dbh, $engine, @queries);
	$dbh->commit();
	$this->close_db_handle();
}

sub make_new_collection {
	my $this       = shift;
	my $title      = shift;
	my $description= shift;
	my $mc         = $this->get_smash_conf_value("collection_prefix");
	my $prefix     = $mc;

	my $st1        = 'SELECT collection_id FROM metagenome_collection';
	my $st2        = 'INSERT INTO metagenome_collection(collection_id, title, description) VALUES(?, ?, ?)';
	my $collection = $this->safe_make_new_entry($prefix, $st1, $st2, [], [$title, $description]);

	# make the data_dir

	mkpath($this->collection_dir($collection), {mode => $this->file_perm});

	# init the database

	$this->init_collection_database($collection);
	return $collection;
}

sub remove_collection {
	my $this       = shift;
	my $collection = shift;
	my $dbh        = $this->get_smashdb_handle;
	my $engine     = $this->get_conf_value("SmashDB", "database_engine");

	# Remove entry from SmashDB table

	{
		my $sth  = $dbh->prepare('DELETE FROM metagenome_collection WHERE collection_id=?');
		$sth->execute($collection);
		$sth->finish();
	}

	# Remove the collection database

	if ($engine eq "mysql") {
		$dbh->func("dropdb", $collection, 'admin');
	} elsif ($engine eq "sqlite3") {
		unlink $this->get_collection_sqlite_file($collection);
	}

	$dbh->commit() unless $engine eq "mysql";
	$this->close_smashdb_handle();

	# Remove all the files

	rmtree($this->collection_dir($collection));
}

sub make_new_metagenome {
	my $this       = shift;
	my $collection = shift;
	my $title      = shift;
	my $description= shift;
	my $mg         = $this->get_smash_conf_value("metagenome_prefix");
	my $prefix     = "$collection.$mg";

	my $st1        = 'SELECT metagenome_id FROM metagenome WHERE collection_id=?';
	my $st2        = 'INSERT INTO metagenome(metagenome_id, collection_id, title, description) VALUES(?, ?, ?, ?)';
	my $metagenome = $this->safe_make_new_entry($prefix, $st1, $st2, [$collection], [$collection, $title, $description]);

	return $metagenome;
}

sub remove_metagenome_by_name {
	my $this = shift;
	my $name = shift;
	my $dbh  = $this->get_smashdb_handle;
	{
		my $sth  = $dbh->prepare('DELETE FROM metagenome WHERE metagenome_id=?');
		$sth->execute($name);
		$sth->finish();
	}
	$dbh->commit();
	$this->close_smashdb_handle();
}

sub make_new_assembly {
	my $this       = shift;
	my $metagenome = shift;
	my $assembler  = shift;
	my $as         = $this->get_smash_conf_value("assembly_prefix");
	my $prefix     = "$metagenome.$as";

	my $st1        = 'SELECT external_id FROM assembly WHERE metagenome_id=?';
	my $st2        = 'INSERT INTO assembly(external_id, metagenome_id, assembler) VALUES(?, ?, ?)';
	my $assembly   = $this->safe_make_new_entry($prefix, $st1, $st2, [$metagenome], [$metagenome, $assembler]);

	return $assembly;
}

sub remove_assembly_by_name {
	my $this = shift;
	my $name = shift;
	my $dbh  = $this->get_smashdb_handle;
	{
		my $sth  = $dbh->prepare('DELETE FROM assembly WHERE external_id=?');
		$sth->execute($name);
	}
	$dbh->commit();
	$this->close_smashdb_handle();
}

=item B<make_new_genepred>

Makes a new gene prediction. Any instance of subclasses of L<Smash::Analyses::GenePredictor> should call this function to get a new gene prediction id from Smash. It creates an entry in the C<gene_prediction> table with the C<program_id> of this instance.

=item B<remove_genepred_by_name>

Removes the gene prediction from the C<gene_prediction> table given the external id.

=back

=cut

sub make_new_genepred {
	my $this       = shift;
	my $assembly   = shift;
	my $predictor  = shift;
	my $assembly_id= $this->get_id_by_name("assembly", $assembly);
	my $gp         = $this->get_smash_conf_value("genepred_prefix");
	my $prefix     = "$assembly.$gp";

	my $st1        = 'SELECT external_id FROM gene_prediction WHERE assembly_id = ?';
	my $st2        = 'INSERT INTO gene_prediction(external_id, assembly_id, predictor) VALUES(?, ?, ?)';
	my $genepred   = $this->safe_make_new_entry($prefix, $st1, $st2, [$assembly_id], [$assembly_id, $predictor]);

	return $genepred;
}

sub remove_genepred_by_name {
	my $this = shift;
	my $name = shift;
	my $dbh  = $this->get_smashdb_handle;
	{
		my $sth  = $dbh->prepare('DELETE FROM gene_prediction WHERE external_id=?');
		$sth->execute($name);
	}
	$dbh->commit();
	$this->close_smashdb_handle();
}

=head2 Other methods

=cut

#########################
# Fasta utility functions
#########################

=over 4

=item C<filter_input_file($input_file, $output_file, $field_idx, $filter_hash)>

Selects a subset of input file based on field-level filtering. You can filter a file using
a defined list of accepted values in a specified field.
Reads $input_file and checks the field at field position C<$field_idx> (zero-based white-space delimited).
If this field has a key in C<%$filter_hash>, then that line is printed to $output_file.

e.g.,

	$smash->filter_input_file($in, $out, 2, {seq1=>1, seq2=>1, seq5=>2})

will print each line in C<$in> that has C<seq1> or C<seq2> or C<seq5> in 3rd column (0-based field index 
2 refers to the 3rd column)

=cut

=item B<get_gc_percent>

Returns GC percent of the given string, within the bounds. It is called as:

	$smash->get_gc_percent($string, $lower_bound, $upper_bound);

For example,

	$smash->get_gc_percent("ggggcccccgggcgcgcgcgacggcgcgcc");         
		# returns 98 (49 out of 50)
	$smash->get_gc_percent("ggggcccccgggcgcgcgcgacggcgcgcc", 20, 90); 
		# returns 90 (since 90 is the upper bound)

=item B<pretty_fasta>

Returns a string broken into fixed number of characters per line. Useful to format Fasta sequences.
By default it breaks them into 80 characters per line, but this can be specified in functional call.

	$smash->pretty_fasta($string);     
		# returns a string with newline every 80 characters
	$smash->pretty_fasta($string, 50); 
		# returns a string with newline every 50 characters

=item B<pretty_qual>

Returns a string broken into fixed number of quality values per line. Useful to format quality sequences.
By default it breaks them into 17 characters per line, but this can be specified in functional call.
Input string is a space delimited set of integers.

	$smash->pretty_qual($string);     
		# returns a string with newline every 17 characters
	$smash->pretty_qual($string, 10); 
		# returns a string with newline every 10 characters


=back

=cut

sub get_gc_percent {
	my $this  = shift;
	my $seq   = shift;
	my $lower = shift || 0;
	my $upper = shift || 100;

	if (!$seq || $seq eq "") {
		return $lower;
	}

	my $at  = uc($seq);
	my $gc  = uc($seq);
	my $gcp = $lower;

	# Get the GC percent
	$at      =~ s/[^AT]//g;
	$gc      =~ s/[^GC]//g;
	$gcp     = 100 * length($gc)/(length($gc)+length($at));
	if ($gcp < $lower) {
		$gcp = $lower;
	} elsif ($gcp > $upper) {
		$gcp = $upper;
	}
	return $gcp;
}

sub pretty_fasta {
	my $this = shift;
	my $seq  = shift;
	my $line = shift || 80;
	my $skip = 0;
	my $length = length($seq);
	my $pretty = "";
	while ($skip < $length) {
		my $sub = substr($seq, $skip, $line);
		$pretty .= ($sub."\n");
		$skip += $line;
	}
	return $pretty;
}

sub pretty_qual {
	my $this = shift;
	my $seq  = shift;
	my $line = shift || 17;
	my @val  = split(/\s+/, $seq);
	my $skip = 0;
	my $length = scalar(@val);
	my $pretty = "";
	while ($skip < $length) {
		my $end = $skip+$line-1;
		if ($end > $#val) {
			$end = $#val;
		}
		$pretty .= (join(" ", @val[$skip..$end]));
		$pretty .= "\n";
		$skip += $line;
	}
	return $pretty;
}

sub revcomp {
	my $this = shift;
	my $seq = shift;
	$seq = reverse($seq);
	$seq =~ tr/atgcATGC/tacgTACG/;
	return $seq;
}

# Get the full path of a script

sub get_script_path {
	my $this   = shift;
	my $script = shift;
	if (-f "$SMASH_SCRIPT_LOCATION/$script") {
		return "$SMASH_SCRIPT_LOCATION/$script";
	}

	# Look for the PERL5LIB

	my $smash_location = $INC{"Smash/Global.pm"};
	$smash_location =~ s^Smash/Global\.pm^^;
	if (-f "$smash_location/bin/$script") {
		return "$smash_location/bin/$script";
	}

	die "Cannot find $script in your Smash installation!\n";
}

############################################
# Get the first word of the fasta header
############################################

sub process_fasta_header {
	shift; # Don't need the object!
	my $def = shift;
	$def =~ s/>//;
	$def =~ s/\s.*//;
#	if ($def =~ /gi\|\d+\|gb\|([^\|]+)\|/) {
#		$def = $1;
#	}
	return $def;
}

############################################
# General utilities
############################################

# Print a progress bar

sub progress_bar {
	use POSIX qw(ceil);
	my $this    = shift;
	my $number  = shift;
	my $dotsize = shift || 10;

	# zero gets 0 printed

	return "0" if $number == 0;

	# See if this needs to be printed

	return "" if ($number%$dotsize != 0);

	# otherwise print $number/$dotsize

	my $n = int($number/$dotsize);       # number to display
	my $ten = 10*ceil($n/10);            # it's next 10er
	my $str = sprintf("%d", $ten);       # next 10er as string
	my $ord = length($str);              # length of the next 10er string
	my $p;
	if ($ten-$n < $ord) {                # Are we within the length of next 10er string?
		$p = substr($str, $ord-$ten+$n-1, 1);
		if ($n%100 == 0) {
			$p .= "\n";
		}
	} else {
		$p = ".";
	}
	return $p;
}

sub filter_input_file {
	my $this = shift;
	my ($input_file, $output_file, $field_idx, $filter_hash) = @_;

	# Read the text file
	open(INPUT, "<$input_file") || die "Cannot open $input_file: $!";
	open(OUTPUT, ">$output_file") || die "Cannot open $output_file: $!";
	HSP:while (<INPUT>) {
		next HSP if (m/^#/);
		my ($field) = (split(/\s+/))[$field_idx];
		if (defined($filter_hash->{$field})) {
			print OUTPUT $_;
		}
	}
	close(OUTPUT);
	close(INPUT);
	return 0;
}

sub safe_open_file {
	use Smash::Utils::File;
	my $file_name = shift;
	my $gzip = 0;

	# already REF glob?

	if (ref($file_name) =~ /GLOB/) {
		die "Cannot pass a file reference to safe_open_file!";
	}

	# compressed?

	if ($file_name =~ /\.gz$/) {
		$gzip = 1;
	}
	my $file = new Smash::Utils::File(NAME => $file_name, GZIP => $gzip);
	return $file;
}

sub safe_close_file {
	my $file = shift;
	$file->close();
}

sub suicide {
	my $usage = shift;
	my $msg   = shift;
	print STDERR "Error: $msg\n";
	print $usage;
	exit(1);
}

sub get_next_uncommented_line {
	my $this = shift;
	my $fh = shift;
	my $line;

	return 0 if eof $fh;

	# get the next non-comment line

	do {
		$line = <$fh>;
	} while ($line && $line =~ /^#/);

	return 0 unless $line;
	return 0 if ($line =~ /^#/);
	chomp($line);
	return $line;
}

=head2 Non-object-oriented functions

=over 4

=item C<get_median(@list)>

returns the median value from the given list

=back

=cut

sub max2 {
	my ($a, $b) = @_;
	return $a if ($a > $b);
	return $b;
}

sub get_median {
	my $list = shift;
	if (ref($list) eq "HASH") {
		my @x = values(%$list);
		$list = \@x;
	} elsif (ref($list) eq "ARRAY") {
	} else {
		die "get_median() requires reference to array or hash\n";
	}
	my $n = scalar(@$list);
	my $mid = int($n/2);
	my $median;
	if ($n%2 == 0) {
		$median = ($list->[$mid-1] + $list->[$mid])/2;
	} else {
		$median = $list->[$mid];
	}
	return $median;
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

	#print "#$command\n";
	$status = system($command);
	if ($status != 0) {
		warn "Error executing command:\n\t$command\nError Status:\t$status\n";
	}
	return $status;
}

1;

=head1 DEPENDENCIES

Smash requires certain perl modules to be installed on your local system. The most notable ones are:

L<FAlite|lib::FAlite>, L<FQlite|lib::FQlite>, L<XML::Parser>, L<DBI>

Depending on the database you use, L<DBD::mysql> or L<DBD::SQLite> must also be installed.

L<FAlite|lib::FAlite> (courtesy of Ian Korf, ifkorf@ucdavis.edu) and L<FQlite|lib::FQlite> are included in the Smash distribution.

=cut

