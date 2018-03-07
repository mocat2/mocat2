############################################
#  Generic sequence assembly program handler
#  The SUPER class
############################################

package Smash::Analyses::Assembler;

#define parent
use base "Smash::Analyses";

# define children
use Smash::Analyses::Assembler::Arachne;
use Smash::Analyses::Assembler::Celera;
use Smash::Analyses::Assembler::Forge;
use strict;

############################################
# Constructor
############################################

sub new {
	my $class  = shift;
	my %params = @_;

	if (!defined($params{'ASSEMBLER'})) {
		$params{'ASSEMBLER'} = "Generic";
	}
	if (!defined($params{'TYPE'})) {
		$params{'TYPE'} = "Assembler";
	}

	bless {
		%params
	}, $class;
}

############################################
############################################
##    Member variables through methods    ##
############################################
############################################

sub _overload    {shift->{ASSEMBLER}}
sub name         {shift->{NAME}}

sub single_genome{shift->{SINGLE_GENOME}}
sub extra_options{shift->{EXTRA_OPTIONS}}
sub assembly_dir {shift->{ASSEMBLY_DIR}}

############################################
# This is the output from the assembler
############################################

sub assembly_fasta {
	return shift->{ASSEMBLY_FASTA};
}

sub set_assembly_fasta {
	my $this = shift;
	$this->{ASSEMBLY_FASTA} = shift;

}

############################################
# This is the final output from Smash
############################################

sub output_fasta {
	return shift->{OUTPUT_FASTA};
}

sub set_output_fasta {
	my $this = shift;
	$this->{OUTPUT_FASTA} = shift;
}

sub genome_size {
	return shift->{GENOME_SIZE};
}

sub set_genome_size {
	my $this = shift;
	$this->{GENOME_SIZE} = shift;
}

sub workspace {
	my $this      = shift;
	if ($this->{WORKSPACE}) {
		return $this->{WORKSPACE};
	}
	my $workspace = sprintf("%s/%s/%s/%s", $this->data_dir, $this->get_smash_conf_value("workspace_dir"), $this->type, $this->software_name);
	return $workspace;
}

############################################
############################################
##    Pipeline functions                  ##
############################################
############################################

#########################################
# Init the object
# Right now it does:
#	. SUPER::init()
#########################################

sub init {
	use File::Path;
	my $this = shift;

	$this->parse_config();

	# parse collection id from metagenome id
	# you need collection to be set up before Core::init() is called

	my ($collection, $metagenome);
	if ($this->assembly) {
		($collection, $metagenome) = $this->parse_assembly_id($this->assembly);
		$this->{METAGENOME} = $metagenome;
	} else {
		$metagenome = $this->metagenome;
		$collection = $this->parse_metagenome_id($metagenome);
	}
	$this->{COLLECTION} = $collection;

	$this->SUPER::init();

	# Init software details
	
	$this->init_software_details();
}

sub prepare {
	my $this = shift;

	# Did the caller specify an assembly_id?
	# If not, make a new assembly id for this run and use it

	my $metagenome        = $this->metagenome;
	my $collection        = $this->parse_metagenome_id($metagenome);
	my $assembly          = $this->assembly;
	if ($assembly) {
		my $assembly_id = $this->get_id_by_name("assembly", $assembly);
		if (!$assembly_id) {
			warn "WARNING: You have requested Smash to use $assembly as an assembly id, but this does not exist!\n";
			$this->abort();
		}
	} else {
		my $prog_id   = $this->get_program_id("$this", $this->version, $this->options);
		$assembly     = $this->make_new_assembly($metagenome, $prog_id);
	}
	$this->{NAME} = $assembly;

	# figure out input and output dirs

	my $data_dir          = $this->data_dir;
	my $assembly_dir      = "$data_dir/metagenome_collections/$collection/metagenomes/$metagenome/assemblies/$assembly";

	mkpath $assembly_dir, {mode => $this->file_perm};
	$this->{ASSEMBLY_DIR} = $assembly_dir;
}

############################################
# If you run into error, remove the assembly id you just added to the database.
# There might still be cases where errors are not traceable, in which case the
# assemblies remain in the database!
############################################

sub abort {
	my $this = shift;
	$this->remove_assembly_by_name($this->name);
	$this->finish();
	die "ERROR  : Aborting assembly - see message above.\n";
}

############################################
# Assemble!
############################################

sub run {
	my $this = shift;

	$this->prepare();

	# If it is only finishing up, then skip the validation and assembly
	if ($this->{FINISH} != 1) {
		$this->validate();
		my $status = $this->assemble();
		if ($status == 31) {
			return;
		} elsif ($status != 1) { # not success
			$this->abort();
		} 
	}

	$this->post_assembly();
	print "<output>".$this->name."</output>\n";
	print "********************************************************\n";
	print "  Assembly id assigned for this assembly: ".$this->name."\n";
	print "********************************************************\n";
}

############################################
# Run the assembler program
############################################

sub assemble {
	my $this    = shift;
	my $command = $this->get_command_line();
	return $this->execute($command);
}

# pretty much copy the files sent in the list to the assembly directory
sub copy_assembly_files {
	use File::Copy;
	my $this = shift;
	my @files = @_;
	
	foreach my $file (@files) {
		copy($file, $this->assembly_dir) || die "Cannot copy $file: $!";
	}
}

1;

############################################
#  Generic sequence assembly program handler
############################################

package Smash::Analyses::Assembler::external;
use base "Smash::Analyses::Assembler";
use strict;

sub run {
	die "Smash::Analyses::Assembler::external is a placeholder only."
}

sub init_software_details {
}

sub options {return shift->{OPTIONS}}

sub _overload {return shift->{ASSEMBLER} || "external";}

1;
