#########################################################################################
#########################################################################################
## Generic cluster management
#########################################################################################
#########################################################################################

package Smash::Utils::Cluster;
use strict;
use warnings;
our @ISA = qw(Smash::Utils);

=head1 NAME

Smash::Utils::Cluster - General Cluster Management library

=head1 SYNOPSIS

	use Smash::Utils::Cluster;

	my $cluster = new Smash::Utils::Cluster::sigma();

=cut

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

use overload '""' => 'name';

sub name      { shift->{NAME}}
sub queue     { shift->{QUEUE}}
sub type      { shift->{TYPE}}
sub qsub      { shift->{QSUB}}
sub qstat     { shift->{QSTAT}}
sub job_id    {shift->{JOB_ID}}
sub sub_task_id { shift->{SUB_TASK_ID}}
sub preamble  { shift->{PREAMBLE}}

1;

#########################################################################################
#########################################################################################
## PBS family
#########################################################################################
#########################################################################################

package Smash::Utils::Cluster::PBS;
use strict;
use warnings;
our @ISA = qw(Smash::Utils::Cluster);

############################################
# Constructor
############################################

sub new {
	my $class      = shift;
	my %params     = @_;
	my $self       = $class->SUPER::new(@_);
	$self->{TYPE}  = "PBS";
	$self->{QSUB}  = "qsub"; 
	$self->{QSTAT} = "qstat -f";
	$self->{SUB_TASK_ID} = '$PBS_ARRAY_INDEX';
	$self->{JOB_ID} = '$PBS_JOBID';
	return $self
}

sub get_array_job_options {
	my $this       = shift;
	my $array_size = shift;
	if ($array_size == 1) {
		return "";
	}
	return "-J 1-$array_size";
}

sub get_hardware_options {
	my $this = shift;
	my %requests = @_;
	my $mem = $requests{MEMORY};
	   $mem = sprintf("%d", $mem*1.1);
	my $cpus = $requests{CPUS};

	return "-l select=1:ncpus=$cpus:mem=${mem}mb";
}

sub set_defaults {
	my $this = shift;
	$this->{QSUB}     = "/usr/pbs/bin/qsub";
	$this->{QSTAT}    = "/usr/pbs/bin/qstat -f";
	$this->{PREAMBLE} = "source /usr/local/bin/pbs.sh;";
}

sub qsub {
	my $this = shift;
	my $qsub = $this->{QSUB};
	if ($this->queue) {
		$qsub .= (" -q ".$this->queue);
	} else {
		warn "PBS: No queue specified. Using default queue.\n";
	}
	return $qsub;
}

###########################################
# Check status of a job in the cluster
#
# returns:
#	0 - job is done
#	1 - job is still running
###########################################
sub get_job_status {
	my $this         = shift;
	my $job_id       = shift;
	my $qstat        = $this->qstat;
	my $output       = $this->execute_capture("$qstat $job_id");

	# Empty output usually means job doesnt exist
	if ($output eq "") { # job is finished!
		return 0;
	}

	# Output must contain at least a line with "job_state = "
	my @lines = split(/\n/, $output);
	my @state = grep {$_ =~ /job_state = /} @lines;

	# If "job_state = " was not found, may be server was unreachable
	if (scalar(@state) != 1) {
		warn "Error while checking status for job $job_id. Perhaps PBS master server was unreachable.";
	} else {
		my $real_state = $state[0];
		$real_state =~ s/job_state = //;
		$real_state =~ s/\s*$//;
		$real_state =~ s/^\s*//;
		if ($real_state eq "X") {
			return 0;
		}
	}
	return 1;
}

1;

#########################################################################################
#########################################################################################
## SGE family
#########################################################################################
#########################################################################################

package Smash::Utils::Cluster::SGE;
use strict;
use warnings;
our @ISA = qw(Smash::Utils::Cluster);

sub new {
	my $class       = shift;
	my %params      = @_;
	my $self        = $class->SUPER::new(@_);
	$self->{TYPE}   = "SGE";
	$self->{QSUB}   = "qsub";
	$self->{QSTAT}  = "qstat -j";
	$self->{JOB_ID} = '$JOB_ID';
	$self->{SUB_TASK_ID} = '$SGE_TASK_ID';
	return $self;
}

sub get_array_job_options {
	my $this       = shift;
	my $array_size = shift;
	if ($array_size == 1) {
		return "";
	}
	return "-t 1-$array_size";
}

sub get_hardware_options {
	my $this = shift;
	my %requests = @_;
	my $mem = $requests{MEMORY};
	   $mem = sprintf("%d", $mem*1.1);
	my $cpus = $requests{CPUS};
	my $cpu_options = "";
	if ($cpus && $cpus > 1) {
		$cpu_options = $this->get_cpu_options($cpus);
	}

	return "-l h_vmem=${mem}M $cpu_options";
}

=begin DISABLE_INSIDE_EMBL

sub get_cpu_options {
	my $this = shift;
	die "SGE parallel environment is queue specific. Please implement get_cpu_options for $this\n";
}

=cut

sub get_cpu_options {
	my $this = shift;
	my $cpus = shift;
	return "-pe smp $cpus";
}

###########################################
# Check status of a job in the cluster
#
# returns:
#	0 - job is done
#	1 - job is still running
###########################################
sub get_job_status {
	my $this         = shift;
	my $job_id       = shift;
	my $qstat        = $this->qstat;
	my $output       = $this->execute_capture("$qstat $job_id");
	if ($output eq "") { # job is finished!
		return 0;
	}
	my @lines = split(/\n/, $output);
	my @state = grep {$_ =~ /job_state = /} @lines;
	if (scalar(@state) > 1) {
		die "Error while checking status for job $job_id";
	}
	my $real_state = $state[0];
	$real_state =~ s/job_state = //;
	$real_state =~ s/\s*$//;
	$real_state =~ s/^\s*//;
	if ($real_state eq "X") {
		return 0;
	}
	return 1;
}

1;

#############################################################
#############################################################
##             EMBL site specific clusters                 ##
#############################################################
#############################################################

#############################################################
##                  otto                                   ##
#############################################################

package Smash::Utils::Cluster::otto;
use strict;
use warnings;
our @ISA = qw(Smash::Utils::Cluster::SGE);

sub new {
	my $class  = shift;
	my $self   = $class->SUPER::new(@_);
	$self->{QSUB}  = "/net/c1/local/sge/bin/lx24-amd64/qsub";
	$self->{QSTAT} = "/net/c1/local/sge/bin/lx24-amd64/qstat -j";
	return $self;
}

sub preamble {
	return "source /etc/profile.d/sge_settings.sh;";
}

1;

#############################################################
##                  sigma                                  ##
#############################################################

package Smash::Utils::Cluster::sigma;
use strict;
use warnings;
our @ISA = qw(Smash::Utils::Cluster::SGE);

sub new {
	my $class  = shift;
	my $self   = $class->SUPER::new(@_);
	$self->{QSUB}  = "/usr/sgeadmin/bin/lx24-amd64/qsub";
	$self->{QSTAT} = "/usr/sgeadmin/bin/lx24-amd64/qstat -j";
	return $self;
}

sub preamble {
	return "source /etc/profile.d/settings.sh;";
}

sub get_cpu_options {
	my $this = shift;
	my $cpus = shift;
	return "-pe smp $cpus";
}

1;

#############################################################
##                  epsilon                                ##
#############################################################

package Smash::Utils::Cluster::epsilon;
use strict;
use warnings;
our @ISA = qw(Smash::Utils::Cluster::SGE);

sub new {
	my $class  = shift;
	my $self   = $class->SUPER::new(@_);
	$self->{QSUB}  = "/opt/gridengine/bin/lx26-amd64/qsub";
	$self->{QSTAT} = "/opt/gridengine/bin/lx26-amd64/qstat -j";
	return $self;
}

sub preamble {
	return "source /etc/profile.d/sge-binaries.sh;";
}

sub get_cpu_options {
	my $this = shift;
	my $cpus = shift;
	return "-pe smp $cpus";
}

1;
