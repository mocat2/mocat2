#########################################################################################
#########################################################################################
## Generic cluster job management
#########################################################################################
#########################################################################################

package Smash::Utils::ClusterJob;
use base("Smash::Utils");
use Smash::Utils::Cluster;
use Cwd;

=head1 NAME

Smash::Utils::ClusterJob - General cluster job management module

=head1 SYNOPSIS

	my %options = ( 
			NAME        => "testjob",
			TYPE        => "SGE",
			MEMORY      => 8000,      #default 2000
			WORKING_DIR => "mydir",   #default cwd
			EODIR       => "eofiles", #default cwd
			QUEUE       => "all.q",
			POLLING     => 90,        #default 60
			CPUS        => 4,         #default 1
			GROUP       => 2,         #default 1
			LAMMPI     => 1,         #default 0
		);
	my $job = new Smash::Utils::ClusterJob(%options);
	$job->submit_commands("echo first_command", "echo second_command", "echo last_command");

This will submit an array job with the three commands in the argument to 
C<submit_commands()> to the queue called "all.q" in the SGE batch queuing
system that is visible to the current host 
with the following options.

	hard memory limit: 8000 MB
	working directory: mydir/
	error/output to  : eofiles/
	cpus allotted    : 4

GROUP specifies if the commands should be grouped together at all. This is useful for very short jobs, where the cluster management overhead time is more than the actual run time of the job. In this example, commands are grouped in sets of two. Let's assume this was submitted as job array 100. Task 1 of Job 100 will execute the two lines:

	echo first_command
	echo second_command

Task 2 will execute the following line:

	echo last_command

=cut

sub new {
	use Env qw[$HOME];
	my $class  = shift;
	my %params = @_;

	my $system = $params{TYPE}  || die "Needs TYPE for a job to be run";
	my $queue  = $params{QUEUE};
	my $instance = "Smash::Utils::Cluster::$system"->new(QUEUE => $queue);
	$params{INSTANCE} = $instance;

	if (!defined($params{WORKING_DIR})) {
		$params{WORKING_DIR} = Cwd::cwd;
	}
	if (!defined($params{MEMORY})) {
		$params{MEMORY} = 2000;
	}
	if (!defined($params{EODIR})) {
		$params{EODIR} = Cwd::cwd;
	}
	if (!defined($params{SHELL})) {
		$params{SHELL} = "/bin/bash";
	}
	if (!defined($params{CPUS})) {
		$params{CPUS} = 1;
	}
	if (!defined($params{GROUP})) {
		$params{GROUP} = 1;
	}
	if (!defined($params{EXTRA_ARGS})) {
		$params{EXTRA_ARGS} = "";
	}
	if (!defined($params{POLLING})) {
		$params{POLLING} = 60;
	}

	# single and group dont go together
	if ($params{SINGLE}) {
		$params{GROUP} = 1;
	}

	bless {
		%params
	}, $class;
}

sub poll_pause {shift->{POLLING}}
sub mem_limit  {shift->{MEMORY}}
sub job_name   {shift->{NAME}}
sub working_dir{shift->{WORKING_DIR}}
sub queue      {shift->{QUEUE}}
sub instance   {shift->{INSTANCE}}
sub shell      {shift->{SHELL}}
sub eo_dir     {shift->{EODIR}}
sub group      {shift->{GROUP}}
sub single     {shift->{SINGLE}}
sub extra_args {shift->{EXTRA_ARGS}}
sub cpus       {shift->{CPUS}}
sub lammpi    {shift->{LAMMPI}}

sub lam_prolog {
	use File::Basename;
	my $this = shift;
	my (undef, $binpath, undef) = fileparse($0);
	return "$binpath/lamboot.sh ".$this->instance->type;
}

sub lam_epilog {
	use File::Basename;
	my $this = shift;
	my (undef, $binpath, undef) = fileparse($0);
	return "$binpath/lamhalt.sh ".$this->instance->type;
}

########
# submit the jobs in the array and return the job_id
########
sub submit_commands {
	use File::Basename;
	use File::Temp;
	use Cwd;
	use Fcntl;
	use POSIX qw(ceil);

	my $this         = shift;
	my (@commands)   = @_;

	# job properties

	my $working_dir  = $this->working_dir;
	my $eo_dir       = $this->eo_dir;
	my $job_name     = $this->job_name;
	my $group        = $this->group;
	my $shell        = $this->shell;
	my $xtra_args    = $this->extra_args;

	if ($this->single) {
		$group = scalar(@commands);
	}

	# cluster properties useful here

	my $instance     = $this->instance;
	my $preamble     = $instance->preamble;
	my $qsub         = $instance->qsub;
	my $SUB_TASK_ID  = $instance->sub_task_id;
	my $JOB_ID       = $instance->job_id;

	# Set up configurable job properties

	my $mem_limit    = sprintf("%d", $this->mem_limit*1024*1.1);
	my $hardware_options = $instance->get_hardware_options(MEMORY => $this->mem_limit, CPUS => $this->cpus);

	my $array_size   = ceil(@commands/$group);
	my $array_options= $instance->get_array_job_options($array_size);

	# LAM/MPI?
	# If this is LAM/MPI, add the lamboot and lamhalt steps. 
	# Gets complicated since PBS has it integrated and SGE doesnt.
	# I have shell scripts called lamboot.sh and lamhalt.sh in Smash's bin directory that 
	# sets things up and calls lamboot/lamhalt.

	if ($this->lammpi == 1) {
		@commands = map {$this->lam_prolog." && $_; ".$this->lam_epilog} @commands;  # && because if lamboot failed we want to abort!
	}

	# Job script

	my $job_file     = new File::Temp(TEMPLATE => "smashjobXXXXXXXX", DIR => Cwd::cwd(), SUFFIX => '.queue', UNLINK => 1);
				# Can be used as name in string context and handle in handle context
				# We don't have to unlink it - the object does it when it quits.

	# Generate the job file
	
	####
	# Prologue
	####
	my $sub_task_id_lvalue = $SUB_TASK_ID;
	$sub_task_id_lvalue =~ s/^\$//;
	print $job_file <<EOF;
echo "Job $JOB_ID started  on \$HOSTNAME at \`date\`";
source \$HOME/.bashrc;
#ulimit -v $mem_limit; # commented out since Yan is fixing things with h_vmem as well
if [ -z "$SUB_TASK_ID" -o "$SUB_TASK_ID" == "undefined" ]; then
	$sub_task_id_lvalue=1;
fi
case $SUB_TASK_ID in
EOF

	####
	# Job itself
	####
	my $job_array_size = 0;
	for (my $i = 0; $i <= $#commands; $i+=$group) {
		$job_array_size++;
		print $job_file <<EOF;
$job_array_size)
	echo '<cwd>$working_dir</cwd>'
EOF
		for (my $j = $i; $j < $i+$group; $j++) {
			if ($j <= $#commands) {
				print $job_file <<EOF;
	echo '<Command>'
	cat <<EOL
	$commands[$j]
EOL
	echo '</Command>'
## ACTUAL JOB
	cd $working_dir;
	$commands[$j]
EOF
			}
		}
		print $job_file "	ERR=\$?\n";
		print $job_file "	;;\n";
	}
	print $job_file "esac\n";

	####
	# Epilogue
	####
	print $job_file <<EOF;
if [ "\$ERR" != "0" ]; then 
	echo "Your program exited with error \$ERR. It could be because you exceeded your memory request ($mem_limit MB) or your job sucks." 1>&2; 
fi
echo "Job $JOB_ID finished on \$HOSTNAME at \`date\`";
EOF
	close($job_file);

	# submit the job to the queue
	my $real_job_name = substr($job_name, 0, 15);
	my $command       = "$preamble $qsub $xtra_args -e $eo_dir -o $eo_dir -N $real_job_name -S $shell -V -v PATH $hardware_options $array_options $job_file"; # string context of $job_file
	my ($job_array)   = $this->execute_capture($command);
	#my $job_array = "Something";
	print "#$command\n";

	if (!defined($job_array)) {
		warn "Job not submitted\n";
		return -1;
	}

	chomp($job_array);
	if ($job_array =~ /^Your job (\d+) /) {
		$job_array = $1;
	} elsif ($job_array =~ /^Your job-array (\d+)\./) {
		$job_array = $1;
	} elsif ($job_array =~ /^(\d+\[?\]?)\./) {
		$job_array = $1;
	} else {
		warn "Cannot parse output. Potential error\n";
		warn "Please check configuration of your queueing environment\n";
		return -1;
	}
	print "Job submitted as $job_array\n";
	return ($job_array, $array_size);
}

########
# wait for a given job_id
########
sub wait_for_job {
	my $this         = shift;
	my $job_id       = shift;
	my $instance      = $this->instance;
	my $qstat        = $instance->qstat;

	# keep checking for the running job
	my $output = $this->execute_capture("$qstat $job_id");
	chomp($output);
	while ($output ne "") {
		sleep($this->poll_pause);
		$output = $this->execute_capture("$qstat $job_id");
		chomp($output);
	}
}

########
# submit job and wait for completion
########
sub run_commands {
	my $this        = shift;
	my (@commands)  = @_;

	#submit the job and wait for completion
	$this->wait_for_job($this->submit_commands(@commands));

}

1;
