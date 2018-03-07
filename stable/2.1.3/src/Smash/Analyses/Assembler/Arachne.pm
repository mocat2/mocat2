############################################
#  Arachne sequence assembly program handler
############################################

package Smash::Analyses::Assembler::Arachne;
use base "Smash::Analyses::Assembler";
use strict;

=head1 NAME

Smash::Analyses::Assembler::Arachne - Implementation of Arachne assembly software pipeline

=head1 SYNOPSIS

=head1 DESCRIPTION

This module performs iterative (or superiterative) assemblies
using Arachne assembler.

=head2 Default options

Arachne uses two executables: C<Assemblez> and C<Assemble>. For metagenomic assembly,
we have compared a few assemblies and chosen C<Assemble> to be the default since it
assembles more reads although the longest scaffold and contig N50 are shorter. The
default parameters used by Smash for C<Assemble> are:

			    ACE=True
		   one_ace_file=True
	  aggressive_correction=False
		    min_overlap=10
		 REINDEX_SUPERS=True
	 ignore_version_warning=True
		  FORCE_VERSION=True
		ENLARGE_CONTIGS=True
		 IMPROVE_SUPERS=True
		     PATCH_GAPS=True
		    k_for_merge=12
		   check_plates=True
		       maxcliq1=500;
		       maxcliq2=500;
	
C<Assemblez> provides better assemblies of single genomes, and is recommended by
the developers of Arachne. Here are the defaults used by Smash for C<Assemblez>.

			    ACE=True
		   one_ace_file=True
	  aggressive_correction=False
		    min_overlap=10
		       FAST_RUN=True
		       maxcliq1=500;
		       maxcliq2=500;
	    recycle_bad_contigs=True;
		    SW_GAP_STEP=True;
		       FAST_RUN=False;
		 mc_min_overlap=30;

=cut

sub new {
	my $class  = shift;
	my %params = @_;
	my $self   = $class->SUPER::new(%params);
	return $self;
}

############################################
############################################
##    Member variables through methods    ##
############################################
############################################

sub arachne_pre {shift->{ARACHNE_PRE}}
sub data        {shift->{DATA}}
sub iteration   {shift->{ITERATION}}
sub iterative   {shift->{ITERATIVE}} # is this assembly iterative?
sub mode        {shift->{MODE}}
sub software_name {"Arachne"};

sub exclusion_file{shift->{EXCLUSION_FILE}}

############################################
############################################
##  Non-Member variables through methods  ##
############################################
############################################

############################################
# Extra options for Arachne
############################################

############################################
# Construct the output directory from the 
# parameters
############################################

sub output_dir {
	my $this        = shift;
	my $arachne_pre = $this->arachne_pre;
	my $data        = $this->data;
	my $name        = $this->name;
	my $iteration   = $this->iteration;
	return "$arachne_pre/$data/${name}_run${iteration}";
}

sub assembly_fasta {
	my $this = shift;
	return $this->output_dir."/assembly.bases.gz";
}


############################################
############################################
##    Pipeline functions                  ##
############################################
############################################

############################################
# Initialize the object. Right now it does:
#	. Init the parent class
#	. Set the pkg directory
############################################

sub init {
	my $this = shift;
	$this->SUPER::init();
}

=head1 FUNCTIONS

=over 4

=item C<prepare()>

Prepares the assembly files in the Smash working directory: specifically, 
the reads, quals and xml files in 
F<fasta>, F<qual>, F<traceinfo> directories, respectively. It also copies the
F<reads_config.xml> file that is required by Arachne.

=cut

sub prepare {
	use File::Path;
	use File::Basename;

	my $this = shift;
	$this->SUPER::prepare();

	my $metagenome  = $this->metagenome;
	my $assembly    = $this->name;
	my $data        = $assembly;
	my $arachne_pre = $this->workspace."/Arachne_data";

	############
	# Set Arachne specific vars
	############

	$this->{DATA}        = $data;
	$this->{ARACHNE_PRE} = $arachne_pre;
	if ($this->single_genome) {
		$this->{ITERATIVE} = 0;
		$this->{MODE}      = "Assemblez";
	} else {
		$this->{ITERATIVE} = 1;
		$this->{MODE}      = "Assemble";
	}

	############
	# Prepare the assembly directory
	############

	my $base        = "$arachne_pre/$data";

	############
	# link the input files in a new directory in ARACHNE_PRE
	############

	foreach my $dir ("$base", "$base/fasta", "$base/qual", "$base/traceinfo") {
		mkpath $dir, {mode => $this->file_perm};
	}
	foreach my $file ($this->fasta_files($metagenome)) {
		symlink "$file", "$base/fasta/@{[basename($file)]}";
	}
	foreach my $file ($this->qual_files($metagenome)) {
		symlink "$file", "$base/qual/@{[basename($file)]}";
	}
	foreach my $file ($this->xml_files($metagenome)) {
		symlink "$file", "$base/traceinfo/@{[basename($file)]}";
	}

	my $exclusion_file = "$base/reads.to_exclude"; # This is the default file
	my $xmlconfig_file = "$base/reads_config.xml"; # This is the default file

	unlink $xmlconfig_file; # remove xmlconfig file from old runs, if any
	symlink "$arachne_pre/reads_config.xml", $xmlconfig_file;

	unlink $exclusion_file; # remove exclude file from old runs, if any
	$this->{EXCLUSION_FILE} = $exclusion_file;
}

############################################
# Validate if everything is set up right
# for Arachne
############################################

sub validate {
	my $this = shift;

	# Check for executable

	my $pkg_dir = $this->pkg_dir;
	my $exe     = $this->exe;
	my $status  = $this->execute("cd $pkg_dir && $exe 1>/dev/null");
	if ($status != 0) {
		warn "Cannot execute $pkg_dir/$exe. Please check your set-up!\n";
		$this->abort();
	}
}

############################################
# prepare the genome.size file for Arachne
############################################

=item C<prepare_genome_size()>

Arachne uses a file to specify the genome size. This creates this file at every iteration
so the right genome size is used by that iteration.

=cut

sub prepare_genome_size {
	my $this = shift;

	my $genome_size      = $this->genome_size;
	my $genome_size_file = sprintf("%s/%s/genome.size", $this->arachne_pre, $this->data);

	# make genome.size
	open(GENOMESIZE, ">$genome_size_file") || die "Cannot open $genome_size_file: $!";
	print GENOMESIZE "$genome_size\n";
	close(GENOMESIZE);
}

############################################
# Construct Arachne command line
############################################

=item C<get_command_line()>

Creates the command to run from all the arguments and options.

=cut

sub get_command_line {
	my $this       = shift;
	my $pkg_dir    = $this->pkg_dir;
	my $arachne_pre= $this->arachne_pre;
	my $data       = $this->data;
	my $name       = $this->name;
	my $iteration  = $this->iteration;
	my $exe        = $this->exe;
	my $output_dir = $this->output_dir;
	my $options    = $this->options;
	$name          = "${name}_run${iteration}";

	my $cmd_line = "cd $pkg_dir && ./$exe PRE=$arachne_pre DATA=$data RUN=$name $options 1>$output_dir/$name.log";
	return $cmd_line;
}

=item C<superiterate($max_iterations)>

Specific to Arachne. If C<--superiterate> was selected in 
L<doAssembly.pl|doAssembly>, this sub runs the superiteration.
Here is the pseudocode of the superiteration:

	init: genome_size, max=$max_iterations, k=12, k_max=20, m=10, m_max=19

	i=1
	while (genome_size > 1Mb)
		option="aggressive_correction=True k_for_merge=<k> min_overlap=<m>"
		assemble()
		if (i == max)
			genome_size = genome_size / 2
			m = min(m+1, m_max)
			k = min(k+1, k_max)
		i = (i+1)%max

	genome_size = 800Kb
	i=1
	while (genome_size > 100Kb)
		option="aggressive_correction=True k_for_merge=<k> min_overlap=<m> IMPROVE_SUPERS=False MERGE_SUPERCONTIGS=False"
		assemble()
		if (i == max)
			genome_size = genome_size / 2
			m = min(m+1, m_max)
			k = min(k+1, k_max)
		i = (i+1)%max

=cut

sub superiterate {
	my $this              = shift;
	my $max_iterations    = shift;
	my @modes             = @_;
	my $k_for_merge       = 12;
	my $k_for_merge_limit = 20;
	my $min_overlap       = 10; # Assemble(z) default
	my $min_overlap_limit = 19; # Other steps in Arachne use min_overlap=20, so keep it less

	$this->prepare();
	$this->validate();

	# until genome size is more than 1Mb
	my $count = 0;
	my $genome_size = $this->genome_size;
	while ($genome_size >= 1000000) {
		$this->{GENOME_SIZE}   = $genome_size;
		$this->{OPTIONS} = "aggressive_correction=True k_for_merge=$k_for_merge min_overlap=$min_overlap";
		foreach my $mode (@modes) {
			$this->{MODE} = $mode;
			my $status = $this->assemble();
			if ($status != 1) { # not success
				$this->abort();
			} 
		}
		if ($count == ($max_iterations-1)) {
			$genome_size /= 2;
			if ($min_overlap < $min_overlap_limit) {
				$min_overlap += 1;
			}
			if ($k_for_merge < $k_for_merge_limit) {
				$k_for_merge += 1;
			}
		}
		$count = ($count+1)%$max_iterations;
	}

	$genome_size = 800000;
	$count = 0;
	while ($genome_size >= 100000) {
		$this->{GENOME_SIZE}   = $genome_size;
		$this->{OPTIONS} = "aggressive_correction=True k_for_merge=$k_for_merge min_overlap=$min_overlap IMPROVE_SUPERS=False MERGE_SUPERCONTIGS=False";
		foreach my $mode (@modes) {
			$this->{MODE} = $mode;
			my $status = $this->assemble();
			if ($status != 1) { # not success
				$this->abort();
			} 
		}
		if ($count == ($max_iterations-1)) {
			$genome_size /= 2;
			if ($min_overlap < $min_overlap_limit) {
				$min_overlap += 1;
			}
			if ($k_for_merge < $k_for_merge_limit) {
				$k_for_merge += 1;
			}
		}
		$count = ($count+1)%$max_iterations;
	}
	$this->post_assembly();
	print "********************************************************\n";
	print "   Assembly id assigned for this assembly: ".$this->name."\n";
	print "********************************************************\n";
}

############################################
# Over-rides the parent method since Arachne 
# needs iterative assembly.
############################################

=item C<assemble()>

Arachne overrides C<assemble()> from L<Analyses::Assembler> since
it uses iterative assembly. Iterative assembly forces the assembler to 
keep assembling as long it can. In simpler terms, after every assembly, 
it takes the unassembled reads and reassembles them using the same
parameters, until it can no longer assemble. Each iteration runs as
a new assembly with its own data directory, called something like
F<MC20.MG1.AS1_run1>, F<MC20.MG1.AS1_run2> and so on. The last 
iteration is stored as a global variable C<$this->{ITERATION}>
so that the next call to C<assemble()> knows where it left off. 

=cut

sub assemble {
	my $this = shift;
	my $iteration = 1;

	$this->prepare_genome_size();

	ITERATION: while (1) {

		$this->{ITERATION} = $iteration;

		############
		# make output directory
		############

		my $output_dir     = $this->output_dir; # output for this iteration
		my $assembly_fasta = $this->assembly_fasta;

		mkpath $output_dir, {mode => $this->file_perm};

		############
		# run Arachne
		############

		my $command = $this->get_command_line($iteration);

		if (-f $assembly_fasta) {
			printf STDERR "Skipping iteration $iteration:\n\t%s exists\n", $assembly_fasta;
			printf STDERR "To redo this iteration, cleanup this assembly directory:\n\t";
			printf STDERR $output_dir."\n";
		} else {
			my $status = $this->execute($command);

			# Some error!
			# If this is a non-iterative assembly, then
			# an error is a problem. So report this.
			# If it is an iterative assembly, then an error
			# just means that this is the last iteration. 
			# So don't worry about it.

			if ($status != 0 && !$this->iterative) {
				return -1;
			}
		}

		############
		# set things up for next run
		############

		if (!-f $assembly_fasta) {last ITERATION;}

		# Make it look like $iteration=2 was attempted and finished with no assembly
		# thus non-iterative will run the iteration only once

		if (!$this->iterative) {
			$this->{ITERATION} = $iteration+1;
			last ITERATION;
		}

		$this->local_post_assembly();
		$iteration++;
	}
	return 1;
}

=item C<local_post_assembly()>

This is the post-assembly step that summarizes the assembly and makes the contig-to-read maps
and the contig fasta files after each iteration inside C<assemble()>.

=cut

sub local_post_assembly {
	my $this           = shift;
	my $output_dir     = $this->output_dir;
	my $exclusion_file = $this->exclusion_file;
	my @exclude        = ();

	# Check the assembly.reads file for readinfo
	if (open(ARACHNE, "<$output_dir/assembly.reads")) {
		while (<ARACHNE>) {
			chomp;
			if ( $_ =~ /^#/) { next;}
			if ( $_ eq "") { next;}
			my $name = (split(/\t/))[0];
			push(@exclude, $name);
		}
		close(ARACHNE);
	}

	# Check the assembly.unplaced file for low_quality or vector_or_host or previously excluded reads
	if (open(UNPLACED, "<$output_dir/assembly.unplaced")) {
		while (<UNPLACED>) {
			if ($_ =~ /deliberate/ || $_ =~ /low_quality/ || $_ =~ /vector_or_host/ || $_ =~ /mitochondrial/) {
				push(@exclude, (split(/\t/, $_))[0]);
			}
		}
		close(UNPLACED);
	}

	# Put them in exclusion_file
	my %hash = map { $_, 1 } @exclude;
	@exclude = keys %hash;

	open(EXCLUDE, ">$exclusion_file") || die "Cannot open file $exclusion_file: $!";
	foreach (@exclude) {
		print EXCLUDE "$_\n";
	}
	close(EXCLUDE);
}

=item C<post_assembly()>

This is the global post-assembly step at the end of the assembly after all the 
iterative/superiterative assemblies are done. It generates:

	1. Contig fasta file
	2. Contig-to-read mapping file in GFF format
	3. Scaffold-to-contig mapping file in GFF format

=back

=cut

sub post_assembly {
	use File::Temp;
	use File::Spec::Functions qw/tmpdir/;
	my $this        = shift;
	my $total_runs  = $this->iteration;

	# Output files
	my $contig2read;
	my $contig_fa;

	# Corresponding handles
	my $CONTIG_FA;
	my $CONTIG2READ;
	my $SCAF2CONTIG;

	############################################
	# Combine contigs from multiple iterations
	# and unassembled reads from the last 
	# successful iteration in iteration=1
	############################################

	$this->{ITERATION} = 1;
	my $contig_fa    = sprintf("%s/%s.contigs.fa", $this->output_dir, $this->name);
	my $contig2read  = sprintf("%s/%s.contig2read.gff", $this->output_dir, $this->name);
	my $scaf2contig  = sprintf("%s/%s.scaf2contig.gff", $this->output_dir, $this->name);
	open($CONTIG_FA,   ">$contig_fa")   || die "Cannot open $contig_fa: $!";
	open($CONTIG2READ, ">$contig2read") || die "Cannot open $contig2read: $!";
	open($SCAF2CONTIG, ">$scaf2contig") || die "Cannot open $scaf2contig: $!";

	for (my $i=1; $i<$total_runs; $i++) {
		####
		# Remap the contig headers and write to a new fasta file
		# Write the contig2read map for contigs
		####

		$this->{ITERATION} = $i;
		my $assembly_fasta = $this->assembly_fasta;
		open(ASSEMBLY, "gunzip -c $assembly_fasta |") || die "Cannot open gunzip pipe for $assembly_fasta: $!";
		my $fasta = new FAlite(\*ASSEMBLY);
		while (my $entry = $fasta->nextEntry) {
			my $def = $entry->def;
			my $seq = $entry->seq;
			$def =~ s/>contig_//;
			$def++; # Arachne starts contigs with 0
			$def = sprintf("%s.I%d.C%d", $this->name, $this->iteration, $def);
			print $CONTIG_FA ">$def\n";
			print $CONTIG_FA $this->pretty_fasta($seq);
		}
		close(ASSEMBLY);

		$this->contig2read($CONTIG2READ);
		$this->scaf2contig($SCAF2CONTIG);
	}

	####
	# Make a list of unplaced reads
	####

	my @unplaced_reads = ();
	if ($total_runs > 1) {
		$this->{ITERATION} = $total_runs-1;
		print STDERR "Generating list of unplaced reads\n";
		open(UNPLACED, "<@{[$this->output_dir]}/assembly.unplaced") || die "Cannot open unplaced file:$!";
		while (<UNPLACED>) {
			if ($_ !~ /deliberate/ && $_ !~ /low_quality/ && $_ !~ /vector_or_host/) {
				push(@unplaced_reads, (split(/\t/, $_))[0]);
			}
		}
		close(UNPLACED);
	}

	####
	# Write the unplaced reads to the new fasta file
	# and write the contig2read maps for these reads
	####

	$this->process_unplaced_reads(\@unplaced_reads, $CONTIG_FA, $SCAF2CONTIG, $CONTIG2READ, $this->name.".I".$this->iteration, $total_runs == 1);

	close($CONTIG_FA);
	close($CONTIG2READ);
	close($SCAF2CONTIG);

	$this->copy_assembly_files($contig_fa, $contig2read, $scaf2contig);
}

# Get sequences given by a list from a file containing a superset of sequences
# Gets handle as arguments, not the file names. Since these are from File::Temp, do not try to write to the files directly using the name
sub process_unplaced_reads {
	my $this = shift;
	my ($unplaced_reads, $fasta_handle, $scaf2contig, $contig2read, $prefix, $NO_ASSEMBLY) = @_;
	my $read_index = 1;

	# Get the list of definitions to be extracted
	my %List = map {$_ => 1} @$unplaced_reads;

	# Get the list of input fasta files
	my $input_dir = sprintf("%s/%s/fasta", $this->arachne_pre, $this->data);
	opendir(DIR, "$input_dir") || die "Cannot open dir $input_dir: $!";
	my @files = grep { -f "$input_dir/$_" } readdir(DIR);
	close(DIR);
	@files = map { "$input_dir/$_" } grep { $_ =~ /fasta/ } @files;

	# Get the relevant entries from each of those files
	$, = "\t";
	foreach my $file (@files) {
		open(FASTA, "<$file") || die "Cannot open $file: $!";
		my $fasta = new FAlite(\*FASTA);
		while(my $entry = $fasta->nextEntry) {
			my $def;
			my $name;
			my @words;
			$def   = $entry->def;
			$def   =~ s/^>//;
			@words = split(/\s/, $def);
			$name  = $words[$#words];
			$def   = $words[0];
			if ($NO_ASSEMBLY == 1 || defined($List{$name})) { # contig and scaffold ids are the same here
				my $length = length($entry->seq);
				my $contig = "$prefix.R$read_index";
				print $contig2read $contig, "$this", "read", 1, $length, $length, "+", ".", "read \"$def\";\n";
				print $scaf2contig $contig, "$this", "contig", 1, $length, $length, "+", ".", "contig \"$contig\";\n";
				print $fasta_handle ">$prefix.R$read_index\n";
				print $fasta_handle $this->pretty_fasta($entry->seq);
				$read_index++;
			}
		}
		close(FASTA);
	}
}

sub scaf2contig {
	my $this         = shift;
	my $OUTPUT_FH    = shift;
	my $links_file   = sprintf("%s/assembly.links", $this->output_dir);

	my $start        = 0;
	my $prev_scaf_id = -1;

	# Process scaf2contig output in assembly.links
	open(LINKS, "<$links_file") || die "Cannot open $links_file: $!";
	$, = "\t";
	while (<LINKS>) {
		chomp;
		if ( $_ =~ /^#/) { next;}
		if ( $_ eq "") { next;}
		my ($scaffold_id, $scaf_length, $contig_count, $contig_idx, $contig_id, $contig_length, $gap_before, $gap_after) = split(/\t/);
		###
		# Arachne starts scaffolds and contigs at 0
		###

		$scaffold_id++;
		$contig_id++;

		if ($prev_scaf_id != $scaffold_id) {
			$start = 0;
		}
		$start += $gap_before;
		my $contig = sprintf("%s.I%d.C%d", $this->name, $this->iteration, $contig_id);
		my $scaffold =sprintf("%s.I%d.S%d", $this->name, $this->iteration, $scaffold_id); 
		print $OUTPUT_FH $scaffold, "$this", "contig", $start+1, $start+$contig_length, $scaf_length, "+", ".", "contig \"$contig\";\n";
		$start += $contig_length;
		$prev_scaf_id = $scaffold_id;
	}
	close(LINKS);
}

sub contig2read {
	my $this         = shift;
	my $OUTPUT_FH    = shift;
	my $input_dir    = sprintf("%s/%s/fasta", $this->arachne_pre, $this->data);
	my $reads_file   = sprintf("%s/assembly.reads", $this->output_dir);

	# Generate a map of read name to fasta header since arachne has a weird convention

	# Process arachne output
	open(READS, "<$reads_file") || die "Cannot open $reads_file: $!";
	$, = "\t";
	while (<READS>) {
		chomp;
		if ( $_ =~ /^#/) { next;}
		if ( $_ eq "") { next;}
		my ($name, $status, $read_length, $left_trim, $trim_read_length, $contig, $contig_length, $read_start, $read_end, $strand, $partner, $partner_status, $partner_contig, $insert_size) = split(/\t/);

		###
		# Arachne starts scaffolds and contigs at 0
		###

		$contig++;
		if ($partner_contig ne "") {
			$partner_contig++;
		}

		my $contig_id = sprintf("%s.I%d.C%d", $this->name, $this->iteration, $contig);
		print $OUTPUT_FH $contig_id, "$this", "read", $read_start+1, $read_end+1, $contig_length, $strand, ".", "read \"$name\";";
		if ($partner_contig ne "") {
			printf $OUTPUT_FH (" mate_pair \"%s\"; contig \"%s.I%d.C%d\";", $partner, $this->name, $this->iteration, $partner_contig);
			if ($insert_size ne "") {
				print $OUTPUT_FH " insert_size \"$insert_size\";";
			}
		}
		print $OUTPUT_FH "\n";
	}
	close(READS);
}

sub exe     {shift->{MODE}}

sub options {
	my $this  = shift;
	my $mode  = $this->mode;
	my $extra = $this->extra_options; # from the wrapper script
	my %parameters = ();
	if ($mode eq "Assemble") {
		%parameters = ( 
						    ACE => "True",
					   one_ace_file => "True",
				  aggressive_correction => "False",
					    min_overlap => "10",
					 REINDEX_SUPERS => "True",
				 ignore_version_warning => "True",
					  FORCE_VERSION => "True",
					ENLARGE_CONTIGS => "True",
					 IMPROVE_SUPERS => "True",
					     PATCH_GAPS => "True",
					    k_for_merge => "12",
					   check_plates => "True"
				 );
		if (!$this->single_genome) {
			$parameters{maxcliq1} = "500";
			$parameters{maxcliq2} = "500";
		}
	} elsif ($mode eq "Assemblez") {
		%parameters = ( 
						    ACE => "True",
					   one_ace_file => "True",
				  aggressive_correction => "False",
					    min_overlap => "10",
					       FAST_RUN => "True"
				 );
		if (!$this->single_genome) {
			$parameters{maxcliq1}            = "500";
			$parameters{maxcliq2}            = "500";
			$parameters{recycle_bad_contigs} = "True";
			$parameters{SW_GAP_STEP}         = "True";
			$parameters{FAST_RUN}            = "False";
			#$parameters{n_haplotypes}       = "4";
			$parameters{mc_min_overlap}      = "30";
		}
	}

	# options set inside this module

	if ($this->{OPTIONS}) {
		map {
			my ($key, $value) = split("="); 
			if ($parameters{$key}) {
				$parameters{$key} = $value;
			}
		} split(" ", $this->{OPTIONS});
	}

	# if any of the set options need to be changed through EXTRA_OPTIONS, do it now
	# anything that has not been set so far will be ignored (as a failsafe against unknown option failure)

	if ($extra) {
		map {
			my ($key, $value) = split("="); 
			if ($parameters{$key}) {
				$parameters{$key} = $value;
			}
		} split(" ", $extra);
	}

	# return the newly constructed options string

	return join(" ", 
		map {
			sprintf "%s=%s", $_, $parameters{$_}
		} sort {$a cmp $b} keys %parameters
	);
}

1;
