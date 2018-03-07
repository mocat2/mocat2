############################################
#  Forge sequence assembly program handler
############################################

package Smash::Analyses::Assembler::Forge;
use base "Smash::Analyses::Assembler";
use strict;

############################################
############################################
##    Member variables through methods    ##
############################################
############################################

############################################
############################################
##  Non-Member variables through methods  ##
############################################
############################################

############################################
# Output directory
############################################

sub output_dir {
	return shift->workspace;
}

############################################
# Extra options for Forge
############################################

sub options {
	return "--metagenomic";
}

sub cpus {
	return "2";
}

sub tech    {shift->{TECH}}
sub refseq  {shift->{REFSEQ}}

sub workspace {
	my $this = shift;
	my $workspace = $this->SUPER::workspace();
	return $workspace."/".$this->name;
}

############################################
# Construct Forge command line
############################################

sub get_command_line {
	my $this        = shift;
	my $name        = $this->name;
	my $pkg_dir     = $this->pkg_dir;
	my $genome_size = $this->genome_size;
	my $tech        = $this->tech;
	my $options     = $this->options;
	my $cpus        = $this->cpus;
	my $cluster     = $this->cluster;
	my $refseq      = $this->refseq;
	my $workspace   = $this->workspace;
	my $cmd_line    = "python $pkg_dir/scripts/runforge.py $name $genome_size --refseq=$refseq --forge=$pkg_dir --tech=$tech --cpu=$cpus $options 1>$workspace/$name.out 2>$workspace/$name.err";
	#print $cmd_line."\n";
	return $cmd_line;
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
	my $name = $this->name;

	$this->{REFSEQ}    = $this->pkg_dir."/C_elegans.fa";

	$this->set_assembly_fasta(sprintf("%s/%s.keep.fasta", $this->workspace, $name));
	$this->set_output_fasta(sprintf("%s/%s.contigs.fa", $this->workspace, $name));
}

############################################
# Validate if everything is set up right
# for Forge
############################################

sub validate {
	use Env qw[$HOME];
	use File::Path;
	use File::Copy;

	my $this       = shift;
	my $name       = $this->name;
	my $metagenome = $this->metagenome;
	my $fasta      = "$metagenome.fasta";
	my $qual       = "$metagenome.qual";
	my $pkg_dir    = $this->pkg_dir;
	my $read_dir   = $this->read_dir($metagenome);
	my $workspace  = $this->workspace;
	my @techs      = ();

	# Check for fasta and quality files
	if (! -f "$read_dir/$fasta") {
		die "Cannot find input Fasta file: $read_dir/$fasta";
	}
	if (! -f "$read_dir/$qual") {
		die "Cannot find input Quality file: $read_dir/$qual";
	}

	mkpath $workspace, {mode => $this->file_perm};
	symlink "$read_dir/$fasta", "$workspace/$name.fasta";
	symlink "$read_dir/$qual", "$workspace/$name.qual";

	# Check for technologies
	{
		my $dbh = $this->get_db_handle;
		my $sth = $dbh->prepare('SELECT DISTINCT type FROM library l JOIN sample s ON l.sample_id=s.sample_id WHERE s.metagenome_id=? ORDER BY type ASC');
		$sth->execute($metagenome);
		while (my ($type) = $sth->fetchrow_array()) {
			push(@techs, $type);
		}
	}

	# Check for .pairs, .seqBin and _sim.args

	# Make .pairs for sanger reads
	if (scalar(grep {$_ eq "sanger"} @techs) > 0) {
		if (! -f "$workspace/$name.pairs") {
			#my $now = time;
			#utime $now, $now, "$workspace/$name.pairs";
			$this->execute("python $pkg_dir/xml/parseTraceInfo.py $workspace $name $read_dir/$metagenome.xml");
			rename "$workspace/${name}pairs.dat", "$workspace/${name}.pairs";
			unlink "$workspace/${name}diffs.dat";
			unlink "$workspace/${name}singletons.dat";
		}
	}

	#make seqbin file
	{
		my $dbh = $this->get_db_handle;
		my $sth = $dbh->prepare('SELECT library_id FROM library l JOIN sample s ON l.sample_id=s.sample_id WHERE s.metagenome_id=? and type=?');
		my $sth2 = $dbh->prepare('SELECT count(read_id) FROM readinfo WHERE library_id=?');
		open(SEQBIN, ">$workspace/$name.seqbin") || die "Cannot open $name.seqBin: $!";

		################################
		# Sanger/454
		# For each tech
		# Get number of reads
		# 1. Get the libraries
		# 2. Get count for each
		################################

		my $reads_so_far = 0;
		foreach my $tech (@techs) {
			my $tech_count = 0;
			$sth->execute($metagenome, $tech);
			while (my ($library_id) = $sth->fetchrow_array()) {
				$sth2->execute($library_id);
				my ($count) = $sth2->fetchrow_array();
				$tech_count += $count;
				$sth2->fetchrow_array();
			}
			for (my $i = $reads_so_far; $i < ($reads_so_far + $tech_count); $i++) {
				print SEQBIN "$i\t$tech\t0\n";
			}
			$reads_so_far = $tech_count;
		}

		$this->{TECH} = join("+", @techs);

		close(SEQBIN);
		$sth->finish();
		$sth2->finish();
	}
	if (! -f "$workspace/${name}_sim.args") {
		$this->execute("echo makeHash -e 5000000 > $workspace/${name}_sim.args");
	}

	copy("$pkg_dir/forgeG", "$HOME/.forgeG") || die "Cannot copy Forge config file: $!";
}

############################################
# Override the SUPER::assemble method
# since Forge runs on MPI
############################################

sub assemble {
	use Smash::Utils::ClusterJob;
	my $this       = shift;
	my $output_dir = $this->workspace;
	my $name       = $this->name;
	my $clust_name = $this->cluster;
	my $job        = Smash::Utils::ClusterJob->new( JOB_NAME    => "${name}_forge",
							MEM         => 8000,
							WORKING_DIR => $output_dir,
							EODIR       => $output_dir,
							CPUS        => 2*$this->cpus,
							LAM_MPI     => 1,
							CLUSTER     => $clust_name);
	my $command    = $this->get_command_line();
	my $job_id     = $job->submit_commands($command);
	print "$this assembly of $name submitted as $job_id in $clust_name.\n";
	print "Please check $output_dir/$name.{out,err} when it finishes.\n";
	print "If it looks good, then run $0 with the same options and '--finish' to finish up\n";
	return 31;
}

sub post_assembly {
	use File::Temp;
	use File::Spec::Functions qw/tmpdir/;
	my $this        = shift;
	my $workspace   = $this->workspace;
	my $contig2read = sprintf("%s/%s.contig2read.gff", $workspace, $this->name);
	my $scaf2contig = sprintf("%s/%s.scaf2contig.gff", $workspace, $this->name);
	my $output_fasta= $this->output_fasta;
	my $last_contig;
	my $CONTIG_FA;
	my $CONTIG2READ;
	my $SCAF2CONTIG;

	############################################
	# Get clean contigs
	############################################
	my $command = sprintf("cd $workspace && python %s/scripts/scaffCleanup.py --cwd %s %d 0 0 > /dev/null", $this->pkg_dir, $this->name, $this->genome_size);
	print STDERR "Generating clean contigs...";
	if ($this->execute($command) != 0) {
		die "Error calling scafCleanup: $command";
	}
	print STDERR "done\n";

	############################################
	# Combine contigs and unassembled reads.
	############################################

	open($CONTIG_FA,   ">$output_fasta") || die "Cannot open $output_fasta: $!";
	open($CONTIG2READ, ">$contig2read") || die "Cannot open $contig2read: $!";
	open($SCAF2CONTIG, ">$scaf2contig") || die "Cannot open $scaf2contig: $!";

	####
	# Remap the contig headers and write to a new fasta file
	# Write the contig2read map for contigs
	####

	print STDERR "Renaming contigs...";
	open(ASSEMBLY, "<".$this->assembly_fasta) || die "Cannot open assembly file @{[$this->assembly_fasta]}: $!";
	my $fasta = new FAlite(\*ASSEMBLY);
	while (my $entry = $fasta->nextEntry) {
		my $def = $entry->def;
		my $seq = $entry->seq;
		$def =~ s/\s.*//;
		$def =~ s/>scaffold_//;
		$last_contig = $def;
		$def = sprintf("%s.C%d", $this->name, $def);
		print $CONTIG_FA ">$def\n";
		print $CONTIG_FA $this->pretty_fasta($seq);
	}
	close(ASSEMBLY);
	print STDERR "done\n";

	print STDERR "Mapping reads to contigs...";
	$this->map_reads_to_scaffold($this->name, $CONTIG2READ, $SCAF2CONTIG);
	print STDERR "done\n";

	close($CONTIG_FA);
	close($CONTIG2READ);
	close($SCAF2CONTIG);

	$this->copy_assembly_files($output_fasta, $contig2read, $scaf2contig);
}

sub map_reads_to_scaffold {
	use FAlite;
	use FQlite;

	my $this        = shift;
	my $project     = shift;
	my $CONTIG2READ = shift;
	my $SCAF2CONTIG = shift;
	my $workspace   = $this->workspace;
	my @files       = <$workspace/$project.ml.scaf.*>;

	$, = "\t";
	foreach my $file (@files) {
		my $fasta;
		my $qual;
		my ($ML_SCAF_HANDLE, $ML_READ_HANDLE, $ML_CONS_HANDLE, $ML_CONQ_HANDLE);
		my $part = $file;
		$part =~ s/.*($project)\.ml.scaf\.//; # get just the part 
		open($ML_SCAF_HANDLE, "<$file") || die "Cannot open $file: $!";
		open($ML_READ_HANDLE, "<$workspace/$project.ml.read.$part") || die "Cannot open $project.ml.read.$part: $!";
		open($ML_CONS_HANDLE, "<$workspace/$project.ml.consensus.$part") || die "Cannot open $project.ml.consensus.$part: $!";
		open($ML_CONQ_HANDLE, "<$workspace/$project.ml.conqual.$part") || die "Cannot open $project.ml.conqual.$part: $!";
		$fasta = new FAlite($ML_CONS_HANDLE);
		$qual  = new FQlite($ML_CONQ_HANDLE);
		while(<$ML_SCAF_HANDLE>) {
			chomp;

			my $entry;
			my $sequence;
			my $contig_length;
			my @bases;
			my %Map;
			my @read_locations;
			my $qscores;
			my @quals;
			my @clean_quals;
			my @AFLINES;
			my @NONAFLINES;

			my ($scaffold, $num_reads, $num_bases, $num_pad, $num_gap, $frac_gap, $depth, $gc) = split(/\t/);
			$entry = $fasta->nextEntry();
			if ($entry->def ne ">scaffold_$scaffold") {
				die "def mismatch: >scaffold_$scaffold vs $entry->def";
			}
			$sequence = $entry->seq;
			$sequence =~ s/\-//g;
			$contig_length = length($sequence);

			$sequence = $entry->seq;
			$sequence =~ s/\-/*/g;
			@bases    = split(//, $sequence);
			if ($num_bases - $num_pad - $num_gap < 1) {
				for (my $i=0; $i<$num_reads; $i++) {
					my $line = <$ML_READ_HANDLE>;
				}
				next;
			}
			#print STDERR "Processing scaffold_$scaffold:";
			
			# Parse the reads

			for (my $i=0; $i<$num_reads; $i++) {
				my $line = <$ML_READ_HANDLE>;
				chomp($line);
				my @words = split(/\t/, $line);
				my ($def, $left, $right, $direction, $read) = ($words[2], $words[3], $words[4], $words[5], $words[6]);
				if ($words[0] != $scaffold) {
					die "scaffold mismatch: $words[0] vs $scaffold";
				}
				push(@read_locations, "$def:$left:");
				$Map{$def}{'strand'} = ($direction eq 'f')?"+":"-";
				$Map{$def}{'read'}   = $read;
			}
			# Stored "EWCK25H01CH331:149"

			#print STDERR "($#bases bases)...";
			for (my $i=0, my $count=0; $i <= $#bases; $i++) {

				# Get the reads starting at this position
			
				#if    ($i%1000 == 0) {print STDERR "$i";}
				#elsif ($i%100 == 0) { print STDERR ".";}
				my @hits = grep { $_ =~ /:($i):/} @read_locations;
				foreach my $hit (@hits) {
					my ($def, $start) = split(':', $hit);
					my $strand = $Map{$def}{'strand'};

					my $read        = $Map{$def}{'read'};
					my $read_length = length($read);
					my $line        = "";
					my $offset      = 0;

					# Write ace format with the padded read

					$read =~ s/-/*/g;

					# remove the pads

					my @read_bases  = split(//, $read);
					for (my $j = $start, my $bcount=0; $bcount < $read_length; $j++, $bcount++) {
						if ($bases[$j] eq '*') {
							$read_bases[$bcount] = '#';
						}
					}

					# if read starts at a padded position, 
					# it will automatically go to the next unpadded position
					# since $count will not be incremented until the next unpadded position

					# reassemble the read

					$read      = join('', @read_bases);
					$read      =~ s/#//g;
					$read_length = length($read);

					print $CONTIG2READ "$project.C$scaffold", "$this", "read", $count+1,  $count+$read_length, $contig_length, $strand, ".", "read \"$def\";\n";

				}
				if ($bases[$i] ne '*') {
					$count++;
				}
			}
			print $SCAF2CONTIG "$project.S$scaffold", "$this", "contig", 1,  $contig_length, $contig_length, "+", ".", "contig \"$project.C$scaffold\";\n";
			#print STDERR ".done\n";
		}
		close($ML_SCAF_HANDLE);
		close($ML_READ_HANDLE);
		close($ML_CONS_HANDLE);
		close($ML_CONQ_HANDLE);
	}
}

1;
