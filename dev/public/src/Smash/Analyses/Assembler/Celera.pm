############################################
#  Celera sequence assembly program handler
############################################

package Smash::Analyses::Assembler::Celera;

use strict;
use warnings;
use Smash::Global qw(:all);

our @ISA = qw(Smash::Analyses::Assembler);

=head1 NAME

Smash::Analyses::Assembler::Celera - Implementation of Celera assembly software pipeline

=head1 SYNOPSIS

=head1 DESCRIPTION

This module performs assemblies
using Celera assembler.

=head2 Default options

Here are the default options used by Smash for a metagenomic assembly:

	utgErrorRate         = 0.12
	ovlErrorRate         = 0.14
	cnsErrorRate         = 0.14
	cgwErrorRate         = 0.14
	merSize              = 14
	overlapper           = ovl
	doFragmentCorrection = 0
	doExtendClearRanges  = 1
	utgBubblePopping     = 0
	doOverlapBasedTrimming = 1
	unitigger            = bog
	merOverlapperSeedBatchSize = 100000
	utgGenomeSize        = $this->genome_size

	# The following numbers are optimized for a 1GB memory limit running on a laptop.
	# This can only run very small assemblies.
	# Must override for larger assemblies run on clusters

	ovlMemory                  = 1GB
	ovlCorrBatchSize           = 100000
	frgCorrBatchSize           = 100000


And here are the defaults used for a single genome assembly:


	utgErrorRate         = 0.03
	ovlErrorRate         = 0.06
	cnsErrorRate         = 0.06
	cgwErrorRate         = 0.10
	merSize              = 22
	overlapper           = mer
	doFragmentCorrection = 1
	doExtendClearRanges  = 2
	utgBubblePopping     = 1
	doOverlapBasedTrimming = 1
	unitigger            = bog
	merOverlapperSeedBatchSize = 100000
	utgGenomeSize        = $this->genome_size

	# The following numbers are optimized for a 1GB memory limit running on a laptop.
	# This can only run very small assemblies.
	# Must override for larger assemblies run on clusters

	ovlMemory                  = 1GB
	ovlCorrBatchSize           = 100000
	frgCorrBatchSize           = 100000

As mentioned inside the options, C<ovlMemory>, C<ovlCorrBatchSize> and C<frgCorrBatchSize>
must be changed for running large assemblies. We recommend C<ovlMemory=8GB> if you have such
a server. 

=head2 Machine specific options

=head2 Overriding default options

You can override any of the options mentioned earlier, or send other options understood by
Celera assembler using the C<extra_options> parameter as follows:

	doAssembly.pl ... --assembler=Celera \
		--extra_options="ovlMemory=8GB merSize=19 overlapper=mer"

=cut

my $PROGRESS;

############################################
############################################
##    Member variables through methods    ##
############################################
############################################

sub installation_name {'wgs'}

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

sub first_assembly_dir {
	my $this = shift;
	my $base = $this->workspace."/".$this->name;
	$base .= "/9-terminator";
	return $base;
}

sub toggle_assembly_dir {
	my $this = shift;
	my $base = $this->workspace."/".$this->name;
	$base .= "/10-toggledAsm";
	$base .= "/9-terminator";
	return $base;
}

sub results_dir {
	my $this = shift;
	return ($this->single_genome)?($this->first_assembly_dir):($this->toggle_assembly_dir);
}

=head1 FUNCTIONS

=over 4

=item C<assemble()>

Complicated. To do!

=item C<post_assembly()>

This is the global post-assembly step at the end of the assembly after all the 
assemblies are done. It generates:

	1. Contig fasta file
	2. Contig-to-read mapping file in GFF format
	3. Scaffold-to-contig mapping file in GFF format

=back

=cut

############################################
# Extra options for Celera
############################################

sub options_string_to_hash {
	my $this   = shift;
	my $string = shift;
	my $delim  = shift || " ";
	my $hash   = {};
	map {
		my ($key, $value) = split("="); 
		$hash->{$key} = $value;
	} split($delim, $string);

	return $hash;
}

sub options_hash_to_string {
	my $this  = shift;
	my $hash  = shift;
	my $delim = shift || " ";
	return join($delim, map {"$_=$hash->{$_}"} sort {$a cmp $b} keys %$hash);
}

sub merge_and_override_options {
	my $this = shift;
	my ($options, $override) = @_;
	while (my ($key, $value) = each %$override) {$options->{$key} = $value} 
}

sub options_hash {
	my $this = shift;

	# defaults from Celera itself, good for single genomes
	my $options =  {
			cleanup              => "none",
			utgErrorRate         => "0.03",
			ovlErrorRate         => "0.06",
			cnsErrorRate         => "0.06",
			cgwErrorRate         => "0.10",
			doOverlapBasedTrimming => "1",
			unitigger            => "bog",
			overlapper           => "mer",
			utgBubblePopping     => "1",
			doFragmentCorrection => "1",
			doExtendClearRanges  => "2",
			merSize              => "22",
			merOverlapperSeedBatchSize => "100000",
			utgGenomeSize        => $this->genome_size,

			# the following numbers are optimized for a 2GB memory limit

			ovlHashBits                => "23",
			ovlHashBlockLength         => "30000000",

			ovlCorrBatchSize           => "100000",
			frgCorrBatchSize           => "100000",
			};

	# merge the options specific to metagenome assembly mode

	if (!$this->single_genome) {
		my $special =  {
				utgErrorRate         => "0.12",
				ovlErrorRate         => "0.14",
				cnsErrorRate         => "0.14",
				cgwErrorRate         => "0.14",
				merSize              => "14",
				overlapper           => "ovl",
				doFragmentCorrection => "0",
				doExtendClearRanges  => "1",
				utgBubblePopping     => "0",
				};
		$this->merge_and_override_options($options, $special);
	}

	# merge the options specified in system wide config file if running under cluster
	# if nothing is specified, send empty string instead of undef
	# delimiter between options is space, delimiter between key and value is =

	if (my $cluster = $this->cluster) {
		my $data_dir    = $this->data_dir;
		my $config_dir  = $this->get_smash_conf_value("config_dir");
		my $config_file = "$data_dir/$config_dir/celera.@{[$this->version]}.$cluster.spec";
		if (-f $config_file) {
			my $parser      = new Smash::Config::ConfigParser($config_file);
			my $config      = $parser->parse("=");
			$this->merge_and_override_options($options, $config->{global});
		}
	}


	# merge the options specified in command line
	# if nothing is specified, send empty string instead of undef
	# delimiter between options is space, delimiter between key and value is =

	my $extra = $this->options_string_to_hash($this->extra_options || "", " ");
	$this->merge_and_override_options($options, $extra);

	return $options;
}

# Get a flat string of options, delimited by space

sub options {
	my $this = shift;
	$this->options_hash_to_string($this->options_hash, " ");
}

sub workspace {
	my $this = shift;
	my $workspace = $this->SUPER::workspace();
	return $workspace."/".$this->name;
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
	$PROGRESS = \*STDERR;
	select($PROGRESS); $| = 1; select(STDOUT);

	# machine specific search

	my $target = `uname`;
	chomp($target);

	# only Linux is accepted

	if ($target ne "Linux") {
		warn "WARNING: Only Linux platform is supported by Smash\n";
		$this->abort();
	}

	# now get more specific

	$target = `uname -m`;
	chomp($target);
	my $subdir;
	if ($target eq "x86_64") {
		$subdir="Linux-amd64";
	} elsif ($target eq "i686") {
		$subdir="Linux-i686";
	} else {
		die "Unknown Linux architecture found: $target";
	}
	$this->{PKG_DIR} .= "/$subdir";
}

sub prepare {
	my $this = shift;
	$this->SUPER::prepare();

	#my $name = $this->name;
	#$this->set_output_fasta(sprintf("%s/%s.contigs.fa", $this->workspace, $name));
}

############################################
# Validate if everything is set up right
# for Celera
############################################

sub validate {
	my $this = shift;

	# Check for executable

	my $pkg_dir = $this->pkg_dir;
	my $status  = $this->execute("perl $pkg_dir/bin/runCA 1>/dev/null");
	if ($status != 0) {
		warn "WARNING: Cannot execute perl $pkg_dir/bin/runCA. Please check your set-up!\n";
		$this->abort();
	}
}


############################################
# Override the SUPER::assemble method
# since Celera runs on cluster
############################################

sub rename_fasta_entries {
	my $this = shift;
	my ($infile, $outfile, $Remap) = @_;
	open(FASTA, "<$infile") || die "Cannot open $infile: $!";
	open(OUT, ">$outfile") || die "Cannot open $outfile: $!";
	my $fasta = new FAlite(\*FASTA);
	while (my $entry = $fasta->nextEntry) {
		my $def = substr($entry->def, 1);
		$def = $Remap->{$def};
		print OUT ">$def\n";
		print OUT $this->pretty_fasta($entry->seq);
	}
	close(OUT);
	close(FASTA);
}

sub rename_qual_entries {
	my $this = shift;
	my ($infile, $outfile, $Remap) = @_;
	open(FASTA, "<$infile") || die "Cannot open $infile: $!";
	open(OUT, ">$outfile") || die "Cannot open $outfile: $!";
	my $qual = new FQlite(\*FASTA);
	while (my $entry = $qual->nextEntry) {
		my $def = substr($entry->def, 1);
		$def = $Remap->{$def};
		print OUT ">$def\n";
		print OUT $this->pretty_qual($entry->seq);
	}
	close(OUT);
	close(FASTA);
}

sub assemble {
	use Env qw[$HOME];
	use File::Path;
	use File::Copy;
	use Smash::Utils::ClusterJob;

	my $this       = shift;
	my $name       = $this->name;
	my $metagenome = $this->metagenome;
	my $read_dir   = $this->read_dir($metagenome);
	my @fasta      = $this->fasta_files($metagenome);
	my @qual       = $this->qual_files($metagenome);
	my @xml        = $this->xml_files($metagenome);
	my $pkg_dir    = $this->pkg_dir;
	my $workspace  = $this->workspace;
	my $maui_dir   = $this->maui_dir;
	my @frg_files;

	# Get Sanger and 454 sequences separately

	mkpath $workspace, {mode => $this->file_perm};

	# Sanger

	print $PROGRESS "Retrieving Sanger read information ...";
	my @xml_files = ();
	my $dbh = $this->get_db_handle;
	{
		my $sth_samples_for_tech = $dbh->prepare('SELECT DISTINCT sample_id FROM library l INNER JOIN sample s USING (sample_id) WHERE s.metagenome_id=? AND type=?');
		my $sth_reads_for_sample = $dbh->prepare('SELECT read_id, defline FROM readinfo INNER JOIN library USING (library_id) WHERE sample_id=?');
		$sth_samples_for_tech->execute($metagenome, "sanger");
		while (my ($sample_id) = $sth_samples_for_tech->fetchrow_array()) {
			my $frg_file = "$metagenome.$sample_id.sanger.frg";
			if (! -f "$read_dir/$frg_file") { # if frg file doesnt exist, add that xml file in the "to be processed" list
				my $fasta_file = "$metagenome.$sample_id.sanger.fasta";
				my $qual_file  = "$metagenome.$sample_id.sanger.qual";
				my $xml_file   = "$metagenome.$sample_id.sanger.xml";
				if (! -f "$read_dir/$fasta_file" && ! -f "$read_dir/$qual_file") { # if fasta/qual files are not separated by sample id's already, separate them yourself
					die "Missing $read_dir/$fasta_file and $read_dir/$qual_file. Celera needs fasta/qual/xml files separated by samples.\n";
				} else {
					my ($read_id, $defline);
					my $Def = {};
					$sth_reads_for_sample->execute($sample_id);
					$sth_reads_for_sample->bind_columns(\$read_id, \$defline);
					while ($sth_reads_for_sample->fetch()) {
						$Def->{$read_id} = $defline;
					}
					$this->rename_fasta_entries("$read_dir/$fasta_file", "$workspace/fasta.sanger.$sample_id", $Def);
					$this->rename_qual_entries ("$read_dir/$qual_file",  "$workspace/qual.sanger.$sample_id", $Def);
					symlink "$read_dir/$xml_file",   "$workspace/xml.sanger.$sample_id";
					push(@xml_files, "xml.sanger.$sample_id");
					push(@frg_files, "sanger.2.$sample_id.frg.bz2");
				}
			} else {
				symlink "$read_dir/$frg_file", "$workspace/$frg_file";
				push(@frg_files, $frg_file);
			}
		}
	}
	$this->close_db_handle();
	print $PROGRESS " done\n";

	# make frg files from the "to be processed" xml file list

	if (@xml_files) {
		print $PROGRESS "Generating FRG files for Celera ...\n";
		my $status;
		foreach my $xml (@xml_files) {
			$status = $this->execute("cd $workspace && perl $pkg_dir/bin/tracedb-to-frg.pl -xml $xml");
			if ($status != 0) {
				warn "WARNING: tracedb-to-frg.pl -xml failed!";
				$this->abort();
			}
		}
		$status = $this->execute("cd $workspace && perl $pkg_dir/bin/tracedb-to-frg.pl -lib ".join(" ", @xml_files));
		if ($status != 0) {
			warn "WARNING: tracedb-to-frg.pl -xml failed!";
			$this->abort();
		}
		foreach my $xml (@xml_files) {
			$status = $this->execute("cd $workspace && perl $pkg_dir/bin/tracedb-to-frg.pl -frg $xml");
			if ($status != 0) {
				warn "WARNING: tracedb-to-frg.pl -xml failed!";
				$this->abort();
			}
		}

		# the steps above make these frg files as well

		@frg_files = ("sanger.1.lib.frg", @frg_files, "sanger.3.lkg.frg");

		print $PROGRESS " done\n";

		# remove the fasta and qual files you just made.
		# remember that @xml_files doesnt have full path

		print $PROGRESS "Deleting intermediate files ...";
		foreach my $file (@xml_files) {
			my $remove = $file;
			$remove =~ s/^xml\.//;
			unlink "$workspace/fasta.$remove";
			unlink "$workspace/qual.$remove";
		}
		print $PROGRESS "done\n";

		# zip the file

		map {my $f = $_; $f =~ s/\.bz2$//; print $PROGRESS "bzipping $f ..."; $this->execute("cd $workspace && bzip2 -9 $f") if (-f "$workspace/$f"); print $PROGRESS " done\n";} grep {m/^sanger\.2\..*\.frg\.bz2$/} @frg_files;

	}

	# 454 Titanium

	print $PROGRESS "Retrieving 454 read information ...";
	$dbh = $this->get_db_handle;
	{
		my $sth_libs_for_tech    = $dbh->prepare('SELECT sample_id, library_id, external_id FROM library l INNER JOIN sample s USING (sample_id) WHERE s.metagenome_id=? AND type=?');
		my $sth_reads_for_lib    = $dbh->prepare('SELECT read_id, defline FROM readinfo WHERE library_id=?');
		$sth_libs_for_tech->execute($metagenome, "454");
		while (my ($sample_id, $library_id, $library_name) = $sth_libs_for_tech->fetchrow_array()) {
			my $frg_file = "$metagenome.$sample_id.454.$library_id.frg";
			my ($read_id, $defline);
			my $LIST = new File::Temp();

			# Get readnames for this library

			my $Def = {};
			$sth_reads_for_lib->execute($library_id);
			$sth_reads_for_lib->bind_columns(\$read_id, \$defline);
			while ($sth_reads_for_lib->fetch()) {
				print $LIST "$read_id\n";
				$Def->{$read_id} = $defline;
			}

			if (! -f "$read_dir/$frg_file") { # if frg file doesnt exist, make them
				my $fasta_file = "$metagenome.$sample_id.454.$library_id.fasta";
				my $qual_file  = "$metagenome.$sample_id.454.$library_id.qual";
				if (! -f "$read_dir/$fasta_file" && ! -f "$read_dir/$qual_file") { # if fasta/qual files are not separated by library id's already, separate them yourself
					my $tmp1 = new File::Temp();
					my $tmp2 = new File::Temp();

					# Get the reads for this library
					# and rename the headers

					$this->execute("cat @{[join(' ', @fasta)]} > $tmp1");
					$this->execute("$maui_dir/filterFasta --input=$tmp1 --list=$LIST --output=$tmp2");
					$this->rename_fasta_entries("$tmp2", "$workspace/$fasta_file", $Def);

					# Get the quals for this library
					# and rename the headers

					$this->execute("cat @{[join(' ', @qual)]}  > $tmp1");
					$this->execute("$maui_dir/filterQual  $tmp1  $LIST $tmp2");
					$this->rename_qual_entries ("$tmp2", "$workspace/$qual_file",  $Def);
				} else {
					symlink "$read_dir/$fasta_file", "$workspace/$fasta_file";
					symlink "$read_dir/$qual_file",  "$workspace/$qual_file";
				}
				$this->execute("perl $pkg_dir/bin/convert-fasta-to-v2.pl -l $library_name -454 -s $workspace/$fasta_file -q $workspace/$qual_file | bzip2 -9c > $workspace/$frg_file.bz2");
				push(@frg_files, "$frg_file.bz2");
				unlink("$workspace/$fasta_file", "$workspace/$qual_file");
			} else {

				# if frg file exists, then it doesn't have the MC1.MG1.readname format, since it was straight from sffToCA.
				# so we dont need to rename the entries as above or in the case of Sanger reads.
				# However, the read id's in the read fasta file are in MC1.MG1.readname format and so is the read_id in the
				# readinfo table. So there is no problem in remapping the names in creating the contig2read entries in
				# post_assembly().

				symlink "$read_dir/$frg_file", "$workspace/$frg_file";
				push(@frg_files, $frg_file);
			}
		}
	}
	$this->close_db_handle();
	print $PROGRESS " done\n";

	# round 1

	my $options_hash = $this->options_hash;
	my $options      = $this->options_hash_to_string($options_hash, "\n");
	my $spec_file    = "$workspace/$name.spec";
	open(SPEC, ">$spec_file") || die "Cannot open $spec_file: $!";
	print SPEC "$options\n";
	close(SPEC);

	my ($r1_command, $r2_command, $cleanup_command);
	$r1_command  = "cd $workspace && perl $pkg_dir/bin/runCA -s $spec_file -p $name -d $name ";
	$r1_command .= join(" ", @frg_files);

	# metagenome specific:
	# If the first round finished successfully, then try toggling

	my $toggle = {
				   doToggle => "1",
			doExtendClearRanges => "0",
			 toggleUnitigLength => "250",
			 toggleNumInstances => "0"
		     };
	$this->merge_and_override_options($options_hash, $toggle);
	$options   = $this->options_hash_to_string($options_hash, "\n");
	$spec_file = "$workspace/$name.toggle.spec";
	open(SPEC, ">$spec_file") || die "Cannot open $spec_file: $!";
	print SPEC "$options\n";
	close(SPEC);

	$r2_command  = "cd $workspace && perl $pkg_dir/bin/runCA -s $spec_file -p $name -d $name ";
	$r2_command .= join(" ", @frg_files);

	# prepare the cleanup that runs after assembly is done

	my $cleanup = {cleanup => "aggressive"};
	$this->merge_and_override_options($options_hash, $cleanup);
	$options   = $this->options_hash_to_string($options_hash, "\n");
	$spec_file    = "$workspace/$name.cleanup.spec";
	open(SPEC, ">$spec_file") || die "Cannot open $spec_file: $!";
	print SPEC "$options\n";
	close(SPEC);
	$cleanup_command  = "cd $workspace && perl $pkg_dir/bin/runCA -s $spec_file -p $name -d $name ";
	$cleanup_command .= join(" ", @frg_files);

	# On cluster or not?

	my $clust_name = $this->cluster;
	my $cpus       = $this->{CPUS} || 4;
	if ($clust_name) {
		my $job        = Smash::Utils::ClusterJob->new( NAME        => "r1_${name}",
								MEM         => 8000,
								WORKING_DIR => $workspace,
								EODIR       => $workspace,
								CPUS        => $cpus,
								TYPE        => "SGE",
								CLUSTER     => $clust_name);
		my ($job_id, undef) = $job->submit_commands($r1_command);
		return -1 if ($job_id == -1);
		if (!$this->single_genome) { # need to submit the toggle job
			my $tjob = Smash::Utils::ClusterJob->new(       NAME        => "r2_${name}",
									MEM         => 8000,
									WORKING_DIR => $workspace,
									EODIR       => $workspace,
									CPUS       => $cpus,
									TYPE        => "SGE",
									EXTRA_ARGS  => "-h",
									CLUSTER     => $clust_name);
			my ($tjob_id, undef) = $tjob->submit_commands($r2_command);
print <<EOF;
Step 1 of $this assembly of $name submitted as job $job_id in cluster $clust_name.
Step 2 of $this assembly of $name submitted as $tjob_id with a 'hold'.
1. Wait till all the Step 1 jobs finish. Usually these contain $name in the job-name.
2. Check @{[$this->first_assembly_dir]}/$name.qc.
3. If it looks good, then release the hold on Step 2 using
	qrls $tjob_id
4. Wait till all the Step 2 jobs finish. These also contain $name in the job-name.
5. Check @{[$this->toggle_assembly_dir]}/$name.qc.
6. If it looks good, then run $SMASH_SCRIPT_NAME as follows to finish up
	$SMASH_PERL $0 --assembler=$this --version=@{[$this->version]} --finish --assembly=$name
EOF
		} else {
			print "$this assembly of $name submitted as $job_id in $clust_name.\n";
			print "Please check @{[$this->first_assembly_dir]}/$name.qc when it finishes.\n";
			print "If it looks good, then run $SMASH_SCRIPT_NAME with the same options and '--finish' to finish up\n";
		}
		return 31;
	} else {
		# We are running in non-cluster mode. We need to run two steps and track them properly.

		my $r1_success  = 0;
		my $r2_required = 1; # default is required. single genome will uncheck this.
		my $r2_success  = 0;

		if ($this->single_genome) {
			$r2_required = 0;
		}

		print $PROGRESS "<Assembly command=\"$r1_command\">\n";
		$r1_success = ($this->execute($r1_command) == 0);

		# Round 1 finished. Was it good?

		if ($r1_success == 0) {
			return 0;
		} else {

			print $PROGRESS "</Assembly>\n";

			# Do I have to do the toggle step? It's a good idea to do it for metagenomes
			# if the degenerate reads are more than say 5%. Otherwise, there wont be any
			# new contigs found, and the toggle step will fail!
			# The info we need is in the qc file in the [Reads] section.

			# [Reads]
			# TotalReadsInput=NA
			# TotalUsableReads=534423
			# AvgClearRange=324
			# ContigReads=79482(14.87%)
			# BigContigReads=9923(1.86%)
			# SmallContigReads=69559(13.02%)
			# DegenContigReads=144442(27.03%)
			# SurrogateReads=0(0.00%)
			# PlacedSurrogateReads=0(0.00%)
			# SingletonReads=310499(58.10%)
			# ChaffReads=310270(58.06%)

			my $qc_file = $this->first_assembly_dir."/$name.qc";
			my $parser  = new Smash::Config::ConfigParser($qc_file);
			my $config  = $parser->parse("=");
			if ($config->{Reads}) {
				my $degen = $config->{Reads}->{DegenContigReads};
				my $surr  = $config->{Reads}->{SurrogateReads};
				if ($degen =~ /\d+\(([\d\.]+)%\)/ ) {
					$degen = $1;
				} else {
					$degen = 0;
				}
				if ($surr =~ /\d+\(([\d\.]+)%\)/ ) {
					$surr = $1;
				} else {
					$surr = 0;
				}
				if ($degen + $surr < 1) {
					printf $PROGRESS "NOTE: Toggled assembly, if requested, will be skipped since \n      degenerate/surrogate contigs account for %d%% of reads\n", $degen+$surr;
					$r2_required = 0;
				}
			}

			# Now $r2_required is updated according to the r1 assembly status.

			if ($r2_required) {
				print $PROGRESS "<Assembly command=\"$r2_command\">\n";
				$r2_success = ($this->execute($r2_command) == 0);
				if ($r2_success == 1) {
					print $PROGRESS "</Assembly>\n";
				} else {
					return 0;
				}
			}

			# Smash cleaner has arrived!
			# Cleaner should onl cleanup if r1 was successful. 
			# If r2 was required and not successful, then you can still cleanup since r2 is not guaranteed anyway!

			print $PROGRESS "Smash cleaner has arrived!\n";
			print $PROGRESS "<Assembly command=\"$cleanup_command\">\n";
			$this->execute($cleanup_command);
			print $PROGRESS "</Assembly>\n";
			return 1;
		}
	}
}

sub post_assembly {
	use File::Temp;
	use File::Spec::Functions qw/tmpdir/;
	my $this        = shift;
	my $metagenome  = $this->metagenome;
	my $name        = $this->name;
	my $workspace   = $this->workspace;
	my $results_dir = $this->results_dir;
	my $maui_dir    = $this->maui_dir;

	my $contig_fasta= sprintf("%s/%s.contigs.fa",      $workspace, $name);
	my $contig2read = sprintf("%s/%s.contig2read.gff", $workspace, $name);
	my $scaf_fasta  = sprintf("%s/%s.scaffolds.fa",    $workspace, $name);
	my $scaf2contig = sprintf("%s/%s.scaf2contig.gff", $workspace, $name);
	my $tmp         = new File::Temp();
	my $contig_reads= "$tmp.cont";
	my $CONTIG_FA;
	my $CONTIG2READ;
	my $SCAF_FA;
	my $SCAF2CONTIG;

	# Did toggle really happen?

	if ($results_dir eq $this->toggle_assembly_dir) {
		if (! -d $results_dir) {
			print $PROGRESS "WARNING: $results_dir does not exist.\n         Toggled assembly did not complete successfully!\n";
			$results_dir = $this->first_assembly_dir;
			print $PROGRESS "         Checking original assembly results in $results_dir\n";
		}
	}
	if (! -d $results_dir) {
		print $PROGRESS "ERROR: $results_dir does not exist. Assembly did not complete successfully!\n";
		$this->abort();
	}

	############################################
	# Get contig lengths and min contig id
	############################################

	print $PROGRESS "Finding contig/scaffold lengths and min id ...";
	my %Min;
	my %Lengths;
	foreach my $type qw(ctg scf) {
		my $min = 10**20;
		open(FASTA, "<$results_dir/$name.$type.fasta");
		my $fasta = new FAlite(\*FASTA);
		while (my $entry = $fasta->nextEntry) {
			my $l = length($entry->seq);
			my $d = $entry->def;
			$d =~ s/^>${type}//;
			$Lengths{$d} = $l;
			if ($d < $min) {
				$min = $d;
			}
		}
		close(FASTA);
		$Min{$type} = $min-1; # so that they start at 1
	}
	print $PROGRESS " done\n";

	print $PROGRESS "Retrieving read name maps ...";

	my $ReadId2Def = {};
	my $Def2ReadId = {};
	my $dbh        = $this->get_db_handle;
	{
		my $sth_reads_for_metagenome = $dbh->prepare('SELECT read_id, defline FROM readinfo INNER JOIN library USING (library_id) INNER JOIN sample USING (sample_id) WHERE metagenome_id=?');
		my ($read_id, $defline);
		$sth_reads_for_metagenome->execute($metagenome);
		$sth_reads_for_metagenome->bind_columns(\$read_id, \$defline);
		while ($sth_reads_for_metagenome->fetch()) {
			if ($defline =~ m/ti\|(\S+)\s/) { # this is exactly what tracedb-to-frg.pl does!
				$defline = $1;
			}
			$Def2ReadId->{$defline} = $read_id;
			$ReadId2Def->{$read_id} = $defline;
		}
	}
	$this->close_db_handle();
	print $PROGRESS " done\n";

	############################################
	# Open file handles
	############################################

	open($CONTIG_FA,   ">$contig_fasta") || die "Cannot open $contig_fasta: $!";
	open($SCAF_FA,     ">$scaf_fasta")   || die "Cannot open $scaf_fasta: $!";
	open($CONTIG2READ, ">$contig2read")  || die "Cannot open $contig2read: $!";
	open($SCAF2CONTIG, ">$scaf2contig")  || die "Cannot open $scaf2contig: $!";
	open(ASM_READS,    ">$contig_reads") || die "Cannot open $contig_reads: $!";

	print $PROGRESS "Parsing position maps:\n";

	############################################
	# Get contig posmaps
	############################################

	print $PROGRESS "\tcontig2read ...";
	open(IN_MAP, "<$results_dir/$name.posmap.frgctg");
	my $contig_id_offset = $Min{"ctg"};
	while (<IN_MAP>) {
		chomp();
		my ($frg, $ctg, $start, $end, $dir) = split(/\t/);

		# If we used SFF to import, then the SFF file contains original names in defline field and the read_id wont match the frg in assembly.
		# If we imported fasta files, then read_id will match the frg.
		# Since we dont know how the reads were imported, we just check if either defline or read_id matches the frg.

		my $read = $Def2ReadId->{$frg};
		if (!$read) {
			$read = $frg if $ReadId2Def->{$frg};
		}
		die "Unknown read frg:$frg in assembly:$name" unless $read;

		my $length = $Lengths{$ctg};
		my $strand = ($dir eq "f")?"+":"-";
		$ctg =~ s/^ctg//;
		$ctg -= $contig_id_offset;
		$start++;
		print $CONTIG2READ "$name.C$ctg	Celera	read	$start	$end	$length	$strand	.	read \"$read\";\n";
		print ASM_READS "$read\n";
	}
	close(IN_MAP);
	print $PROGRESS " done\n";

	############################################
	# Get scaffold posmaps
	############################################

	print $PROGRESS "\tscaf2contig ...";
	open(IN_MAP, "<$results_dir/$name.posmap.ctgscf");
	my $scaffold_id_offset = $Min{"scf"};
	while (<IN_MAP>) {
		chomp();
		my ($ctg, $scf, $start, $end, $dir) = split(/\t/);
		my $length = $Lengths{$scf};
		my $strand = ($dir eq "f")?"+":"-";
		$ctg =~ s/^ctg//;
		$ctg -= $contig_id_offset;
		$scf =~ s/^scf//;
		$scf -= $scaffold_id_offset;
		$start++;
		print $SCAF2CONTIG "$name.S$scf	Celera	contig	$start	$end	$length	$strand	.	contig \"$name.C$ctg\";\n";
	}
	close(IN_MAP);
	print $PROGRESS " done\n";

	############################################
	# Combine contigs and unassembled reads.
	############################################

	####
	# Remap the contig headers and write to a new fasta file
	####

	print $PROGRESS "Renaming contigs ...";
	open(ASSEMBLY, "<$results_dir/$name.ctg.fasta") || die "Cannot open assembly file $results_dir/$name.ctg.fasta: $!";
	my $fasta = new FAlite(\*ASSEMBLY);
	while (my $entry = $fasta->nextEntry) {
		my $def = $entry->def;
		my $seq = $entry->seq;
		$def =~ s/^>ctg//;
		$def = sprintf("%s.C%d", $name, $def-$contig_id_offset);
		print $CONTIG_FA ">$def\n";
		print $CONTIG_FA $this->pretty_fasta($seq);
	}
	close(ASSEMBLY);
	print $PROGRESS " done\n";

	####
	# Remap the scaffold headers and write to a new fasta file
	####

	print $PROGRESS "Renaming scaffolds ...";
	open(ASSEMBLY, "<$results_dir/$name.scf.fasta") || die "Cannot open assembly file $results_dir/$name.scf.fasta: $!";
	$fasta = new FAlite(\*ASSEMBLY);
	while (my $entry = $fasta->nextEntry) {
		my $def = $entry->def;
		my $seq = $entry->seq;
		$def =~ s/^>scf//;
		$def = sprintf("%s.S%d", $name, $def-$scaffold_id_offset);
		print $SCAF_FA ">$def\n";
		print $SCAF_FA $this->pretty_fasta($seq);
	}
	close(ASSEMBLY);
	print $PROGRESS " done\n";

	####
	# Read the singleton fasta file and write to the contig fasta file
	####

	print $PROGRESS "Parsing singletons ...";
	my $singleton_fasta = "$tmp.sngl";
	my $reads           = "$tmp.read";
	$this->execute("cat @{[join(' ', $this->fasta_files($metagenome))]} > $reads");
	$this->execute("$maui_dir/filterFasta --exclude --input=$reads --list=$contig_reads --output=$singleton_fasta");
	open(SINGLETON, "<$singleton_fasta") || die "Cannot open singleton file $singleton_fasta: $!";
	$fasta = new FAlite(\*SINGLETON);
	my $singleton = 0;
	while (my $entry = $fasta->nextEntry) {
		$singleton++;
		my $def    = $entry->def;
		my $seq    = $entry->seq;
		my $length = length($seq);
		my $ctg    = sprintf("%s.R%d", $name, $singleton);
		$def =~ s/\s.*//;
		$def =~ s/^>//;
		print $CONTIG_FA ">$ctg\n";
		print $CONTIG_FA $this->pretty_fasta($seq);
		print $SCAF_FA   ">$ctg\n";
		print $SCAF_FA   $this->pretty_fasta($seq);
		print $CONTIG2READ "$ctg	Celera	read	1	$length	$length	+	.	read \"$def\";\n";
		print $SCAF2CONTIG "$ctg	Celera	contig	1	$length	$length	+	.	contig \"$ctg\";\n";
	}
	close(SINGLETON);
	unlink $singleton_fasta;
	unlink $contig_reads;
	unlink $reads;
	print $PROGRESS " done\n";

	close($CONTIG_FA);
	close($SCAF_FA);
	close($CONTIG2READ);
	close($SCAF2CONTIG);

	print $PROGRESS "Copying assembly files ...";
	$this->copy_assembly_files($contig_fasta, $contig2read, $scaf_fasta, $scaf2contig);
	print $PROGRESS " done\n";
}

1;
