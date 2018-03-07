############################################
#  Generic gene prediction program handler
#  The SUPER class
############################################

package Smash::Analyses::GenePredictorWrapper;

use strict;
use warnings;

use base qw(Exporter);
our @EXPORT_OK = qw(run_predictor);

use Cwd;
use Smash::Global qw(:all);
use Smash::Analyses::GenePredictor::GeneMark;
use Smash::Analyses::GenePredictor::MetaGene;
use Smash::Utils::ClusterJob;


############################################
# Constructor
############################################

sub new {
	my $class  = shift;
	my %params = @_;

	if (!defined($params{'PREDICTOR'})) {
		$params{'PREDICTOR'} = "Generic";
	}
	if (!defined($params{'TYPE'})) {
		$params{'TYPE'} = "GenePredictor";
	}

	bless {
		%params
	}, $class;
}

sub object {shift->{OBJECT}}

############################################
# parallelization of the event
############################################

sub parallelize_on_assembly {
	use File::Path;
	use File::Basename;
	use Config;

	# figure out input and output dirs

	my $this            = shift;
	my $version         = $this->{VERSION};
	my $cluster         = $this->{CLUSTER};

	# object contained in the wrapper

	my $object          = $this->object;
	my $genepred        = $object->name || die "parallelize needs genepred id assigned ahead!";
	my $predictor       = $object->predictor;
	my $self_train      = $object->self_train || 0;
	my $assembly        = $object->assembly;
	my $assembly_id     = $object->get_id_by_name("assembly", $assembly);
	my $assembly_fa     = sprintf("%s/%s.contigs.fa", $object->assembly_dir($assembly), $assembly);
	my $output_dir      = $object->genepred_dir($genepred);
	my $train_min       = $object->train_min;
	my $maui_dir        = $object->maui_dir;
	my $PARTITION_MAX   = 10000000;

	my %ScafInfo;
	my %ScafLen;

	mkpath $output_dir, {mode => $object->file_perm};


	# parse through scaffolds in reverse order of length
	# break if scaffold bigger than $train_min or current file has 1Mb sequence

	my $dbh = $object->get_db_handle();
	{
		my $sth = $dbh->prepare("SELECT scaffold_id, c.external_id, s.length, c.length FROM scaffold s LEFT JOIN scaffold2contig sc USING (scaffold_id) LEFT JOIN contig c USING (contig_id) WHERE c.assembly_id=?");
		$sth->execute($assembly_id);
		my ($scaf_id, $contig_id, $scaf_length, $contig_length);
		$sth->bind_columns( \($scaf_id, $contig_id, $scaf_length, $contig_length) );
		while ($sth->fetchrow_arrayref()) {
			$ScafInfo{$scaf_id}{$contig_id} = $contig_length;
			$ScafLen{$scaf_id} = $scaf_length;
		}
	}
	$object->close_db_handle();

	my $pred_job_file = "$output_dir/$genepred.pred.sh";
	my $merge_job_file = "$output_dir/$genepred.merge.sh";
	open(JOB_FILE, ">$pred_job_file") || die "Cannot open $pred_job_file: $!";

	# decide on the scripts etc

	my (undef, $script_location) = fileparse($0);
	$script_location =~ s|/$||;

	# scaffolds bigger than $train_min, self-trainable

	my $file_counter = 0;
	my @scafs;

	@scafs = sort {$a <=> $b} grep {$ScafLen{$_} > $train_min} keys %ScafLen;
	for (my $i=0; $i < @scafs; $i++) {
		my $scaf = $scafs[$i];
		my @contig_list = sort {$a cmp $b} keys %{$ScafInfo{$scaf}}; # get the contigs from this scaf

		$file_counter++;

		my $part          = "scaffold_${file_counter}";
		my $local_outdir  = "$output_dir/$part";
		my $scaffold_file = "$local_outdir/$part.fasta"; # It's stupid, but we have to do this for GeneMark!
		my $list_file     = "$local_outdir/$part.list";

		mkpath $local_outdir, {mode => $object->file_perm};
		unlink($scaffold_file);
		open(LIST, ">$list_file");
		print LIST join("\n", @contig_list), "\n";
		close(LIST);
		$object->execute(sprintf("$maui_dir/filterFasta --input=%s --list=%s --output=%s", $assembly_fa, $list_file, $scaffold_file));

		print JOB_FILE "$SMASH_PERL $0 --predictor=$predictor --genepred=$genepred --label=$part --fasta=$scaffold_file --output_dir=$local_outdir";
		if ($version) {
			printf JOB_FILE " --version=%s", $version;
		}
		if ($self_train == 1) {
			print JOB_FILE " --self_train";
		}
		print JOB_FILE "\n";
	}

	# small scaffolds, no self-training
	# split them into files so that they are all ~$PARTITION_MAX long
	# make a new part when running length-sum goes over $PARTITION_MAX. 

	my $part_size   = 0;
	my @contig_list = ();
	@scafs = sort {$a <=> $b} grep {$ScafLen{$_} <= $train_min} keys %ScafLen;
	for (my $i=0; $i < @scafs; $i++) {
		my $scaf = $scafs[$i];
		my @list = sort {$a cmp $b} keys %{$ScafInfo{$scaf}}; # get the contigs from this scaf
		push(@contig_list, @list);
		$part_size += $ScafLen{$scaf};

		if ($part_size > $PARTITION_MAX || $i == $#scafs) {  # if over size limit, or the last scaf
			$file_counter++;
			my $part          = "scaffold_${file_counter}";
			my $local_outdir  = "$output_dir/$part";
			my $scaffold_file = "$local_outdir/$part.fasta"; # It's stupid, but we have to do this for GeneMark!
			my $list_file     = "$local_outdir/$part.list";

			mkpath $local_outdir, {mode => $object->file_perm};
			unlink($scaffold_file);
			open(LIST, ">$list_file");
			print LIST join("\n", @contig_list), "\n";
			close(LIST);
			$object->execute(sprintf("$maui_dir/filterFasta --input=%s --list=%s --output=%s", $assembly_fa, $list_file, $scaffold_file)); 

			print JOB_FILE "$SMASH_PERL $0 --predictor=$predictor --genepred=$genepred --label=$part --fasta=$scaffold_file --output_dir=$local_outdir";
			if ($version) {
				printf JOB_FILE " --version=%s", $version;
			}
			print JOB_FILE "\n";

			$part_size = 0;
			@contig_list = ();
		}
	}

	close(JOB_FILE);

	# Abort when no contig sequences found!

	if ($file_counter == 0) {
		warn "WARNING: No sequences found for $assembly. Please verify that contigs are loaded into the database.\n";
		unlink($pred_job_file);
		$object->abort();
	}

	# Merging files

	open(MERGE_FILE, ">$merge_job_file") || die "Cannot open $merge_job_file: $!";
	print MERGE_FILE "$SMASH_PERL $SMASH_SCRIPT_LOCATION/doGenePrediction.pl --predictor=$object --genepred=$genepred --merge --pieces=$file_counter\n";
	close(MERGE_FILE);

	if ($cluster) {
		my $job;
		my $job_id;
		my @commands;

		# Prediction job

		$job = Smash::Utils::ClusterJob->new( 
				NAME        => "pred_$genepred",
				TYPE        => $cluster,
				MEMORY      => 4000,
				);

		open(JOB, "<$pred_job_file") || die "Cannot open $pred_job_file: $!";
		@commands = <JOB>;
		close(JOB);
		map {chomp($_)} @commands;

		$job_id = $job->submit_commands(@commands);
		print "Prediction job submitted as $job_id\n";

		# Loading job

		$job->{NAME} = "merge_$genepred";
		$job->{EXTRA_ARGS} = "-hold_jid $job_id";

		open(JOB, "<$merge_job_file") || die "Cannot open $merge_job_file: $!";
		@commands = <JOB>;
		close(JOB);
		map {chomp($_)} @commands;
		$job_id = $job->submit_commands(@commands);
		print "Merging job submitted as $job_id\n";
	} else {
		print "Prediction script: $pred_job_file\n";
		print "Merging script: $merge_job_file\n";

		print "Please run $pred_job_file to generate gene prediction\n";
		print "Afterwards, please run $merge_job_file to merge genes\n";
	}
}

sub parallelize_on_file {
	use File::Path;
	use File::Basename;
	use Config;
	use FAlite;

	# figure out input and output dirs

	my $this            = shift;
	my $version         = $this->{VERSION};
	my $cluster         = $this->{CLUSTER};
	my $fasta_file      = $this->{FASTA_FILE};

	# object contained in the wrapper

	my $object          = $this->object;
	my $genepred        = $object->name || die "parallelize needs genepred id assigned ahead!";
	my $predictor       = $object->predictor;
	my $self_train      = $object->self_train || 0;
	my $assembly        = $object->assembly;
	my $assembly_id     = $object->get_id_by_name("assembly", $assembly);
	my $assembly_fa     = sprintf("%s/%s.contigs.fa", $object->assembly_dir($assembly), $assembly);
	my $output_dir      = $object->genepred_dir($genepred);
	my $train_min       = $object->train_min;
	my $maui_dir        = $object->maui_dir;

	mkpath $output_dir, {mode => $object->file_perm};

	open(FASTA, "<$fasta_file") || die "Cannot open $fasta_file: $!";
	my $fasta = new FAlite(\*FASTA);

	# parse through scaffolds 
	# break if scaffold bigger than $train_min or current file has 1Mb sequence

	my $pred_job_file = "$output_dir/$genepred.pred.sh";
	my $merge_job_file = "$output_dir/$genepred.merge.sh";
	open(JOB_FILE, ">$pred_job_file") || die "Cannot open $pred_job_file: $!";

	# scaffolds bigger than $train_min, self-trainable

	# small scaffolds, no self-training
	# split them into files so that they are all ~1Mb long
	# make a new part when running length-sum goes over 1Mb

	my $PARTITION_MAX = 10000000;
	my $part_size     = 0;
	my $file_counter  = 1;
	my $part          = "scaffold_${file_counter}";
	my $local_outdir  = "$output_dir/$part";
	my $scaffold_file = "$local_outdir/$part.fasta"; # It's stupid, but we have to do this for GeneMark!
	my @non_trainable = ($file_counter); # keeping track of the non-trainable multifasta files

	mkpath $local_outdir, {mode => $object->file_perm};
	open(NONTRAIN, ">$scaffold_file") || die "Cannot open $scaffold_file: $!";

	ENTRY:while (my $entry = $fasta->nextEntry) {
		my $def    = $entry->def;
		my $length = length($entry->seq);

		if ($length > $train_min) {
			$file_counter++;
			my $part          = "scaffold_${file_counter}";
			my $local_outdir  = "$output_dir/$part";
			my $scaffold_file = "$local_outdir/$part.fasta"; # It's stupid, but we have to do this for GeneMark!

			mkpath $local_outdir, {mode => $object->file_perm};
			open(TRAIN, ">$scaffold_file") || die "Cannot open $scaffold_file: $!";
			print TRAIN "$entry";
			close(TRAIN);

			print JOB_FILE "$SMASH_PERL $0 --predictor=$predictor --genepred=$genepred --label=$part --fasta=$scaffold_file --output_dir=$local_outdir";
			if ($version) {
				printf JOB_FILE " --version=%s", $version;
			}
			if ($self_train == 1) {
				print JOB_FILE " --self_train";
			}
			print JOB_FILE "\n";
		} else {

			print NONTRAIN "$entry";
			$part_size += $length;

			if ($part_size > $PARTITION_MAX) {  # if over size limit, or the last scaf
				close(NONTRAIN);
				$file_counter++;
				push(@non_trainable, $file_counter);

				my $part          = "scaffold_${file_counter}";
				my $local_outdir  = "$output_dir/$part";
				my $scaffold_file = "$local_outdir/$part.fasta"; # It's stupid, but we have to do this for GeneMark!

				mkpath $local_outdir, {mode => $object->file_perm};
				open(NONTRAIN, ">$scaffold_file") || die "Cannot open $scaffold_file: $!";


				$part_size = 0;
			}
		}
	}

	# Nothing written to the last open file, so get rid of it
	if ($part_size == 0) {
		pop(@non_trainable);
	}

	foreach my $item (@non_trainable) {
		my $part          = "scaffold_${item}";
		my $local_outdir  = "$output_dir/$part";
		my $scaffold_file = "$local_outdir/$part.fasta"; # It's stupid, but we have to do this for GeneMark!
		print JOB_FILE "$SMASH_PERL $0 --predictor=$predictor --genepred=$genepred --label=$part --fasta=$scaffold_file --output_dir=$local_outdir";
		if ($version) {
			printf JOB_FILE " --version=%s", $version;
		}
		print JOB_FILE "\n";
	}

	close(NONTRAIN);
	close(JOB_FILE);

	# Abort when no contig sequences found!

	if ($file_counter == 0) {
		warn "No sequences found for $assembly. Please verify that contigs are loaded into the database.\n";
		unlink($pred_job_file);
		$object->abort();
	}

	# Merging files

	open(MERGE_FILE, ">$merge_job_file") || die "Cannot open $merge_job_file: $!";
	print MERGE_FILE "$SMASH_PERL $SMASH_SCRIPT_LOCATION/doGenePrediction.pl --genepred=$genepred --predictor=$object --merge --pieces=$file_counter\n";
	close(MERGE_FILE);

	print "Prediction script: $pred_job_file\n";
	print "Merging script: $merge_job_file\n";

	if ($cluster) {
		my $job;
		my $job_id;
		my @commands;

		# Prediction job

		$job = Smash::Utils::ClusterJob->new( 
				NAME        => "pred_$genepred",
				TYPE        => $cluster,
				MEMORY      => 4000,
				);

		open(JOB, "<$pred_job_file") || die "Cannot open $pred_job_file: $!";
		@commands = <JOB>;
		close(JOB);
		map {chomp($_)} @commands;
		$job_id = $job->submit_commands(@commands);
		print "Prediction job submitted as $job_id\n";

		# Loading job

		$job->{JOB_NAME} = "merge_$genepred";
		$job->{EXTRA_ARGS} = "-hold_jid $job_id";

		open(JOB, "<$merge_job_file") || die "Cannot open $merge_job_file: $!";
		@commands = <JOB>;
		close(JOB);
		map {chomp($_)} @commands;
		$job_id = $job->submit_commands(@commands);
		print "Merging job submitted as $job_id\n";
	} else {
		print "Please run $pred_job_file to generate gene prediction\n";
		print "Afterwards, please run $merge_job_file to merge genes\n";
	}
}

sub merge_pieces {
	my $this     = shift;
	my $genepred = shift;
	my $pieces   = shift;
	my $output_dir = $this->object->genepred_dir($genepred);
	foreach my $extension ("", qw(.gene.fa .protein.fa .contig2gene.gff)) {
		my $file = "$output_dir/${genepred}${extension}";
		open(FILE, ">$file") || die "Cannot open $file: $!";
		for (my $i=1; $i <= $pieces; $i++) {
			my $part = "scaffold_$i";
			my $in   = "$output_dir/$part/${part}${extension}";
			open(IN, "<$in") || die "Cannot open $in: $!";
			while(<IN>) {
				print FILE $_;
			}
			close(IN);
		}
		close(FILE);
	}
}

############################################
# Actual execution of the object
############################################

# incoming hash args are lower case
# outgoing hash args (to the object constructors) are upper case
# if we need to set something here for the outgoing hash, use upper case
# so that we can differentiate between values set before we got here, and values set here
# e.g., the check for $options{fasta_file} will be true only if it came through %options.
# if we set it based on $options{assembly}, it would not be true

sub run_predictor {
	my %options = @_;

	# Get the wrapper instance

	my $wrapper = new Smash::Analyses::GenePredictorWrapper(map {uc($_) => $options{$_}} keys %options);

	# label is a convenient term to use commonly between doGenePrediction.pl and loadGenePrediction.pl
	# GenePredictor only understands NAME, so populate NAME with value from LABEL
	if ($options{label}) {
		$options{name} = $options{label};
		delete $options{label};
	}

	# Get the predictor instance
	# This is where the new genepred id is assigned, in the init() call

	my $class     = "Smash::Analyses::GenePredictor::".$options{predictor};
	my $predictor = "$class"->new(map {uc($_) => $options{$_}} keys %options);
	$predictor->init();
	$wrapper->{OBJECT} = $predictor;

	if ($options{merge}) {
		$wrapper->merge_pieces($options{genepred}, $options{pieces});
		$predictor->finish();
		return;
	} elsif ($options{assembly}) {
		my $assembly         = $options{assembly};
		$predictor->{FASTA_FILE} = sprintf("%s/%s.contigs.fa", $predictor->assembly_dir($assembly), $assembly);
	}

	# TO BE IMPLEMENTED:
	# If self training is requested, then sequences > train_min needs to be separated.
	# Parallelization does it automatically, so we can reuse it with a HUGE partition size.
	# This way, there will be one piece for each trainable set, and one huge piece for the rest.
	# But self_train is checked first, so that if self_train and parallelize are both selected,
	# you still get 1Mb.

	if ($options{parallelize}) {
		if ($options{fasta_file}) {
			$wrapper->parallelize_on_file();
		} elsif ($options{assembly}) {
			$wrapper->parallelize_on_assembly();
		}
	} else {
		$predictor->run();
	}
	print "<output>".$predictor->genepred."</output>\n";
	print "********************************************************\n";
	print "  Prediction id assigned to this gene set: ".$predictor->genepred."\n";
	print "********************************************************\n";
	$predictor->finish();
}

1;
