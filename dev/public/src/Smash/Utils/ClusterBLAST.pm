package Smash::Utils::ClusterBLAST;
use strict;
use warnings;
use Smash::Global qw($SMASH_SCRIPT_LOCATION);
use Smash::Core;
use Smash::Utils::BLAST;
use base("Smash::Utils");

sub run_blast_cluster {
	use File::Basename;
	use File::Path;
	use File::Spec;
	use File::Glob;

	my $this          = shift;
	my %params        = @_;
	my $MyNameIs      = "ClusterBLAST";

	# From arguments

	map {
		die "$_ undefined in $MyNameIs" unless defined($params{$_});
	} qw(FLAVOR BLAST DATABASE QUERY LABEL OUTDIR EVALUE TABULAR CPUS SUBJECTS LOCALDISK);

	my $label         = $params{'LABEL'};
	my $pieces        = $params{'PIECES'};
	my $query_size    = $params{'SPLITSIZE'};
	my $bit_threshold = $params{'BITSCORE'};
	my $post_processor = $params{'POSTPROCESSOR'} || "";
	my $pre_processor = $params{'PREPROCESSOR'} || "";

	# Write the variable definitions for this run

	my $blast_obj                = new Smash::Utils::BLAST(%params);
	my $query                    = $blast_obj->query;
	my $remote_db_location       = (fileparse(File::Spec->rel2abs($blast_obj->database)))[1];
	my $database                 = (fileparse($blast_obj->database))[0];
	my $local_db_location        = $params{'LOCALDISK'};

	my $output_dir               = $params{'OUTDIR'};
	if (!-d $output_dir) {
		die "Directory $output_dir does not exist";
	}

	my $remote_output_location   = File::Spec->rel2abs($output_dir);
	my $local_output_location    = "/tmp/".(getpwuid($<)||"smash")."/out";

	   $blast_obj->{DATABASE}    = "$local_db_location/$database";
	my $extensions               = $blast_obj->get_blast_database_extension();
	my $flavor                   = $blast_obj->flavor;
	my $blast                    = $blast_obj->blast;
	my $blast_exe                = $blast_obj->get_blast_exe();
	my $blast_database_arg       = $blast_obj->get_database_arg();
	my $blast_query_arg          = $blast_obj->get_query_arg();
	my $blast_parameters         = $blast_obj->get_blast_parameters();
	my $cpus                     = $blast_obj->cpus;
	my $cpu_arg                  = $blast_obj->get_cpu_arg();

	# Locally decided

	my $shell_script             = "$label.$database.$blast.defs.sh";

	my (undef, $script_location) = fileparse($0);
	my $cleanup_flag  = 0; # Just start the jobs and quit. Can combine later. Leave the temp folder alone
	my $query_name    = fileparse($query, qr/\.[^.]*/);

	my @files;
	my $file_count;
	my @commands;

	my $MAX_CONCURRENT_NFS_COPIES = 10;
	my $nfs_lock_counter          = int(rand($MAX_CONCURRENT_NFS_COPIES)); # 0 <= $nfs_lock_counter < $MAX_CONCURRENT_NFS_COPIES

	# Split query file into smaller files for cluster
	if ($pieces > 1) {
		my $smash         = new Smash::Core();
		   $smash->init();
		my $maui_dir      = $smash->maui_dir;
		my $split_suffix  = "fa";
		my $split_dir     = "$output_dir/$label.$database.split";
		if (-d $split_dir) {
			die "ERROR: Directory $split_dir (to place split queries) already exists. \n       Please remove it, or change the label of the current BLAST\n";
		} else {
			mkpath ($split_dir, mode => 0770);
		}
		$this->execute("$maui_dir/splitFasta --mode filesize --input $query --prefix $split_dir/ --size $query_size") == 0 || die "Cannot run /$maui_dir/splitFasta successfully";
		opendir(DIR, "$split_dir") || die "Cannot open $split_dir: $!";
		@files = map { "$split_dir/$_" } grep { -f "$split_dir/$_" } readdir(DIR);
		close(DIR);
		$file_count = scalar(@files);

		my ($split_subdir) = fileparse($split_dir);
		$local_output_location .= "/$split_subdir";
		$remote_output_location = File::Spec->rel2abs($split_dir);

		####
		# JOBFILE has one command per line to run, preferably by a batch queue submission system.
		# It invokes the shell script generated earlier with each file.
		####

		####
		# Generate BLAST commands
		####

		my $random_wait_span = 60;
		if ($file_count > $random_wait_span) {
			$random_wait_span = $file_count;
		}
		my $queue_file = "$output_dir/$label.$database.$blast.queue";
		open(JOBFILE, ">$queue_file") || die "Cannot open $queue_file: $!";
		foreach my $input_file (@files) {
			my ($filename)  = fileparse($input_file, (".$split_suffix"));
			$blast_obj->{QUERY} = "\$remote_output_location/\$query.$split_suffix";
			$blast_query_arg = $blast_obj->get_query_arg();
			my $random_wait_seconds = int(rand($random_wait_span));

			my $script_dir = File::Spec->abs2rel($output_dir, $output_dir); # to help copying them to $HOME in clusters
			my $command  = "query='$filename'";
			   $command .= " bash $SMASH_SCRIPT_LOCATION/run_blast.sh $script_dir/$shell_script $nfs_lock_counter $random_wait_seconds";

			print JOBFILE "$command\n";
			$nfs_lock_counter++;
			$nfs_lock_counter %= $MAX_CONCURRENT_NFS_COPIES;
		}
		close(JOBFILE);

		####
		# Summarize
		####
		my $sge_extra_args = "";
		my $pbs_extra_args = "";
		if ($cpus > 1) {
			$sge_extra_args = "--cpus $cpus";
			$pbs_extra_args = "--cpus $cpus";
		}
		print <<EOS;
PROGRESS:
---------
	Split $query into $file_count files in $split_dir
	Wrote definitions to $shell_script
	Wrote commands to $queue_file

SUCCESS

NEXT:
-----
	Submit $queue_file to the cluster, e.g., 
	cd $output_dir && submitJobs.pl --type PBS --name $label --memory 4000 $pbs_extra_args @{[basename($queue_file)]}
	cd $output_dir && submitJobs.pl --type SGE --name $label --memory 4000 $sge_extra_args @{[basename($queue_file)]}

FINAL:
------
	cat $split_dir/*.$blast.filtered > $output_dir/$label.$database.$blast
	cat $split_dir/*.$blast > $output_dir/$label.$database.$blast.full
IN CASE OF ERROR:
-----------------
	If any individual piece needs to be reBLASTed, run the following command:
		query='007' bash $SMASH_SCRIPT_LOCATION/run_blast.sh $output_dir/$shell_script $nfs_lock_counter 0
EOS

	} else {
		$blast_obj->{QUERY} = File::Spec->rel2abs($query);
		$blast_query_arg = $blast_obj->get_query_arg();
		my $random_wait_seconds = int(rand(60));
		my $script_dir = File::Spec->abs2rel($output_dir); # to help copying them to $HOME in clusters
		my $command .= "query='$query_name'";
		   $command .= " bash $SMASH_SCRIPT_LOCATION/run_blast.sh $script_dir/$shell_script $nfs_lock_counter $random_wait_seconds";

		####
		# Summarize
		####
		print <<EOS;
PROGRESS:
---------
	Wrote definitions to $shell_script

SUCCESS

NEXT:
-----
	Run it as follows:
$command
EOS
	}

	open(SCRIPT, ">$output_dir/$shell_script") || die "Cannot open $output_dir/$shell_script: $!";

	if ($pieces == 1) {
		print SCRIPT "label=$label\n";
	}

	print SCRIPT <<EOF;
remote_db_location=$remote_db_location;
remote_output_location=$remote_output_location;
local_db_location=$local_db_location;
local_output_location=$local_output_location;
database=$database;
extensions=$extensions;
flavor=$flavor;
blast=$blast;
blast_exe='$blast_exe';
blast_database_arg='$blast_database_arg';
blast_query_arg='$blast_query_arg';
blast_parameters='$blast_parameters';
cpu_arg='$cpu_arg';
pre_processor='$pre_processor';
post_processor='$post_processor';
EOF
	print SCRIPT "bit_threshold=$bit_threshold\n" if $bit_threshold;
	close(SCRIPT);
}

1;
