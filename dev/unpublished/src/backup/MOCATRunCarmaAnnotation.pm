package MOCATRunCarmaAnnotation;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

sub create_job {
	### DEFINE VARIABLES AND OPEN OUTPUT FILE ###
	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	open JOB, '>', "$cwd/MOCATJob_$job\_$date" or die "ERROR & EXIT: Cannot write $cwd/MOCATJob_$job\_$date. Do you have permission to write in $cwd?\n";
	
	
	unless ( -e $carma_file_to_annotate ) {
		die "ERROR & EXIT: Missing input file $carma_file_to_annotate\n";
	}
	unless ( $run_carma_annotate_host1 ) {
		die "ERROR & EXIT: Missing -host1 \n";
	}
	unless ( $run_carma_annotate_host2 ) {
		die "ERROR & EXIT: Missing -host2 \n";
	}
	
	if ($carma_file_to_annotate =~ m/\//) {
		unless ($carma_file_to_annotate =~ m/^\//) {
			die "ERROR & EXIT: Please specify the carma file with a full path.\n";
		}
	} else {
		die "ERROR & EXIT: Please specify the carma file with a full path.\n";
	}

	chomp( my $user     = `whoami` );
	chomp( my $basename = `basename $carma_file_to_annotate` );
	my $pathname = "annotate.$basename";

	# Varibles
	my $carma_config        = $conf{rca_carma_cfg};
	my $smash_parallelblast = $conf{rca_parallel_blast};
	my $smash_submit        = $conf{rca_submit_jobs};
	my $smash_config        = $conf{rca_smash_config};
	my $NR_path             = $conf{rca_nr_db_path};
	my $NR_path2            = $conf{rca_nr_copy_temp_path};
	my $NR_path3            = "$NR_path2\\/nr";

	for my $c ($carma_config, $smash_parallelblast, $smash_submit, $smash_config) {
		unless (-e $c) {
			die "ERROR & EXIT: Missing $c\n";
		}
	}

	# Other definitions
	my $sample = "annotate.$basename";
	@samples = "";
	$samples[0] = "annotate.$basename";
	my $LOG  = " >> $cwd/logs/$job/samples/MOCATJob_$job.$sample.log 2>> $cwd/logs/$job/samples/MOCATJob_$job.$sample.log ";
	my $cwd2 = $cwd;
	$cwd2 =~ s/\//\\\//g;

	print localtime() . ": Creating $job jobs...";

	# Prepare
	print JOB "rm -rf $cwd/$pathname/carma.annotation && ";
	print JOB "mkdir -p $cwd/$pathname/carma.annotation && ";
	print JOB "cd $cwd/$pathname/carma.annotation && cp $smash_config . && ";
	
	# Make soft link. Never used, but nice for the user to know which file was used
	print JOB "ln -fs $carma_file_to_annotate $cwd/$pathname/$basename && ";
	
	# Make blast jobs
	print JOB "$smash_parallelblast --flavor=NCBI --blast=blastp --database=$NR_path --query=$carma_file_to_annotate --label=$basename --outdir=$cwd/$pathname/carma.annotation  --splitsize=2500000 --evalue=10.0 --subjects=300 --stage=tmp --cpus=$processors $LOG && ";

	# carma doesn't do this, so we skip it too, it says, that with -F T, it may crash
	# print JOB "sed -i 's/-F F/-F T/' $basename.nr.blastp.defs.sh && ";
	# also, carma has the option -C 0 set, but after discussions with Kristoffer, we have it set to default -C 2

	# Edit blast config file
	print JOB "sed -i 's/-m 8/-m 9/' $basename.nr.blastp.defs.sh && ";
	print JOB "sed -i 's/local_db_location=.*/local_db_location=$NR_path2;/' $basename.nr.blastp.defs.sh && ";
	print JOB "sed -i \"s/blast_database_arg=.*/blast_database_arg=\\'-d $NR_path3\\';/\" $basename.nr.blastp.defs.sh && ";
	print JOB "sed -i \"s/local_output_location=.*/local_output_location=$cwd2\\/$pathname\\/temp;/\" $basename.nr.blastp.defs.sh && ";

	# Run blast
	print JOB "ssh $run_carma_annotate_host1 \"cd $cwd/$pathname/carma.annotation && $smash_submit --type SGE --name N$basename --cpus $processors --memory 4000 --extra_args '-sync y' $basename.nr.blastp.queue $LOG\" && ";

	# Prepare carma
	print JOB "cp $carma_config $cwd/$pathname/carma.annotation/$basename.nr.split && cd $cwd/$pathname/carma.annotation/$basename.nr.split && ";

	# Create carma queue jobs
	print JOB "rm -f $cwd/$pathname/carma.annotation/$basename.nr.split/execute.carma.sh && for f in $cwd/$pathname/carma.annotation/$basename.nr.split/*.blastp; do g=`echo \$f | sed 's/.blastp/.fa/'`; echo \"$bin_dir/carma --classify-blast --type p --database p --input \$f --fasta-input \$g --output \$f.carma\" >> $cwd/$pathname/carma.annotation/$basename.nr.split/execute.carma.sh; done && ";

	# Copy carma files to memory
	print JOB "ssh $run_carma_annotate_host2 \"for f in `qstat -f | grep \'BIP\' | cut -f 1 -d\" \" | cut -f 2 -d\"@\"`; do ssh \$f  \"rm -fr /dev/shm/carma; mkdir -p /dev/shm/carma; cp `which blastall` /dev/shm/carma; cp `which fastacmd` /dev/shm/carma; cp `which formatdb` /dev/shm/carma\"; done \" && ";

	# Edit carma config file
	print JOB "sed -i 's/blastall_script =.*/blastall_script = \\/dev\\/shm\\/carma\\/blastall/' $cwd/$pathname/carma.annotation/$basename.nr.split/carma.cfg && ";
	print JOB "sed -i 's/fastacmd_script =.*/fastacmd_script = \\/dev\\/shm\\/carma\\/fastacmd/' $cwd/$pathname/carma.annotation/$basename.nr.split/carma.cfg && ";
	print JOB "sed -i 's/formatdb_script =.*/formatdb_script = \\/dev\\/shm\\/carma\\/formatdb/' $cwd/$pathname/carma.annotation/$basename.nr.split/carma.cfg && ";
	print JOB "sed -i 's/cluster_tmp_dir = .*/cluster_tmp_dir = $cwd2\\/$pathname\\/carma.annotation\\/temp/' $cwd/$pathname/carma.annotation/$basename.nr.split/carma.cfg && ";
	print JOB "sed -i 's/blast_nr_database = .*/blast_nr_database = $NR_path3/' $cwd/$pathname/carma.annotation/$basename.nr.split/carma.cfg && ";

	# Submit carma jobs
	print JOB "ssh $run_carma_annotate_host2 \"cd $cwd/$pathname/carma.annotation/$basename.nr.split/ && $smash_submit --type SGE --name N$basename --cpus 2 --memory 8000 --extra_args '-sync y' $cwd/$pathname/carma.annotation/$basename.nr.split/execute.carma.sh\" && ";

	# Cat output files, and sort
	print JOB "cat $cwd/$pathname/carma.annotation/$basename.nr.split/*.blastp.carma > $cwd/$pathname/carma.annotation/$basename.carma && ";
	print JOB "ln -sf $cwd/$pathname/carma.annotation/$basename.carma $cwd/$pathname/$basename.carma && ";
	print JOB "sort $cwd/$pathname/$basename.carma > $cwd/$pathname/$basename.carma.sorted && ";

	# Remove files from memory
	print JOB "ssh $run_carma_annotate_host2 \"for f in `qstat -f | grep \'BIP\' | cut -f 1 -d\" \" | cut -f 2 -d\"@\"`; do ssh \$f  \"rm -rf dev/shm/carma\"; done\" \n";

	close JOB;
	print " OK!\n";
}

1;
