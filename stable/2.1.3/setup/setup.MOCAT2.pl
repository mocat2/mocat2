#!/usr/bin/perl
use strict;
use warnings;
use Term::ANSIColor;

my $system = `uname -s`;
chomp($system);

my $agree;
my $version = "v 2.1.3";
my $sep     = "=================================================================================";
my $space = " ";
sub welcome {
	print "$sep\n";
	print color 'bold';
	print "                  MOCAT";
	print color 'reset';
	print " - Metagenomics Analysis Toolkit                $space";
	print color 'bold';
	print "$version\n";
	print color 'reset';
	print " by Jens Roat Kultima, Luis Pedro Coelho, Shinichi Sunagawa @ Bork Group, EMBL
$sep\n";
}
welcome();

unless ( -e "./src/MOCAT.pl" ) {
	die "ERROR & EXIT: This file must be executed from the MOCAT base directory\nIn this directory there should exist the /src, /bin, /data, /lib directories for MOCAT.\n";
}

print "\n";
print "SOFTWARE AGREEMENTS INFORMATION\n";
print "To proceed you have to agree with the licence agreements of all included software\n";

print "DOWNLOAD SIZE AND INSTALLATION TIME\n";
print "The Human Genome 19 (hg19) database can be downloaded and built.\n";
print "The download size is 6.5GB. The largest sample dataset to\n";
print "download is the simulated dataset, which requires a 3.5GB download\n\n";

print "Please type 'yes' and press Enter, if you agree with the licence agreements\n";
print "for each software used by MOCAT, and wish to proceed (otherwise type 'no'):\n";
print "(yes/no)\n";

chomp( $agree = <> );
unless ( $agree eq 'yes' || $agree eq 'YES' || $agree eq 'Yes' || $agree eq 'y' || $agree eq 'Y' ) {
	die "You did not agree with the terms. Thank you for considering MOCAT!\n";
}

print "\n";
print "SETUP THE CONFIG FILE\n";
print "What queuing system do you use on your system?\n";
print "This can be changed in the config file later.\n";
print "Others may be supported than listed below. Check the config file after installation.\n";
print "If you don't know, type 'none'\n";
print "('SGE', 'PBS', 'LSF' or 'none')\n";

my $samples = 0;
my $queue;
chomp( $agree = <> );
if ( $agree eq 'SGE' || $agree eq 'sge' ) {
	$queue = "SGE";
}
elsif ( $agree eq 'pbs' || $agree eq 'PBS' ) {
	$queue = "PBS";
}
elsif ( $agree eq 'lsf' || $agree eq 'LSF' ) {
	$queue = "LSF";
}
else {
	$queue = "none";
}

my @pes;
my $PE = "make";
if ( $queue eq 'SGE' ) {
	print "SGE SPECIFIC QUESTIONS\n";
	print "Please specify a parellel environment to use, choose one of the ones listed below.\n";
	print "This can be manually changed in the config file later. To list these envornments,\n";
	print "run 'qconf -spl' in a Terminal Window.\n";
	print "\n";
	print "Please enter one of these envornments to use:\n";
	@pes = `qconf -spl`;
	chomp(@pes);

	foreach my $pe (@pes) {
		print "- $pe\n";
	}
	print "\n";
	chomp( $PE = <> );

	my $exist = 0;
	foreach my $pe (@pes) {
		if ( $pe eq $PE ) {
			$exist = 1;
		}
	}
	unless ( $exist == 1 ) {
		print "\n";
		print "WARNING! THE PARALLELL ENVIRONMENT YOU SPECIFIED DID NOT EXIST IN THE LIST OF POSSIBLE ENVIRONMENTS!\n";
		print "To resolve this the environment has been set to 'make'\n";
		$PE = 'make';
	}

}

print "\n";
print "DOWNLOAD EXTERNAL DATABASES\n";
print "Would you like to download and install the hg19 reference database (approx 6.5GB)? (yes/no)\n";
print "This database is provided as a database to screen reads against (e.g. remove human contaminants).\n";
print "(yes/no):\n";

my $hg19 = 0;
chomp( $agree = <> );
if ( $agree eq 'yes' || $agree eq 'YES' || $agree eq 'Yes' || $agree eq 'y' || $agree eq 'Y' ) {
	$hg19 = 1;
}

print "\n";
print "DOWNLOAD ARTICLE DATASETS AND EXAMPLE FILES\n";
print "There are two different sample files that can be downloaded:\n";
print "- The simulated metagenome (100 strains, used in MOCAT paper)\n";
print "- The simulated metagenome (87 species, used in mOTU paper)\n";
print "- the HMP mock community data\n";

print "\n";
print "Would you like to download THE SIMULATED METAGENOME (100 strains) used inthe MOCAT paper (approx 3.5GB)? (yes/no)\n";
my $simulated = 0;
chomp( $agree = <> );
if ( $agree eq 'yes' || $agree eq 'YES' || $agree eq 'Yes' || $agree eq 'y' || $agree eq 'Y' ) {
	$simulated = 1;
}

print "\n";
print "Would you like to download THE SIMULATED METAGENOME (87 species) used in the mOTU paper (approx 2.3GB)? (yes/no)\n";
my $simulated2 = 0;
chomp( $agree = <> );
if ( $agree eq 'yes' || $agree eq 'YES' || $agree eq 'Yes' || $agree eq 'y' || $agree eq 'Y' ) {
	$simulated2 = 1;
}

print "\n";
print "Would you like to download THE MOCK COMMUNITY (approx 500 MB)? (yes/no)\n";
my $mock = 0;
chomp( $agree = <> );
if ( $agree eq 'yes' || $agree eq 'YES' || $agree eq 'Yes' || $agree eq 'y' || $agree eq 'Y' ) {
	$mock = 1;
}

$samples = 0;

print "\n";
print "Will you be running the SIMULATED METAGENOME or MOCK COMMUNITY on a system with less than 8 cores?\n";
print "(yes/no)\n";
my $NCORES = "na";
chomp( $agree = <> );
if ( $agree eq 'yes' || $agree eq 'YES' || $agree eq 'Yes' || $agree eq 'y' || $agree eq 'Y' ) {
	print "\nHow many cores?\n";
	chomp( $NCORES = <> );
}

print "\n";
print "DOWNLOADING AND SETTING UP DATABSES & DATASETS...\n";

system "mkdir -p data";
if ( $system =~ m/Darwin/ ) {
	system "curl -# http://vm-lux.embl.de/~kultima/share/MOCAT/data/1506MG.tar.gz > data/1506MG.tar.gz";
}
else {
	system "wget -t 5 -O data/1506MG.tar.gz \"http://vm-lux.embl.de/~kultima/share/MOCAT/data/1506MG.tar.gz\"";
}

if ( $hg19 == 1 ) {
	if ( $system =~ m/Darwin/ ) {
		system "curl -# http://vm-lux.embl.de/~kultima/share/MOCAT/data/hg19.tar.gz > data/hg19.tar.gz";
	}
	else {
		system "wget -t 5 -O data/hg19.tar.gz \"http://vm-lux.embl.de/~kultima/share/MOCAT/data/hg19.tar.gz\"";
	}
}

my @db = ( 'no', 'no' );
print "Extracting 1506MG DB...\n";
system "cd data && tar -zvxf 1506MG.tar.gz";
unless ( -s "data/1506MG" ) {
	print " FAILED! Perhaps file data/1506MG.tar.gz was not correctly downloaded from http://vm-lux.embl.de/~kultima/share/MOCAT/data/1506MG.tar.gz? Try downloading manually and re-run setup.\n";
}
else {
	$db[0] = "yes";
}
if ( $hg19 == 1 ) {
	print "Extracting hg19 DB...\n";
	system "cd data && tar -zvxf hg19.tar.gz";
	unless ( -s "data/hg19" ) {
		print " FAILED! Perhaps file data/hg19.tar.gz was not correctly downloaded from http://vm-lux.embl.de/~kultima/share/MOCAT/data/hg19.tar.gz? Try downloading manually and re-run setup.\n";
	}
	else {
		$db[1] = "yes";
	}
}

if ( $simulated == 1 ) {
	print "\n\nDOWNLOADING SIMULATED METAGENOME (100 strains):\n";
	system "mkdir -p article_datasets/simulated_metagenome_100strains/100strains";
	if ( $system =~ m/Darwin/ ) {
		system "curl -# http://vm-lux.embl.de/~kultima/share/MOCAT/data/illumina_100strains.1.fq.gz > article_datasets/simulated_metagenome_100strains/100strains/illumina_100strains.1.fq.gz";
		system "curl -# http://vm-lux.embl.de/~kultima/share/MOCAT/data/illumina_100strains.2.fq.gz > article_datasets/simulated_metagenome_100strains/100strains/illumina_100strains.2.fq.gz";
		system "curl -# http://vm-lux.embl.de/~kultima/share/MOCAT/data/100_strains_references.gz > article_datasets/simulated_metagenome_100strains/100_strains_references.gz";
		print "Extracting simulated metagenome reference DB...\n";
		system "gzcat -dc article_datasets/simulated_metagenome_100strains/100_strains_references.gz > article_datasets/simulated_metagenome_100strains/100_strains_references";
	}
	else {
		system "wget -t 5 -O article_datasets/simulated_metagenome_100strains/100strains/illumina_100strains.1.fq.gz \"http://vm-lux.embl.de/~kultima/share/MOCAT/data/illumina_100strains.1.fq.gz\"";
		system "wget -t 5 -O article_datasets/simulated_metagenome_100strains/100strains/illumina_100strains.2.fq.gz \"http://vm-lux.embl.de/~kultima/share/MOCAT/data/illumina_100strains.2.fq.gz\"";
		system "wget -t 5 -O article_datasets/simulated_metagenome_100strains/100_strains_references.gz \"http://vm-lux.embl.de/~kultima/share/MOCAT/data/100_strains_references.gz\"";
		print "Extracting simulated metagenome reference DB...\n";
		system "gzip -d article_datasets/simulated_metagenome_100strains/100_strains_references.gz";
	}
	print "\nCOMPLETED DOWNLOAD.\n\n";
}

if ( $simulated2 == 1 ) {
	print "\n\nDOWNLOADING SIMULATED METAGENOME (87 species):\n";
	system "mkdir -p article_datasets/simulated_metagenome_87/87_species";
	if ( $system =~ m/Darwin/ ) {
		system "curl -# http://vm-lux.embl.de/~kultima/share/MOCAT/data/simulated.87species.tar.gz > article_datasets/simulated_metagenome_87/87_species/simulated.87species.tar.gz";
	}
	else {
		system "wget -t 5 -O article_datasets/simulated_metagenome_87/87_species/simulated.87species.tar.gz \"http://vm-lux.embl.de/~kultima/share/MOCAT/data/simulated.87species.tar.gz\"";
	}
	print "Extracting simulated metagenome reference DB...\n";
	system "cd article_datasets/simulated_metagenome_87/87_species && tar -zvxf simulated.87species.tar.gz";
	print "\nCOMPLETED DOWNLOAD.\n\n";
}

if ( $mock == 1 ) {
	print "\n\nDOWNLOADING MOCK COMMUNITY:\n";
	system "mkdir -p article_datasets/mock_community/";
	system "mkdir -p article_datasets/mock_community/even_sample";
	if ( $system =~ m/Darwin/ ) {
		system "curl -# http://vm-lux.embl.de/~kultima/share/MOCAT/data/SRR172902.fq.gz > article_datasets/mock_community/even_sample/SRR172902.fq.gz";
		system "curl -# http://vm-lux.embl.de/~kultima/share/MOCAT/data/mock_community_ref.gz > article_datasets/mock_community/mock_community_ref.gz";
		print "Extracting mock community reference DB...\n";
		system "gzcat -dc article_datasets/mock_community/mock_community_ref.gz > article_datasets/mock_community/mock_community_ref";
	}
	else {
		system "wget -t 5 -O article_datasets/mock_community/even_sample/SRR172902.fq.gz \"http://vm-lux.embl.de/~kultima/share/MOCAT/data/SRR172902.fq.gz\"";
		system "wget -t 5 -O article_datasets/mock_community/mock_community_ref.gz \"http://vm-lux.embl.de/~kultima/share/MOCAT/data/mock_community_ref.gz\"";
		print "Extracting mock community reference DB...\n";
		system "gzip -d article_datasets/mock_community/mock_community_ref.gz";
	}
	print "\nCOMPLETED DOWNLOAD.\n\n";
}

chomp( my $cwd = `pwd` );
my $cwd2 = $cwd;
$cwd2 =~ s/\//\\\//g;

if ( $system =~ m/Darwin/ ) {
	system "sed -i '' 's/MOCAT_dir.*/MOCAT_dir                : $cwd2\\//' GETTING_STARTED/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_dir.*/MOCAT_dir                : $cwd2\\//' article_datasets/mock_community/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_dir.*/MOCAT_dir                : $cwd2\\//' article_datasets/simulated_metagenome_100strains/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_dir.*/MOCAT_dir                : $cwd2\\//' article_datasets/simulated_metagenome_87/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_dir.*/MOCAT_dir                : $cwd2\\//' article_datasets/make_and_annotate_gene_catalog/MOCAT.cfg";
	

	system "sed -i '' 's/MOCAT_qsub_system.*/MOCAT_qsub_system        : $queue/' GETTING_STARTED/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_qsub_system.*/MOCAT_qsub_system        : $queue/' article_datasets/mock_community/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_qsub_system.*/MOCAT_qsub_system        : $queue/' article_datasets/simulated_metagenome_100strains/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_qsub_system.*/MOCAT_qsub_system        : $queue/' article_datasets/simulated_metagenome_87/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_qsub_system.*/MOCAT_qsub_system        : $queue/' article_datasets/make_and_annotate_gene_catalog/MOCAT.cfg";

	if ( $queue eq 'SGE' ) {
		system "sed -i '' 's/MOCAT_SGE_parallell_env.*/MOCAT_SGE_parallell_env  : $PE [" . join( ",", @pes ) . "]/' GETTING_STARTED/MOCAT.cfg";
		system "sed -i '' 's/MOCAT_SGE_parallell_env.*/MOCAT_SGE_parallell_env  : $PE [" . join( ",", @pes ) . "]/' article_datasets/mock_community/MOCAT.cfg";
		system "sed -i '' 's/MOCAT_SGE_parallell_env.*/MOCAT_SGE_parallell_env  : $PE [" . join( ",", @pes ) . "]/' article_datasets/simulated_metagenome_100strains/MOCAT.cfg";
		system "sed -i '' 's/MOCAT_SGE_parallell_env.*/MOCAT_SGE_parallell_env  : $PE [" . join( ",", @pes ) . "]/' article_datasets/simulated_metagenome_87/MOCAT.cfg";
		system "sed -i '' 's/MOCAT_SGE_parallell_env.*/MOCAT_SGE_parallell_env  : $PE [" . join( ",", @pes ) . "]/' article_datasets/make_and_annotate_gene_catalog/MOCAT.cfg";
	}

}
else {
	system "sed -i 's/MOCAT_dir.*/MOCAT_dir                : $cwd2\\//' GETTING_STARTED/MOCAT.cfg";
	system "sed -i 's/MOCAT_dir.*/MOCAT_dir                : $cwd2\\//' article_datasets/mock_community/MOCAT.cfg";
	system "sed -i 's/MOCAT_dir.*/MOCAT_dir                : $cwd2\\//' article_datasets/simulated_metagenome_100strains/MOCAT.cfg";
	system "sed -i 's/MOCAT_dir.*/MOCAT_dir                : $cwd2\\//' article_datasets/simulated_metagenome_87/MOCAT.cfg";
	system "sed -i 's/MOCAT_dir.*/MOCAT_dir                : $cwd2\\//' article_datasets/make_and_annotate_gene_catalog/MOCAT.cfg";

	system "sed -i 's/MOCAT_qsub_system.*/MOCAT_qsub_system        : $queue [SGE,PBS,LSF,none]/' GETTING_STARTED/MOCAT.cfg";
	system "sed -i 's/MOCAT_qsub_system.*/MOCAT_qsub_system        : $queue [SGE,PBS,LSF,none]/' article_datasets/mock_community/MOCAT.cfg";
	system "sed -i 's/MOCAT_qsub_system.*/MOCAT_qsub_system        : $queue [SGE,PBS,LSF,none]/' article_datasets/simulated_metagenome_100strains/MOCAT.cfg";
	system "sed -i 's/MOCAT_qsub_system.*/MOCAT_qsub_system        : $queue [SGE,PBS,LSF,none]/' article_datasets/simulated_metagenome_87/MOCAT.cfg";
	system "sed -i 's/MOCAT_qsub_system.*/MOCAT_qsub_system        : $queue [SGE,PBS,LSF,none]/' article_datasets/make_and_annotate_gene_catalog/MOCAT.cfg";

	if ( $queue eq 'SGE' ) {
		system "sed -i 's/MOCAT_SGE_parallell_env.*/MOCAT_SGE_parallell_env  : $PE [" . join( ",", @pes ) . "]/' GETTING_STARTED/MOCAT.cfg";
		system "sed -i 's/MOCAT_SGE_parallell_env.*/MOCAT_SGE_parallell_env  : $PE [" . join( ",", @pes ) . "]/' article_datasets/mock_community/MOCAT.cfg";
		system "sed -i 's/MOCAT_SGE_parallell_env.*/MOCAT_SGE_parallell_env  : $PE [" . join( ",", @pes ) . "]/' article_datasets/simulated_metagenome_100strains/MOCAT.cfg";
		system "sed -i 's/MOCAT_SGE_parallell_env.*/MOCAT_SGE_parallell_env  : $PE [" . join( ",", @pes ) . "]/' article_datasets/simulated_metagenome_87/MOCAT.cfg";
		system "sed -i 's/MOCAT_SGE_parallell_env.*/MOCAT_SGE_parallell_env  : $PE [" . join( ",", @pes ) . "]/' article_datasets/make_and_annotate_gene_catalog/MOCAT.cfg";
	}
}

unless ( $NCORES eq 'na' ) {
	if ( $system =~ m/Darwin/ ) {
		system "sed -i '' 's/\\&\\&/-p $NCORES \\&\\&/' article_datasets/mock_community/RUN_PIPELINE.sh";
		system "sed -i '' 's/\\&\\&/-p $NCORES \\&\\&/' article_datasets/simulated_metagenome_100strains/RUN_PIPELINE.sh";
		
	}
	else {
		system "sed -i 's/\\&\\&/-p $NCORES \\&\\&/' article_datasets/mock_community/RUN_PIPELINE.sh";
		system "sed -i 's/\\&\\&/-p $NCORES \\&\\&/' article_datasets/simulated_metagenome_100strains/RUN_PIPELINE.sh";
	}
}

if ( $system =~ m/Darwin/ ) {
	system "sed -i '' 's/MOCAT_zip_program.*/MOCAT_zip_program        : gzip [gzip]/' GETTING_STARTED/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_zip_program.*/MOCAT_zip_program        : gzip [gzip]/' article_datasets/mock_community/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_zip_program.*/MOCAT_zip_program        : gzip [gzip]/' article_datasets/simulated_metagenome_100strains/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_zip_program.*/MOCAT_zip_program        : gzip [gzip]/' article_datasets/simulated_metagenome_87/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_zip_program.*/MOCAT_zip_program        : gzip [gzip]/' article_datasets/make_and_annotate_gene_catalog/MOCAT.cfg";

	system "sed -i '' 's/MOCAT_mapping_mode.*/MOCAT_mapping_mode       : random [allbest *NOT SUPPORTED ON OSX*,random,unique]/' GETTING_STARTED/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_mapping_mode.*/MOCAT_mapping_mode       : random [allbest *NOT SUPPORTED ON OSX*,random,unique]/' article_datasets/mock_community/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_mapping_mode.*/MOCAT_mapping_mode       : random [allbest *NOT SUPPORTED ON OSX*,random,unique]/' article_datasets/simulated_metagenome_100strains/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_mapping_mode.*/MOCAT_mapping_mode       : random [allbest *NOT SUPPORTED ON OSX*,random,unique]/' article_datasets/simulated_metagenome_87/MOCAT.cfg";
	system "sed -i '' 's/MOCAT_mapping_mode.*/MOCAT_mapping_mode       : random [allbest *NOT SUPPORTED ON OSX*,random,unique]/' article_datasets/make_and_annotate_gene_catalog/MOCAT.cfg";

}

system "echo \"\n### EXPORTED BY MOCAT ###\nexport PERL5LIB=\\\$PERL5LIB:$cwd/src\nPATH=\\\$PATH:$cwd/src\n### EXPORTED BY MOCAT ###\n\" >> ~/.bashrc";
system "echo \"\n### EXPORTED BY MOCAT ###\nexport PERL5LIB=\\\$PERL5LIB:$cwd/src\nPATH=\\\$PATH:$cwd/src\n### EXPORTED BY MOCAT ###\n\" >> ~/.bash_profile";

print "CONFIG FILE INFORMATION\n";
print "Changed MOCAT.cfg for this system.\n";
print "Setup set the flag 'MOCAT_dir' to $cwd\n";

print "\n";
print "UNIX / OSX STARTUP FILE CHANGES\n";
print "Exported $cwd/src to \$PERL5LIB and \$PATH variable\n";
print "to ~/.bashrc and ~/.bash_profile\n";

system "chmod -R uga+rwx bin/*";

if ( scalar $db[0] eq "yes" && ( $db[1] eq "yes" || $hg19 == 0 ) ) {
	print "\n$sep\nSUCCESS! MOCAT has now been setup and is ready to be used!\nYou start it by running 'MOCAT.pl' from any directory.\n$sep\n\n";
	print "LAST THING TO DO BEFORE RUNNING MOCAT.pl\n";
	print "Please execute commands below, to ensure Perl libraries are correctly loaded\n";
	print "EXECUTE: source ~/.bashrc; source ~/.bash_profile\n";
	print "and then you can run MOCAT using 'MOCAT.pl'\n";
	exit 0;
}
else {
	if ( $db[0] eq "no" ) {
		print "\nNOTE! The 1506MG database was not correctly setup. Perhaps download from http://vm-lux.embl.de/~kultima/share/MOCAT/data/1506MG.tar.gz and re-run this setup?\n";
	}
	if ( $hg19 == 1 ) {
		if ( $db[1] eq "no" ) {
			print "\nNOTE! The hg19 database was not correctly setup. Perhaps download from http://vm-lux.embl.de/~kultima/share/MOCAT/data/hg19.tar.gz and re-run this setup?\n";
		}
	}
	die "MOCAT setup exited with errors.\n";
}

exit 0;
