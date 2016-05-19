package MOCATRedoScaftig;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;
use Scalar::Util qw(looks_like_number);

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

sub create_job {
	### DEFINE VARIABLES AND OPEN OUTPUT FILE ###

	my $job        = $_[0];
	my $processors = $_[1];
	$ZIP =~ s/pigz.*/pigz -p $processors/;
	
	my $db   = $_[3];
	my $num  = looks_like_number($db);
	my $scaftig_date = `date +\%Y\%b\%d_\%H\%M\%S`;
	chomp($scaftig_date);
	if ( $num == 0 ) {
		die "ERROR & EXIT: Missing scaftig length setting or is not numeric.\n";
	}

	### JOB SPECIFIC ###
	open JOB, '>', "$cwd/MOCATJob_$job";
	print localtime() . ": Creating $job jobs...";
	foreach my $sample (@samples) {
		( my $max, my $avg, my $kmer ) = MOCATCore::get_kmer( $sample, $reads, "-scaf" );
		my $folder = "$cwd/$sample/$assembly.$reads.$conf{MOCAT_data_type}.K$kmer";
		my $file   = "$sample.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer";

		if ( $assembly eq 'assembly' ) {
			unless ( -e "$folder/$file.scafSeq" ) {
				die "\nERROR & EXIT: Missing scafSeq file\n";
			}
			print JOB "mv -f $folder/$file.scaftig $folder/$file.scaftig.$scaftig_date && $scr_dir/MOCATAssembly_getScaf.pl $db $folder/$file.scafSeq $sample > $folder/$file.scaftig\n";
		}
		if ( $assembly eq 'assembly.revised' ) {
			unless ( -e "$folder/$file.scaftig" ) {
				die "\nERROR & EXIT: Missing scaftig file\n";
			}
			print JOB "mv $folder/$file.scaftig $folder/$file.scaftig.$scaftig_date && $scr_dir/MOCATRedoScaftig_f2t.pl $folder/$file.scaftig.$scaftig_date | awk '{split(\$2,n,\"=\"); if (n[2] >= $db) {print \$0}}' | sed 's/\\t/\\n/' > $folder/$file.scaftig\n";
		}
	}
	print " OK!\n";
	close JOB;
}

1;
