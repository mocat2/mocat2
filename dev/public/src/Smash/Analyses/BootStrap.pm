package Smash::Analyses::BootStrap;

use strict;
use warnings;
use FAlite;
use File::Path;
use Smash::Analyses::Resampler qw(:all);

sub bootstrap_dir {"bootstrap"};

# $array_ref is a ref to array containing metagenome ids
# $depth_file is a file with tab-delimited depth information for each item
sub bootstrap_list {
	my ($type, $label, $REPLICATES, $PROGRESS, $array_ref, $replacement, $depth_file) = @_;
	my @items = @$array_ref;
	my $depth_hash_ref;
	if (defined($depth_file)) {
		print $PROGRESS "Reading depth information ... ";
		open(DEPTH, "<$depth_file") || die "Cannot open $depth_file: $!";
		while (<DEPTH>) {
			chomp();
			my ($entry, $depth) = split(/\t+/);
			$depth_hash_ref->{$entry} = $depth;
		}
		close(DEPTH);
		print $PROGRESS "done\n";
	}
	my $smash = new Smash::Core(COLLECTION => 'MC2');
	$smash->init();
	print $PROGRESS "Making $REPLICATES virtual bootstrap replicates of $type ... \n";
	ITEM: foreach my $item (@items) {
		my (undef, $metagenome) = $smash->parse_concat_id($item);
		my $boot_dir = sprintf("%s/%s", $smash->bootstrap_dir($metagenome), $label);
		my $fasta;
		if ($type eq "reads") {
			my $read_dir = $smash->read_dir($item);
			$fasta = "$read_dir/$item.fasta";
		} elsif ($type eq "contigs") {
			my $assembly_dir = $smash->assembly_dir($item);
			$fasta = "$assembly_dir/$item.contig.fa";
		} elsif ($type eq "genes") {
			my $genepred_dir = $smash->genepred_dir($item);
			$fasta = "$genepred_dir/$item.gene.fa";
		}
		printf $PROGRESS "\t%-12s ", $item;
		bootstrap_fasta($fasta, $boot_dir, $REPLICATES, $PROGRESS, $replacement, $depth_hash_ref);
		print $PROGRESS "\n";
	}
	$smash->finish();
}

sub bootstrap_fasta {
		my ($fasta_file, $boot_dir, $REPLICATES, $PROGRESS, $replacement, $depth_ref) = @_;
		my $must_simulate = 0;
		REPLICATE: for (my $i=0; $i<$REPLICATES; $i++) {
			my $list_file = "$boot_dir/replicate.$i";
			if (! -f $list_file) {
				$must_simulate = 1;
				last REPLICATE;
			}
		}

		if ($must_simulate == 0) {
			print $PROGRESS "in tact";
			return;
		}

		my %Depth;
		if ($depth_ref) {
			%Depth = %$depth_ref;
		}

		# Read all fasta names 
		my @names = ();
		open(FASTA, "<$fasta_file") || die "Cannot open $fasta_file: $!";
		my $fasta = new FAlite(\*FASTA);
		while(my $entry = $fasta->nextEntry) {
			my $def = $entry->def;
			$def =~ s/\s+.*//;
			$def =~ s/^>//;
			for (my $i=0; $i<($Depth{$def}||1); $i++) { # If there is no depth information, then add it once. depth cannot be 0
				push(@names, $def);
			}
		}
		close(FASTA);

		# Generate only the replicates that do not exist.
		# It is ok if some replicates are made at a different
		# point of time. It's actually better that way!
		# Sample the same number of names with replacement

		mkpath "$boot_dir";
		for (my $i=0; $i<$REPLICATES; $i++) {
			my $list_file = "$boot_dir/replicate.$i";
			if (! -f $list_file) {
				open(LIST, ">$list_file") || die "Cannot open $list_file: $!";
				my $list;
				if ($replacement==1) {
					$list = sample_with_replacement(\@names);
				} else {
					$list = sample_without_replacement(\@names);
				}
				for (my $j=0; $j<@$list; $j++) {
					print LIST $list->[$j];
					print LIST "\n";
				}
				close(LIST);
			}
			print $PROGRESS ".";
		}
}

1;
