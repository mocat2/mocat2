package MOCATSampleStatus_dev;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012
# This code is released under GNU GPL v3.

sub run {
	print localtime() . ": Calculating full status";

	my %hash;
	my %max;
	my %sample_max;
	my $folder_max = 0;
	open OUT, ">$sample_file.status.full";
	foreach my $sample (@samples) {
		my @folders = `find $cwd/$sample -type d | grep -v '\/temp' | grep -v '\/stats'`;
		foreach my $folder (@folders) {
			chomp($folder);
			my @files;
			if ($folder =~ m/$sample$/) {
				@files = `ls -1 $folder | grep '.fq'`;
			} else {
				@files = `ls -1 $folder | grep -v 'qual_stats' | grep -v '\\.log\$'`;
			}
			$folder =~ s/$cwd//;
			$folder =~ s/.solexaqa//;
			$folder =~ s/.fastx//;
			$folder =~ s/.K\d+//;
			$folder =~ s/^\/$sample//;
			$folder =~ s/^\///;
			
			if ($folder eq "") {
				$folder = "A raw lanes";
			}
			
			if ($folder =~ m/^reads.(extracted|screened|mapped|filtered)\.(.*)/) {
				$folder = "C $2 ($1)";
			}
			if ($folder =~ m/^(base|insert)\.coverage\.(.*)/) {
				$folder = "F $2 ($1 coverage)";
			}
			if ($folder =~ m/^(assembly|assembly.revised)\.(.*)/) {
				$folder = "D $2 (assembly)";
			}
			if ($folder =~ m/^(reads.processed)/) {
				$folder = "B $1";
			}
			if ($folder =~ m/^gene.prediction.(assembly.revised|assembly).(.*)/) {
				$folder = "E $2 (GP $1)";
			}
			
			
			
			unless ( $folder eq "" ) {

				$hash{$folder}{$sample} = "X" x scalar @files;
				unless ( $max{$folder} ) {
					$max{$folder} = scalar @files;
				}
				if ( scalar @files > $max{$folder} ) {
					$max{$folder} = scalar @files;
				}
				unless ( $sample_max{$sample} ) {
					$sample_max{$sample} = scalar @files;
				}
				if ( scalar @files > $sample_max{$sample} ) {
					$sample_max{$sample} = scalar @files;
				}

				if ( scalar @files > length( $max{$folder} ) ) {
					$max{$folder} = scalar @files;
				}
				if ( length($folder) > $folder_max ) {
					$folder_max = length($folder);
				}
			}
			print ".";
		}
	}
	foreach my $sample (@samples) {
		unless ( $sample_max{$sample} ) {
			$sample_max{$sample} = length($sample);
		}
		if ( $sample_max{$sample} < length($sample) ) {
			$sample_max{$sample} = length($sample);
		}
	}
	my $total_length = $folder_max - 1;
	print OUT "SAMPLE" . " " x ( $folder_max - 8 ) . "|";
	foreach my $sample ( sort @samples ) {
		print OUT "$sample" . " " x ( $sample_max{$sample} - length($sample) ) . "|";
		$total_length = $total_length + $sample_max{$sample} + 1;
	}
	print OUT "\n";
	my $old_break = "";
	my $new_break;
	foreach my $folder ( sort keys %hash ) {
		$new_break = substr ($folder, 0, 1);
		if ($new_break ne $old_break) {
			print OUT "-" x $total_length . "\n";
		}
		my $nice_folder = $folder;
		if ($folder =~ m/\w \S+/) {
			$nice_folder = substr($folder, 2, length($folder)-1)
		}
		$old_break = $new_break;
		print OUT "$nice_folder" . " " x ( $folder_max - length($folder) ) . "|";
		foreach my $sample ( sort @samples ) {
			unless ( $hash{$folder}{$sample} ) {
				print OUT " " x ( $sample_max{$sample} ) . "|";
			}
			else {
				print OUT "$hash{$folder}{$sample}" . " " x ( $sample_max{$sample} - length( $hash{$folder}{$sample} ) ) . "|";
			}
		}
		print OUT "\n";
	}
	print OUT "-" x $total_length . "\n";
	close OUT;
	print " OK!\n";
}

1;
