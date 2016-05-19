package MOCATSampleStatus_dev;
use strict;
use warnings;
use MOCATCore;
use MOCATVariables;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

sub system_ {
	my $cmd = shift;
	system($cmd) == 0 or die("\nERROR & EXIT: system($cmd) failed: $!\n");
}

sub run {
	print localtime() . ": Calculating full status";

	my %hash;
	my %star;
	my %max;
	my %sample_max;
	my $folder_max = 0;

	my $original_sample_file = $sample_file;
	if ($original_sample_file =~ m/\//) {
		my @a = split '/', $original_sample_file;
		$original_sample_file = $a[-1];
	}
	my $sample_file          = "$cwd/SUMMARIES/$original_sample_file/$original_sample_file";
	my $sample_file_folder   = "$cwd/SUMMARIES/$original_sample_file/";
	system_(" mkdir -p $cwd/SUMMARIES/$original_sample_file");

	open OUT, ">$sample_file.status.full" or die "ERROR & EXIT: Cannot write to $sample_file.status.full";
	foreach my $sample (@samples) {
		my @folders = `find -L $cwd/$sample -type d -o -name '*.zip' | grep -v '\/temp' | grep -v '\/stats' | grep -v 'fastqc'`;
		foreach my $folder (@folders) {
			my $total_nonzero_files = 0;
			chomp($folder);
			my @files;
			if ( $folder =~ m/$sample$/ ) {
				@files = `ls -1 $folder | grep '.fq'`;

			}
			elsif ( $folder =~ m/\.zip$/ ) {
				$files[0] = "";
			}
			else {
				@files = `ls -1 $folder | grep -v 'qual_stats' | grep -v '\\.log\$'`;
			}

			my $folder_org = $folder;
			$folder =~ s/$cwd//;
			$folder =~ s/.solexaqa//;
			$folder =~ s/.fastx//;
			$folder =~ s/.K\d+//;
			$folder =~ s/^\/$sample//;
			$folder =~ s/^\///;

			if ( $folder eq "" ) {
				$folder = "A raw lanes";
			}

			my $type = "";
			my $star = " ";
			if ( $folder =~ m/^reads.(extracted|screened)\.(.*)/ ) {
				$folder = "D $2 ($1)";
				$type   = "screened";
			}
			if ( $folder =~ m/^reads.(mapped)\.(.*)/ ) {
				$folder = "C $2 ($1)";
				$type   = 'mapped';
			}
			if ( $folder =~ m/^reads.(filtered)\.(.*)/ ) {
				$folder = "E $2 ($1)";
				$type   = 'filtered';
			}
			if ( $folder =~ m/^(base|insert)\.coverage\.(.*)/ ) {
				$folder = "F $2 ($1 coverage)";
			}
			if ( $folder =~ m/^(assembly|assembly.revised)\.(.*)/ ) {
				$folder = "G $2 (assembly)";
				$type   = 'assembly';
			}
			if ( $folder =~ m/^(reads.processed)/ ) {
				$folder = "B $1";
			}
			if ( $folder =~ m/^gene.prediction.(assembly.revised|assembly).(.*)/ ) {
				$folder = "H $2 (GP_$1)";
				$type   = 'gene prediction';
			}
			if ( $folder =~ m/^$sample.(\S*)\.profile.(.*)/ ) {
				$folder = "I $2 ($1\_profile)";
				$type   = 'profile';
			}
			if ( $folder =~ m/^genes.fetched.MGs.(.*)/ ) {
				$folder = "H $1 (fetched_MGs)";
				$type   = 'fetch MG';
			}



			foreach my $file (@files) {
				chomp $file;
				my $file_size = -s "$folder_org/$file";
				my $toCheck   = "$folder_org/$file";
				if ( $file eq '' ) {
					$file_size = -s "$folder_org";
					$toCheck   = "$folder_org";
				}
				if ($file_size) {
					if ( $file_size > 0 ) {
						$total_nonzero_files++;
					}
				}
				if ( $type eq 'mapped' ) {
					if ( $toCheck =~ m/.soap.gz$/ ) {
						$star = "*";
					}
				}
				if ( $type eq 'filtered' ) {
					if ( $toCheck =~ m/.bam$/ ) {
						$star = "*";
					}
				}
				if ( $type eq 'assembly' ) {
					if ( $toCheck =~ m/.scaftig.gz$/ ) {
						$star = "*";
					}
				}
				if ( $type eq 'gene prediction' ) {
					if ( $toCheck =~ m/.tab$/ ) {
						$star = "*";
					}
				}
				if ( $type eq 'fetch MG' ) {
					if ( $toCheck =~ m/marker_genes_scores.table$/ ) {
						$star = "*";
					}
				}
				if ( $type eq 'profile' ) {
					if ( $toCheck =~ m/.zip$/ ) {
						$star = "*";
					}
				}

			}
			if ( $type eq 'screened' ) {
				my @files2;
				foreach my $f (@files) {
					if ($f =~ m/(single|pair).*\.fq/) {
						push @files2, $f;
					}
					if ($f =~ m/(single|pair).*\.fq.gz/) {
						push @files2, $f;
					}
				}				
				if ( scalar @files2 % 3 == 0 ) {
					$star = "*";
				} else {
					$star = " ";
				}
			}

			unless ( $folder eq "" ) {

				$hash{$folder}{$sample} = $total_nonzero_files;
				$star{$folder}{$sample} = $star;
				unless ( $max{$folder} ) {
					$max{$folder} = $total_nonzero_files;
				}
				if ( $total_nonzero_files > $max{$folder} ) {
					$max{$folder} = $total_nonzero_files;
				}

				#unless ( $sample_max{$sample} ) {
				#	$sample_max{$sample} = scalar @files;
				#}
				#if ( scalar @files > $sample_max{$sample} ) {
				#	$sample_max{$sample} = scalar @files;
				#}
				if ( $total_nonzero_files > length( $max{$folder} ) ) {
					$max{$folder} = $total_nonzero_files;
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
		$new_break = substr( $folder, 0, 1 );
		if ( $new_break ne $old_break ) {
			print OUT "-" x $total_length . "\n";
		}
		my $nice_folder = $folder;
		if ( $folder =~ m/\w \S+/ ) {
			$nice_folder = substr( $folder, 2, length($folder) - 1 );
		}
		$old_break = $new_break;
		print OUT "$nice_folder" . " " x ( $folder_max - length($folder) ) . "|";
		foreach my $sample ( sort @samples ) {
			unless ( $hash{$folder}{$sample} ) {
				$hash{$folder}{$sample} = "." x ( $sample_max{$sample} - 1 );
				$star{$folder}{$sample} = " ";

				#print OUT " " x ( $sample_max{$sample} ) . "|";
			}

			#else {
			print OUT "$hash{$folder}{$sample}$star{$folder}{$sample}" . " " x ( $sample_max{$sample} - length( $hash{$folder}{$sample} ) - 1 ) . "|";

			#}
		}
		print OUT "\n";
	}
	print OUT "-" x $total_length . "\n";
	close OUT;

	system_(
		"sed \'s/|/\\t/g\' $sample_file.status.full | sed \'s/ //g\' | grep -v \'^-\' | awk -F\"\\t\" \'
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = \$i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str\"\\t\"a[i,j];
        }
        print str
    }
}\' > $sample_file.status.full.reversed"
	);

	print " OK!\n";
}

1;
