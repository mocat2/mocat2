#!/usr/bin/env perl

package Smash::Analyses::Comparative;
use strict;
use warnings;
use File::Path;
use Smash::Core;
use Smash::Analyses::Comparative::Functional;
use Smash::Analyses::Comparative::RefGenome;
use Smash::Analyses::Comparative::16S;
use Smash::Utils::Taxonomy qw($BergeyTree $NCBITree %Smash2RDP get_taxonomic_rank);
use Smash::Utils::MatrixIO qw(:all);

my $REPLICATES = 100;
my $PRINT_PRECISION = 8;
my $UNIFORM_PRIOR = 0.000001;
my $FUNC2NAME;
my $Tree;

sub new {
	my $class = shift;
	my %param = @_;
	my $this  = bless {%param}, $class;
	select(STDERR); $| = 1; select(STDOUT);
	$this->{PROGRESS} = \*STDERR;
	return $this;
}

sub init {
	my $this     = shift;
	$this->{FEATURE_COUNT} = {};
	$this->{LABELS} = [];
}

sub finish {
	my $this = shift;
	if ($this->smash) {
		$this->smash->finish();
	}
}

sub attach_sample {
	my $this     = shift;
	my $collection = shift || die "attach_sample needs collection";
	my $smash    = $this->smash;

	$smash->finish() if $smash;
	$smash = new Smash::Core(COLLECTION => $collection);
	$smash->init();
	$this->{SMASH} = $smash;
}

sub attach_genome {
	my $this     = shift;
	my $genome   = shift || die "attach_genome needs genome";
	my $label    = shift || die "attach_genome needs label";
	my $smash    = $this->smash;

	$smash->close_refgenomedb_handle() if $smash;
	$smash->finish() if $smash;
	$smash = new Smash::Core();
	$smash->init();
	$this->{SMASH} = $smash;

	my $dbh   = $smash->get_refgenomedb_handle();
	my $sth   = $dbh->prepare_cached("SELECT gene_count FROM gene_stats WHERE taxonomy_id=?");
	$sth->execute($genome);
	my ($count) = $sth->fetchrow_array();
	$sth->fetchrow_array();
	$this->{LABEL} = $label;
	$this->{DATA_COUNT}->{$label} = $count;
	push(@{$this->{LABELS}}, $label);
}

# passed to the constructor

sub name      {shift->{NAME}}
sub type      {shift->{TYPE}}
sub experiment{shift->{EXPERIMENT}}
sub refdb     {shift->{REF_DB}}
sub replicate {shift->{REPLICATE}}
sub label     {shift->{LABEL}}
sub identity  {shift->{IDENTITY}}
sub multilevel{shift->{MULTILEVEL}}
sub normalize {shift->{NORMALIZE}}
sub level     {shift->{LEVEL}}
sub is_local  {shift->{IS_LOCAL}}
sub filterlow {shift->{FILTERLOW}}
sub tree      {shift->{TREE}}

sub feature_count {shift->{FEATURE_COUNT}}
sub data_count    {shift->{DATA_COUNT}}
sub smash         {shift->{SMASH}}
sub labels        {shift->{LABELS}}
sub progress      {shift->{PROGRESS}}

# convert the tax_id's in the input feature_vector to names
# by applying &funcp_tax_to_name on the input field1

sub convert_id_to_name {
	my $this = shift;
	my $feature_count = shift;
	return apply_to_field1($feature_count, $this->funcp_id_to_name);
}

# return one of the functional pointers, depending on the tree.

sub funcp_id_to_name {
	my $this = shift;
	my $type = $this->type;
	# 16S doesnt init Taxonomy, so do it here - but tree specific
	if ($type eq "refgenome" || $type eq "16S") {
		if ($this->tree =~ /bergey/i) {
			Smash::Utils::Taxonomy::init("BERGEY");
			Smash::Utils::Taxonomy::init("NCBI.DIV");
			$Tree = $BergeyTree;
		} else {
			Smash::Utils::Taxonomy::init("NCBI.DIV");
			$Tree = $NCBITree;
		}
		return \&tax_to_name;
	} elsif ($type eq "functional") {
		$FUNC2NAME = $this->{uc($this->level)."2NAME"};
		return \&func_to_name;
	}
}

# pointer to a function that converts NCBI/Bergey tax_id to name

sub tax_to_name {
	my $tax  = shift;
	return $tax if $tax == -1;
	my $node = $Tree->nodes->{$tax};
	if (!$node) {
		die "No node for $tax!\n";
	}
	my $name =  $Tree->nodes->{$tax}->name || $tax;
	   $name =~ s/\s/_/g;
	return $name;
}

# pointer to a function that converts OG id to description

sub func_to_name {
	my $query = shift;
	my $name = $FUNC2NAME->{$query};
	return $query unless $name;
	$name =~ s/\s/_/g;
	$name = "${query}__${name}";
	return $name;
}

sub get_popup_info {};

sub merge_feature_vectors {
	my $this = shift;
	my $feature_count = shift;
	my $merge_level   = shift;

	my $MergedFeatureCount;
	my $labels = get_column_labels($feature_count);

	my @ids = grep {$_ >= 0} keys %$feature_count; # skip unmapped = -1 here

	if ($this->tree =~ /NCBI/i) {
		foreach my $genome_id (@ids) {
			# the following line is obsolete since refgenome freeze4
			# my $ncbi_tax_id = $this->get_ncbi_tax_id($tax_id);
			my ($tax_id) = split(/\./, $genome_id);
			my $rank = get_taxonomic_rank($NCBITree, $tax_id, $merge_level);
			foreach my $label (@$labels) {
				$MergedFeatureCount->{$rank}->{$label} += ($feature_count->{$genome_id}{$label}||0);
			}
		}
	} elsif ($this->tree =~ /bergey/i) {

		# if this was refgenome mapping, then remap to RDP

		TAX_ID:foreach my $genome_id (@ids) {
			my $rdp_tax_id = $genome_id;
			if ($this->type eq "refgenome") {
				$rdp_tax_id = $Smash2RDP{$genome_id};
				if (!defined($rdp_tax_id)) {
					warn "$genome_id has no match\n";
					next TAX_ID;
				}
			}
			my $rank = get_taxonomic_rank($BergeyTree, $rdp_tax_id, $merge_level);
			foreach my $label (@$labels) {
				$MergedFeatureCount->{$rank}->{$label} += ($feature_count->{$genome_id}{$label}||0);
			}
		}
	}

	# handle unmapped
	foreach my $label (@$labels) {
		$MergedFeatureCount->{-1}->{$label} += ($feature_count->{-1}{$label}||0);
	}

	return $MergedFeatureCount;
}

sub bergey2ncbi_feature {
	my $this   = shift;
	my $fc     = shift;
	my $new_fc = {};

	Smash::Utils::Taxonomy::init("NCBI.DIV");

	my @missing; # id's that will be missing in the other tree

	my $labels = get_column_labels($fc);

	foreach my $feature (keys %{$fc}) {
		my $new_feature;
		if ($feature == -1) {
			$new_feature = $feature;
		} else {
			$new_feature = Smash::Utils::Taxonomy::bergey2ncbi($feature);
		}
		if ($new_feature) {
			foreach my $label (@$labels) {
				$new_fc->{$new_feature}->{$label} = $fc->{$feature}->{$label};
			}
		} else {
			push(@missing, $feature);
		}
	}

	if (@missing) {
		warn "WARNING: The following taxon from RDP Taxonomy could not be ported to NCBI Taxonomy:\n";
		foreach my $f (@missing) {
			$BergeyTree->nodes->{$f}->print_node(1, \*STDERR);
		}
	}

	return $new_fc;
}

sub ncbi2bergey {
	my $this   = shift;
	my $fc     = shift;
}

# only implemented by the subclasses

sub convert_to_abundance {shift}

sub bootstrap_elements {
	my $this          = shift;
	my $run_replicate = $this->replicate;
	my $feature_count = $this->feature_count;
	my @labels        = @{$this->labels};
	print {$this->progress} "Bootstrapping distance matrices ";

	my $output_dir    = "bootstrap_files";
	mkpath($output_dir, {mode => $this->smash->file_perm});

	my ($min, $max);
	if (defined($run_replicate) && $run_replicate >= 0) {
		$min = $run_replicate;
		$max = $run_replicate;
	} else {
		$min = 0;
		$max = $REPLICATES - 1;
	}

	foreach my $replicate ($min..$max) {
		my $resampled_features = {};
		foreach my $label (@labels) {

			# make a concatenated range of numbers
			# and generate a random number within the range
			# depending on the position of the random number, a feature is sampled
			# e.g., (a => 2, b => 4, c => 8) will result in [0,2), [2,6), [6,14)
			# and 14 samples will be drawn in [0,14) using rand(14). They then go to
			# the right feature depending on the range that they fall in.

			# make the list

			my $sample_count = 0;
			my @features =  sort keys %{$feature_count};
			my $lower_bounds = [0];
			foreach my $feature (@features) {
				my $value = $feature_count->{$feature}{$label} || 0;
				$sample_count += $value;
				push(@$lower_bounds, $sample_count);
			}

			# initialize resampling to 0

			map {$resampled_features->{$_}{$label} = 0} @features;

			# make $sample_count random numbers and assign them to the right class using probs.
			# this is done using a binary search with only lower bounds, so be careful if you decide to change this.
			# check binary_search for more details.

			foreach my $i (1..$sample_count) {
				my $f = $features[binary_search($lower_bounds, rand($sample_count))];
				$resampled_features->{$f}{$label}++;
			} 
		}

		# from counts of KO(functional) or templates(refgenome), we get an abundance that takes care of
		# size of module/pathway (functional) or genome(refgenome)

		$resampled_features = $this->convert_to_abundance($resampled_features) if $this->{GENOMESIZE_NORMALIZE};

		# merge things that are subgroups of a given level: e.g., all species under a genus will be summed

		if ($this->level) {
			$resampled_features = $this->merge_feature_vectors($resampled_features, $this->level); # now %FeatureCount is merged
		}

		zero_strip_matrix($resampled_features);

		# make distance matrices

		$this->make_distance_matrices($output_dir, $resampled_features, $replicate);

		print {$this->progress} Smash::Core->progress_bar($replicate+1);
	}
	print {$this->progress} " done\n";
}

sub binary_search {
	my ($array, $x) = @_;
	my $L = 0;
	my $R = $#$array;
	my ($p, $M);
	do {
		$p = int(($R-$L)/2);
		$M = $L+$p;
		if ($x < $array->[$M]) {
			$R = $M;
		} elsif ($x >= $array->[$M]) {
			if ($x < $array->[$M+1]) {
				return $M;
			} else {
				$L = $M;
			}
		}
	} while ($p > 0);
	die "Couldn't find $x in ".join(",", @$array);
}

sub normalize_data_by_col_sum {
	my $this          = shift;
	my $feature_count = shift;
	my $samples       = shift;
	my $new_feature_count = {};
	my @features = keys %$feature_count;
	foreach my $label (@$samples) {
		my $sum = 0;
		foreach my $feature (@features) {
			$sum += ($feature_count->{$feature}->{$label} || 0);
		}
		foreach my $feature (@features) {
			$new_feature_count->{$feature}->{$label} = ($feature_count->{$feature}->{$label}||0) / $sum;
		}
	}
	return $new_feature_count;
}

# 12.11.2009. The input data will be fraction already. They MUST sum to 1.
# This frees this script from the burden of tracking data counts for normalization.
# Also, tracking unassigned as {-1} is also not necessary, since it MUST be done
# before we get here.
# priors are in terms of percentage, so no need to get the sum over all samples

sub add_uniform_prior_and_renormalize {

	my $this        = shift;
	my $prior       = shift;
	my $counts_ref  = shift;
	my $samples_ref = shift;
	my %Counts      = %$counts_ref;
	my @samples     = @$samples_ref;

	if (!defined($prior)) {
		die "prior should be specified";
	}

	my %RemappedCounts;
	my @features = keys %Counts;
	foreach my $sample (@samples) {
		foreach my $feature (@features) {
			####
			# 1. add prior
			# 2. renormalize so it sums up to 1
			####
			my $frac = (($Counts{$feature}{$sample}||0) + $prior);
			$frac /= (1+@features*$prior);
			$RemappedCounts{$feature}{$sample} = $frac;
		}
	}
	return \%RemappedCounts;
}

sub make_distance_matrices {
	use Smash::Utils::MatrixIO qw(normalize_cols_by_sum get_column_labels);
	my $this          = shift;
	my $output_dir    = shift;
	my $feature_count = shift;
	my $replicate     = shift;
	my $labels        = get_column_labels($feature_count); # has to be consistent across header/data
	my $remapped_counts;

	# make jensen-shannon and kullback-leibler distances
	# add uniform prior before hand, so we dont get into division by zero

	if ($this->normalize eq "hits") {
		delete $feature_count->{-1};
	}

	$remapped_counts = normalize_cols_by_sum($feature_count);
	$remapped_counts = $this->add_uniform_prior_and_renormalize($UNIFORM_PRIOR, $remapped_counts, $labels);
	#delete $feature_count->{-1};
	
	$this->make_distance_matrix($output_dir, "kld", $remapped_counts, $labels, $replicate);
	$this->make_distance_matrix($output_dir, "jsd", $remapped_counts, $labels, $replicate);

	# make euclidean distance after standard-normalizing
	# also, remove unmapped

	$remapped_counts = copy_matrix($feature_count);
	if (defined($remapped_counts->{-1})) {
		delete $remapped_counts->{-1};
	}
	$remapped_counts = normalize_cols_by_sum($remapped_counts);
	if ($this->filterlow) {
		$this->filter_low_occurrences($this->filterlow, $remapped_counts);
	}
	$this->make_distance_matrix($output_dir, "euclidean", $remapped_counts, $labels, $replicate);
}

sub make_distance_matrix {

	my $this        = shift;
	my $output_dir  = shift;
	my $metric      = shift;
	my $counts_ref  = shift;
	my $samples_ref = shift;
	my $replicate   = shift;
	my %Counts      = %$counts_ref;
	my @samples     = @$samples_ref;

	my $print_precision = $PRINT_PRECISION;
	if ($metric eq "rawcount") {
		$print_precision = 6;
	}

	# choose filenames

	my ($distance_file, $feature_file) = map {sprintf("$output_dir/%s.%s.%s.%s.%s.%s", $this->name, $_, $metric, $this->type, $this->level, $this->refdb)} qw(distance feature_vector);
	if (defined($replicate)) {
		$distance_file .= ".$replicate";
		$feature_file  .= ".$replicate";
	}

	# make distance matrix
	# ignore feature labels (first column) and keep sample labels (first row)
	# thus cut -f2-
	# then pass on to distance_matrix
	# noise = less than 0.01% of the total of the matrix will be removed by distance_matrix, if --noise
	# columns are normalized by total data count by this perl script
	# rows will be normalized by distance_matrix, if --rownorm

	open(FEATURES, ">$feature_file");

	# Hack for comparing Jeroen's matrix
	#my @labels = qw(AM-AD-1 JP-AD-8 JP-AD-9 JP-IN-1 JP-IN-2 JP-IN-3 JP-IN-4 IT-AD-1 IT-AD-2 IT-AD-3 IT-AD-4 AM-AD-2 IT-AD-5 IT-AD-6 FR-AD-1 FR-AD-2 FR-AD-3 FR-AD-4 FR-AD-5 DA-AD-1 DA-AD-2 JP-AD-1 JP-AD-2 JP-AD-3 JP-AD-4 JP-AD-5 JP-AD-6 JP-AD-7);

	print FEATURES join("\t", @samples);
	print FEATURES "\n";

	foreach my $feature (sort keys %Counts) {
		my $name;
		if ($this->type =~ /^refgenome$|^16S$/i) {
			if ($this->tree =~ /bergey/i) {
				$name = $BergeyTree->nodes->{$feature}->name;
			} else {
				$name = $NCBITree->nodes->{$feature}->name;
			}
			if ($name) {
				$name =~ s/\s/_/g;
			}
		}
		$name = $feature unless $name;
		print FEATURES ($name);
		foreach my $sample (@samples) {
			printf FEATURES "\t%.${print_precision}f", ($Counts{$feature}{$sample} || 0);
		}
		print FEATURES "\n";
	}
	close(FEATURES);

	if ($metric ne "jsd" && $metric ne "euclidean" && $metric ne "kld") {
		warn "$this does not understand distance metric $metric";
		return;
	}

	# make the distance matrix

	my $maui_dir = $this->smash->maui_dir;
	my $command = "(echo -n 'dummy	'; cat $feature_file) | cut -f2- | $maui_dir/distanceMatrix --distance=$metric --features=@{[scalar(keys %Counts)]} --samples=@{[scalar(@samples)]} --bootstrap_count=1 --bootstrap_type=rows";
	if ($metric eq "euclidean") {
		$command .= " --rownorm --noise";
	}
	system("$command - > $distance_file");

	return;
}

sub filter_low_occurrences {
	my $this      = shift;
	my $threshold = shift;
	my $Counts    = shift;
	foreach my $feature (keys %$Counts) {
		my $delete = 1;
		foreach my $label (sort {$a cmp $b} keys %{$Counts->{$feature}}) {
			my $frac = ($Counts->{$feature}{$label} || 0);
			if ($frac > $threshold) {
				$delete = 0;
			}
		}
		if ($delete == 1) {
			delete $Counts->{$feature};
		}
	}
}

sub sum_array {
	my $sum = 0;
	foreach (@_) {
		if ($_) {
			$sum += $_;
		}
	}
	return $sum;
}

1;
