package Smash::Analyses::Resampler;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(sample_with_replacement sample_without_replacement downsample_reads fisher_yates_shuffle fisher_yates_shuffle_in_place);
our %EXPORT_TAGS = ('all' => [@EXPORT_OK]);

=head1 NAME

Smash::Analyses::Resampler - Module for resampling data in a list

=head1 DESCRIPTION

Smash::Analyses::Resampler provides functions to resample a list. This is useful
in generating bootstrap replicates for many analyses.

=head1 FUNCTIONS

=over 4

=item C<downsample_reads>

=cut

=item C<sample_with_replacement>

given a reference to a list with C<N> elements, returns reference to a list of C<N> elements 
after sampling the input with replacement.

=item C<sample_without_replacement>

given a reference to a list with C<N> elements, returns reference to a list of C<N> elements 
after shuffling the elements using Fisher Yates shuffling method.

=back

=cut

sub sample_with_replacement {
	my $in_list = shift;
	my $n    = scalar(@$in_list);
	my $list = [map {$in_list->[int(rand($n))]} 0..($n-1)];
	return $list;
}

sub sample_without_replacement {
	return fisher_yates_shuffle(shift);
}

sub fisher_yates_shuffle {
	my $in_list = shift;
	my $n    = scalar(@$in_list);
	my $list = [map {$_} @$in_list];
	while ($n > 0) {
		$n--;
		my $k = int(rand($n+1)); #(0<=k<=n)
		if ($n != $k) {
			($list->[$k], $list->[$n]) = ($list->[$n], $list->[$k]);
		}
	}
	return $list;
}

sub fisher_yates_shuffle_in_place {
	my $list = shift;
	my $n    = scalar(@$list);
	while ($n > 0) {
		my $k = int(rand($n)); #(0<=k<n), which is good since this is 0-based
		$n--;
		if ($n != $k) {
			($list->[$k], $list->[$n]) = ($list->[$n], $list->[$k]);
		}
	}
}

1;
