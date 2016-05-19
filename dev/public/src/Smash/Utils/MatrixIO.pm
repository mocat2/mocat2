package Smash::Utils::MatrixIO;

use strict;
use warnings;
use POSIX;
use Math::Round qw(nearest);

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = (	
		# simple hash

			"read_two_column_hash",
			"write_two_column_hash",
			"write_two_column_hash_sorted",
			"invert_hash",
			"make_distribution_of_keys",
			"make_binning_of_keys",
			"make_distribution_of_values",
			"make_binning_of_values",

		# phylip matrix - symmetric, since it is a self relationship for the same set

			"read_phylip_matrix",
			"write_phylip_matrix",

		# 2-D hash, as R format tab-delimited file, or three column file

			"read_R_matrix",
			"write_R_matrix",
			"write_R_matrix_custom",
			"write_R_matrix_fixed_length",
			"read_three_column_matrix",
			"write_three_column_matrix",
			"write_three_column_matrix_sorted",
			"write_array",

		# arbitrary dimensional hash, as multi column file

			"read_multi_column_matrix",
			"write_multi_column_matrix",

		# 2-D matrix related: $h->{v}->{w} = $x;

			"get_column_labels",
			"t",
			"transpose_matrix",
			"copy_matrix",
			"apply_to_matrix",
			"apply_to_matrix_keys",
			"apply_to_field1",
			"apply_to_field2",
			"multiply_matrices",
			"scalar_multiply_matrix",
			"merge_matrices",
			"fill_matrix",
			"strip_matrix",
			"zero_fill_matrix",
			"zero_strip_matrix",
			"filter_matrix",
			"normalize_cols_by_sum",
			"normalize_rows_by_sum",
			"rowSums",
			"colSums",

		# Arbitry dimensional hash

			"copy_hash",
			"sum_hash",
			"count_hash",
			"average_hash",
			"min_hash",
			"max_hash"
		   );

our %EXPORT_TAGS = ('all' => [@EXPORT_OK]); 

=head1 NAME

Smash::Utils::MatrixIO - Utilities for I/O of matrices stored in hashes.

=head1 SYNOPSIS

	use Smash::Utils::MatrixIO qw(:all);
	my $hash = read_R_matrix("matrix_file.txt");
	print "Matrix[$a,$b] = ".$hash->{$a}->{$b};
	write_R_matrix("new_matrix_file.txt", $hash);

=head1 DESCRIPTION

Smash:Utils::MatrixIO provides functions to read and write matrices stored
as hashes.

=head1 FUNCTIONS

=head2 Matrix/hash read/write methods.

=over 4

=item B<read_R_matrix>

read an R format matrix file into a hash

=item B<write_R_matrix>

write a hash into an R format matrix file

=item B<read_phylip_matrix>

read a phylip format matrix file into a hash

=item B<write_phylip_matrix>

write a hash into a phylip format matrix file

=cut

my %print_format = ('int' => "%d", 'float' => "%.8f", 'double' => "%.8f", 'string' => "%s");

sub read_R_matrix {
	my $file = shift;
	my $FH;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, "<$file") || die "Cannot open $file: $!";
	}
	my $Matrix = {};
	my $header = <$FH>;
	chomp($header);
	my @col_labels = split(/\t/, $header);
	while(<$FH>) {
		chomp();
		my @words = split(/\t/);
		my $row_label = shift(@words);
		map {$Matrix->{$row_label}->{$col_labels[$_]} = $words[$_]} 0..$#words;
	}
	if (ref($file) !~ /GLOB/) {
		close($FH);
	}
	return $Matrix;
}

sub write_R_matrix_custom {
	my ($file, $Matrix, $func) = @_;
	my $FH;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, ">$file") || die "Cannot open $file: $!";
	}
	my @rows = sort keys %$Matrix;

	my %Col;
	foreach my $row (@rows) {
		map {$Col{$_} = 1} keys %{$Matrix->{$row}};
	}
	my @cols = sort keys %Col;

	print $FH join("\t", @cols)."\n";
	foreach my $row (@rows) {
		print $FH $row;
		foreach my $col (@cols) {
			printf $FH "\t%s", $func->($Matrix->{$row}->{$col});
		}
		print $FH "\n";
	}
	if (ref($file) !~ /GLOB/) {
		close($FH);
	}
}

sub write_R_matrix {
	my ($file, $Matrix) = @_;
	write_R_matrix_custom($file, $Matrix, \&to_string);
}

sub write_R_matrix_fixed_length {
	my ($file, $Matrix, $length) = @_;
	my $FH;
	$length = 32 unless $length;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, ">$file") || die "Cannot open $file: $!";
	}
	my @rows = sort keys %$Matrix;

	my %Col;
	foreach my $row (@rows) {
		map {$Col{$_} = 1} keys %{$Matrix->{$row}};
	}
	my @cols = sort keys %Col;

	print $FH join("\t", @cols)."\n";
	foreach my $row (@rows) {
		printf $FH "%-${length}s", substr($row, 0, $length);
		foreach my $col (@cols) {
			printf $FH "\t%s", to_string($Matrix->{$row}->{$col});
		}
		print $FH "\n";
	}
	if (ref($file) !~ /GLOB/) {
		close($FH);
	}
}

# NOTE: phylip matrices are always square matrices!
# So you might see @rows and @cols used interchangably in scalar context
# They are also technically symmetric matrices

sub read_phylip_matrix {
	my $file = shift;
	my $FH;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, "<$file") || die "Cannot open $file: $!";
	}
	my $Matrix = {};
	my $count = <$FH>;
	chomp($count);

	my $lines_read = 0;
	my @col_labels;
	while(<$FH>) {
		chomp();
		my @words = split(/\s+/);
		my $row_label = shift(@words);
		$col_labels[$lines_read] = $row_label;
		map {$Matrix->{$row_label}->{$col_labels[$_]} = $words[$_]} 0..($lines_read-1);
		map {$Matrix->{$col_labels[$_]}->{$row_label} = $words[$_]} 0..($lines_read-1);
		$lines_read++;
	}
	foreach my $label (@col_labels) {
		$Matrix->{$label}->{$label} = 0;
	}
	if (ref($file) !~ /GLOB/) {
		close($FH);
	}
	return $Matrix;
}

sub write_phylip_matrix {
	my ($file, $Matrix) = @_;
	my $FH;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, ">$file") || die "Cannot open $file: $!";
	}
	my @rows = sort keys %$Matrix;

	my %Col;
	foreach my $row (@rows) {
		map {$Col{$_} = 1} keys %{$Matrix->{$row}};
	}
	my @cols = sort keys %Col;

	print $FH scalar(@cols)."\n";
	foreach my $row (@rows) {
		printf $FH "%-10s", $row;
		foreach my $col (@cols) {
			printf $FH " %f", defined($Matrix->{$row}->{$col})?$Matrix->{$row}->{$col}:"0";
		}
		print $FH "\n";
	}
	if (ref($file) !~ /GLOB/) {
		close($FH);
	}
}

=item B<read_three_column_matrix>

reads a three column text file representing a matrix into a hash.
e.g.,

	#sample   species   abundance
	JP-AD-1   Bifidobacterium_longum   0.1234
	JP-AD-1   Bacteroides_caccae       0.0546

will set

	$hash->{JP-AD-1}->{Bifidobacterium_longum} = 0.1234;
	$hash->{JP-AD-1}->{Bacteroides_caccae}     = 0.0546;

=item B<write_three_column_matrix>

writes a matrix in three columns as shown above.

=item B<read_two_column_hash>

reads a two column text file representing a simple hash.

=item B<C<write_two_column_hash>>

writes a simple hash as a two column text file.

=back

=cut

sub read_three_column_matrix {
	my $file = shift;
	my $FH;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, "<$file") || die "Cannot open $file: $!";
	}
	my $Matrix = {};
	while(<$FH>) {
		chomp();
		if (!m/^#/) {
			my ($row_label, $col_label, $value) = split(/\t/);
			$Matrix->{$row_label}->{$col_label} = $value;
		}
	}
	if (ref($file) !~ /GLOB/) {
		close($FH);
	}
	return $Matrix;
}

sub write_three_column_matrix {
	my ($file, $Matrix) = @_;
	my $FH;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, ">$file") || die "Cannot open $file: $!";
	}
	my @rows = sort keys %$Matrix;
	foreach my $row (@rows) {
		while (my ($col, $value) = each %{$Matrix->{$row}}) {
			printf $FH ("%s\t%s\t%s\n", $row, $col, to_string($value)) if $value;
		}
	}
	if (ref($file) !~ /GLOB/) {
		close($FH);
	}
}

sub write_three_column_matrix_sorted {
	my ($file, $Matrix) = @_;
	my $FH;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, ">$file") || die "Cannot open $file: $!";
	}
	my @rows = sort {$a <=> $b} keys %$Matrix;
	foreach my $row (@rows) {
		while (my ($col, $value) = each %{$Matrix->{$row}}) {
			printf $FH ("%s\t%s\t%s\n", $row, $col, to_string($value)) if $value;
		}
	}
	if (ref($file) !~ /GLOB/) {
		close($FH);
	}
}

=head2 Multilevel hash read/write and manipulation methods

=over 4

=item B<read_multi_column_matrix>

reads a multiple column text file into a multi-level hash/matrix. For example, a line
such as:

	Root	Bacteria	Bacteroidetes	Bacteroides	1

will execute something to the effect of:

	$h->{Root}->{Bacteria}->{Bacteroidetes}->{Bacteroides} = 1;


=item B<write_multi_column_matrix>

writes a multi-level hash as multi-column text file as shown above.

=item B<print_hash($FH, $hash)>

writes a multi-level hash to $FH as a multi-column text file.

=item B<count_hash>

count the total number of elements (leaf nodes if you will) in the
multi-level hash.

=item B<sum_hash>

sums all the numbers at the leaf nodes.

=item B<average_hash>

gets the average of all the numbers at the leaf nodes.

=item B<max_hash>

gets the maximum value stored at the leaf nodes.

=item B<min_hash>

gets the minimum value stored at the leaf nodes.

=back

=cut

sub read_multi_column_matrix {
	my $file = shift;
	my $FH;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, "<$file") || die "Cannot open $file: $!";
	}
	my $Matrix = {};
	while(<$FH>) {
		chomp();
		if (!m/^#/) {

			# Beware of the complexity of this piece of code.
			# Since we dont know how many levels there are, we use a loop
			# to construct the assignment code fragment and then
			# evaluate it using eval. The first part is in single quotes,
			# since the name of the variable should be literal. But the later
			# parts are double quoted, since keys and values must be evaluated
			# to their values.

			my @words = split(/\t/);
			my $value = pop(@words);
			my $code = '$Matrix';
			foreach my $key (@words) {
				$code .= "->{\"$key\"}";
			}
			$code .= " = \"$value\";";
			eval $code;
			#eval 'print $code;';
		}
	}
	if (ref($file) !~ /GLOB/) {
		close($FH);
	}
	return $Matrix;
}

sub write_multi_column_matrix {
	my ($file, $Matrix) = @_;
	my $FH;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, ">$file") || die "Cannot open $file: $!";
	}
	print_hash($FH, $Matrix);
	if (ref($file) !~ /GLOB/) {
		close($FH);
	}
}

sub print_hash {
	my ($FH, $hash) = @_;
	recursive_print_hash($FH, "", $hash);
}

# this is a recursive function that traverses through a hash and prints a tab-delimited multi-column
# text file representing the hash.
# this must be called from the top level as:
# recursive_print_hash($FH, "", $hash);
# which is what print_hash does. So use print_hash and not recursive_print_hash directly.

sub recursive_print_hash {
	my ($FH, $prefix, $hash) = @_;
	if (ref($hash) =~ /HASH/) {
		while (my ($key, $val) = each %$hash) {
			recursive_print_hash($FH, "$prefix$key\t", $val);
		}
	} else {
		printf $FH "%s%s\n", $prefix, to_string($hash);
	}
}

sub average_hash {
	my $hash = shift;
	return sum_hash($hash)/count_hash($hash);
}

sub sum_hash {
	my $hash = shift;
	return recursive_sum_hash(0, $hash);
}

sub recursive_sum_hash {
	my ($sum, $hash) = @_;
	if (ref($hash) =~ /HASH/) {
		while (my ($key, $val) = each %$hash) {
			$sum = recursive_sum_hash($sum, $val);
		}
		return $sum;
	} else {
		return $sum + $hash;
	}
}

sub count_hash {
	my $hash = shift;
	return recursive_count_hash(0, $hash);
}

sub recursive_count_hash {
	my ($count, $hash) = @_;
	if (ref($hash) =~ /HASH/) {
		while (my ($key, $val) = each %$hash) {
			$count = recursive_count_hash($count, $val);
		}
		return $count;
	} elsif (defined $hash) {
		return $count + 1;
	} else {
		return $count;
	}
}

sub min2 {
	my ($x, $y) = @_;
	[$x => $y]->[$x >= $y];
}

sub max2 {
	my ($x, $y) = @_;
	[$x => $y]->[$x <= $y];
}

sub min_hash {
	my $hash = shift;
	return recursive_func_hash(INT_MAX, $hash, \&min2, sub {shift});
}

sub max_hash {
	my $hash = shift;
	return recursive_func_hash(INT_MIN, $hash, \&max2, sub {shift});
}

sub recursive_func_hash {
	my ($result, $hash, $func_leaf_defined, $func_leaf_undef) = @_;
	if (ref($hash) =~ /HASH/) {
		while (my ($key, $val) = each %$hash) {
			$result = recursive_func_hash($result, $val, $func_leaf_defined, $func_leaf_undef);
		}
		return $result;
	} elsif (defined $hash) {
		return $func_leaf_defined->($result, $hash);
	} else {
		return $func_leaf_undef->($result, $hash);
	}
}

# Simple hash

sub read_two_column_hash {
	my $file = shift;
	my $FH;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, "<$file") || die "Cannot open $file: $!";
	}
	my $hash = {};
	while(<$FH>) {
		chomp();
		if (!m/^#/) {
			my ($key, $value) = split(/\t/);
			$hash->{$key} = $value || undef;
		}
	}
	if (ref($file) !~ /GLOB/) {
		close($FH);
	}
	return $hash;
}

sub write_two_column_hash_fixed_length {
	my ($file, $hash, $length) = @_;
	my $FH;
	$length = 32 unless $length;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, ">$file") || die "Cannot open $file: $!";
	}
	while (my ($key, $value) = each %$hash) {
		printf $FH "%-${length}s\t%s\n", substr($key, 0, $length), to_string($value);
	}
	if (ref($file) !~ /GLOB/) {
		close($FH);
	}
}

sub write_two_column_hash {
	my ($file, $hash) = @_;
	my $FH;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, ">$file") || die "Cannot open $file: $!";
	}
	while (my ($key, $value) = each %$hash) {
		printf $FH "%s\t%s\n", $key, to_string($value);
	}
	if (ref($file) !~ /GLOB/) {
		close($FH);
	}
}

sub write_two_column_hash_sorted {
	my ($file, $hash) = @_;
	my $FH;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, ">$file") || die "Cannot open $file: $!";
	}
	foreach my $key (sort {$a <=> $b} keys %$hash) {
		my $value = $hash->{$key};
		printf $FH "%s\t%s\n", $key, to_string($value);
	}
	if (ref($file) !~ /GLOB/) {
		close($FH);
	}
}

sub invert_hash {
	my $hash = shift;
	my $inv  = {};
	while (my ($key, $value) = each %$hash) {
		$inv->{$value} = $key;
	}
	return $inv;
}

sub make_binning_of_keys {
	my $hash = shift;
	my $precision = shift || 1;
	my $Distribution = {};
	my $count = 0;
	foreach my $key (keys %$hash) {
		my $nearest = nearest($precision, $key);
		$Distribution->{$nearest}++;
		$count++;
	}
}

sub make_distribution_of_keys {
	my $hash = shift;
	my $precision = shift || 1;
	my $Distribution = {};
	my $count = 0;
	foreach my $key (keys %$hash) {
		my $nearest = nearest($precision, $key);
		$Distribution->{$nearest}++;
		$count++;
	}

	# normalize

	my @keys = keys %$Distribution;
	foreach my $key (@keys) {
		$Distribution->{$key} /= $count;
	}
	
	return $Distribution;
}

sub make_binning_of_values {
	my $hash = shift;
	my $precision = shift || 1;
	my $result = {};
	recursive_binning_of_values($result, $hash, $precision);
	return $result;
}

sub recursive_binning_of_values {
	my ($result, $hash, $precision) = @_;
	if (ref($hash) =~ /HASH/) {
		while (my ($key, $val) = each %$hash) {
			recursive_binning_of_values($result, $val, $precision);
		}
		return;
	} else {
		my $nearest = nearest($precision, $hash);
		$result->{$nearest}++;
	}
}

sub make_binning_of_values_old {
	my $hash = shift;
	my $precision = shift || 1;
	my $Binning = {};
	foreach my $value (values %$hash) {
		my $nearest = nearest($precision, $value);
		$Binning->{$nearest}++;
	}

	return $Binning;
}

sub make_distribution_of_values {
	my $hash = shift;
	my $precision = shift || 1;
	my $Binning = make_binning_of_values($hash, $precision);

	# normalize

	my @keys = keys %$Binning;
	my $count = sum_hash($Binning);

	my $Distribution = {};
	foreach my $key (@keys) {
		$Distribution->{$key} = $Binning->{$key}/$count;
	}
	
	return $Distribution;
}

=head2 Simple array manipulation methods

=over 4

=item B<write_array($array_ref)>

=back

=cut

sub write_array {
	my ($file, $array) = @_;
	my $FH;
	if (ref($file) =~ /GLOB/) {
		$FH = $file;
	} else {
		open($FH, ">$file") || die "Cannot open $file: $!";
	}
	foreach my $e (@$array) {
		print $FH "$e\n";
	}
	if (ref($file) !~ /GLOB/) {
		close($FH);
	}
}

=head2 2-D Matrix manipulation methods

=over 4

=item B<transpose_matrix($m)>

returns the transpose matrix of $m

=item B<t($m)>

short-cut to C<transpose_matrix($m)>

=item B<multiply_matrices($a, $b)>

returns the matrix product of $a and $b.

=item B<scalar_multiply_matrix($h, $scalar)>

multiplies the matrix (passed by ref) by the given scalar

=item B<merge_matrices($a, $b)>

merges the matrices $a and $b. Rows and columns in either one are present in the 
merged matrix. 

=item B<zero_fill_matrix($m)>

fill the undefined cells with 0.

=item B<zero_strip_matrix($m)>

removes all cells of the matrix with value 0. If a full row (or column) has 0's
all over and they were all removed, the row (or column) will also be 
removed.

=item B<filter_matrix($m, $min_value)>

removes all cells of the matrix B<BELOW> C<$min_value>. If a full row (or column)
has been removed, the row (or column) will also be removed.

=item B<normalize_cols_by_sum($m)>

normalizes each value by the column sum.

=item B<normalize_rows_by_sum($m)>

normalizes each value by the row sum.

=item B<merge_matrices($a, $b, ...)>

merges all the matrices and returns a new matrix with values from all.

=back

=cut

# Apply the given function to every element in the 2D matrix

sub apply_to_matrix {
	my $hash = shift;
	my $func = shift;
	my $copy  = {};
	foreach my $key1 (keys %$hash) {
		foreach my $key2 (keys %{$hash->{$key1}}) {
			$copy->{$key1}->{$key2} = $func->($hash->{$key1}->{$key2});
		}
	}
	return $copy;
}

# Apply the given function to every element in the 2D matrix
# But also pass $row, $col along with the value

sub apply_to_matrix_keys {
	my $hash = shift;
	my $func = shift;
	my $copy  = {};
	foreach my $key1 (keys %$hash) {
		foreach my $key2 (keys %{$hash->{$key1}}) {
			$copy->{$key1}->{$key2} = $func->($key1, $key2, $hash->{$key1}->{$key2});
		}
	}
	return $copy;
}

# Apply the given function to key1 in the 2D matrix

sub apply_to_field1 {
	my $hash = shift;
	my $func = shift;
	my $copy  = {};
	foreach my $key1 (keys %$hash) {
		my $newkey1 = $func->($key1);
		foreach my $key2 (keys %{$hash->{$key1}}) {
			$copy->{$newkey1}->{$key2} = $hash->{$key1}->{$key2};
		}
	}
	return $copy;
}

# Apply the given function to key2 in the 2D matrix

sub apply_to_field2 {
	my $hash = shift;
	my $func = shift;
	my $tra  = transpose_matrix($hash);
	apply_to_field1($tra, $func);
	return transpose_matrix($tra);
}

# copy the matrix and leave the original untouched

sub copy_matrix {
	my $hash = shift;
	my $copy  = {};
	foreach my $key1 (keys %$hash) {
		foreach my $key2 (keys %{$hash->{$key1}}) {
			$copy->{$key1}->{$key2} = $hash->{$key1}->{$key2};
		}
	}
	return $copy;
}

sub t {
	return transpose_matrix(shift);
}

sub transpose_matrix {
	my $hash = shift;
	my $tra  = {};
	foreach my $key1 (keys %$hash) {
		foreach my $key2 (keys %{$hash->{$key1}}) {
			$tra->{$key2}->{$key1} = $hash->{$key1}->{$key2};
		}
	}
	return $tra;
}

sub multiply_matrices {
	my ($M1, $M2) = @_;
	my $product = {};

	my @common = sort keys %$M2;
	my @rows   = sort keys %$M1;
	my %Col;
	foreach my $row (@common) {
		map {$Col{$_} = 1} keys %{$M2->{$row}};
	}
	my @cols = sort keys %Col;

	foreach my $row (@rows) {
		foreach my $col (@cols) {
			my $sum = 0;
			foreach my $x (@common) {
die "M1:$row, $x!" unless defined $M1->{$row}->{$x};
die "M2:$x, $col!" unless defined $M2->{$x}->{$col};

				$sum += ($M1->{$row}->{$x} * $M2->{$x}->{$col});
			}
			$product->{$row}->{$col} = $sum;
		}
	}
	return $product;
}

sub merge_matrices {
	my $merged = {};
	foreach my $M (@_) {
		foreach my $row (keys %$M) {
			foreach my $col (keys %{$M->{$row}}) {
				$merged->{$row}->{$col} = $M->{$row}->{$col};
			}
		}
	}
	return $merged;
}
 
sub fill_matrix {
	my $Matrix = shift;
	my $value  = shift;

	my @rows = sort keys %$Matrix;
	my %Col;
	foreach my $row (@rows) {
		map {$Col{$_} = 1} keys %{$Matrix->{$row}};
	}
	my @cols = sort keys %Col;

	foreach my $row (@rows) {
		foreach my $col (@cols) {
			if (!defined($Matrix->{$row}->{$col})) {
				$Matrix->{$row}->{$col} = $value;
			}
		}
	}
}
 
sub zero_fill_matrix {
	my $Matrix = shift;

	my @rows = sort keys %$Matrix;
	my %Col;
	foreach my $row (@rows) {
		map {$Col{$_} = 1} keys %{$Matrix->{$row}};
	}
	my @cols = sort keys %Col;

	foreach my $row (@rows) {
		foreach my $col (@cols) {
			if (!defined($Matrix->{$row}->{$col})) {
				$Matrix->{$row}->{$col} = 0;
			}
		}
	}
}

sub strip_matrix {
	my $Matrix = shift;
	my $item   = shift;
	strip_matrix_row($Matrix, $item);
}

sub strip_matrix_row {
	my $Matrix = shift;
	my $item   = shift;
	my @rows = sort keys %$Matrix;
	my %Col;
	foreach my $row (@rows) {
		map {$Col{$_} = 1} keys %{$Matrix->{$row}};
	}
	my @cols = sort keys %Col;

	foreach my $row (@rows) {
		my $keep_row = 0;
		foreach my $col (@cols) {
			if (defined($Matrix->{$row}->{$col}) && $Matrix->{$row}->{$col} eq $item) {
				delete $Matrix->{$row}->{$col};
			} else {
				$keep_row = 1;
			}
		}
		delete $Matrix->{$row} unless $keep_row;
	}
}

sub zero_strip_matrix {
	my $Matrix = shift;
	zero_strip_matrix_row($Matrix);
}

sub zero_strip_matrix_row {
	my $Matrix = shift;
	my @rows = sort keys %$Matrix;
	my %Col;
	foreach my $row (@rows) {
		map {$Col{$_} = 1} keys %{$Matrix->{$row}};
	}
	my @cols = sort keys %Col;

	foreach my $row (@rows) {
		my $keep_row = 0;
		foreach my $col (@cols) {
			if (defined($Matrix->{$row}->{$col}) && $Matrix->{$row}->{$col} == 0) {
				delete $Matrix->{$row}->{$col};
			} else {
				$keep_row = 1;
			}
		}
		delete $Matrix->{$row} unless $keep_row;
	}
}

sub filter_matrix {
	my $Matrix = shift;
	my $min    = shift;

	my @rows = sort keys %$Matrix;
	my %Col;
	foreach my $row (@rows) {
		map {$Col{$_} = 1} keys %{$Matrix->{$row}};
	}
	my @cols = sort keys %Col;

	foreach my $row (@rows) {
		my $keep_row = 0;
		foreach my $col (@cols) {
			if (defined($Matrix->{$row}->{$col}) && $Matrix->{$row}->{$col} < $min) {
				delete $Matrix->{$row}->{$col};
			} else {
				$keep_row = 1;
			}
		}
		delete $Matrix->{$row} unless $keep_row;
	}
}

sub scalar_multiply_matrix {
	my $Matrix = shift;
	my $scalar = shift;

	my @rows = sort keys %$Matrix;
	my %Col;
	foreach my $row (@rows) {
		map {$Col{$_} = 1} keys %{$Matrix->{$row}};
	}
	my @cols = sort keys %Col;

	foreach my $row (@rows) {
		foreach my $col (@cols) {
			if (defined($Matrix->{$row}->{$col})) {
				$Matrix->{$row}->{$col} *= $scalar;
			}
		}
	}
	return 0;
}

sub colSums {
	my $Matrix = shift;

	my @rows = sort keys %$Matrix;
	my %Col;
	foreach my $row (@rows) {
		map {$Col{$_} = 1} keys %{$Matrix->{$row}};
	}
	my @cols = sort keys %Col;

	my $Sums = {};
	foreach my $col (@cols) {
		my $sum = 0;
		foreach my $row (@rows) {
			$sum += ($Matrix->{$row}->{$col} || 0);
		}
		$Sums->{$col} = $sum;
	}
	return $Sums;
}

sub rowSums {
	my $Matrix = shift;

	my @rows = sort keys %$Matrix;
	my %Col;
	foreach my $row (@rows) {
		map {$Col{$_} = 1} keys %{$Matrix->{$row}};
	}
	my @cols = sort keys %Col;

	my $Sums = {};
	foreach my $row (@rows) {
		my $sum = 0;
		foreach my $col (@cols) {
			$sum += ($Matrix->{$row}->{$col} || 0);
		}
		$Sums->{$row} = $sum;
	}
	return $Sums;
}

sub normalize_cols_by_sum {
	my $Matrix = shift;

	my @rows = sort keys %$Matrix;
	my %Col;
	foreach my $row (@rows) {
		map {$Col{$_} = 1} keys %{$Matrix->{$row}};
	}
	my @cols = sort keys %Col;

	my $new_matrix = {};
	foreach my $col (@cols) {
		my $sum = 0;
		foreach my $row (@rows) {
			$sum += ($Matrix->{$row}->{$col} || 0);
		}
		foreach my $row (@rows) {
			if ($Matrix->{$row}->{$col}) {
				$new_matrix->{$row}->{$col} = $Matrix->{$row}->{$col} / $sum;
			} else {
				$new_matrix->{$row}->{$col} = 0;
			}
		}
	}
	return $new_matrix;
}

sub normalize_rows_by_sum {
	my $Matrix = shift;

	my @rows = sort keys %$Matrix;
	my %Col;
	foreach my $row (@rows) {
		map {$Col{$_} = 1} keys %{$Matrix->{$row}};
	}
	my @cols = sort keys %Col;

	my $new_matrix = {};
	foreach my $row (@rows) {
		my $sum = 0;
		foreach my $col (@cols) {
			$sum += ($Matrix->{$row}->{$col} || 0);
		}
		foreach my $col (@cols) {
			$new_matrix->{$row}->{$col} = ($Matrix->{$row}->{$col} || 0) / $sum;
		}
	}
	return $new_matrix;
}

=head2 Utility methods

=over 4

=item B<get_column_labels>

=item B<to_string>

=item B<typeof>

=back

=cut

sub get_column_labels {
	my $Matrix = shift;
	my @rows = sort keys %$Matrix;
	my %Col;
	foreach my $row (@rows) {
		map {$Col{$_} = 1} keys %{$Matrix->{$row}};
	}
	return [sort keys %Col];
}

sub to_string {
	my $x = shift;
	my $type = typeof($x);
	if ($type eq "int") {
		return "$x";
	} elsif ($type eq "string") {
		return $x;
	} elsif ($type eq "float" || $type eq "double") {
		return sprintf("%.8f", $x);
	} else {
		return $x;
	}
}

sub typeof {
	my $val = shift;
	if (!defined($val)) {
		return 'null';
	} elsif (!ref($val)) {
		if ($val =~ /^-?\d+$/) {
			return 'int';
		} elsif ($val =~ /^-?\d+(\.\d+)?$/) {
			return 'float';
		} elsif ($val =~ /^-?\d(\.\d+)?e[\-\+]\d+$/) {
			return 'double';
		} else {
			return 'string';
		}
	} else {
		my $type = ref($val);
		if ($type eq 'HASH' || $type eq 'ARRAY') {
			return lc($type);
		} elsif ( $type eq 'CODE' || $type eq 'REF' || $type eq 'GLOB' || $type eq 'LVALUE' ) {
			return $type;
		} else {
			return 'object';
		}
	}
}

1;
