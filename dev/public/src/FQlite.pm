package FQlite;
use strict;
sub new {
	my ($class, $fh) = @_;
	if (ref $fh !~ /GLOB/)
		{die ref $fh, "\n", "FQlite ERROR: expect a GLOB reference\n"}
	my $this = bless {};
	$this->{FH} = $fh;
	while(<$fh>) {last if $_ =~ /\S/} # not supposed to have blanks, but...
	my $firstline = $_;
	if (not defined $firstline) {warn "FQlite: Empty\n"; return $this}
	if ($firstline !~ /^>/) {warn "FQlite: Not FASTA formatted\n"; return $this}
	$this->{LASTLINE} = $firstline;
	chomp $this->{LASTLINE};
	return $this;
}
sub nextEntry {
	my ($this) = @_;
	return 0 if not defined $this->{LASTLINE};
	my $fh = $this->{FH};
	my $def = $this->{LASTLINE};
	my @seq;
	my $lines_read = 0;
	while(<$fh>) {
		$lines_read++;
		if ($_ =~ /^>/) {
			$this->{LASTLINE} = $_;
			chomp $this->{LASTLINE};
			last;
		}
		chomp $_;
		push @seq, $_;
	}
	return 0 if $lines_read == 0;
	chomp @seq;
	my $entry = FQlite::Entry::new($def, \@seq);
	return $entry;
}

package FQlite::Entry;
use overload '""' => 'all';
sub new {
	my ($def, $seqarry) = @_;
	my $this = bless {};
	$this->{DEF} = $def;
	$this->{SEQ} = join(" ", @$seqarry);
	$this->{SEQ} =~ s/\s+/ /g; # just in case more spaces
	$this->{SEQ} =~ s/^\s+//g; # spaces in the beginning
	$this->{SEQ} =~ s/\s+$//g; # spaces in the end
	return $this;
}
sub def {shift->{DEF}}
sub seq {shift->{SEQ}}
sub all {my $e = shift; return $e->{DEF}."\n".$e->{SEQ}."\n"}

1;

__END__

=head1 NAME

FQlite;

=head1 SYNOPSIS

 use FQlite;
 my $fasta = new FQlite(\*STDIN);
 while(my $entry = $fasta->nextEntry) {
     $entry->def;
     $entry->seq;
 }

=head1 DESCRIPTION

FQlite is a package for parsing quality files and databases. The FASTA format is
widely used in bioinformatics. It consists of a definition line followed by
sequence with an arbitrary number of lines and line lengths.

A quality file looks like this:

 >identifier descriptive text
 8 12 16 17 28 34 35 36 35

A FASTA database looks like this:

 >identifier1 some text describing this entry
 8 12 16 17 28 34 35 36 35
45 47 49 42 40 48 45 44 43
 >identifier2 some text describing this entry
 8 19 21 23 29 34 38 39 45
45 48 39 48 42 47 48 43 46

=head2 Object

FQlite has two kinds of objects, the file and the entry.

 my $fasta_file = new FQlite(\*STDIN); # or any other filehandle
 $entry = $fasta_file->nextEntry; # single fasta fle
 while(my $entry = $fasta_file->nextEntry) {
     # canonical form of use for fasta database
 }

The entry has two attributes (def and seq).

 $entry->def; # access the def line
 $entry->seq; # access the sequence
 "$entry";    # overload to fasta file ($entry->def . "\n" . $entry->seq)

=head1 AUTHOR

Mani Arumugam (mozhiyan@cse.wustl.edu)

=head1 BASED ON

FAlite.pm from:
Ian Korf (ikorf@sapiens.wustl.edu, http://sapiens.wustl.edu/~ikorf)

=head1 ACKNOWLEDGEMENTS

FQlite was developed at the Laboratory of Computational Genomics at
Washington University, St. Louis, MO

FAlite was developed at the Genome Sequencing Center at Washington
University, St. Louis, MO.

=head1 COPYRIGHT

Copyright (C) 2006 Mani Arumugam. All Rights Reserved.

FAlite Copyright (C) 1999 Ian Korf. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut
