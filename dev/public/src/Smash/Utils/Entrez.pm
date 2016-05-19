#!/usr/bin/env perl

package Smash::Utils::Entrez;

use strict;
use warnings;
use Smash::Utils::HTML qw(:all);

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(retrieve_single_gbff_for_query retrieve_gbff_for_query retrieve_contents_for_query retrieve_gbff_for_lsof write_gbff_for_id get_details_for_id get_ids_for_query);

# ===========================================================================
#
# Based on the example script available at NCBI with the following notice
#
# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
#
# Author:  Oleg Khovayko
#
# File Description: eSearch/eFetch calling example
#  
# ===========================================================================

my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";

=head1 NAME

Smash::Utils::Entrez - Utility for communicating with the NCBI Entrez server
using Entrez Programming Utilities and Entrez Eutils.

=head1 SYNOPSIS

	use Smash::Utils::Entrez qw(:all);
	my $database = "genome";
	my $query    = "Bacteria[organism]+AND+WGS";
	my $format   = "gbff"; # gbff, fasta
	retrieve_contents_for_query($query, $database, $format);

=head1 DESCRIPTION

Smash::Utils::Entrez provides several useful functions to communicate with the
NCBI Entrez server. Most of these are queries to specific databases using
search terms or specific ids. Here are some examples:

=head2 Fetching the search results for a specific search term

Suppose you want to see all the whole-genome shotgun (WGS) projects for 
bacterial species. You should search the NCBI B<genome> database using 
the search term B<"bacteria[organism] AND WGS">. To do this programatically,
you would say:

	my @ids = get_ids_for_query("bacteria[organism] AND WGS", "genome");

This would retrieve a list of primary id's of the results in the genome 
database. Note that this id is not the B<Genome Project ID> assigned by
GenBank. That's an entirely different id altogether. If this is what you 
wanted, you are done. If you needed more information, read on.

=head2 Fetching the details of a given entry in the database

Suppose you want to see the details for one (or all) of the results of the search
you performed above. You can get the details for the entry 6994 in the genome
database by saying:

	my $details = get_details_for_id(6994, "genome");

This will retrieve a reference to a hash containing key-value pairs of the details.
If you want specific details, then you should ask for the right key. For example:

	my $tax_id = $details->{TaxId};
	my $accession = $details->{Caption};
	my $update_date = $details->{UpdateDate}

=head2 Retrieving a whole record in GenBank format

Suppose you want to retrieve the whole GenBank record for genome entry 6994. You
would then say:

	write_gbff_for_id(6994, "genome");

This will write the GenBank record to a file called F<679926.NC_014507.gbff>. Here
C<679926> is the NCBI taxonomy id for this genome and C<NC_014597> is the accession
number of this record.

=head2 Retrieving all records for a search term

Given the above examples, it should be trivial to get all the records for a 
search term and write them into individual GenBank files:

	my @ids = get_ids_for_query("bacteria[organism] AND WGS", "genome");
	foreach my $id (@ids) {
		my $filename = write_gbff_for_id($id, "genome");
		print "Wrote $id to $filename\n" if $filename;
	}

=head1 FUNCTIONS

=over 4

=item C<retrieve_gbff_for_lsof($lsof, $db)>

retrieves the GenBank formatted files for the ids in the given list
and writes them all locally

=item C<retrieve_gbff_for_query($query, $db)>

retrieves the GenBank formatted files for the ids that are results for
the given query and writes them all locally

=item C<write_gbff_for_id($id, $db)>

writes the GenBank formatted file for the given id in the given database in
F<E<lt>tax_idE<gt>.E<lt>accessionE<gt>.gbff> and returns the filename.

=item C<get_details_for_id($id, $db)>

gets the following details for the given id in the given database:

	my ($tax_id, $accession, $update_date) = get_details_for_id(6994, "genome");

=item C<get_ids_for_query($query, $db)>

gets the list of ids that are the results of searching the given database
with the given query (B<NOTE:> this has a fixed limit of 10000 id's being 
returned; if your query would retrieve more hits, then you should change the
limit in the code or restrict your query, or use a different method to get
your results).

=cut

sub retrieve_gbff_for_lsof {
	my ($lsof, $db) = @_;
	open(LIST, "<$lsof") || die "Cannot open file $lsof: $!";
	my @ids = <LIST>;
	map {chomp()} @ids;
	close(LIST);
	my @files = ();
	foreach my $id (@ids) {
		my $file = write_gbff_for_id($id, $db);
		push(@files, $file) if $file;
	}
	return @files;
}

sub retrieve_single_gbff_for_query {
	my ($query, $db, $filename) = @_;
	my @ids = get_ids_for_query($query, $db);
	if (@ids) {
		write_single_gbff_for_ids(\@ids, $db, $filename);
	} else {
		return 0;
	}
}

sub write_single_gbff_for_ids {
	use File::Temp;
	my $ids = shift;
	my $db  = shift;
	my $filename = shift;

	my $limit = 400;
	my $count = scalar(@$ids);
	my $term = join(",", @$ids);
	my $written = 0;

	my $tmp = File::Temp->new();
	open(FILE, ">>$filename") || die "Cannot open $filename: $!";
	do {
		my @list = splice(@$ids, 0, $limit);
		$term = join(",", @list);
		my $efetch = "$utils/efetch.fcgi?rettype=gbwithparts&retmode=text&retmax=$limit&db=$db&id=$term";
		if (retrieve_remote_file_persistent($efetch, $tmp->filename)) {
			seek $tmp, 0, 0;
			while(<$tmp>) {
				print FILE $_;
				chomp();
				if (m|^//$|) {
					$written++;
				}
			}
		} else {
			die "Failed to retrieve list into one file!"
		}
	} while(@$ids);
	if ($written != $count) {
		print STDERR "ERROR: $filename: Wrote $written out of $count records\n";
	}
	close(FILE);
}

sub retrieve_contents_for_query {
	my ($query, $db, $format) = @_;
	my @ids = get_ids_for_query($query, $db);
	my @files = ();
	foreach my $id (@ids) {
		my $file = write_contents_for_id($id, $db, $format);
		push(@files, $file) if $file;
	}
	return @files;
}

sub retrieve_gbff_for_query {
	my ($query, $db) = @_;
	my @ids = get_ids_for_query($query, $db);
	my @files = ();
	foreach my $id (@ids) {
		my $file = write_gbff_for_id($id, $db);
		push(@files, $file) if $file;
	}
	return @files;
}

sub get_ids_for_accession {
	my $acc    = shift;
	my $query  = $acc."[ACCN]";
	my $db      = "Nucleotide";
	my $esearch = "$utils/esearch.fcgi?db=$db&term=";

	my $esearch_result = get($esearch . $query);
	$esearch_result =~ m|<Count>(\d+)</Count>|s;
	my $Count    = $1;

	$esearch_result =~ m|<IdList>(.*?)</IdList>|s;
	my $ids    = $1;
	my (@ids) = $ids =~ m|<Id>(.*?)</Id>|s;

	print STDERR "Query: $esearch$query\nNo. results: $Count\nResults: @{[join(',', @ids)]}\n";

	return @ids;
}

sub get_ids_for_query {
	my $query   = shift;
	my $db      = shift;

	$query =~ s/\s+/+/g;

	my $esearch = "$utils/esearch.fcgi?db=$db&retmax=10000&term=";
	my $esearch_result = get_html_page($esearch . $query);
	$esearch_result =~ m|<Count>(\d+)</Count>|s;
	my $Count    = $1;

	$esearch_result =~ m|<IdList>(.*?)</IdList>|s;
	my $ids    = $1;
	my (@ids) = $ids =~ m|<Id>(.*?)</Id>|gs;

	print STDERR "Query: $esearch$query\nNo. results: $Count\nResults: @{[join(',', @ids)]}\n";

	return @ids;
}

sub get_details_for_id {
	my $id    = shift;
	my $db    = shift;
	my $esummary = "$utils/esummary.fcgi?db=$db&id=$id";

	my $details = undef;
	my $Details = {};

	do {
		my $esummary_result = get_html_page($esummary);
		$esummary_result =~ m|<DocSum>(.*?)</DocSum>|s;
		$details    = $1;
	} while (!defined($details));

	my (@items) = $details =~ m|(<Item .*?</Item>)|gs;

	my ($accession, $update_date, $tax_id);
	foreach my $item (@items) {
		my ($key, $value) = $item =~ m|<Item Name="([^"]+)" Type="[^"]+">(.*)</Item>|;
		$Details->{$key} = $value;
	}
	return $Details;
}

sub write_gbff_for_id {
	my $id = shift;
	my $db = shift;

	my $efetch = "$utils/efetch.fcgi?rettype=gbwithparts&retmode=text&retmax=1&db=$db&id=$id";
	print STDERR "EFETCH_QUERY=$efetch\n";     
	my $details = get_details_for_id($id, $db);
	my ($tax_id, $accession, $update_date) = map {$details->{$_}} qw(TaxId Caption UpdateDate);
	my $file_name = "$tax_id.$accession.gbff";
	print "$file_name\t$tax_id\t$accession\t$update_date\n";
	if (retrieve_remote_file_persistent($efetch, "$file_name")) {
		return $file_name;
	} else {
		return;
	}
}

sub write_contents_for_id {
	my $id = shift;
	my $db = shift;
	my $format = shift;

	my $extension = $format;
	my $rettype = $format;

	if ($rettype eq "gbff" || $rettype eq "genbank") {
		$rettype   = "gbwithparts";
		$extension = "gbff";
	}

	my $efetch = "$utils/efetch.fcgi?rettype=$rettype&retmode=text&retmax=1&db=$db&id=$id";
	print STDERR "EFETCH_QUERY=$efetch\n";     
	my $details = get_details_for_id($id, $db);
	my ($tax_id, $accession, $update_date) = map {$details->{$_}} qw(TaxId Caption UpdateDate);
	my $file_name = "$tax_id.$accession.$extension";
	print "$file_name\t$tax_id\t$accession\t$update_date\n";
	if (retrieve_remote_file_persistent($efetch, "$file_name")) {
		return $file_name;
	} else {
		return;
	}
}

sub get_gbff_for_id {
	my $id = shift;
	my $db = shift;

	# ---------------------------------------------------------------------------
	# this area defines a loop which will display $retmax citation results from 
	# Efetch each time the the Enter Key is pressed, after a prompt.

	my $efetch = "$utils/efetch.fcgi?" .
		"rettype=gbwithparts&retmode=text&retmax=1&" .
		"db=$db&id=$id";
	#print STDERR "EF_QUERY=$efetch\n";     

	return get($efetch);
}

sub get_fasta_for_id {
	my $id = shift;
	my $db = "Nucleotide";

	my $efetch = "$utils/efetch.fcgi?" .
		"rettype=fasta&retmode=text&retmax=1&" .
		"db=$db&id=$id";
	#print STDERR "EF_QUERY=$efetch\n";     

	return get($efetch);
}

sub get_contents_for_id {
	my $id     = shift;
	my $db     = shift;
	my $format = shift;

	my $rettype = $format;
	$rettype = "gbwithparts" if $rettype eq "gbff";

	my $efetch = "$utils/efetch.fcgi?" .
		"rettype=$rettype&retmode=text&retmax=1&" .
		"db=$db&id=$id";
	#print STDERR "EF_QUERY=$efetch\n";     

	return get($efetch);
}

1;
