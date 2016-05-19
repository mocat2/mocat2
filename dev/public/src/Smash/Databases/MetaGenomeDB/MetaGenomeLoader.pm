package Smash::Databases::MetaGenomeDB::MetaGenomeLoader;

use strict;
use warnings;
use File::Path;
use Smash::Core;

our @ISA = qw(Smash::Databases::MetaGenomeDB::Loader);

my %FastaInfo;      # Key for %FastaInfo is the first word for regular, last word for weird_fasta
my %ReverseLookup;  # A way to look up the changes. Useful to lookup original data from modified data
                    # Normally, the modified name is, 
                    # $FastaInfo{$name}{ID} = "$metagenome.$name", so
                    # $ReverseLookup{$FastaInfo{$name}{ID}} = $name
my %LibraryInfo;
my %Trace;          # Used only in sanger data
my %ReadRemains;
my $weird_fasta_header; # Last word is the actual def in these weird cases
my $chop_fasta_header;
my $PROGRESS;
my $quality_trim;   # since XML parsers cannot get quality trim from the $this object, I store it here.

############################################
############################################
##    Member variables through methods    ##
############################################
############################################

sub type          {shift->{TYPE}}
sub quality_trim  {shift->{QUALITY_TRIM}}

############################################
############################################
##    Pipeline functions                  ##
############################################
############################################

#########################################
# Init the object
# Right now it does:
#	. SUPER::init()
#	. USER_OPTIONS = ""
#########################################

sub init {
	my $this = shift;
	$this->parse_config();
	#my $metagenome = $this->metagenome;
	$this->SUPER::init();
	$chop_fasta_header = 1; # always take the first word of the fasta header
	$weird_fasta_header = ($this->{WEIRD_FASTA_HEADER} || 0);
	$PROGRESS = \*STDERR;
	select($PROGRESS); $| = 1; select(STDOUT);
	$quality_trim = $this->quality_trim;
}

############################################
# Load!
############################################

sub readFasta {
	my $this = shift;
	my @fasta_files = @_;
	foreach my $file (@fasta_files) {
		open(FASTA, "<$file") || die "Cannot open $file: $!";
		my $fasta = new FAlite(\*FASTA);
		while (my $entry = $fasta->nextEntry) {
			my $def = $entry->def;
			$def =~ s/^>//;
			my $length = length($entry->seq);
			my @words = split(/\s/, $def);
			my $name = $words[0];
			if ($weird_fasta_header) {
				$name = $words[$#words];
			}
			$FastaInfo{$name}{DEF} = $def;
			$FastaInfo{$name}{LEN} = $length;
		}
	}
	return 1;
}

# Trace archive XML handling
{
	my $this; # reference to the object
	my %trace;
	my $tag;
	my $value;
	my $count = 0;
	my $in_data = 0;
	my $curr_data = "";
	my $metagenome;
	my $sample_id;
	my $remap_xml_fh;
	my $SmashLibraryIds = {};

	sub readXML {
		use XML::Parser;
		$this = shift;
		$sample_id = shift; # sample_id is set here, and stays
		my @xml_files = @_;
		$metagenome = $this->metagenome;

		my $parser = new XML::Parser;
		$parser->setHandlers(
					XMLDecl => \&xmlDeclaration,
					Start   => \&startElement, 
					End     => \&endElement, 
					Default => \&default
				);
		$count = 0;
		foreach my $file (@xml_files) {
			$parser->parsefile($file);
		}
		return 1;
	}

	sub remapXML {
		use XML::Parser;
		$this         = shift;
		$remap_xml_fh = shift;
		my @xml_files = @_;
		$metagenome   = $this->metagenome;

		my $parser = new XML::Parser;
		$parser->setHandlers(
					XMLDecl => \&xmlDeclaration,
					Start   => \&startElement, 
					End     => \&endElementRemap, 
					Default => \&default
				);
		$count = 0;
		foreach my $file (@xml_files) {
			$parser->parsefile($file);
		}
		return 1;
	}

	sub filterXML {
		use XML::Parser;
		$this         = shift;
		$remap_xml_fh = shift;
		my @xml_files = @_;
		$metagenome   = $this->metagenome;

		my $parser = new XML::Parser;
		$parser->setHandlers(
					XMLDecl => \&xmlDeclaration,
					Start   => \&startElement, 
					End     => \&endElementFilter, 
					Default => \&default
				);
		$count = 0;
		foreach my $file (@xml_files) {
			$parser->parsefile($file);
		}
		return 1;
	}

	sub xmlDeclaration { }

	sub startElement {
		my ($parseinst, $element, %attrs ) = @_;
		my $lelement = lc($element);
		if ($lelement eq "trace") {
			%trace = ();
		} elsif ($lelement eq "trace_volume") {
		} else {
			$tag = $element;
		}
	}

	sub endElement {
		my ($parseinst, $element ) = @_;
		my $lelement = lc($element);
		if ($lelement eq "trace") {
			$count++;
			if (($count+1)%10000 == 0) {
				print $PROGRESS $this->progress_bar($count+1, 10000);
			}

			my $trace_name  = $trace{"TRACE_NAME"};
			if (defined($FastaInfo{$trace_name})) {
				my ($library_id, $template_id, $trace_end, $clip_left, $clip_right, $insert_size, $insert_stdev);
				($library_id = $trace{"LIBRARY_ID"} || $trace{"SEQ_LIB_ID"}) || die "LIBRARY_ID not defined";
				$insert_size = $trace{"INSERT_SIZE"} || -1;
				$insert_stdev= $trace{"INSERT_STDEV"} || -1;
				$template_id = $trace{"TEMPLATE_ID"};
				$trace_end   = $trace{"TRACE_END"};
				if ($trace_end && $trace_end =~ /(f|fwd|forward)/i) {
					$trace_end = "F";
				} else {
					$trace_end = "R";
				}

				##########################
				# Process LIBRARY_ID
				##########################

				my $smash_library_id = $SmashLibraryIds->{$library_id}->{$sample_id}->{'sanger'}->{$insert_size}->{$insert_stdev};
				if (!$smash_library_id) {
					my ($smash_insert_size, $smash_insert_stdev) = ($insert_size, $insert_stdev);
					$smash_insert_size  = undef if ($smash_insert_size  == -1);
					$smash_insert_stdev = undef if ($smash_insert_stdev == -1);
					$smash_library_id = $this->get_library_id($library_id, $sample_id, 'sanger', $smash_insert_size, $smash_insert_stdev);
					$SmashLibraryIds->{$library_id}->{$sample_id}->{'sanger'}->{$insert_size}->{$insert_stdev} = $smash_library_id;
				}

				($insert_size != -1)  && ($LibraryInfo{$smash_library_id}{"INSERT_SIZE"}  = $insert_size);
				($insert_stdev != -1) && ($LibraryInfo{$smash_library_id}{"INSERT_STDEV"} = $insert_stdev);

				##########################
				# Process clipping coordinates
				##########################

				my ($val1, $val2);
				$val1       = $trace{"CLIP_VECTOR_LEFT"};
				$val2       = $trace{"CLIP_QUALITY_LEFT"};
				$clip_left  = $trace{"CLIP_LEFT"}  || ($val1 && $val2 && (($val1>$val2)?$val1:$val2)) || $val1 || $val2;
				$val1       = $trace{"CLIP_VECTOR_RIGHT"};
				$val2       = $trace{"CLIP_QUALITY_RIGHT"};
				$clip_right = $trace{"CLIP_RIGHT"} || ($val1 && $val2 && (($val1<$val2)?$val1:$val2)) || $val1 || $val2;

				##########################
				# Store Trace info
				##########################

				$Trace{$trace_name}{TEMPLATE_ID} = $template_id;
				$Trace{$trace_name}{TRACE_END}   = $trace_end;
				$Trace{$trace_name}{LIB}         = $smash_library_id;
				$Trace{$trace_name}{LEFT}        = $clip_left;
				$Trace{$trace_name}{RIGHT}       = $clip_right;
				$Trace{$trace_name}{CLIP_VECTOR_LEFT}  = $trace{"CLIP_VECTOR_LEFT"};
				$Trace{$trace_name}{CLIP_VECTOR_RIGHT} = $trace{"CLIP_VECTOR_RIGHT"};
			}

		} elsif ($lelement eq "trace_volume") {
			print $PROGRESS " done (Processed $count traces.)";
		} else {
			$trace{uc($tag)} = $value;
		}
		$in_data = 0;
	}

	sub endElementRemap {
		my ($parseinst, $element ) = @_;
		my $lelement = lc($element);
		if ($lelement eq "trace") {
			my $name = $trace{TRACE_NAME};
			if (defined($FastaInfo{$name})) {
				if ($quality_trim) {
					delete $trace{CLIP_LEFT};
					delete $trace{CLIP_VECTOR_LEFT};
					delete $trace{CLIP_QUALITY_LEFT};
					delete $trace{CLIP_RIGHT};
					delete $trace{CLIP_VECTOR_RIGHT};
					delete $trace{CLIP_QUALITY_RIGHT};
				}
				delete $trace{SEQ_LIB_ID};
				$trace{LIBRARY_ID} = $Trace{$name}{LIB};
				$trace{TRACE_NAME} = $FastaInfo{$name}{ID};
				write_trace($remap_xml_fh, %trace);
			}
		} elsif ($lelement eq "trace_volume") {
		} else {
			$trace{uc($tag)} = $value;
		}
		$in_data = 0;
	}

	sub endElementFilter {
		my ($parseinst, $element ) = @_;
		my $lelement = lc($element);
		if ($lelement eq "trace") {
			my $name = $trace{TRACE_NAME};
			delete $trace{SEQ_LIB_ID};
			if ($quality_trim) {
				# All the clipping information goes, since we have trimmed the
				# sequences and these coordinates are not valid anymore
				delete $trace{CLIP_LEFT};
				delete $trace{CLIP_QUALITY_LEFT};
				delete $trace{CLIP_VECTOR_LEFT};
				delete $trace{CLIP_RIGHT};
				delete $trace{CLIP_QUALITY_RIGHT};
				delete $trace{CLIP_VECTOR_RIGHT};
			}
			if (defined($ReadRemains{$ReverseLookup{$name}})) {
				write_trace($remap_xml_fh, %trace);
			}
		} elsif ($lelement eq "trace_volume") {
		} else {
			$trace{uc($tag)} = $value;
		}
		$in_data = 0;
	}

	sub default {
		my ($parseint, $data) = @_;
		if ($data !~ /\S/) { return;}
		if ($in_data == 1) {
			$curr_data .= $data;
		} else {
			$curr_data  = $data;
		}
		$value = $curr_data;
		$in_data = 1;
	}

	sub write_trace {
		my ($fh, %trace) = @_;
		print $fh "\t<TRACE>\n";
		foreach my $tag (sort keys %trace) {
			if (defined($trace{$tag})) {
				printf $fh "\t\t<%s>%s</%s>\n", $tag, $trace{$tag}, $tag;
			}
		}
		print $fh "\t</TRACE>\n";
	}
}

############################################
# Load reads and xml information into db
############################################

sub run {
	use File::Path;
	use FAlite;
	use FQlite;
	my $this = shift;

	my $metagenome  = $this->metagenome;

	if (!$this->get_metagenome_label($metagenome)) {
		die "Metagenome $metagenome does not exist!\nPlease add $metagenome to the master table.\n";
	}

	my $read_dir    = $this->read_dir($metagenome);
	my $type        = $this->type;

	mkpath $read_dir, {mode => $this->file_perm};

	# Local vars

	my @fasta_files = @{$this->{READS}};
	my @qual_files  = @{$this->{QUALS}};
	my @xml_files   = @{$this->{XMLS}};
	my $sample      = $this->{SAMPLE} || $metagenome; # If no sample specified, use default
	my $sample_id   = $this->get_sample_id($metagenome, $sample);

	if ($type eq "sanger") {

		print $PROGRESS "Reading fasta files ...";
		if ($this->readFasta(@fasta_files) != 1) {
			die "Encountered error in readFasta";
		}
		print $PROGRESS " done\n";

		# Generate new fasta, qual and xml files with id's as names
		# Changes:
		# xml: trace_name, library_id
		# fasta, qual: fasta header

		my $prefix = "$read_dir/$metagenome.$sample_id.$type";
		print $PROGRESS "Reading traceinfo xml files: processing traces (counts in 10000) ";
		if ($this->readXML($sample_id, @xml_files) != 1) {
			die "Encountered error in readXML";
		}
		print $PROGRESS "\n";
		if (-f "$prefix.fasta") {
			warn "WARNING: $prefix.fasta exists and will be overwritten. Please add Sanger data belonging to the same sample in one step.\n";
		}
		$this->write_remapped_fasta_file(IN => \@fasta_files, OUT => "$prefix.fasta");

		if (scalar(@qual_files)>0) {
			if (-f "$prefix.qual") {
				warn "WARNING: $prefix.qual exists and will be overwritten. Please add Sanger data belonging to the same sample in one step.\n";
			}
			$this->write_remapped_quality_file(IN => \@qual_files, OUT => "$prefix.qual");
		}
		print $PROGRESS "Generating remapped xml file ...";
		my $first_line  = "<?xml version=\"1.0\"?>";
		my $fh;
			if (-f "$prefix.xml") {
				warn "WARNING: $prefix.xml exists and will be overwritten. Please add Sanger data belonging to the same sample in one step.\n";
			}
		open($fh, ">$prefix.xml") || die "Cannot open $prefix.xml: $!";
		print $fh "$first_line\n";
		print $fh "<TRACE_VOLUME>\n";
		$this->remapXML($fh, @xml_files);
		print $fh "</TRACE_VOLUME>\n";
		close($fh);
		print $PROGRESS " done\n";

		# Do quality trimming
		if ($this->quality_trim) {
			$this->trim_sequences($prefix);
		}

		print $PROGRESS "Loading read information into database ...";

		my $dbh  = $this->get_db_handle;
		my $sth  = $dbh->prepare('INSERT INTO readinfo(read_id, defline, template_id, direction, library_id, clip_left, clip_right, length) VALUES(?, ?, ?, ?, ?, ?, ?, ?);');

		foreach my $name (sort keys %FastaInfo) {
			$sth->execute($FastaInfo{$name}{ID}, $FastaInfo{$name}{DEF}, $Trace{$name}{TEMPLATE_ID}, $Trace{$name}{TRACE_END}, $Trace{$name}{LIB}, $Trace{$name}{LEFT}, $Trace{$name}{RIGHT}, $FastaInfo{$name}{LEN});
		}
		print $PROGRESS " done\n";

		$sth->finish();
		$dbh->commit();
		$this->close_db_handle();
	} elsif ($type eq "454") {
		my @sff_files       = @{$this->{SFFS}};
		my $library_name    = $this->{LIBRARY} || ($this->metagenome.".454");
		my $pe_insert_size  = $this->{INSERT_SIZE};
		if ($pe_insert_size && !$this->{INSERT_STDEV}) {
			$this->{INSERT_STDEV} = int(0.5 + $pe_insert_size / 10);
		}
		my $pe_insert_stdev = $this->{INSERT_STDEV};
		my $library_id      = $this->get_library_id($library_name, $sample_id, $type, $pe_insert_size, $pe_insert_stdev);
		my $prefix          = "$read_dir/$metagenome.$sample_id.$type.$library_id";
		my $frg;
		if (@sff_files) {
			my ($fasta, $qual);
			($frg, $fasta, $qual) = $this->convert_sff_to_fasta(@sff_files);
			@fasta_files = ($fasta);
			@qual_files  = ($qual);
		}

		print $PROGRESS "Reading fasta files ...";
		if ($this->readFasta(@fasta_files) != 1) {
			die "Encountered error in readFasta";
		}
		print $PROGRESS " done\n";

		if (-f "$prefix.fasta") {
			warn "WARNING: $prefix.fasta exists and will be overwritten. Please add 454 data belonging to the same library in one step.\n";
		}
		$this->write_remapped_fasta_file(IN => \@fasta_files, OUT => "$prefix.fasta");

		if (scalar(@qual_files)>0) {
			$this->write_remapped_quality_file(IN => \@qual_files, OUT => "$prefix.qual");
		}

		if (@sff_files) {

			unlink @fasta_files;
			unlink @qual_files;

			# we dont remap frg files anymore, because it is easier to use the original read names in Celera.
			# we keep the frg files the same, but the fasta and qual files will be remapped.
			# see corresponding note in Celera.pm

			#$this->write_remapped_frg_file(IN => [$frg], OUT => "$prefix.frg");
			#unlink $frg;
			move($frg, "$prefix.frg");
		} else {
			# Do quality trimming
			if ($this->quality_trim) {
				$this->trim_sequences($prefix);
			}
		}

		print $PROGRESS "Loading read information into database ...";

		my $dbh  = $this->get_db_handle;
		my $sth  = $dbh->prepare('INSERT INTO readinfo(read_id, defline, template_id, direction, library_id, clip_left, clip_right, length) VALUES(?, ?, ?, ?, ?, ?, ?, ?);');

		foreach my $name (sort keys %FastaInfo) {
			if (!$FastaInfo{$name}{ID}) {
				print "$name\n";
			}
			$sth->execute($FastaInfo{$name}{ID}, $FastaInfo{$name}{DEF}, $FastaInfo{$name}{ID}, undef, $library_id, undef, undef, $FastaInfo{$name}{LEN});
		}
		print $PROGRESS " done\n";

		$sth->finish();
		$dbh->commit();
		$this->close_db_handle();
	} elsif ($type eq "external") {
		my $library_name    = $this->{LIBRARY} || ($this->metagenome.".external");
		my $library_id      = $this->get_library_id($library_name, $sample_id, $type, undef, undef);
		my $prefix          = "$read_dir/$metagenome.$sample_id.$type.$library_id";

		print $PROGRESS "Reading fasta files ...";
		if ($this->readFasta(@fasta_files) != 1) {
			die "Encountered error in readFasta";
		}
		print $PROGRESS " done\n";

		if (-f "$prefix.fasta") {
			warn "WARNING: $prefix.fasta exists and will be overwritten. Please add external data belonging to the same library in one step.\n";
		}
		$this->write_remapped_fasta_file(IN => \@fasta_files, OUT => "$prefix.fasta");

		if (scalar(@qual_files)>0) {
			$this->write_remapped_quality_file(IN => \@qual_files, OUT => "$prefix.qual");
		}

		print $PROGRESS "Loading read information into database ...";

		my $dbh  = $this->get_db_handle;
		my $sth  = $dbh->prepare('INSERT INTO readinfo(read_id, defline, template_id, direction, library_id, clip_left, clip_right, length) VALUES(?, ?, ?, ?, ?, ?, ?, ?);');

		foreach my $name (sort keys %FastaInfo) {
			if (!$FastaInfo{$name}{ID}) {
				print "$name\n";
			}
			$sth->execute($FastaInfo{$name}{ID}, $FastaInfo{$name}{DEF}, $FastaInfo{$name}{ID}, undef, $library_id, undef, undef, $FastaInfo{$name}{LEN});
		}
		print $PROGRESS " done\n";

		$sth->finish();
		$dbh->commit();
		$this->close_db_handle();
	}
}

sub convert_sff_to_fasta {
	use File::Temp;
	my $this = shift;
	my $library_name = $this->{LIBRARY};
	my @sff_files = @_;

	my $celera = new Smash::Analyses::Assembler::Celera(METAGENOME => $this->metagenome, ASSEMBLER => "Celera", NAME => "sffToFasta", CLUSTER => $this->cluster);
	$celera->init();
	my $exe = $celera->pkg_dir."/bin/sffToCA";
	$celera->finish();
	if (! -f $exe) {
		warn "$exe does not exist. Loading SFF files requires Celera assembler to be installed\n";
		$this->abort();
	}

	print $PROGRESS "Converting SFF to FASTA/QUAL/FRG ...";
	my $tmpfile = new File::Temp();
	my $frg = "$tmpfile.frg";
	my $lib_options = "-libraryname $library_name";
	if ($this->{INSERT_SIZE} && $this->{INSERT_STDEV}) {
		$lib_options = join(" ", $lib_options, "-insertsize", $this->{INSERT_SIZE}, $this->{INSERT_STDEV}, "-linker", $this->{TECH});
	}
	my $command = "$exe -clear 454 -trim chop $lib_options -output $frg ".join(" ", @sff_files);
	my $status = $this->execute("$command 2>$tmpfile.err");
	if ($status != 0) {
		warn "sffToCA failed with exit code $status\n";
		warn "Here's what I captured:\n";
		open(E, "<$tmpfile.err");
		while (<E>) {
			warn $_;
		}
		close(E);
		unlink "$tmpfile.err";
	}

	my ($fasta_file, $qual_file) = ("$tmpfile.fasta", "$tmpfile.qual");
	open(FASTA, ">$fasta_file") || die "Cannot open $fasta_file: $!";
	open(QUAL,  ">$qual_file" ) || die "Cannot open $qual_file: $!";

	my $acc;
	open(FRG, "<$frg") || die "Cannot open $frg";
	while (<FRG>) {
		chomp();
		# { to cancel the extra one opened inside the regexp in next line
		if (m/^}$/) { 
			$acc = undef;
		} elsif (m/^seq:$/) {
			my $seq = "";
			SEQ:while (<FRG>) {
				chomp();
				if (m/^\.$/) {
					last SEQ;
				}
				$seq .= $_;
			}
			print FASTA ">$acc\n";
			print FASTA $this->pretty_fasta($seq);
		} elsif (m/^qlt:$/) {
			my $qual = "";
			QUAL:while (<FRG>) {
				chomp();
				if (m/^\.$/) {
					last QUAL;
				}
				$qual .= $_;
			}
			my $zero = ord('0');
			print QUAL  ">$acc\n";
			print QUAL  $this->pretty_qual(join(" ", map {ord($_)-$zero} split(//, $qual)));
		} elsif (m/^acc:(.+)/) {
			$acc = $1;
		}
	}

	close(FRG);
	close(FASTA);
	close(QUAL);

	print $PROGRESS " done\n";

	return ($frg, $fasta_file, $qual_file);
}

sub write_remapped_fasta_file {
	my $this        = shift;
	my %args        = @_;
	my $metagenome  = $this->metagenome;
	my $filename    = ${args{OUT}};
	my @fasta_files = @{$args{IN}};
	my $fh;
	print $PROGRESS "Generating remapped fasta file ...";
	open($fh, ">$filename") || die "Cannot open $filename: $!";
	foreach my $file (@fasta_files) {
		open(FASTA, "<$file") || die "Cannot open $file: $!";
		my $fasta = new FAlite(\*FASTA);
		while (my $entry = $fasta->nextEntry) {
			my $def = $entry->def;
			my $seq = $entry->seq;
			$def =~ s/^>//;
			if ($chop_fasta_header == 1) {
				my @words = split(/\s/, $def);
				if ($weird_fasta_header) {
					$def = $words[$#words];
				} else {
					$def = $words[0];
				}
			}
			my $name                       = $def;                             # name is the processed def
			my $remapped_name              = "$metagenome.$name";     # This is where fasta def is modified to include MC2.MG1 etc
			$FastaInfo{$name}{ID}          = $remapped_name;
			$ReverseLookup{$remapped_name} = $name;
			print $fh ">$metagenome.$name\n";
			print $fh $this->pretty_fasta($seq);
		}
	}
	close($fh);
	print $PROGRESS " done\n";
}

sub write_remapped_quality_file {
	my $this        = shift;
	my %args        = @_;
	my $filename    = ${args{OUT}};
	my @qual_files  = @{$args{IN}};
	my $fh;

	my $skipped = 0;
	print $PROGRESS "Generating remapped qual file ...";
	open($fh, ">$filename") || die "Cannot open $filename: $!";
	foreach my $file (@qual_files) {
		open(QUAL, "<$file") || die "Cannot open $file: $!";
		my $qual = new FQlite(\*QUAL);
		while (my $entry = $qual->nextEntry) {
			my $def = $entry->def;
			my $seq = $entry->seq;
			$def =~ s/^>//;
			if ($chop_fasta_header == 1) {
				my @words = split(/\s/, $def);
				if ($weird_fasta_header) {
					$def = $words[$#words];
				} else {
					$def = $words[0];
				}
			}
			my $name = $def;                # name is the processed def
			if (defined($FastaInfo{$name})) {
				$name = $FastaInfo{$name}{ID};  # get the modified fasta def from the hash
				print $fh ">$name\n";
				print $fh $this->pretty_qual($seq);
			} else {
				$skipped++;
			}
		}
	}
	close($fh);
	if ($skipped > 0) {
		print $PROGRESS " ... (skipped $skipped entries in quality file) ...";
	}
	print $PROGRESS " done\n";
}

sub write_remapped_frg_file {
	my $this        = shift;
	my %args        = @_;
	my $filename    = ${args{OUT}};
	my @frg_files   = @{$args{IN}};
	my $fh;
	print $PROGRESS "Generating remapped frg file ...";
	open($fh, ">$filename") || die "Cannot open $filename: $!";
	foreach my $file (@frg_files) {
		open(FRG, "<$file") || die "Cannot open $file: $!";
		while (<FRG>) {
			if (m/^acc:(.+)/) {
				my $acc = $1;
				if (defined($FastaInfo{$acc})) {
					$acc = $FastaInfo{$acc}{ID};  # get the modified fasta def from the hash
				}
				print $fh "acc:$acc\n";
			} else {
				print $fh $_;
			}
		}
		close(FRG);
	}
	close($fh);
	print $PROGRESS " done\n";
}

sub trim_sequences {
	use File::Copy;
	use File::Basename;
	use File::Temp;
	use Smash::Utils;
	use Smash::Analyses::Assembler::Forge;

	my $this         = shift;
	my $prefix       = shift;
	my $type         = $this->type;
	my $metagenome   = $this->metagenome;
	my $read_dir     = $this->read_dir($metagenome);
	my $trimming_dir = "$read_dir/originals";
	my $maui_dir     = $this->maui_dir;
	my @reads        = $this->fasta_files($metagenome);
	my $quality_trim = $this->quality_trim;

	mkpath "$trimming_dir";

	# trim sequences

	my $read = "$prefix.fasta";
	my $qual = "$prefix.qual";
	my $xml  = "$prefix.xml";

	# Make a list of files to backup and then process

	my @list_of_files = ($read, $qual);
	if (lc($type) eq "sanger") {
		push(@list_of_files, $xml);
	}

	# Do the clipping
	# Otherwise, do full fledged trimming

	if (defined($quality_trim)) {
		print $PROGRESS "Trimming reads for low quality sequences ...";

		# Backup

		foreach my $file2bak (@list_of_files) {
			$this->execute("gzip -c < $file2bak > $trimming_dir/@{[basename($file2bak)]}.bak.gz");
		}

		print $PROGRESS ".";

		# Copy the files to working directory

		my $tmp_fasta = "$trimming_dir/tmp.fasta";
		my $tmp_qual  = "$trimming_dir/tmp.qual";
		copy $read, $tmp_fasta;
		copy $qual, $tmp_qual;
		print $PROGRESS ".";

		# Clip quality using forge/lucy

		if ($quality_trim eq "forge") {

			my $forge     = new Smash::Analyses::Assembler::Forge(ASSEMBLER => "Forge", HOSTNAME => ($this->cluster || $this->host));
			$forge->parse_config();
			$forge->init_software_details();
			my $forge_dir = $forge->pkg_dir();

			# prepare host file for lamboot
			my $host_file = new File::Temp(UNLINK => 1);
			print $host_file "localhost cpu=2\n";
			close($host_file);

			execute_and_report($PROGRESS, "cd $trimming_dir && lamboot $host_file 1>/dev/null && mpirun -np 2 $forge_dir/makeClip -t tmp 1>$trimming_dir/tmp.out && lamhalt 1>/dev/null");
			Smash::Utils->forge2maui("$trimming_dir/tmp.clip", "$trimming_dir/tmp.trimInfo", "both");
			execute_and_report($PROGRESS, "$maui_dir/trimFastaHash $tmp_fasta $trimming_dir/tmp.trimInfo > $read");
			execute_and_report($PROGRESS, "$maui_dir/trimQualHash $tmp_qual $trimming_dir/tmp.trimInfo > $qual");

			foreach my $ext qw(clip clipReport trimInfo out) {
				unlink "$trimming_dir/tmp.$ext";
			}

		} elsif ($quality_trim eq "lucy") {

			my $smash = new Smash::Core();
			$smash->init();
			my ($lucy_dir) = $smash->software_dir("lucy", "current");

			my ($f, $q);
			$f = "$tmp_fasta.lucy";
			$q = "$tmp_qual.lucy";

			execute_and_report($PROGRESS, "$lucy_dir/lucy -quiet -minimum 100 -output $f $q $tmp_fasta $tmp_qual"); # consider using -xtra N for multithreading

			# Reformat them

			# Open FAlite and FQlite objects of tmp files

			open(FASTA, "<$f") || die "Cannot open $f: $!";
			open(QUAL,  "<$q") || die "Cannot open $q: $!";
			my $falite = new FAlite(\*FASTA);
			my $fqlite = new FQlite(\*QUAL);

			# Open write handles to the read/qual files

			open(FOUT, ">$read") || die "Cannot open $read to write: $!";
			open(QOUT, ">$qual") || die "Cannot open $qual to write: $!";

			# Trim defs and write them

			my ($qentry, $fentry);
			while (($fentry = $falite->nextEntry) && ($qentry = $fqlite->nextEntry)) {
				my $def = $fentry->def;
				$def =~ s/\s.*//;
				if ($def ne $qentry->def) {
					die "$f and $q are not synchronized";
				}
				print FOUT $def."\n";
				print FOUT $this->pretty_fasta($fentry->seq);
				print QOUT $def."\n";
				print QOUT $this->pretty_qual($qentry->seq);
			}
			close(FOUT);
			close(QOUT);
			close(FASTA);
			close(QUAL);
			foreach my $ext qw(fasta.lucy qual.lucy) {
				unlink "$trimming_dir/tmp.$ext";
			}
		} elsif ($quality_trim eq "xml" && lc($type) eq "sanger") {

			# Clip CLIP_VECTOR_LEFT entries, if any

			# Open FAlite and FQlite objects of tmp files

			open(FASTA, "<$tmp_fasta") || die "Cannot open $tmp_fasta: $!";
			open(QUAL,  "<$tmp_qual")  || die "Cannot open $tmp_qual: $!";
			my $falite = new FAlite(\*FASTA);
			my $fqlite = new FQlite(\*QUAL);

			# Open write handles to the read/qual files

			open(FOUT, ">$read") || die "Cannot open $read to write: $!";
			open(QOUT, ">$qual") || die "Cannot open $qual to write: $!";

			# Trim reads and write them

			my ($qentry, $fentry);
			while (($fentry = $falite->nextEntry) && ($qentry = $fqlite->nextEntry)) {
				if ($fentry->def ne $qentry->def) {
					die "$tmp_fasta and $tmp_qual are not synchronized";
				}
				my $def      = $fentry->def;
				my $sequence = $fentry->seq;
				my $length   = length($sequence);
				$def =~ s/^>//;
				my @words = split(/\s/, $def);
				$def = $words[$#words];
				my $name = $ReverseLookup{$def};
				my $vector_left = $Trace{$name}{LEFT};
				my $vector_right = $Trace{$name}{RIGHT};
				if (defined($vector_left)) {
					$vector_left--;
				} else {
					$vector_left = 0;
				}
				if (defined($vector_right)) {
					$vector_right--;
				} else {
					$vector_right = $length-1;
				}
				if ($vector_left != 0 || $vector_right != $length-1) {
					$sequence = substr($sequence, $vector_left, $vector_right - $vector_left + 1);
					my @quality  = split(" ", $qentry->seq);
					@quality     = @quality[$vector_left..$vector_right];
					print FOUT $fentry->def."\n";
					print FOUT Smash::Core->pretty_fasta($sequence);
					print QOUT $qentry->def."\n";
					print QOUT Smash::Core->pretty_qual(join(" ", @quality));
				} else {
					print FOUT $fentry->def."\n";
					print FOUT Smash::Core->pretty_fasta($sequence);
					print QOUT $qentry->def."\n";
					print QOUT Smash::Core->pretty_qual($qentry->seq);
				}
			}

			close(QOUT);
			close(FOUT);
			close(QUAL);
			close(FASTA);
		} else {
			Smash::Core::suicide("Quality trimming mode $quality_trim is not implemented");
		}
		unlink $tmp_fasta, $tmp_qual;
		print $PROGRESS " done\n";
	}

	# Make a list of reads that remain after clipping
	print $PROGRESS "Making list of reads that remain: ";
	open(FASTA, "<$read") || die "Cannot open $read: $!";
	my $fasta = new FAlite(\*FASTA);
	while (my $entry = $fasta->nextEntry) {

		# By the time we are here, the fasta headers are already processed.
		# So no need to process them again. They can serve as $name in the hash
		# Change it if that assumption changes.

		my $def = $entry->def;
		$def =~ s/^>//;
		my $name = $ReverseLookup{$def};
		$FastaInfo{$name}{LEN} = length($entry->seq);
		if (lc($type) eq "sanger") {
			delete $Trace{$name}{LEFT};
			delete $Trace{$name}{RIGHT};
		}
		$ReadRemains{$name} = 1;
	}
	close(FASTA);
	printf $PROGRESS "%d of %d reads remain!\n", scalar(keys %ReadRemains), scalar(keys %FastaInfo);

	print $PROGRESS "Removing completely trimmed reads ...";
	my @keys = sort keys %FastaInfo;
	foreach my $name (@keys) {
		if (!defined($ReadRemains{$name})) {
			delete $FastaInfo{$name};
		}
	}
	print $PROGRESS " done\n";

	# Write the new xml file only for sequences that exist in the fasta files

	if (lc($type) eq "sanger") {
		print $PROGRESS "Generating new xml file for the trimmed reads ...";
		my $first_line  = "<?xml version=\"1.0\"?>";
		my $new_xml     = "$trimming_dir/tmp.xml";
		my $fh;
		open($fh, ">$new_xml") || die "Cannot open $new_xml: $!";
		print $fh "$first_line\n";
		print $fh "<TRACE_VOLUME>\n";
		$this->filterXML($fh, $xml);
		print $fh "</TRACE_VOLUME>\n";
		close($fh);
		move($new_xml, $xml);
		print $PROGRESS " done\n";
	}
}

sub unload_db {
	my $this       = shift;

	my $metagenome = $this->metagenome || die "metagenome_id not defined. Remove aborted!";
	my $dbh;
	
	# Remove from MC* database

	$dbh        = $this->get_db_handle;

	print $PROGRESS "Removing entries from: ";

	# Remove reads from readinfo and library
	print $PROGRESS "readinfo, library";
	{
		my @libs;
		my $sth1 = $dbh->prepare('SELECT library_id FROM library LEFT JOIN sample USING (sample_id) WHERE metagenome_id=?');
		my $sth2 = $dbh->prepare('DELETE FROM readinfo WHERE library_id=?');
		my $sth3 = $dbh->prepare('DELETE FROM library WHERE library_id=?');
		$sth1->execute($metagenome);
		while (my ($library_id) = $sth1->fetchrow_array()) {
			push(@libs, $library_id);
		}
		foreach my $library_id (@libs) {
			$sth2->execute($library_id);
			$sth3->execute($library_id);
		}
	}

	# Remove entries from sample
	print $PROGRESS ", sample";
	{
		my $sth = $dbh->prepare('DELETE FROM sample WHERE metagenome_id=?');
		$sth->execute($metagenome);
	}
	print $PROGRESS ".\n";

	$dbh->commit();
	$this->close_db_handle();

	# Remove reads

	my $read_dir = $this->read_dir($metagenome);
	print $PROGRESS "Removing files from $read_dir\n";
	rmtree($read_dir);
}

sub wipe_out {
	my $this       = shift;
	my $metagenome = $this->metagenome || die "metagenome_id not defined. Remove aborted!";

	print $PROGRESS "Removing entries from SmashDB.metagenome\n";
	$this->remove_metagenome_by_name($metagenome);
}

sub execute_and_report {
	my $FH      = shift;
	my $command = shift;
	my $status  = system($command);
	print $FH ".";
	return $status;
}

1;
