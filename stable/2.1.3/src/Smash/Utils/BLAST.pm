package Smash::Utils::BLAST;
use strict;
use warnings;
use File::Path;
use Smash::Core;
use base("Smash::Utils");

sub new {
	my $class = shift;
	my %param = @_;
	my $this  = bless {%param}, $class;
	my $smash = new Smash::Core();
	   $smash->init();
	$this->{SMASH} = $smash;
	$this->set_database;
	$this->set_query;
	return $this;
}

sub set_database {
	my $this = shift;
	my $external_dir = $this->smash->get_smash_conf_value("data_dir")."/external";
	my $db_type = $this->db_type;
	my $status = 0;
	my $given_database = $this->database;

	# First try the file in DATABASE

	$status = $this->check_blast_database;

	# Second try database in external directory

	if ($status == 0) {
		$this->{DATABASE} = "$external_dir/$given_database";
		$status = $this->check_blast_database;
	}

	# Third try making database from a file in the external directory, if --makedb

	if ($status == 0 && $this->{MAKEDB}) {
		my $db_input = "$external_dir/$given_database.$db_type.fa";
		if (-f $db_input) {
			$this->make_blast_database_from_file($db_input);
			$status = 1;
		} elsif (-f "$db_input.gz") {
			$this->make_blast_database_from_file("$db_input.gz");
			$status = 1;
		} else {
			die "ERROR: $given_database is not a $db_type blast database that SMASH can make!\n";
		}
	}

	if ($status == 0) {
		die "ERROR: $given_database is not a $db_type blast database in SMASH!\n";
	}
}

##############
# Set the query
##############

# If QUERY is already a file, dont bother figuring it out. Just use it.
# Otherwise, see if it is a genepred/assembly/metagenome, and choose the right file.

sub set_query {
	my $this = shift;
	my $given_query = $this->query;
	my $smash = $this->smash;
	if (!-f $given_query) {
		my $query;
		my $output_dir;
		my ($collection, $metagenome, $assembly, $genepred) = $smash->parse_concat_id($given_query);
		if ($genepred) {
			my $genepred_dir = $smash->genepred_dir($genepred);
			$query           = "$genepred_dir/$genepred.protein.fa";
			$output_dir      = $genepred_dir;
		} elsif ($assembly) {
			my $assembly_dir = $smash->assembly_dir($assembly);
			$query           = "$assembly_dir/$assembly.contigs.fa";
			$output_dir      = $assembly_dir;
		} elsif ($metagenome) {
			my @files        = $smash->fasta_files($metagenome);
			if (@files == 1) { # Just one file, so just use it
				$query   = $files[0];
			} else {           # Multiple files, so make a concatenated file
				$query   = "$metagenome.reads.fasta";
				$smash->execute("cat ".join(" ", @files)." > $query");
			}
			$output_dir      = $smash->analyses_dir($metagenome);
			if (! -d "$output_dir") {
				mkpath($output_dir, {mode => $smash->file_perm});
			}
		} 

		if (! -f $query) {
			die "ERROR: $given_query is neither a fasta file nor a valid SMASH entity!\n";
		}

		$this->{QUERY} = $query;
		$this->{OUTPUT_DIR} = $output_dir;
	}
}

sub smash      {shift->{SMASH}}
sub flavor     {shift->{FLAVOR}}
sub blast      {shift->{BLAST}}
sub version    {shift->{VERSION} || "current"}
sub database   {shift->{DATABASE}}
sub query      {shift->{QUERY}}
sub output_dir {shift->{OUTPUT_DIR}}
sub subjects   {shift->{SUBJECTS}}
sub tabular    {shift->{TABULAR}}
sub cpus       {shift->{CPUS}}
sub evalue     {shift->{EVALUE}}
sub wordsize   {shift->{WORDSIZE}}
sub extra_args {shift->{EXTRA_ARGS}}
sub match      {shift->{M}}
sub mismatch   {shift->{N}}
sub gap_exist  {shift->{G}}
sub gap_extnd  {shift->{E}}
sub outfmt     {shift->{OUTFMT}}

sub type {
	my $this  = shift;
	my $blast = $this->blast;
	if ($blast eq "blastp") {
		return "protein";
	} elsif ($blast eq "blastn") {
		return "nucleotide";
	} else {
		return "translated";
	}
}

sub db_type {
	my $this  = shift;
	my $blast = $this->blast;
	if (":blastn:tblastx:tblastn:" =~ /:${blast}:/i) {
		return "nucleotide";
	} else {
		return "protein";
	}
}

sub make_blast_database_from_file {
	my $this     = shift;
	my $file     = shift;
	my $database = $this->database;
	my $db_type  = $this->db_type;
	my $flavor   = $this->flavor;
	my $command;
	my $type_arg;
	if ($flavor eq "WU") {
		if ($db_type eq "nucleotide") {
			$type_arg = "-n";
		} else {
			$type_arg = "-p";
		}
		if ($file =~ /\.gz$/) {
			$command = "zcat $file | xdformat $type_arg -o $database -- -";
		} else {
			$command = "xdformat $type_arg -o $database $file";
		}
	} elsif ($flavor eq "NCBI") {
		if ($db_type eq "nucleotide") {
			$type_arg = "-p F";
		} else {
			$type_arg = "-p T";
		}
		my ($blast_dir) = $this->smash->software_dir("ncbi-blast", $this->version);
		if ($file =~ /\.gz$/) {
			$command = "zcat $file | $blast_dir/formatdb $type_arg -n $database -i stdin -v 10000"
		} else {
			$command = "$blast_dir/formatdb $type_arg -n $database -i $file -v 10000"
		}
	} elsif ($flavor eq "NCBI+") {
		if ($db_type eq "nucleotide") {
			$type_arg = "-dbtype nucl";
		} else {
			$type_arg = "-dbtype prot";
		}
		my ($blast_dir) = $this->smash->software_dir("ncbi-blast+", $this->version);
		if ($file =~ /\.gz$/) {
			$command = "zcat $file | $blast_dir/makeblastdb $type_arg -out $database -title $file -in - -max_file_sz 10GB"
		} else {
			$command = "$blast_dir/makeblastdb $type_arg -out $database -in $file -max_file_sz 10GB"
		}
	}
	print STDERR "Making $db_type BLAST database: $command\n";
	system($command);
}

sub get_blast_database_extension {
	my $this = shift;
	my $flavor = $this->flavor;
	my $db_type= $this->db_type;
	my $ext;
	if ($flavor eq "WU") {
		if ($db_type eq "nucleotide") {$ext = "xn?"}
		else                          {$ext = "xp?"}
	} elsif ($flavor eq "NCBI") {
		if ($db_type eq "nucleotide") {$ext = "n??"}
		else                          {$ext = "p??"}
	} elsif ($flavor eq "NCBI+") {
		if ($db_type eq "nucleotide") {$ext = "n??"}
		else                          {$ext = "p??"}
	}
	return $ext;
}

sub check_blast_database {
	my $this = shift;
	my $ext  = $this->get_blast_database_extension;
	my $db   = $this->database;
	my @files = <$db.$ext>;
	if (@files < 3) {
		@files = <$db.[0-9]*.$ext>;
		warn "WARNING: Found a segmented database!\n";
		if (@files < 3) {
			return 0;
		}
	}
	return 1;
}

sub get_blast_exe {
	my $this     = shift;
	my $flavor   = $this->flavor;
	my $blast    = $this->blast;
	my $blast_exe;
	if ($flavor eq "WU") {
		$blast_exe = "$blast";
	} elsif ($flavor eq "NCBI") {
		my ($blast_dir) = $this->smash->software_dir("ncbi-blast", $this->version);
		$blast_exe = "$blast_dir/blastall -p $blast";
	} elsif ($flavor eq "NCBI+") {
		my ($blast_dir) = $this->smash->software_dir("ncbi-blast+", $this->version);
		$blast_exe = "$blast_dir/$blast";
	} else {
		die "Unknown BLAST flavor: $flavor!";
	}
	return $blast_exe;
}

sub get_command_line {
	my $this      = shift;
	my $blast_exe = $this->get_blast_exe();
	my $db_arg    = $this->get_database_arg();
	my $query_arg = $this->get_query_arg();
	my $params    = $this->get_blast_parameters();
	my $cpu_arg   = $this->get_cpu_arg();
	return "$blast_exe $db_arg $query_arg $params $cpu_arg";
}

sub get_query_arg {
	my $this     = shift;
	my $query    = $this->query;
	my $flavor   = $this->flavor;
	my $arg;
	if ($flavor eq "WU") {
		$arg = "$query";
	} elsif ($flavor eq "NCBI") {
		$arg = "-i $query";
	} elsif ($flavor eq "NCBI+") {
		$arg = "-query $query";
	} else {
		die "Unknown BLAST flavor: $flavor!";
	}
	return $arg;
}

sub get_database_arg {
	my $this     = shift;
	my $database = $this->database;
	my $flavor   = $this->flavor;
	my $arg;
	if ($flavor eq "WU") {
		$arg = "$database";
	} elsif ($flavor eq "NCBI") {
		$arg = "-d $database";
	} elsif ($flavor eq "NCBI+") {
		$arg = "-db $database";
	} else {
		die "Unknown BLAST flavor: $flavor!";
	}
	return $arg;
}

sub get_vanilla_command_line {
	my $this     = shift;
	my $blast    = shift;
	my $database = shift;
	my $query    = shift;
	my $flavor   = $this->flavor;
	my $command;

	if ($flavor eq "WU") {
		$command = "$blast $database $query";
	} elsif ($flavor eq "NCBI") {
		$command = "$blast -d $database -i $query";
	} elsif ($flavor eq "NCBI+") {
		my ($blast_dir) = $this->smash->software_dir("ncbi-blast+", $this->version);
		$command = "$blast -db $database -query $query";
	} else {
		die "Unknown BLAST flavor: $flavor!";
	}
	return $command;
}

# This gives you everything except query, database and CPUs.
# Useful for BLAST where these cannot be set in stone!

sub get_cpu_arg {
	my $this   = shift;
	my $cpus   = $this->cpus;
	my $flavor = $this->flavor;
	my $arg;
	if ($flavor eq "WU") {
		$arg = "cpus=$cpus";
	} elsif ($flavor eq "NCBI") {
		$arg = "-a $cpus";
	} elsif ($flavor eq "NCBI+") {
		$arg = "-num_threads $cpus";
	} else {
		die "Unknown BLAST flavor: $flavor!";
	}
	return $arg;
}

sub get_blast_parameters {
	my $this   = shift;
	my $flavor = $this->flavor;
	my $blast  = $this->blast;
	my $tabular= $this->tabular;
	my $type   = $this->type;
	my $db_type= $this->db_type;

	# Set wordsize

	if (!$this->wordsize) {
		if ($type eq "nucleotide") {
			$this->{WORDSIZE} = 11;
		} else {
			$this->{WORDSIZE} = 3;
		}
	}

	# Make command

	my $command;
	if ($flavor eq "WU") {
		$command = sprintf(
				"E=%s B=%s W=%d %s", 
				$this->evalue, 
				$this->subjects, $this->wordsize, $this->extra_args
				);
		# span1 is very good in principle, but it would mess things up when you use mpblast.pl, where span applies to the
		# whole query which is a multiplex of multiple queries. So I disable it. You can still specify span1 in --extra_args
		# $command .= " span1";
		if ($tabular) {
			$command .= " mformat=2";
		}

		# Make them comparable to NCBI blast defaults: Taken from http://blast.advbiocomp.com/doc/cparms.html (03.11.2010)

		if ($blast eq "blastp") {
			$command .= " hitdist=40 T=11 kap s2=41 gaps2=62 x=16 gapx=38";
		} elsif ($blast eq "blastx") {
			$command .= " hitdist=40 T=12 s2=41 gaps2=68 x=16 gapx=38";
		} elsif ($blast eq "tblastn") {
			$command .= " hitdist=40 T=13 s2=41 gaps2=62 x=16 gapx=38";
		} elsif ($blast eq "tblastx") {
			$command .= "";
		} elsif ($blast eq "blastn") {
			#$command .= " kap x=6 gapx=25";
		}
		if ($db_type ne "nucleotide") {
			$command .= " q=12 r=1 gapL=.27 gapK=.047 gapH=.23";
		}

		# disable warnings

		$command .= " warnings notes globalexit";

	} elsif ($flavor eq "NCBI") {
		$command = sprintf(
				"-e %s -b %s -W %d %s", 
				$this->evalue, 
				$this->subjects, $this->wordsize, $this->extra_args
				);
		$command .= " -F F"; # no low complexity filtering
		if ($tabular) {
			$command .= " -m 8";
		}
	} elsif ($flavor eq "NCBI+") {
		$command = sprintf(
				"-evalue %s -max_target_seqs %s -word_size %d %s",
				$this->evalue, 
				$this->subjects, $this->wordsize, $this->extra_args
				);
		#equivalent of span1
		# see comment above for span1 - disabling this
		#$command .= " -culling_limit 1"; 

		# no low complexity filtering

		if ($blast eq "blastn") {
			$command .= " -dust no";
		} else {
			$command .= " -seg no";
		}
		if ($tabular) {
			my $outfmt = "7 std";
			# Next one is to emulate WU-BLAST tabular format, but disabled for now
			#my $outfmt = "7 qseqid sseqid evalue bitscore score length nident positive mismatch pident ppos gapopen gaps qframe qstart qend sframe sstart send";
			$command .= " -outfmt \"$outfmt\"";
		} elsif ($this->outfmt) {
			$command .= " -outfmt '".$this->outfmt."'";
		}
	} else {
		die "Unknown BLAST flavor: $flavor!";
	}
	return $command;
}

# Using MAXSUBJECTS will give all the hits that score as much as the top MAXSUBJECTS'th subject. 
# So it is safe to say MAXSUBJECTS => 1, since it will get all the top hits if they score exactly the same.

sub write_sorted_blast_report {

	my %params = @_;

	my $bit_threshold;
	my $exp_threshold;
	my $max_subjects;
	my $sim_threshold;
	my $pid_threshold;
	my $min_query_aligned_bases;
	my $min_query_aligned_fraction;
	my $query_length_file;
	my $min_subject_aligned_bases;
	my $min_subject_aligned_fraction;
	my $subject_length_file;
	my $top_bits;
	my $prev_query = "BlahBlahBooh";

	my @group          = ();
	my %QueryLengths   = ();
	my %SubjectLengths = ();

	# Mandatory

	my $flavor                    = uc($params{'FLAVOR'})  || die "FLAVOR undefined in BLAST";
	my $blast_file                = $params{'BLAST_FILE'}  || die "BLAST undefined in BLAST";
	my $out_fh                    = $params{'OUT_FH'}      || die "Output file handler required";

	# Optional

	$bit_threshold                = $params{'BITSCORE'}    || 40;
	$top_bits                     = $params{'TOPBITS'};
	$exp_threshold                = $params{'EXP'}         || 10;
	$sim_threshold                = $params{'SIM'}         || 10;
	$pid_threshold                = $params{'IDENTITY'}    || 10;
	$max_subjects                 = $params{'MAXSUBJECTS'};
	$min_query_aligned_bases      = $params{'MINQLEN'}     || 0;
	$min_query_aligned_fraction   = $params{'MINQFRAC'}    || 0;
	$query_length_file            = $params{'QLENFILE'};
	$min_subject_aligned_bases    = $params{'MINSLEN'}     || 0;
	$min_subject_aligned_fraction = $params{'MINSFRAC'}    || 0;
	$subject_length_file          = $params{'SLENFILE'};

	my $BLAST;
	if (ref $blast_file =~ /GLOB/) {    # argument is already a glob
		$BLAST = $blast_file;
	} elsif ($blast_file eq "-") {      # argument is "-", STDIN
		$BLAST = \*STDIN;
	} elsif ($blast_file =~ /\.gz$/) {  # argument is a gzipped file
		open($BLAST, "gunzip -c $blast_file |") || die "Cannot open gunzip pipe for $blast_file: $!";
	} else {                            # agument is a regular file
		open($BLAST, "<$blast_file") || die "Cannot open $blast_file: $!";
	}

	if (defined($query_length_file)) {
		open(FILE, "<$query_length_file") || die "Cannot open $query_length_file: $!";
		while(my $line=<FILE>) {
			chomp($line);
			my ($def, $len) = split(/\t/, $line);
			$QueryLengths{$def} = $len;
		}
		close(FILE);
	}

	if (defined($subject_length_file)) {
		open(FILE, "<$subject_length_file") || die "Cannot open $subject_length_file: $!";
		while(my $line=<FILE>) {
			chomp($line);
			my ($def, $len) = split(/\t/, $line);
			$SubjectLengths{$def} = $len;
		}
		close(FILE);
	}

	if (defined $top_bits) {
		die "TOP_BITS must be between 0 and 1\n" unless ($top_bits >= 0 && $top_bits <= 1);
		$top_bits = 1-$top_bits;
	}

	my @wu_indices   = (0,1,2,4,5,10,11,17,18,20,21);
	my @ncbi_indices = (0,1,10,11,11,2,2,6,7,8,9);
	my $score_field  = 5;
	   $score_field  = 11 if $flavor =~ /^NCBI\+?/;
	my $bits_field  = 4;
	   $bits_field  = 11 if $flavor =~ /^NCBI\+?/;


	if (defined($max_subjects)) {
		LINE: while(my $line=<$BLAST>) {
			next LINE if $line =~ /^\s*#/;
			my $skip = 0;
			chomp($line);
			my @words = split(/\t/, $line);
			my ($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send);

			if ($flavor eq "WU") {
				die "Error: WU BLAST requires at least 22 columns. Check the following line:\n$line\n" unless @words > 21;
				($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send) = @words[@wu_indices];
			} elsif ($flavor =~ /^NCBI\+?/) {
				die "Error: NCBI BLAST requires exactly 12 columns. Check the following line:\n$line\n" unless @words == 12;
				($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send) = @words[@ncbi_indices];
			}

			#if ($query eq "query") {next LINE;}

			if ($bits < $bit_threshold || $E > $exp_threshold || $pcpos < $sim_threshold || $pcid < $pid_threshold) { $skip = 1;}
			
			{
				my $aligned_query_length = abs($qend - $qstart) + 1;
				if (($aligned_query_length < $min_query_aligned_bases) ||
				    ($QueryLengths{$query} && $aligned_query_length / $QueryLengths{$query} < $min_query_aligned_fraction)) { 
					$skip = 1;
				}

				my $aligned_subject_length = abs($send - $sstart) + 1;
				if (($aligned_subject_length < $min_subject_aligned_bases) ||
				    ($SubjectLengths{$subject} && $aligned_subject_length / $SubjectLengths{$subject} < $min_subject_aligned_fraction)) { 
					$skip = 1;
				}
			}

			if (eof || $query ne $prev_query) {

				# process the group of subjects for this query

				# if at end of file, add the line just read before you process the group of subjects

				push(@group, $line) if eof;

				# Sort the lines based on their scores

				my @scores = map {(split(/\t/, $_))[$score_field]} @group;
				@group = @group[ sort {$scores[$b] <=> $scores[$a]} 0..$#group ]; # sort the lines based on the score of those lines
				@scores = (); # end context

				# Dont mess with the sort and map lines unless you know what you are doing
				# Get $max_subjects subjects in the sorted list

				# For each subject, we need to follow its max score in any HSP with this query
				# @group is ordered desc by score, so the first time you see a subject is its highest score!
				
				my %MaxScore = ();
				@scores = ();
				for (@group) {
					my @words = split(/\t/);
					my $sbjct = $words[1];
					my $score = $words[$score_field];
					if (!defined($MaxScore{$sbjct})) {
						$MaxScore{$sbjct} = $score;
						push(@scores, $score);
					}
				}

				# Now we have: $MaxScore{"someone"} defined as max score for "someone"

				# We now get a threshold that is the score of the $max_subjects'th score

				#@scores = sort {$b <=> $a} @scores; # why again? @score should be sorted since @group is sorted by score!
				my $threshold;
				if ($max_subjects <= @scores) {
					$threshold = $scores[$max_subjects-1];
				} else {
					$threshold = $scores[$#scores];
				}

				# find the subjects that have at least one hsp higher than threshold
				# then keep all their hsps
				# we dont care about the rest
				map {delete $MaxScore{$_} if $MaxScore{$_} < $threshold;} keys %MaxScore;
				@group = grep {defined($MaxScore{(split(/\t/))[1]})} @group;

				# Get all the subjects, since we need them to sort
				my @sbjcts = map {(split(/\t/, $_))[1]} @group;
				@group = @group[ sort {$MaxScore{$sbjcts[$b]} <=> $MaxScore{$sbjcts[$a]} || $sbjcts[$b] cmp $sbjcts[$a]} 0..$#group ]; # sort the lines based on the score of those lines and group subjects together as well

				# Print sorted lines
				for (@group) {
					print $out_fh "$_\n";
				}
				%MaxScore = ();

				@group = ();
			}

			# if at the end of file, we have processed the last line already!
			if (!eof && $skip == 0) {
				push(@group, $line);
			}

			$prev_query = $query;
		}
	} elsif (defined($top_bits)) {
		LINE: while(my $line=<$BLAST>) {
			next LINE if $line =~ /^\s*#/;
			my $skip = 0;
			chomp($line);
			my @words = split(/\t/, $line);
			my ($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send);

			if ($flavor eq "WU") {
				die "Error: WU BLAST requires at least 22 columns. Check the following line:\n$line\n" unless @words > 21;
				($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send) = @words[@wu_indices];
			} elsif ($flavor =~ /^NCBI\+?/) {
				die "Error: NCBI BLAST requires exactly 12 columns. Check the following line:\n$line\n" unless @words == 12;
				($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send) = @words[@ncbi_indices];
			}

			#if ($query eq "query") {next LINE;}

			if ($bits < $bit_threshold || $E > $exp_threshold || $pcpos < $sim_threshold || $pcid < $pid_threshold) { $skip = 1;}
			
			{
				my $aligned_query_length = abs($qend - $qstart) + 1;
				if (($aligned_query_length < $min_query_aligned_bases) ||
				    ($QueryLengths{$query} && $aligned_query_length / $QueryLengths{$query} < $min_query_aligned_fraction)) { 
					$skip = 1;
				}

				my $aligned_subject_length = abs($send - $sstart) + 1;
				if (($aligned_subject_length < $min_subject_aligned_bases) ||
				    ($SubjectLengths{$subject} && $aligned_subject_length / $SubjectLengths{$subject} < $min_subject_aligned_fraction)) { 
					$skip = 1;
				}
			}

			if (eof || $query ne $prev_query) {

				# process the group of subjects for this query

				# if at end of file, add the line just read before you process the group of subjects

				push(@group, $line) if eof($BLAST);

				if (@group) {
					# For each subject, we need to follow its max bitscore in any HSP with this query
					# @group is ordered desc by bitscore, so the first time you see a subject is its highest bitscore!
					
					my %MaxScore = ();
					for (@group) {
						my @words = split(/\t/);
						my $sbjct = $words[1];
						my $bits = $words[$bits_field];
						if (!defined($MaxScore{$sbjct}) || $bits > $MaxScore{$sbjct}) {
							$MaxScore{$sbjct} = $bits;
						}
					}

					# Now we have: $MaxScore{"someone"} defined as max bitscore for "someone"
					# Get the maximum of all MaxScores

					my ($top_score) = sort {$b <=> $a} values %MaxScore;
					my $threshold = $top_bits*$top_score;

					# find the subjects that have at least one hsp higher than threshold
					# then keep all their hsps
					# we dont care about the rest
					map {delete $MaxScore{$_} if $MaxScore{$_} < $threshold;} keys %MaxScore;
					@group = grep {defined($MaxScore{(split(/\t/))[1]})} @group;

					# Get all the subjects, since we need them to sort

					# Print sorted lines
					for (@group) {
						print $out_fh "$_\n";
					}
					%MaxScore = ();

					@group = ();
				}
			}

			# if at the end of file, we have processed the last line already!
			if (!eof && $skip == 0) {
				push(@group, $line);
			}
			$prev_query = $query;
		}

	} else {
		LINE: while(my $line=<$BLAST>) {
			next LINE if $line =~ /^\s*#/;
			chomp($line);
			my @words = split(/\t/, $line);

			my ($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send);

			if ($flavor eq "WU") {
				die "Error: WU BLAST requires at least 22 columns. Check the following line:\n$line\n" unless @words > 21;
				($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send) = @words[@wu_indices];
			} elsif ($flavor =~ /^NCBI\+?/) {
				die "Error: NCBI BLAST requires exactly 12 columns. Check the following line:\n$line\n" unless @words == 12;
				($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send) = @words[@ncbi_indices];
			}

			#if ($query eq "query") {next LINE;}

			if ($bits < $bit_threshold || $E > $exp_threshold || $pcpos < $sim_threshold || $pcid < $pid_threshold) { next LINE;}

			if ($min_query_aligned_bases || $min_query_aligned_fraction || $min_subject_aligned_bases || $min_subject_aligned_fraction) {
				my $aligned_subject_length = abs($send - $sstart) + 1;
				if ($aligned_subject_length < $min_subject_aligned_bases) { next LINE;}
				if ($SubjectLengths{$subject} 
					&& $aligned_subject_length / $SubjectLengths{$subject} < $min_subject_aligned_fraction) { 
						next LINE;
				}

				my $aligned_query_length = abs($qend - $qstart) + 1;
				if ($aligned_query_length < $min_query_aligned_bases) { next LINE;}
				if ($QueryLengths{$query} 
					&& $aligned_query_length / $QueryLengths{$query} < $min_query_aligned_fraction 
					) {
						next LINE;
				}

			}

			print $out_fh "$line\n";
		}
	}
	close($BLAST);
}

# Using MAXSUBJECTS will give all the hits that score as much as the top MAXSUBJECTS'th subject. 
# So it is safe to say MAXSUBJECTS => 1, since it will get all the top hits if they score exactly the same.

sub filter_blast_report {

	my %params = @_;

	my $bit_threshold;
	my $exp_threshold;
	my $max_subjects;
	my $sim_threshold;
	my $pid_threshold;
	my $min_query_aligned_bases;
	my $min_query_aligned_fraction;
	my $query_length_file;
	my $min_subject_aligned_bases;
	my $min_subject_aligned_fraction;
	my $subject_length_file;
	my $top_bits;
	my $prev_query = "BlahBlahBooh";

	my @group          = ();
	my %QueryLengths   = ();
	my %SubjectLengths = ();

	# Mandatory

	my $flavor                    = uc($params{'FLAVOR'})  || die "FLAVOR undefined in BLAST";
	my $blast_file                = $params{'BLAST_FILE'}  || die "BLAST undefined in BLAST";
	my $out_fh                    = $params{'OUT_FH'}      || die "Output file handler required";

	# Optional

	$bit_threshold                = $params{'BITSCORE'}    || 40;
	$top_bits                     = $params{'TOPBITS'};
	$exp_threshold                = $params{'EXP'}         || 10;
	$sim_threshold                = $params{'SIM'}         || 10;
	$pid_threshold                = $params{'IDENTITY'}    || 10;
	$max_subjects                 = $params{'MAXSUBJECTS'};
	$min_query_aligned_bases      = $params{'MINQLEN'}     || 0;
	$min_query_aligned_fraction   = $params{'MINQFRAC'}    || 0;
	$query_length_file            = $params{'QLENFILE'};
	$min_subject_aligned_bases    = $params{'MINSLEN'}     || 0;
	$min_subject_aligned_fraction = $params{'MINSFRAC'}    || 0;
	$subject_length_file          = $params{'SLENFILE'};

	my $BLAST;
	if (ref $blast_file =~ /GLOB/) {    # argument is already a glob
		$BLAST = $blast_file;
	} elsif ($blast_file eq "-") {      # argument is "-", STDIN
		$BLAST = \*STDIN;
	} elsif ($blast_file =~ /\.gz$/) {  # argument is a gzipped file
		open($BLAST, "gunzip -c $blast_file |") || die "Cannot open gunzip pipe for $blast_file: $!";
	} else {                            # agument is a regular file
		open($BLAST, "<$blast_file") || die "Cannot open $blast_file: $!";
	}

	if (defined($query_length_file)) {
		open(FILE, "<$query_length_file") || die "Cannot open $query_length_file: $!";
		while(my $line=<FILE>) {
			chomp($line);
			my ($def, $len) = split(/\t/, $line);
			$QueryLengths{$def} = $len;
		}
		close(FILE);
	}

	if (defined($subject_length_file)) {
		open(FILE, "<$subject_length_file") || die "Cannot open $subject_length_file: $!";
		while(my $line=<FILE>) {
			chomp($line);
			my ($def, $len) = split(/\t/, $line);
			$SubjectLengths{$def} = $len;
		}
		close(FILE);
	}

	if (defined $top_bits) {
		die "TOP_BITS must be between 0 and 1\n" unless ($top_bits >= 0 && $top_bits <= 1);
		$top_bits = 1-$top_bits;
	}

	my @wu_indices   = (0,1,2,4,5,10,11,17,18,20,21);
	my @ncbi_indices = (0,1,10,11,11,2,2,6,7,8,9);
	my $score_field  = 5;
	   $score_field  = 11 if $flavor =~ /^NCBI\+?/;
	my $bits_field  = 4;
	   $bits_field  = 11 if $flavor =~ /^NCBI\+?/;


	if (defined($max_subjects)) {
		LINE: while(my $line=<$BLAST>) {
			next LINE if $line =~ /^\s*#/;
			my $skip = 0;
			chomp($line);
			my @words = split(/\t/, $line);
			my ($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send);

			if ($flavor eq "WU") {
				die "Error: WU BLAST requires at least 22 columns. Check the following line:\n$line\n" unless @words > 21;
				($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send) = @words[@wu_indices];
			} elsif ($flavor =~ /^NCBI\+?/) {
				die "Error: NCBI BLAST requires exactly 12 columns. Check the following line:\n$line\n" unless @words == 12;
				($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send) = @words[@ncbi_indices];
			}

			#if ($query eq "query") {next LINE;}

			if ($bits < $bit_threshold || $E > $exp_threshold || $pcpos < $sim_threshold || $pcid < $pid_threshold) { $skip = 1;}
			
			{
				my $aligned_query_length = abs($qend - $qstart) + 1;
				if (($aligned_query_length < $min_query_aligned_bases) ||
				    ($QueryLengths{$query} && $aligned_query_length / $QueryLengths{$query} < $min_query_aligned_fraction)) { 
					$skip = 1;
				}

				my $aligned_subject_length = abs($send - $sstart) + 1;
				if (($aligned_subject_length < $min_subject_aligned_bases) ||
				    ($SubjectLengths{$subject} && $aligned_subject_length / $SubjectLengths{$subject} < $min_subject_aligned_fraction)) { 
					$skip = 1;
				}
			}

			if (eof || $query ne $prev_query) {

				# process the group of subjects for this query

				# if at end of file, add the line just read before you process the group of subjects

				push(@group, $line) if eof;

				# Sort the lines based on their scores

				my @scores = map {(split(/\t/, $_))[$score_field]} @group;
				@group = @group[ sort {$scores[$b] <=> $scores[$a]} 0..$#group ]; # sort the lines based on the score of those lines
				@scores = (); # end context

				# Dont mess with the sort and map lines unless you know what you are doing
				# Get $max_subjects subjects in the sorted list

				# For each subject, we need to follow its max score in any HSP with this query
				# @group is ordered desc by score, so the first time you see a subject is its highest score!
				
				my %MaxScore = ();
				@scores = ();
				for (@group) {
					my @words = split(/\t/);
					my $sbjct = $words[1];
					my $score = $words[$score_field];
					if (!defined($MaxScore{$sbjct})) {
						$MaxScore{$sbjct} = $score;
						push(@scores, $score);
					}
				}

				# Now we have: $MaxScore{"someone"} defined as max score for "someone"

				# We now get a threshold that is the score of the $max_subjects'th score

				#@scores = sort {$b <=> $a} @scores; # why again? @score should be sorted since @group is sorted by score!
				my $threshold;
				if ($max_subjects <= @scores) {
					$threshold = $scores[$max_subjects-1];
				} else {
					$threshold = $scores[$#scores];
				}

				# find the subjects that have at least one hsp higher than threshold
				# then keep all their hsps
				# we dont care about the rest
				map {delete $MaxScore{$_} if $MaxScore{$_} < $threshold;} keys %MaxScore;
				@group = grep {defined($MaxScore{(split(/\t/))[1]})} @group;

				# Get all the subjects, since we need them to sort
				my @sbjcts = map {(split(/\t/, $_))[1]} @group;
				@group = @group[ sort {$MaxScore{$sbjcts[$b]} <=> $MaxScore{$sbjcts[$a]} || $sbjcts[$b] cmp $sbjcts[$a]} 0..$#group ]; # sort the lines based on the score of those lines and group subjects together as well

				# Print sorted lines
				for (@group) {
					print $out_fh "$_\n";
				}
				%MaxScore = ();

				@group = ();
			}

			# if at the end of file, we have processed the last line already!
			if (!eof && $skip == 0) {
				push(@group, $line);
			}

			$prev_query = $query;
		}
	} elsif (defined($top_bits)) {
		LINE: while(my $line=<$BLAST>) {
			next LINE if $line =~ /^\s*#/;
			my $skip = 0;
			chomp($line);
			my @words = split(/\t/, $line);
			my ($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send);

			if ($flavor eq "WU") {
				die "Error: WU BLAST requires at least 22 columns. Check the following line:\n$line\n" unless @words > 21;
				($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send) = @words[@wu_indices];
			} elsif ($flavor =~ /^NCBI\+?/) {
				die "Error: NCBI BLAST requires exactly 12 columns. Check the following line:\n$line\n" unless @words == 12;
				($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send) = @words[@ncbi_indices];
			}

			#if ($query eq "query") {next LINE;}

			if ($bits < $bit_threshold || $E > $exp_threshold || $pcpos < $sim_threshold || $pcid < $pid_threshold) { $skip = 1;}
			
			{
				my $aligned_query_length = abs($qend - $qstart) + 1;
				if (($aligned_query_length < $min_query_aligned_bases) ||
				    ($QueryLengths{$query} && $aligned_query_length / $QueryLengths{$query} < $min_query_aligned_fraction)) { 
					$skip = 1;
				}

				my $aligned_subject_length = abs($send - $sstart) + 1;
				if (($aligned_subject_length < $min_subject_aligned_bases) ||
				    ($SubjectLengths{$subject} && $aligned_subject_length / $SubjectLengths{$subject} < $min_subject_aligned_fraction)) { 
					$skip = 1;
				}
			}

			if (eof || $query ne $prev_query) {

				# process the group of subjects for this query

				# if at end of file, add the line just read before you process the group of subjects

				push(@group, $line) if eof($BLAST);

				if (@group) {
					# For each subject, we need to follow its max bitscore in any HSP with this query
					# @group is ordered desc by bitscore, so the first time you see a subject is its highest bitscore!
					
					my %MaxScore = ();
					for (@group) {
						my @words = split(/\t/);
						my $sbjct = $words[1];
						my $bits = $words[$bits_field];
						if (!defined($MaxScore{$sbjct}) || $bits > $MaxScore{$sbjct}) {
							$MaxScore{$sbjct} = $bits;
						}
					}

					# Now we have: $MaxScore{"someone"} defined as max bitscore for "someone"
					# Get the maximum of all MaxScores

					my ($top_score) = sort {$b <=> $a} values %MaxScore;
					my $threshold = $top_bits*$top_score;

					# find the subjects that have at least one hsp higher than threshold
					# then keep all their hsps
					# we dont care about the rest
					map {delete $MaxScore{$_} if $MaxScore{$_} < $threshold;} keys %MaxScore;
					@group = grep {defined($MaxScore{(split(/\t/))[1]})} @group;

					# Get all the subjects, since we need them to sort

					# Print sorted lines
					for (@group) {
						print $out_fh "$_\n";
					}
					%MaxScore = ();

					@group = ();
				}
			}

			# if at the end of file, we have processed the last line already!
			if (!eof && $skip == 0) {
				push(@group, $line);
			}
			$prev_query = $query;
		}

	} else {
		LINE: while(my $line=<$BLAST>) {
			next LINE if $line =~ /^\s*#/;
			chomp($line);
			my @words = split(/\t/, $line);

			my ($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send);

			if ($flavor eq "WU") {
				die "Error: WU BLAST requires at least 22 columns. Check the following line:\n$line\n" unless @words > 21;
				($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send) = @words[@wu_indices];
			} elsif ($flavor =~ /^NCBI\+?/) {
				die "Error: NCBI BLAST requires exactly 12 columns. Check the following line:\n$line\n" unless @words == 12;
				($query, $subject, $E, $bits, $score, $pcid, $pcpos, $qstart, $qend, $sstart, $send) = @words[@ncbi_indices];
			}

			#if ($query eq "query") {next LINE;}

			if ($bits < $bit_threshold || $E > $exp_threshold || $pcpos < $sim_threshold || $pcid < $pid_threshold) { next LINE;}

			if ($min_query_aligned_bases || $min_query_aligned_fraction || $min_subject_aligned_bases || $min_subject_aligned_fraction) {
				my $aligned_subject_length = abs($send - $sstart) + 1;
				if ($aligned_subject_length < $min_subject_aligned_bases) { next LINE;}
				if ($SubjectLengths{$subject} 
					&& $aligned_subject_length / $SubjectLengths{$subject} < $min_subject_aligned_fraction) { 
						next LINE;
				}

				my $aligned_query_length = abs($qend - $qstart) + 1;
				if ($aligned_query_length < $min_query_aligned_bases) { next LINE;}
				if ($QueryLengths{$query} 
					&& $aligned_query_length / $QueryLengths{$query} < $min_query_aligned_fraction 
					) {
						next LINE;
				}

			}

			print $out_fh "$line\n";
		}
	}
	close($BLAST);
}

# return -1 on failure
sub check_complete_blast {
	use FAlite;
	my ($blast_file, $fasta_file, $flavor) = @_;
	my $expected_words = 0;
	my @words = ();
	my $query;

	if (! -f $blast_file) {
		warn "$blast_file not found.\n";
		return -1;
	} elsif ( -z $blast_file) {
		warn "$blast_file looks corrupt. It has 0 bytes.\n";
		return -1;
	}

	# Check the BLAST file itself for errors.
	open(BLAST, "<$blast_file") || die "Cannot open $blast_file to read: $!";
	my $line;
	while ($line = <BLAST>) {
		last if ($line !~ /^\s*#/);
	}
	chomp($line);
	@words = split(/\t/, $line);
	$expected_words = scalar(@words);
	while (<BLAST>) {
		next if (m/^\s*#/);
		$line = $_;
	}
	close(BLAST);

	# Check number of words in the last line
	@words = split(/\t/, $line);
	if (scalar(@words) != $expected_words) {
		warn "$blast_file looks corrupt. Found ".scalar(@words)." words instead of $expected_words.\n";
		return -1;
	}

	# Check the query in the last line!

	$query = $words[0];

	# See how many queries have been finished
	open(FASTA, "<$fasta_file") || die "Cannot open $fasta_file: $!";
	my $fasta = new FAlite(\*FASTA);
	my $entries = 0;
	my $blast_count = 0;
	while (my $entry = $fasta->nextEntry) {
		$entries++;
		my $def = $entry->def;
		$def =~ s/^>//;
		$def =~ s/\s.*//;
		if ($def eq $query) {
			$blast_count = $entries;
		}
	}
	close(FASTA);
	if ($blast_count / $entries < 0.9) {
		warn "$blast_file looks incomplete. BLAST finished at query $blast_count of $entries.\n";
		return -1;
	}
	return 1;
}

1;
