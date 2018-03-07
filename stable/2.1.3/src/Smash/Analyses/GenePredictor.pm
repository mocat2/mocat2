package Smash::Analyses::GenePredictor;

use strict;
use warnings;

use base qw(Smash::Analyses);

use strict;
#use POSIX;
use File::Temp;
use File::Path;
use Smash::Utils::GFF;
use Smash::Utils::GFFlite;
use FAlite;

=head1 NAME

Smash::Analyses::GenePredictor - extensible gene prediction software pipeline

=head1 DESCRIPTION

Smash::Analyses::GenePredictor is the parent class for all gene prediction software wrappers in
Smash. It provides a working pipeline for performing gene prediction on a DNA fasta file and
generating the following output files:

=over 4

=item *

gff file with gene coordinates,

=item *

gene fasta file,

=item *

protein fasta file.

=back

Smash::Analyses::GenePredictor, hereafter referred to as B<GenePredictor>, provides an easy interface to 
Smash to make gene predictions on DNA fasta files. The minimal steps in the pipeline in GenePredictor are:

	new();
	init();
	predict();
	finish();

=head1 KNOWN SUBCLASSES

L<Smash::Analyses::GenePredictor::GeneMark>, L<Smash::Analyses::GenePredictor::MetaGene>.

=head1 MEMBER VARIABLES

The following member variables are set by B<GenePredictor> at runtime, from the values
passed to it during creation of the object. Any subclass can assume that these will be 
populated.

=over 4

=item B<predictor>

name of the software

=item B<self_train>

numeric value specifying if this instance must self-train (if supported by the program)

=item B<organism>

the organism where the sequence comes from

=item B<translation_table>

translation table corresponding to B<organism>

=item B<fasta_filename>

name of the input fasta file to make gene predictions

=item B<gff>

gff file containing gene coordinates

=item B<output>

concatenated raw output file for the predictor

=item B<gene>

fasta file containing predicted gene transcripts (DNA space)

=item B<protein>

fasta file containing predicted proteins (aminoacid space)

=back

=cut

sub _overload  {shift->{PREDICTOR}}
sub output     {shift->{OUTPUT}}
sub gff        {shift->{GFF}}
sub gene       {shift->{GENE}}
sub protein    {shift->{PROTEIN}}
sub self_train {shift->{SELF_TRAIN}}
sub predictor  {shift->{PREDICTOR}}

sub organism          {shift->{ORGANISM}}
sub translation_table {shift->{TRANSLATION_TABLE}}

sub fasta_filename {shift->{FASTA_FILE}}
sub set_fasta_filename {
	my $this = shift;
	$this->{FASTA_FILE} = shift;
}

=head1 FUNCTIONS

The following functions are implemented in B<GenePredictor>. Any subclass can 
implement its own versions of these functions. I recommend that subclasses perform their
steps, and also call C<SUPER::function()> to let B<GenePredictor> do what it ought
to. If they do not call the function from the superclass, it is upto them to 
additionally implement the functionalities of the superclass.

=over 4

=item B<new>

This is a default constructor for GenePredictor or its subclasses. It blesses the hash argument into
an object of GenePredictor class or its subclass.

=item B<init>

This function initializes the relevant member variables for GenePredictor class or its subclasses. It also
initializes the input and output files that were passed to C<new()>.

=item B<predict> (now obsolete)

This is the workhorse in the pipeline. This function creates filehandles for all output files, and writes to them as
it processes the input file. It reads each sequence in the input fasta file and runs gene prediction on them
individually. (Aside: This was a design decision due to the fact that some gene predictors, like GeneMark, do not
recognize multi-fasta files. For scalability of Smash and to maintain its modularity, I chose to process sequences
one after the other.) It reads a DNA sequence entry, checks if it is long enough for the program, writes it into a file, 
passes it to the subclass through C<get_command_line> to get a command-line to execute, and then executes it.

=item B<make_genes_and_proteins_wrapper>

Makes a fasta file containing gene sequences (see C<gene> above) and another one containing the
protein sequences (C<protein> above). It used the input fasta file (see C<fasta_filename> above) 
and the GFF file containing gene coordinates (see C<gff> above) to do this. Subclasses 
are recommended to leave this as is and stop at making the GFF file. However, if the 
program itself generates the transcripts and the protein sequences, and the subclass
decides to parse it as well, then these should be written to files whose names can be
obtained by calling C<gene> and C<protein>, respectively. If this is the case, then
C<make_genes_and_proteins_wrapper> must be overridden using an empty function, so that 
B<GenePredictor> does not call it again to make transcripts and proteins.

=item B<finish>

Closes all open handles and connections and calls the C<finish> function of Smash.

=cut

############################################
# Constructor
############################################

sub new {
	my $class  = shift;
	my %params = @_;

	if (!defined($params{'PREDICTOR'})) {
		$params{'PREDICTOR'} = "Generic";
	}
	if (!defined($params{'TYPE'})) {
		$params{'TYPE'} = "GenePredictor";
	}

	my $this = bless {%params}, $class;

	# set organism info

	if (!$this->organism) {
		$this->{ORGANISM} = "bacteria";
	}
	if (!$this->translation_table) {
		my %table = (standard => 1, mold_mitochondrial => 4, eukaryote => 1, virus => 1, bacteria => 11); 
		$this->{TRANSLATION_TABLE} = $table{$this->organism};
		if (!defined($this->translation_table)) {
			warn "WARNING: Do not know how to make gene predictions for ".$this->organism."\n";
			$this->abort();
		}
	}

	return $this;
}

############################################
# Init object
############################################

# NEED:
# GENEPRED/ASSEMBLY
# OUTPUT_DIR
# FASTA_FILE
# TRAIN/HEURISTIC

sub init {
	my $this       = shift;

	$this->parse_config();
	$this->SUPER::init();

	# This has to be set first, since program_id depends on these!

	# Did the caller specify genepred id?
	# If not, make a new genepred id for this run and use it

	my $genepred = $this->genepred;
	if ($genepred) {
		my $genepred_id = $this->get_id_by_name("gene_prediction", $genepred);
		if (!$genepred_id) {
			warn "WARNING: You have requested Smash to use $genepred as gene prediction id, but this does not exist!";
			$this->abort();
		}
	} else {
		my $prog_id   = $this->get_program_id("$this", $this->version, $this->parameter_settings);
		$genepred     = $this->make_new_genepred($this->assembly, $prog_id);
		$this->{GENEPRED} = $genepred;
	}

	# Did the caller specify a name to be used for output files?
	# If not, use genepred id as the name

	if (!$this->name) {
		$this->{NAME} = $genepred;
	}


	# Did the caller specify an output dir for the output files?
	# If not, make it the standard Smash genepred location

	if (!$this->output_dir) {
		$this->{OUTPUT_DIR} = $this->genepred_dir($genepred);
	}

	my $output_dir = $this->output_dir();
	mkpath $output_dir, {mode => $this->file_perm};

	# Final output files

	my $name = $this->name;
	$this->{OUTPUT}  = "$output_dir/$name"; # Default program format
	$this->{GFF}     = "$output_dir/$name.contig2gene.gff";
	$this->{GENE}    = "$output_dir/$name.gene.fa";
	$this->{PROTEIN} = "$output_dir/$name.protein.fa";

	if (-f $this->output) {
		unlink $this->output;
	}
}

############################################
# If you run into error, remove the gp id you just added to the database.
# There might still be cases where errors are not traceable, in which case the
# predictions remain in the database!
############################################

=item C<abort()>

Aborts the current run after removing this genepred id from the database.

=cut

sub abort {
	my $this = shift;
	$this->remove_genepred_by_name($this->name);
	$this->finish();
	die "ERROR  : Aborting gene prediction - see message above.\n";
}

############################################
# The real stuff!
# Predicts genes and creates the gene
# gff file for each prediction.
############################################

sub predict {
	my $this  = shift;
	my $model = shift;
	my $fasta;
	my $header;
	my $combo_output;
	#my $gff_handle;
	#my $gff_file    = $this->gff;
	my $fasta_file  = $this->fasta_filename;
	my $min_length  = $this->min_length;
	my $output_ptr  = new File::Temp(UNLINK => 1, DIR => "/tmp", TEMPLATE => "genemarkXXXXXX");
	my $output_file = "$output_ptr.pred";
	my $input_ptr   = new File::Temp(UNLINK => 1, DIR => "/tmp", TEMPLATE => "genemarkXXXXXX");
	my $input_file  = "$input_ptr.fa";

	$combo_output  = $this->output;

	# Open GFF file
	#open($gff_handle, ">$gff_file") || die "Cannot open $gff_file: $!";

	# Read the input fasta file
	open(FASTA, "<$fasta_file") || die "Cannot open $fasta_file: $!";
	$fasta = new FAlite(\*FASTA);

	# Process each sequence
	SEQ:while (my $entry = $fasta->nextEntry) {
		my $def    = $this->process_fasta_header($entry->def);
		my $length = length($entry->seq);

		# Check if minimum length is met

		if ($length < $min_length) {
			next SEQ;
		}

		# Prepare temp fasta file

		open(TMP, ">$input_file") || die "Cannot open $input_file: $!";
		print TMP "$entry";
		close(TMP);

		# Make command line and execute it

		my $command = $this->get_command_line($input_file, $output_file, $model);

		#####
		# If the first attempt fails, try once more. GeneMark does some weird segfaults that cannot be reproduced
		# As of 21.05, trying it twice gets rid of the seg faults.
		#####
		my $attempt = 0;
		ATTEMPT:while ($attempt < 10) {
			if ($this->execute($command) != 0) {
				$attempt++;
				if ($attempt == 10) {die "Fatal Error: $this failed on $def after $attempt attempts!";}
			} else {
				last ATTEMPT;
			}
		}

		# Write genes to the GFF file

		# $this->create_gene_gff($gff_handle, $def, $output_file, $length);

		# Write to the concat raw file

		$header = sprintf("Sequence: %s <%dbp>", $entry->def, $length);
		$this->execute("(echo \"$header\"; cat $output_file; echo) >> $combo_output");

		unlink $input_file, $output_file;
	}
	close(FASTA);
	#close($gff_handle);
	close($input_ptr);
	close($output_ptr);

	$this->create_gene_gff();
}

############################################
# Make genes and proteins from a 
# gene gff file. 
############################################

sub make_genes_and_proteins_wrapper {
	my $this         = shift;
	my $name         = $this->name;
	my $assembly     = $this->assembly;
	my $genepred     = $this->genepred;
	my $genepred_id  = $this->get_id_by_name("gene_prediction", $genepred);

	# Input files

	my $fasta_file   = $this->fasta_filename;
	my $gff_file     = $this->gff;

	# Output files

	my $output_dir   = $this->output_dir;
	my $gene_file    = "$output_dir/$name.gene.fa";
	my $protein_file = "$output_dir/$name.protein.fa";
	
	$this->make_genes_and_proteins(GFF => $gff_file, FASTA => $fasta_file, GENE => $gene_file, PROTEIN => $protein_file);
}

# Dont panic since tell() is called in the def line
# When def line has been read, tell() is already at the
# next line, which is exactly what we want!

sub index_fasta_handle {
	my $this = shift;
	my $fh   = shift;
	my $just_saw_defline = 0;
	my $def;
	my $Index = {};
	while (<$fh>) {
		if (m/^>/) {
			chomp();
			s/^>//;
			s/\s.*//;
			$Index->{$_} = tell($fh);
		}
	}
	return $Index;
}

sub get_sequence {
	my $fh = shift;
	my $sequence = "";
	while (<$fh>) {
		last if m/^>/;
		chomp();
		$sequence .= $_;
	}
	$sequence =~ s/\s//g;
	return $sequence;
}

sub make_genes_and_proteins {

	use Smash::Utils::GFF qw(parse_gff text2flag);

	my $this = shift;
	my %args = @_;
	my ($gff_file, $fasta_file, $gene_file, $protein_file) = map {$args{$_}} qw(GFF FASTA GENE PROTEIN);

	my $gene_handle;
	my $protein_handle;
	my $fasta_handle;
	my $gff_handle;

	# init translator

	my $translator = Translator->new($this->translation_table);

	# Open GTF, protein and gene files

	open ($gene_handle, ">$gene_file") || die "Cannot open $gene_file: $!";
	open ($protein_handle, ">$protein_file") || die "Cannot open $protein_file: $!";
	open ($fasta_handle, "<$fasta_file") || die "Cannot open $fasta_file: $!";
	open ($gff_handle, "<$gff_file") || die "Cannot open $gff_file: $!";

	# index the fasta file

	my $FastaIndex = $this->index_fasta_handle($fasta_handle);

	# open the GFF file

	my $gff = new Smash::Utils::GFFlite($gff_handle);

	####
	# Scan the GFF file, sequence by sequence
	####

	while (my $features = $gff->nextSequence) {
		my $seqname = $features->[0]->seqname;
		my $filepos = $FastaIndex->{$seqname};
		seek($fasta_handle, $filepos, 0);
		my $sequence = get_sequence($fasta_handle);

		####
		# process those genes
		####

		foreach my $gene (@$features) {
			my $external_id = $gene->get_attribute("gene") || $gene->get_attribute("gene_id") || ("$seqname:".$gene->start."-".$gene->end);
			my $strand      = $gene->strand;
			my $location    = $gene->location;

			####
			# Make transcript and protein
			####

			my $transcript = "";
			foreach my $loc (sort {$a->start <=> $b->start} @$location) {
				$transcript .= uc(substr($sequence, $loc->start-1, $loc->end-$loc->start+1));
			}
			if ($strand eq "-") {
				$transcript =~ y/ACTG/TGAC/;
				$transcript = reverse($transcript);
			}

			####
			# Write to protein files, if CDS
			# We do this first, since the transcript might change if the frame is wrong.
			####

			if ($gene->feature eq "CDS") {

				# Fix frame

				my $start_codon = text2flag($gene->get_attribute("start_codon"));
				my $stop_codon  = text2flag($gene->get_attribute("stop_codon"));
				my $length      = length($transcript);
				if (!$start_codon && !$stop_codon && $length%3 != 0) { # incomplete gene, we have to find the right frame
					warn "Determining best frame for $external_id\n";

					my $frame   = $translator->get_best_frame($transcript);

					# Negative Frame indicates in-frame STOP codon

					if ($frame < 0) {
						warn "$external_id has internal STOP codons\n";
						$frame = -$frame;
					}

					$frame--;   # Frame is in [1,3], so convert it to [0,2]

					$transcript = substr($transcript, $frame, ($length-$frame) - ($length-$frame)%3);
				}

				print $protein_handle ">$external_id source:$seqname";
				foreach my $loc (sort {$a->start <=> $b->start} @$location) {
					printf $protein_handle " start:%d end:%d", $loc->start, $loc->end;
				}
				print $protein_handle " strand $strand\n";

				my $protein = $translator->translate($transcript, $start_codon, $external_id);
				print $protein_handle $this->pretty_fasta($protein);
			}

			####
			# Write transcripts for all genes
			# The transcript could have been fixed for frame errors inside the CDS conditional above.
			####

			print $gene_handle ">$external_id source:$seqname";
			foreach my $loc (sort {$a->start <=> $b->start} @$location) {
				printf $gene_handle " start:%d end:%d", $loc->start, $loc->end;
			}
			print $gene_handle " strand $strand\n";
			print $gene_handle $this->pretty_fasta($transcript);

		}
	}
	close($gene_handle);
	close($protein_handle);
	close($fasta_handle);
	close($gff_handle);
}

=item B<run>

Runs the instance of this GenePredictor.

=back

=cut

############################################
# Actual execution of the object
############################################

sub run {
	my $this  = shift;

	my $settings;

	# Can this program handle multifasta, and predict on each sequence independently?
	if ($this->is_multifasta_compatible) {
		# If it is indeed multifasta compatible, it should be able to take in one
		# multifasta and train on it, then take the same file in and predict on it.
		# These should not take more than a few lines of code, so we will let the
		# subclass implement it.
		# NOTE: Subclass must check $this->self_train() and use self-training 
		# procedure if it is available and requested.
		
		$this->predict_multifasta();
	} else {
		if ($this->self_train) {
			$settings = $this->get_self_trained_settings();
		} else {
			$settings = $this->get_generic_settings();
		}
		$this->predict($settings);
	}

	$this->make_genes_and_proteins_wrapper();
}

=head1 FUNCTIONS REQUIRED IN SUBCLASSES

The following functions must be implemented by subclasses, since these will be
queried during every instance they are used.

=over 4

=item B<is_multifasta_compatible()>

Returns C<1> if this program can handle multifasta files correctly (see later).

=item B<is_trainable>

Returns C<1> if this program is trainable.

=item B<parameter_settings>

Returns a unique representation of the parameter settings used in this instance.

=back

Each gene prediction programs behaves differently when given a multifasta file as input.
For example, when given such a file, GeneMark ignores all the fasta header lines and 
treats the input as a single sequence! (This behavior has changed in MetaGeneMark).
But MetaGene processes each sequence in the fasta file as a different sequence, as it
should. To accommodate both kinds of programs, Smash::Analyses::GenePredictor provides
two different types of wrappers. Every subclass should implement a function 
C<is_multifasta_compatible()>. Wrappers for programs that correctly handle multifasta 
files (like MetaGene)  should return C<1> for this call. Wrappers for programs that 
cannot, such as GeneMark, should return C<0>.

=head2 Multifasta Compatible Subclasses

Ideally, you need just one system call to make predictions on the whole file. If 
self-training is involved, that would make it perhaps two steps. Since this is 
straightforward, we require these subclasses to implement just one function. This
function, called C<predict_multifasta()> takes no arguments, but can access the
names of input file, output file, output GFF file from it's parameters. It is 
expected to make gene predictions on the input file, and write these predictions
in GFF format in the output GFF file. Additionally, the subclass should check if 
C<self_train> is set to C<1>, in which case, it should train parameters using the
input file if it has the ability to do so.

=over 4

=item B<predict_multifasta()>

Takes no arguments. Makes gene predictions on the input fasta file C<fasta_filename>, 
and writes the raw output from the program to C<output> and GFF format output to
C<gff>.

=back

=head2 Multifasta Incompatible Subclasses

Such programs are difficult to handle. You have to break the multifasta file into
many fasta files containing one sequence each, and then call the gene predictor on
each such file. Therefore, B<GenePredictor> provides these necessary steps, leaving
the subclass to concentrate on just running the program and parsing the output to
generate GFF files. At a very high level, B<GenePredictor> does the following:

	1. open $INPUT fasta file to read
	2. open $GFF_FH to write
	3. foreach sequence in $INPUT
	4.	decide $output to write the output to
	5. 	write sequence to $file
	6. 	$command = get_command_line($file, $output, $parameters)
	7. 	system($command);
	8. 	parse $output and write GFF format lines to $GFF_FH
	9. end
	
Thus any multifasta incompatible subclass of B<GenePredictor> must implement the following functions:

=over 4

=item B<get_command_line($input_file, $output_file, $parameters)>

Command to run given the name of the input file name containing DNA sequence in FASTA format. If the program needs to
decide the command line based on features of the DNA sequence (e.g., GeneMark uses a heuristic model depending on GC 
content of the sequence), the subclass is free to open the file to obtain the required information. It must close the
file as well after processing. A suggested implementation could be:

	sub get_command_line {
		my $this       = shift;
		my $filename   = shift;
		my $output     = shift;
		my $parameters = shift;
		my $option1    = process_something_in($filename);
		return "/somewhere/BestGenePredictor --option1=$option1 --option2 $parameters --input $filename > $output";
	}

The last argument ($parameters) is used to specify special parameters to get the command line. For example, one could
send different values for self-trained or generic parameter settings.

=item B<create_gene_gff>

Writes GFF format lines to the given filehandle after parsing the raw output of the gene prediction program. 
It is usually called as:

		$this->create_gene_gff($gff_handle, $actual_seqname, $raw_prediction_output);

Subclasses must parse the raw output in C<$raw_prediction_output> and write GFF lines containing gene
coordinates to the filehandle in C<$gff_handle>. The subclass B<must not> close tha handle. Since some gene
predictors do not write the sequence definition in the prediction output, the definition is passed as 
C<$actual_seqname>.

=back

=head2 Conditionally mandatory functions

Multifasta incompatible subclasses of B<GenePredictor> must also implement the following functions.

=over 4

=item B<train_min>

Minimum total length of input DNA that is required to train a model under this program.

=item B<get_self_trained_settings>

This function trains the gene predictor using the input data and returns a token that 
must be recognized by C<get_command_line()> which then will formulate the right 
command line to use the parameters that resulted from this training step.
A typical call to this function would be followed by a call to C<get_command_line>, 
like so:

	my $settings = $this->get_self_trained_settings();
	my $command  = $this->get_command_line($fasta_file, $settings);
	system($command);

In this example, C<get_command_line> should make a command that would use the settings returned by the previous call to 
C<get_self_trained_settings>. It is upto the subclass to implement the finer details of this process. For example,
since GeneMark has two steps in the self-training process (making the new model after self-training, then predicting
using that model), the wrapper module for GeneMark 
implements C<get_self_trained_settings> by training a model and making the command line options that will 
use that model. 

=item B<get_generic_settings>

This is similar to C<get_self_trained_settings> above, but returns the token that must
be recognized by C<get_command_line> to generate the correct commandline for non
self-train mode.

=back

=head2 Suggested functions

Any subclass of Smash::Analyses::GenePredictor is strongly advised to implement the following functions:

=over 4

=item B<new>

If the constructor needs to be more specific than the one provided in GenePredictor, the subclass could implement
its own constructor. But it should leave the B<bless>ing of the object to the parent class' constructor as the last
step in its constructor function. It should call the parent object's B<new> method using C<< $this->SUPER->new(@_) >>.
If the subclass does not implement B<new>, then the constructor for the GenePredictor is called.

=item B<min_length>

Minimum DNA sequence length required by the program to make gene predictions. If this function is not implemented, 
GenePredictor sets it to 0 so that every sequence will be passed to the gene predictor.

=back

=cut

sub min_length {
	return 0;
}

sub is_trainable {
	return 0;
}

sub train_min {
	my $this = shift;
	die "$this does not implement train_min()";
}

sub get_self_trained_settings {
	my $this = shift;
	die "$this does not implement get_self_trained_settings()";
}

sub get_command_line {
	my $this = shift;
	die "$this does not implement get_command_line()";
}

sub create_gene_gff {
	my $this = shift;
	die "$this does not implement create_gene_gff()";
}

sub parameter_settings {return undef;}

sub get_translator {
	my $translation_table = shift;
	my $translator = Translator->new($translation_table);
}

############################################
# Check for the programs
############################################

sub check_executables {
	my $predictor   = shift;
	my $program     = $predictor->program;
	my @executables = $predictor->get_executables;
	for my $exe (@executables) {
			# 127 - UNIX exit code when program not found
		if (system("$exe 1>/dev/null 2>/dev/null") == 127) {
			die "$exe program in $program package not found in path";
		}
	}
}

1;

############################################
#  Generic  gene prediction program handler
############################################

package Smash::Analyses::GenePredictor::external;
our @ISA = qw(Smash::Analyses::GenePredictor);
use strict;
use warnings;

sub run {
	die "Smash::Analyses::GenePredictor::external is a placeholder only";
}

sub _overload { return "external";}

1;

############################################
#  Generic  gene prediction program handler
############################################

package Translator;
use strict;
use warnings;

sub start_codons {shift->{START_CODONS}};
sub codon_table  {shift->{CODON_TABLE}};
sub translation_table  {shift->{TRANSLATION_TABLE}};

sub new {
	my $class = shift;
	my $translation_table = shift || die "Translator needs an explicit translation table!";
	my $self = bless {}, $class;

	if ($translation_table != 11 && $translation_table != 4 && $translation_table != 1) {
		die "Translation table $translation_table is not implemented!";
	}

	my $StartCodonTable = {};
	my $CodonTable  = 
		{ 
		TTT=>"F", TTC=>"F", TTA=>"L", TTG=>"L", TTN=>"X",
		TCT=>"S", TCC=>"S", TCA=>"S", TCG=>"S", TCN=>"S",
		TAT=>"Y", TAC=>"Y", TAA=>"*", TAG=>"*", TAN=>"X",
		TGT=>"C", TGC=>"C", TGA=>"*", TGG=>"W", TGN=>"X",
		CTT=>"L", CTC=>"L", CTA=>"L", CTG=>"L", CTN=>"L",
		CCT=>"P", CCC=>"P", CCA=>"P", CCG=>"P", CCN=>"P",
		CAT=>"H", CAC=>"H", CAA=>"Q", CAG=>"Q", CAN=>"X",
		CGT=>"R", CGC=>"R", CGA=>"R", CGG=>"R", CGN=>"R",
		ATT=>"I", ATC=>"I", ATA=>"I", ATG=>"M", ATN=>"X",
		ACT=>"T", ACC=>"T", ACA=>"T", ACG=>"T", ACN=>"T",
		AAT=>"N", AAC=>"N", AAA=>"K", AAG=>"K", AAN=>"X",
		AGT=>"S", AGC=>"S", AGA=>"R", AGG=>"R", AGN=>"X",
		GTT=>"V", GTC=>"V", GTA=>"V", GTG=>"V", GTN=>"V",
		GCT=>"A", GCC=>"A", GCA=>"A", GCG=>"A", GCN=>"A",
		GAT=>"D", GAC=>"D", GAA=>"E", GAG=>"E", GAN=>"X",
		GGT=>"G", GGC=>"G", GGA=>"G", GGG=>"G", GGN=>"G"
		};

	map {$StartCodonTable->{$_} = $CodonTable->{$_}} keys %$CodonTable;

	if ($translation_table == 1) {
		map {$StartCodonTable->{$_} = "M"} qw(ATG);
	} elsif ($translation_table == 4) {
		map {$StartCodonTable->{$_} = "M"} qw(TTA TTG CTG ATT ATC ATA ATG GTG);
		$CodonTable->{TGA} = "W";
	} elsif ($translation_table == 11) {
		map {$StartCodonTable->{$_} = "M"} qw(TTG CTG ATT ATC ATA ATG GTG);
	}

	$self->{CODON_TABLE} = $CodonTable;
	$self->{START_CODONS} = $StartCodonTable;
	$self->{TRANS_TABLE} = $translation_table;

	return $self;
}

=over 4

=item C<translate($dna, $has_start_codon, $name)>

returns the protein translation of the given DNA sequence ($name is
optional and is only used in reporting errors so that you can debug
which sequence had the error).

=item C<get_best_frame($transcript)>

returns the frame of the longest ORF (if the longest ORF still has
a STOP codon, returns negative frame). B<NOTE:> The range it returns is [1,3]
so that a negative number is visible. If the range was [0,2] then the
in-frame STOP codon diagnosis won't work. If you wanted it in [0,2]
just subtract one after taking the absolute value.

=back

=cut

sub translate {
	my $this = shift;
	my ($dna, $has_start_codon, $name) = @_;
	my $length  = length($dna);
	my $protein = "";
	my $AA = $this->codon_table;
	my $Start = $this->start_codons;

	if ($length%3 != 0) {
		die "Transcript $name length not a multiple of 3\n$dna";
	}

	my $offset = 0;
	if ($has_start_codon) {
		my $codon = uc(substr($dna, 0, 3));
		my $aa = $Start->{$codon};
		if (!defined($aa)) {
			$aa = "X";
		}
		$protein.=$aa;
		$offset += 3;
	}
	while ($offset < $length) {
		my $codon = uc(substr($dna, $offset, 3));
		my $aa    = $AA->{$codon};
		if (!defined($aa)) {
			$aa = "X";
		}
		$protein.=$aa;
		$offset += 3;
	}
	$protein =~ s/\*$//;
	return $protein;
}

sub get_best_frame {
	my $this = shift;
	my $transcript = shift;
	my $length     = length($transcript);
	my $min_stops  = $length;
	my $max_length = 0;
	my $best_frame = 0;
	for (my $i=0; $i<3; $i++) {
		my $tx  = substr($transcript, $i, ($length-$i) - ($length-$i)%3);
		my $ptx = $this->translate($tx, 0);
		my $longest = $this->get_longest_protein($ptx);
		my $in_frame_stops = $ptx =~ s/\*//g;
		if (length($longest) > $max_length) {
			$max_length = length($longest);
			$best_frame = $i;
			$min_stops  = $in_frame_stops;
		}
	}
	$best_frame++;
	if ($min_stops > 0) {
		return -$best_frame;
	} else {
		return $best_frame;
	}
}

sub get_longest_protein {
	my $this    = shift;
	my $protein = shift;
	my @pieces  = split(/\*+/, $protein);
	@pieces = sort {length($b) <=> length($a)} @pieces;
	return $pieces[0];
}

1;
