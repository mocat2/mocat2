############################################
#  GeneMark gene prediction program handler
############################################

package Smash::Analyses::GenePredictor::GeneMark;

use strict;
use warnings;

use base "Smash::Analyses::GenePredictor";

=head1 NAME

Smash::Analyses::GenePredictor::GeneMark - Implementation of GeneMark gene prediction software pipeline

=head1 SYNOPSIS

	my $instance = Smash::Analyses::GenePredictor::GeneMark->new
				(
					'NAME'        => $label,
					'OUTPUT_DIR'  => $output_dir,
					'PREDICTOR'   => 'GeneMark',
					'FASTA_FILE'  => $fasta,
					'SELF_TRAIN'  => 0,
					'VERSION'     => "2.6p"
				);
	$instance->run();

=head1 DESCRIPTION

Smash::Analyses::GenePredictor::GeneMark implements wrapper modules for the GeneMark and MetaGeneMark gene 
prediction software. GeneMark and MetaGeneMark are developed by 
L<http://exon.gatech.edu/GeneMark/|"Mark Borodowsky's group"> at Georgia Tech. You can download your own
copy of the software and the license from their website.

=cut

############################################
# Constructor
############################################

sub new1 {
	my $class  = shift;
	my %params = @_;
	my $self   = $class->SUPER::new(@_);
	$self->{PREDICTOR} = "GeneMark";
	return $self;
}

=head1 FUNCTIONS

=head2 Mandatory and conditionally mandatory functions

=over 4

=item B<is_trainable>

returns C<1>.

=item B<train_min>

returns C<200000>.

=item B<get_self_trained_settings>

Trains a prokaryotic gene model using the input fasta 
file and returns the command line option to use that 
model.

=item B<get_generic_settings>

Chooses a heuritic model with the right translation 
table and the GC content of the input sequence and 
returns the command line option to use that  model.

=item B<get_command_line>

Makes a command line with the model passed.

=item B<parameter_settings>

describes the rule based creation of parameter set
for commandline.

=item B<create_gene_gff>

=back

=cut

sub is_multifasta_compatible {return 0;}
sub is_trainable {return 1;}
sub train_min    {return 200000;}

sub train_combo_parameter {
	my $organism = shift;
	my $option;
	if ($organism eq "bacteria") {
		$option   = "--prok";
	} elsif ($organism eq "virus") {
		$option = "--virus";
	} elsif ($organism eq "eukaryote") {
		$option = "--euk";
	} elsif ($organism eq "mold_mitochondrial") {
		$option = "--gcode 4 --shape circular";
	}
	return $option;
}

sub get_self_trained_settings {
	my $this   = shift;
	my $name   = $this->name;
	my $output_dir = $this->output_dir;
	my $option = train_combo_parameter($this->organism);
	my $command = sprintf("cd $output_dir && %s/gmsn.pl $option --name=$name --maxitr=100 %s", $this->pkg_dir, $this->fasta_filename);
print "$command\n";
	$this->execute($command);
	return "-r -m $output_dir/${name}_hmm_combined.mod";
}

############################################
# Construct the model filename from 
# GC content and translation table
############################################

sub get_generic_settings {
	my $this        = shift;
	my $input_file  = shift;
	my $trans_table = $this->translation_table;
	my $modeldir    = $this->paramdir;

	my $model = "-m $modeldir/heu_${trans_table}_<__GC__>.mod";

	return $model;
}

sub get_command_line {
	my $this        = shift;
	my $input_file  = shift;
	my $output_file = shift;
	my $model_arg   = shift;
	my $pkg_dir     = $this->pkg_dir;
	my $exe         = $this->exe;

	if ($model_arg =~ /<__GC__>/) {
		open(FASTA, "<$input_file") || die "Cannot open $input_file: $!";
		my $fasta    = new FAlite(\*FASTA);
		my $gcp;
		if (my $entry = $fasta->nextEntry) {
			$gcp  = $this->get_gc_percent($entry->seq, 30, 70);
		}
		close(FASTA);
		$gcp   = int($gcp + 0.5); # quick hack for round
		$model_arg =~ s/<__GC__>/$gcp/;
	}


	my $command       = "cd $pkg_dir && ./$exe $model_arg -o $output_file $input_file";

	# Used this for E. coli complete genomes
	#my $command       = "cd $pkg_dir && ./$exe -m $modeldir/ecoli_r.mod -r -o $output_file $input_file";

	return $command;
}

sub parameters {
	my $this        = shift;
	my $trans_table = $this->translation_table;
	my $organism    = $this->organism;
	if ($this->self_train) {
		return "-m self_train_$organism -r";
	} else {
		return "-m heu_${trans_table}_gc";
	}
}

sub parameter_settings {
	my $this        = shift;
	my $trans_table = $this->translation_table;
	my $organism    = $this->organism;
	my $train_min   = $this->train_min;
	my $settings    = "";
	if ($this->self_train) {
		$settings .= "'-m self_train_$organism -r' for length > $train_min; ";
		$settings .= "'-m heu_${trans_table}_gc' otherwise";
	} else {
		$settings .= "'-m heu_${trans_table}_gc'";
	}
	return $settings;
}

############################################
# Process GeneMark's raw output
############################################

sub create_gene_gff_old {
	my $this          = shift;
	my $GFF_FH        = shift;
	my $def           = shift;
	my $pred_file     = shift;
	my $genepred      = $this->genepred;
	my $assembly      = $this->assembly;
	my $IN_FH;

	open($IN_FH, "<$pred_file") || die "Cannot open $pred_file: $!";

	# Skip the header section
	PRESCAN:while (<$IN_FH>) {
		if ($_ =~ /^\s*#\s*Length\s*$/) {
			last PRESCAN;
		}
	}

	# Now to the predictions
	####################
	# FIX THE FRAME
	# GENEMARK always lists complete codons, so frame is always 0
	###################

	$, = "\t";
	LINE:while (<$IN_FH>) {
		my $line = $_;

		# Skip empty predictions
		if ($line =~ /^-+$/) {
			last LINE;
		}

		# Parse the gene lines
		if ($line =~ /^\s*(\d+)\s+(\+|-)\s+(<*)(\d+)\s+(>*)(\d+)\s+/) {
			my ($gene_id, $strand, $langle, $start, $rangle, $end) = ($1, $2, $3, $4, $5, $6);
			my $incomplete5  = 0;
			my $incomplete3  = 0;
			my ($start_codon, $stop_codon);
			if ($langle eq "<") {
				$incomplete5 = 1;
			}
			if ($rangle eq ">") {
				$incomplete3 = 1;
			}
			if ($strand eq "+") {
				$start_codon = ($incomplete5 == 1)?"no":"yes";
				$stop_codon  = ($incomplete3 == 1)?"no":"yes";
			} else {
				$stop_codon  = ($incomplete5 == 1)?"no":"yes";
				$start_codon = ($incomplete3 == 1)?"no":"yes";
			}
			if ($def =~ /${assembly}\.(.*)/) {
				$gene_id = "$genepred.$1.G$gene_id";
			} else {
				$gene_id = "$genepred.$def.G$gene_id";
				#die "Invalid contig id: $def";
			}
			print $GFF_FH $def, "$this", "CDS", $start, $end, ".", $strand, ".", "gene_id \"$gene_id\"; start_codon \"$start_codon\"; stop_codon \"$stop_codon\";\n";
		} else {
			die "Incorrect format found in GeneMark output: $line";
		}
	}
	$, = "";
}

sub create_multifasta_gene_gff {
	shift->create_gene_gff();
}

sub create_gene_gff {
	my $this          = shift;
	my $genepred      = $this->genepred;
	my $assembly      = $this->assembly;
	my $gff_file      = $this->gff;
	my $pred_file     = $this->output;
	#my $gene_file     = $this->gene;
	#my $prot_file     = $this->protein;
	my $IN_FH;
	my $GFF_FH;
	#my $PROT_FH;
	#my $GENE_FH;

	#################
	# Ideally the Smash could save time since GeneMark also
	# generates gene and protein sequences. But it does not translate properly
	# using the codon table 11 properly. Non ATG start codons are not translated
	# as M. Until then, the following are disabled. This will let the standard
	# Smash translation machinery to take care of it more accurately.
	###

	open($IN_FH,  "<$pred_file")   || die "Cannot open $pred_file: $!";
	open($GFF_FH, ">$gff_file")    || die "Cannot open $gff_file: $!";
	#open($PROT_FH, ">$prot_file")  || die "Cannot open $prot_file: $!";
	#open($GENE_FH, ">$gene_file")  || die "Cannot open $gene_file: $!";

	my $parser = new Smash::Analyses::GenePredictor::GeneMarkParser(FH => $IN_FH, ASSEMBLY => $assembly, GENEPRED => $genepred);
	while (my $genes = $parser->parseNextSeq()) {
		foreach my $feature (@$genes) {
			if ($feature) {
				$feature->print_feature_gff($GFF_FH);
				#print $PROT_FH ">".$feature->name."\n";
				#print $PROT_FH Smash::Core->pretty_fasta($feature->get_property("prot_seq"));
				#print $GENE_FH ">".$feature->name."\n";
				#print $GENE_FH Smash::Core->pretty_fasta($feature->get_property("dna_seq"));
			}
		}
	}

	close($IN_FH);
	close($GFF_FH);
	#close($PROT_FH);
	#close($GENE_FH);

}

=head2 Suggested functions

=over 4

=item B<min_length>

returns C<41>.

=back

=cut

sub min_length   {return 41;}

=head2 Local functions

These functions are not known outside of the scope of GeneMark.

=over 4

=item B<check_license_key>

Checks for the GeneMark software license key and copies it to the user's home directory. This is required for GeneMark to
run properly. The license key file should be available as F<B<pkg_dir>/gm_key>. See L</DESCRIPTION> for more details.

=item B<paramdir>

Directory where the parameter files reside under B<pkg_dir> for GeneMark.

=back

=cut

sub exe          {return "gmhmmp";}

############################################
# Initialize pkg variables
############################################

sub init {
	my $this = shift;
	$this->SUPER::init();
	PLACES: foreach my $dir qw(modeldir heuristic_mod) {
		my $candidate = $this->pkg_dir."/$dir";
		if (-d $candidate) {
			$this->set_paramdir($candidate);
			last PLACES;
		}
	}
	$this->check_license_key();
}

############################################
# check for the license key in the user's home directory, or
# changed on 12.04.2009: version/key mismatches are possible.
#                        so dont check, just copy!
# copy the license key from smash directory to the user's home directory
############################################

sub check_license_key {
	use File::Copy;
	use Env qw(HOME);
	my $this = shift;
	my $pkg_dir = $this->pkg_dir;
	my $exe     = $this->exe;
	my $output  = `$pkg_dir/$exe 2>&1`;
	if ($output =~ /License key .* not found/si ||
	    $output =~ /Segmentation fault/si ||
	    $output =~ /Your .* license period has ended/si) {
			copy("$pkg_dir/gm_key", "$HOME/.gm_key") || die "Cannot copy license key from $pkg_dir/gm_key: $!";
	}
}

############################################
# GeneMark parameter file locations
############################################

sub paramdir {shift->{PARAMDIR}}
sub set_paramdir {
	my $this = shift;
	$this->{PARAMDIR} = shift;
}

############################################
# Return the command line to run GeneMark 
# on a sequence
# If $model is not defined, use heuristics
############################################

1;

=head1 LOCAL CLASSES

This module has the following two classes that are contained locally.

=cut

############################################
#  MetaGeneMark gene prediction program handler
############################################

=head1 NAME

Smash::Analyses::GenePredictor::MetaGeneMark - Subclass of 
L<Analyses::GenePredictor::GeneMark> implementing 
MetaGeneMark gene prediction software pipeline.

=head1 SYNOPSIS

	my $instance = Smash::Analyses::GenePredictor::GeneMark->new
				(
					'NAME'        => $label,
					'OUTPUT_DIR'  => $output_dir,
					'PREDICTOR'   => 'MetaGeneMark',
					'FASTA_FILE'  => $fasta,
					'SELF_TRAIN'  => 0,
					'VERSION'     => "current"
				);
	$instance->run();

=head1 DESCRIPTION

=cut

package Smash::Analyses::GenePredictor::MetaGeneMark;

use strict;
use warnings;

use base "Smash::Analyses::GenePredictor::GeneMark";

sub is_multifasta_compatible {return 1;}
sub is_trainable {return 0;}

sub predict_multifasta {
	my $this = shift;

	# instance options

	my $input_fasta = $this->fasta_filename;
	my $output_file = $this->output;
	my $gff_file    = $this->gff;

	# program options

	my $pkg_dir     = $this->pkg_dir;
	my $exe         = $this->exe;

	# Removed '-a -d' since translation is wrong at the start codon. We will do it ourselves until
	# they fix it.

	my $command     = "cd $pkg_dir && ./$exe -m MetaGeneMark_v1.mod -o $output_file $input_fasta";

	# Attempt to run it 10 times, if something goes wrong.
	# Give up after 10 attempts.

	my $attempt = 0;
	ATTEMPT:while ($attempt < 10) {
		if ($this->execute($command) != 0) {
			$attempt++;
			if ($attempt == 10) {die "Fatal Error: $this failed on $input_fasta after $attempt attempts!";}
		} else {
			last ATTEMPT;
		}
	}
	$this->create_multifasta_gene_gff();
}

1;

package Smash::Analyses::GenePredictor::GeneMarkParser;

use strict;
use warnings;
use Smash::Utils::GFF;

=head1 NAME

Smash::Analyses::GenePredictor::GeneMarkParser - Parser for GeneMark/MetaGeneMark output

=head1 SYNOPSIS

	my $parser = new Smash::Analyses::GenePredictor::GeneMarkParser(
			FH => \*STDIN, 
			ASSEMBLY => "MC1.MG1.AS1", 
			GENEPRED => "MC1.MG1.AS1.GP1"
			);
	while (my $genes = $parser->parseNextSeq()) {
		foreach my $feature (@$genes) {

			# print the gene feature as GFF.

			$feature->print_feature_gff(\*STDOUT);

			# print the protein sequence as fasta

			print ">".$feature->name."\n";
			print Smash::Core->pretty_fasta($feature->get_property("prot_seq");

			# print the DNA sequence as fasta

			print ">".$feature->name."\n";
			print Smash::Core->pretty_fasta($feature->get_property("dna_seq");
		}
	}

=head1 FUNCTIONS

=over 4

=cut

sub new {
	my $class  = shift;
	my %params = @_;
	my $this   = bless {%params}, $class;

	$this->getNextLine();
	return $this;
}

sub assembly  {shift->{ASSEMBLY}}
sub genepred  {shift->{GENEPRED}}
sub seqname   {shift->{SEQNAME}}
sub last_line {shift->{LAST_LINE}}

=item C<parseNextSeq>

parses the predicted genes for the next sequence, and returns a reference
to an array of L<Utils::GFF>::Feature's. These features have special
properties, C<prot_seq> and C<dna_seq>, which contain the protein
and DNA sequence of the predicted genes.

=cut

sub parseNextSeq {
	my $this = shift;
	my @genes;
	if ($this->moveToNextSeq) {
		#print "Processing ".$parser->seqname."\n";
		@genes  = $this->getGenes();
		if ($this->last_line =~ /^Predicted proteins:/) {
			my %Seqs = $this->getSequences();
			foreach my $gene (@genes) {
				$gene->set_property("prot_seq", $Seqs{$gene->name});
			}
		}
		if ($this->last_line =~ /^Nucleotide sequence/) {
			my %Seqs = $this->getSequences();
			foreach my $gene (@genes) {
				$gene->set_property("dna_seq", $Seqs{$gene->name});
			}
		}
	} else {
		return 0;
	}
	return \@genes;
}

=item C<moveToNextSeq()>

moves to the next sequence on a multifasta file prediction, or to the end
on a single fasta file prediction.

=cut

sub moveToNextSeq {
	my $this = shift;
	# Get the sequence name

	my $seqname;

	# if we have already moved to the next seq, just parse the line
	# else, move to the next seq and then parse the line

	my $line = $this->last_line;
	if ($line =~ /^Sequence: >(\S+)/) {
		$seqname = $1;
	} else {
		SEQNAME:while (my $line = $this->getNextLine()) {
			if ($line =~ /^Sequence: >(\S+)/) {
				$seqname = $1;
				last SEQNAME;
			}
		}
	}
	return 0 unless $seqname;
	$this->{SEQNAME} = $seqname;

	# Get the buffer filled up

	HEADERS:while (my $line = $this->getNextLine()) {
		last HEADERS if ($line =~ /^\s*#\s*Length/);
	}

	$this->{GeneMark2Smash} = {};

	return 1;
}

=item C<getGenes()>

returns an array containing L<Utils::GFF>::Feature objects, where
each feature is an actual predicted gene.

=cut

sub getGenes {
	my $this = shift;
	my $seqname  = $this->seqname;
	my $genepred = $this->genepred;
	my $assembly = $this->assembly;

	my $gene_count = 0;

	my @genes = ();

	if ($this->last_line !~ /^\s*#\s*Length/) {
		die "getGenes called at the wrong place after:\n".$this->last_line."\n";
	}

	####################
	# GENEMARK always lists complete codons, so frame is always 0
	###################

	while (my $line = $this->getNextLine()) {

		# Skip empty predictions

		return @genes if ($line =~ /^-+$/);

		# Parse the gene lines
		if ($line =~ /^\s*(\d+)\s+(\+|-)\s+(<*)(\d+)\s+(>*)(\d+)\s+/) {
			my ($gene_id, $strand, $langle, $start, $rangle, $end) = ($1, $2, $3, $4, $5, $6);
			   $gene_count = $gene_id;
			my $incomplete5  = 0;
			my $incomplete3  = 0;
			my ($start_codon, $stop_codon);
			if ($langle eq "<") {
				$incomplete5 = 1;
			}
			if ($rangle eq ">") {
				$incomplete3 = 1;
			}
			if ($strand eq "+") {
				$start_codon = ($incomplete5 == 1)?"no":"yes";
				$stop_codon  = ($incomplete3 == 1)?"no":"yes";
			} else {
				$stop_codon  = ($incomplete5 == 1)?"no":"yes";
				$start_codon = ($incomplete3 == 1)?"no":"yes";
			}
			my $smash_gene;
			if ($seqname =~ /${assembly}\.(.*)/) {
				$smash_gene = "$genepred.$1.G$gene_id";
			} else {
				$smash_gene = "$genepred.$seqname.G$gene_id";
				#die "Invalid contig id: $seqname";
			}

			my $location = new Smash::Utils::GFF::Location (START => $start, END => $end);
			my $feature  = new Smash::Utils::GFF::Feature(
						SEQNAME => $seqname,
						SOURCE  => "GeneMark",
						FEATURE => "CDS",
						LOCATION => [$location],
						STRAND  => $strand,
						NAME    => $smash_gene
						);
			$feature->set_attribute("gene_id",     $smash_gene);
			$feature->set_attribute("start_codon", $start_codon);
			$feature->set_attribute("stop_codon",  $stop_codon);
			#print $GFF_FH $seqname, "$this", "CDS", $start, $end, ".", $strand, ".", "gene_id \"$smash_gene\"; start_codon \"$start_codon\"; stop_codon \"$stop_codon\";\n";
			$this->{GeneMark2Smash}->{"gene_$gene_id"} = $smash_gene;
			push(@genes, $feature);
		} elsif ($line =~ /^Predicted proteins:/ || $line =~ /^Nucleotide sequence/ || $line =~ /^Sequence: /) {
			return @genes;
		} else {
			die "Incorrect format found in GeneMark output: $_";
		}
	}
}

=item C<getSequences()>

Parse the protein/nucleotide sequences of predicted genes. The calling
function needs to make sure that the file pointer is at the right 
position. Otherwise this will fail.

=cut

sub getSequences {
	my $this = shift;

	my $seq = "";
	my $prev_header;
	my %Seqs;
	LINE:while (my $line = $this->getNextLine()) {
		if ( $line =~ /^>(gene_\d+)\|/ ||
		     $line =~ /^Predicted proteins:/ || 
		     $line =~ /^Nucleotide sequence/ || 
		     $line =~ /^Sequence:/) {
			my $header = $1;
			$Seqs{$prev_header} = $seq unless $seq eq "";

			# no more genes!

			if ($line !~ /^>gene_/) {
				last LINE;
			}

			$prev_header = $this->{GeneMark2Smash}->{$header};
			$seq = "";
			next LINE;
		}
		$seq .= $line;
	}
	return %Seqs;
}

=item C<getNextLine()>

reads the next line and stores it in C<LAST_LINE>.

=cut

sub getNextLine {
	my $this = shift;
	my $fh   = $this->{FH};

	return 0 if eof $fh;

	# Skip empty lines

	my $line;
	LINE:while ($line = <$fh>) {
		last LINE if ($line !~ /^\s*$/);
	}

	return 0 if eof $fh;

	# here's a non empty line

	chomp($line);
	$this->{LAST_LINE} = $line;
	return $line;
}

=back

=cut

1;
