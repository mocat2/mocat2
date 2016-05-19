############################################
#  MetaGene gene prediction program handler
############################################

package Smash::Analyses::GenePredictor::MetaGene;

use strict;
use warnings;

use base "Smash::Analyses::GenePredictor";

=head1 NAME

Smash::Analyses::GenePredictor::MetaGene - Implementation of MetaGene gene prediction software pipeline

=head1 SYNOPSIS

	my $instance = Smash::Analyses::GenePredictor::MetaGene->new
				(
					'NAME'        => $label,
					'OUTPUT_DIR'  => $output_dir,
					'PREDICTOR'   => 'MetaGene',
					'FASTA_FILE'  => $fasta,
					'SELF_TRAIN'  => 0,
					'VERSION'     => "2008-08-19"
				);
	$instance->run();

=head1 DESCRIPTION

This module implements the wrapper around the MetaGeneAnnotator program.
MGA will treat multiple sequences to be from the same species when C<"-s">
is given. Otherwise, each sequence is treated individually.

=cut

############################################
# Constructor
############################################

sub new {
	my $class  = shift;
	my %params = @_;
	my $self   = $class->SUPER::new(@_);
	return $self
}

sub _overload  {return "MetaGene";}

=head1 FUNCTIONS

Smash has switched to the MetaGeneAnnotator program from the original MetaGene
program.

=head2 MetaGene specific functions

=over 4

=item B<program>

Find the correct version of the software. Smash first gets the architecture of the host
and then looks for the right version. Since most 64bit linux versions can also execute
32bit code, it will fall back to 32bit versions if the 64bit version is not available on
a 64bit machine.

=back

=cut

sub program {
	my $this = shift;
	my $pkg_dir = $this->pkg_dir;

	# machine specific search

	my $target = `uname`;
	chomp($target);

	# only Linux is accepted

	if ($target ne "Linux") {
		warn "WARNING: Only Linux platform is supported by Smash\n";
		$this->abort();
	}

	# now get more specific

	my @files; 
	$target = `uname -m`;
	chomp($target);
	if ($target eq "x86_64") {
		@files = qw(mga_linux_ia64 mga_linux_ia32);
	} elsif ($target eq "i686") {
		@files = qw(mga_linux_ia32);
	} else {
		die "Unknown Linux architecture found: $target";
	}

	foreach my $f (@files) {
		my $program = "$pkg_dir/$f";
		return $program if (-e $program);
	}
	warn "WARNING: Cannot find one of (".join(",", @files).") in $pkg_dir\n";
	$this->abort();
}

=head2 Mandatory and conditionally mandatory functions

=over 4

=item B<is_trainable>

returns C<1>.

=item B<train_min>

returns C<200000>.

=item B<parameter_settings>

describes the rule based creation of parameter set
for commandline.

=item B<create_multifasta_gene_gff>

parses the output from MetaGene and creates a GFF file
containing gene coordinates.

=back

=cut

sub installation_name {'metagene'}
sub is_multifasta_compatible {return 1;}
sub is_trainable {return 1;}
sub train_min    {return 200000;}

sub parameter_settings {
	my $this      = shift;
	my $train_min = $this->train_min;
	my $settings  = "";

	if ($this->self_train) {
		$settings .= "'-s' for length > $train_min; ";
		$settings .= "'-m' otherwise";
	} else {
		$settings .= "'-m'";
	}
	return $settings;
}

############################################
# Initialize pkg variables
############################################

sub init {
	my $this = shift;
	$this->SUPER::init();
}

############################################
# Run MetaGeneAnnotator on a multifasta
############################################

sub predict_multifasta {
	my $this = shift;

	# instance options

	my $input_fasta = $this->fasta_filename;
	my $output_file = $this->output;
	my $gff_file    = $this->gff;

	# program options

	my $pkg_dir     = $this->pkg_dir;
	my $program     = $this->program;
	my $species_option;
	if ($this->self_train) {
		$species_option = "-s";
	} else {
		$species_option = "-m";
	}

	my $command = "$program $species_option $input_fasta > $output_file";

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
	$this->create_multifasta_gene_gff($gff_file, $output_file);
}

############################################
# Process MetaGene's raw output
############################################

sub create_multifasta_gene_gff {
	my $this          = shift;
	my $genepred      = $this->genepred;
	my $assembly      = $this->assembly;
	my $gff_file      = $this->gff;
	my $pred_file     = $this->output;
	my $IN_FH;
	my $GFF_FH;

	open($IN_FH,  "<$pred_file") || die "Cannot open $pred_file: $!";
	open($GFF_FH, ">$gff_file")  || die "Cannot open $gff_file: $!";

	my @headers;
	my $line;

	do {
		# Parse the header section
		# grab all lines starting with #

		HEADER:while ($line = <$IN_FH>) {
			chomp($line);
			if ($line =~ /^\s*#\s*/) {
				push(@headers, $line);
			} else {
				last HEADER;
			}
		}

		# Last line read was not a header, so save it.
		# And reset the header line list

		my $def = $headers[-3];
		$def =~ s/^\s*#\s*//;
		@headers = ();

		# Now to the predictions
		####################
		# FIX THE FRAME
		# METAGENE has detailed information on the frame and completenes
		# so it's now easy to parse that information
		###################

		$, = "\t";
		GENE:for (;$line;$line = <$IN_FH>) { # keep reading while $line defined
			chomp($line);
			if ($line =~ /^\s*#\s*/) {
				push(@headers, $line);
				last GENE;
			}
			# Parse the gene lines

			#my @fields = ($line =~ /^gene_(\d+)\s+(\d+)\s+(\d+)\s+(\+|-)\s+([012])\s+([01]{2})\s+([\d\.]+)\s+([sbap-])\s+(\d+|-)\s+(\d+|-)\s+([\d\.\-]+)$/);

			my @fields = split(/\s+/, $line);
			die "Incorrect format found: $line" unless @fields == 11;

			my ($gene_id, $start, $end, $strand, $frame, $status, $score, $model, $rbs_start, $rbs_end, $rbs_score) = split(/\s+/, $line);

			my $origin;
			$origin = "bacterial" if $model eq "b";
			$origin = "archaeal"  if $model eq "a";
			$origin = "phage"     if $model eq "p";

			my ($start_codon, $stop_codon) = split(//, $status);

			if ($strand eq "+") {
				# if incomplete in 5' of gene, adjust start using frame
				if ($start_codon == 0) {
					$start += $frame;
				}

				# if incomplete in 3' of gene, adjust end using length%3
				if ($stop_codon == 0) {
					$end = ($end - ($end-$start+1)%3);
				}
			} else {
				# if incomplete in 5' of gene, adjust end using frame
				if ($start_codon == 0) {
					$end -= $frame;
				}

				# if incomplete in 3' of gene, adjust start using length%3
				if ($stop_codon == 0) {
					$start = ($start + ($end-$start+1)%3);
				}
			}

			$gene_id =~ s/gene_//;
			if ($def =~ /${assembly}\.(.*)/) {
				$gene_id = "$genepred.$1.G$gene_id";
			} else {
				$gene_id = "$genepred.$def.G$gene_id";
				#die "Invalid contig id: $def";
			}

			my $features = "gene_id \"$gene_id\"; start_codon \"$start_codon\"; stop_codon \"$stop_codon\";";
			if ($origin) {
				$features .= " origin \"$origin\";";
			}
			print $GFF_FH $def, "$this", "CDS", $start, $end, ".", $strand, ".", "$features\n";
		}
	} while ($line);
	$, = "";

	close($GFF_FH);
	close($IN_FH);
}

=head2 Suggested functions

=over 4

=item B<min_length>

returns C<41>.

=back

=cut

sub min_length   {return 41;}

1;
