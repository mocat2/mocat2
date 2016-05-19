#############################################################################
# AQUA: automatic quality improvment for multiple sequence alignment        #
#                                                                           #
# Copyright (C) 2009 Jean Muller (muller@embl.de)                           # 
#############################################################################

Please cite the following article:
----------------------------------
Muller J, Creevey CJ, Thompson JD, Arendt D, Bork P.
AQUA: automated quality improvement for multiple sequence alignments.
Bioinformatics.  2010 Jan 15;26(2):263-5. Epub 2009 Nov 19.

#Installation
#############

1- Download relevant programs

AQUA is Tcl/Tk script running at least with a Tcl/Tk v8.5 installation. Most of Unix and Linux distribution have this by defaults but it can be downloaded (see links below).
AQUA is calling by default 4 programs MUSCLE, MAFFT, NORMD, RASCAL which are available either in the AQUA repository (http://www.bork.embl.de/Docu/AQUA) or at there official repositories (see links below). 

Tcl/Tk can be downloaded from http://www.activestate.com/activetcl/ for any architecture. Make sure you have at least v8.5.
MUSCLE can be downloaded from http://www.drive5.com/muscle/
MAFFT  can be downloaded from http://align.bmr.kyushu-u.ac.jp/mafft/software
RASCAL can be downloaded from ftp://ftp-igbmc.u-strasbg.fr/pub/RASCAL/
NORMD  can be downloaded from ftp://ftp-igbmc.u-strasbg.fr/pub/NORMD/

It is to notice that the underlying protocol of AQUA is not directly linked to any programs and could be setup with other alignment programs, refiners. Tcl/Tk is multi-platform and thus AQUA can run on all platforms where the selected programs to run with can run.

2- Install the programs

2-1 Addiotional programs
Compile if necessary each of the additional programs. Make sure the different programs (MUSCLE, MAFFT, NORMD, RASCAL) are available in your path.
By default, AQUA will try to guess their paths if not specified explicitely in the AQUA.tcl script. To do so, open the AQUA.tcl script and go into the procedure named "Main" and search for the variables named PATHNAMEOFTHEPROGRAM (i.e. PATHMUSCLE), uncomment the line (e.g. remove the "#") and replace the fake path by yours.

Please have a look from line 48 in the "Configure program path" section where you explicitely put your own program paths.

for example if you path to the MUSCLE program is /home/muller/bin/muscle:
#set PATHMUSCLE [file join /DIR/muscle]
should be changed into
#set  PATHMUSCLE [file join /home/muller/bin/muscle]
and finally uncommented
set  PATHMUSCLE [file join /home/muller/bin/muscle]

AQUA is running with the default parameters for each programs. One can easely modify the script to reflect specific parameters.

2-2 AQUA
Make sure that the AQUA.tcl is executable (i.e. on unix chmod a+x)

#Usage
######

1- Input files

The input files for AQUA are mutliple sequences in FASTA format. The name of the input file is used to created the output files.

2- Output files

The output files from the alignment programs are all in FASTA format. The NORMD files are text files with a simple score.

Output files are: 
    <FastaFile>.muscle (for MUSCLE output in fasta format) 
    <FastaFile>.mafft  (for MAFFT output in fasta format)
    <FastaFile>.PROGRAM.rascal (for any refined version in fasta format) 
    <FastaFile>.normd  (for the NORMD score)
    <FastaFile>.best   (for the final MSA in fasta format)

3 - Running AQUA

AQUA is running in command line.

To run AQUA, simply use the following commands:
    ./AQUA.tcl FastaFile OutputDIR (optional)
    <FastaFile> is the sequence input file in fasta format
    <OutputDIR> is the directory where to store the results. Default is setup to the same directory as <FastaFile>. 


