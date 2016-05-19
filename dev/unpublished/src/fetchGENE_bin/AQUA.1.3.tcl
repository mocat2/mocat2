#!/usr/bin/env tclsh 

#############################################################################
# AQUA: automatic quality improvment for multiple sequence alignment        #
#                                                                           #
# This script has been edited by Jens ROat Kultima, 2014                    #
# This specific implementation of AQUA is designed ot be                    #
# a part of the fetchGENE software package                                  #
#                                                                           #
# Copyright (C) 2009 Jean Muller (muller@embl.de)                           # 
#                                                                           #
# Please cite the following article:                                        #
#  Muller J, Creevey CJ, Thompson JD, Arendt D, Bork P.                     #
#  AQUA: automated quality improvement for multiple sequence alignments.    #
#  Bioinformatics.  2010 Jan 15;26(2):263-5. Epub 2009 Nov 19.              #
#                                                                           #
# This is part of AQUA source code.                                         #
#                                                                           #
# This program is free software; you can redistribute it and/or             #
# modify it under the terms of the GNU General Public License               # 
# as published by the Free Software Foundation; either version 2            # 
# of the License, or (at your option) any later version.                    #
#                                                                           #
# This program is distributed in the hope that it will be useful,           # 
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU General Public License for more details.                              #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with this program; if not, write to the Free Software               #
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
#                                                                           #
# Contact information                                                       #
# EMBL - Computational Biology Unit (Bork group)                            #
#                                                                           #
# Meyerhofstrasse 1, D-69117 Heidelberg, Germany                            #
# email: muller@embl.de                                                     #
#                                                                           #
#############################################################################

proc Main {} {

    #This is the main procedure that setup and runs the different steps

    global env
    global argv

    set PATHMUSCLE ""
    set PATHMAFFT  ""
    set PATHCLUSTALO ""
    set PATHRASCAL ""
    set PATHNORMD  ""

    if {$argv=={}} {DisplayUsage;return ""}

    set FastaFile [lindex $argv 0]
    set OutputDIR [lindex $argv 1]
    set fetchGENE_bin [lindex $argv 2]
    set use_clustalo [lindex $argv 3]
    set use_mafft [lindex $argv 4]
    set use_muscle [lindex $argv 5]

    #Configure program path
    #######################
    set PATHMUSCLE [file join $fetchGENE_bin/muscle3.8.31_i86linux64/muscle3.8.31_i86linux64]
    set PATHMAFFT [file join $fetchGENE_bin/mafft-7.127-without-extensions/scripts/mafft]
    set env(MAFFT_BINARIES) [file join $fetchGENE_bin/mafft-7.127-without-extensions/binaries/]
    set PATHRASCAL [file join $fetchGENE_bin/rascal1.34/rascal]
    set PATHNORMD [file join $fetchGENE_bin/norMD1_3/normd]
    set PATHCLUSTALO [file join $fetchGENE_bin/clustalo-1.2.0-Ubuntu-x86_64/clustalo-1.2.0-Ubuntu-x86_64]
    #######################

    #0 Various checking for input and output files
    if {$FastaFile=="" || ![file exists $FastaFile]} {puts "No input file provided.";DisplayUsage;return ""}
    if {[CheckFastaFormat $FastaFile]=="0"} {puts "$FastaFile does not seem to be in fasta format.";DisplayUsage;return ""}

    if {$OutputDIR==""} {set OutputDIR [file dirname $FastaFile]}

    set NameOfFile [file tail $FastaFile]
    #if {[regexp {\.[a-zA-Z]+$} $NameOfFile]} {regsub {\.[a-zA-Z]+$} $NameOfFile "" NameOfFile}
    set OutputFile [file join $OutputDIR ${NameOfFile}.best]

    #1- Setting up program PATHs
    #
    #1-1 MUSCLE
    #
    #if you wish not to use MUSCLE uncomment next line
    if {$use_muscle=="0"} {
	set NoMUSCLE 1
    }
    set Message    ""
    if {$PATHMUSCLE==""} {
	#Otherwise guessing path
	if {[catch {set PATHMUSCLE [exec which muscle]} Message]} {puts "Warning No MUSCLE binary found";set NoMUSCLEFound 1}
    }
    
    #1-2 MAFFT
    #
    #if you wish not to use MAFFT uncomment next line
    if {$use_mafft=="0"} {
	set NoMAFFT 1
    }
    set Message   ""
    if {$PATHMAFFT==""} {
	if {[catch {set PATHMAFFT [exec which mafft]} Message]} {puts "Warning No MAFFT binary found.";set NoMAFFTFound 1}
    }
	
    #1-3 CLUSTALO
    #
    #if you wish not to use CLUSTALO uncomment next line
    if {$use_clustalo=="0"} {
	set NoCLUSTALO 1
    }
    set Message   ""
    if {$PATHCLUSTALO==""} {
	if {[catch {set PATHCLUSTALO [exec which clustalo]} Message]} {puts "Warning No CLUSTALO binary found.";set NoCLUSTALOFound 1}
    }
    
    #1-4 RASCAL
    set Message    ""
    if {$PATHRASCAL==""} {
	if {[catch {set PATHRASCAL [exec which rascal]} Message]} {puts "No RASCAL binary found. Full stop.";return ""}
    }
    
    #1-5 NORMD
    set Message   ""
    if {$PATHNORMD==""} {
	if {[catch {set PATHNORMD [exec which normd]} Message]} {puts "No NORMD binary found. Full stop.";return ""}
    }

    if {[info exists NoMUSCLEFound] && [info exists NoMAFFTFound] && [info exists NoCLUSTALOFound]} {puts "No alignment programs binaries found. Full stop.";return ""}
    
    set L_MSA  {}
    set L_MSAR {}

    #2- Running programs
    #   First the 3 aligners
    #   Second the refiner on all MSA computed
    #   Third assessing the score of each version of MSA (refined and unrefined)
    #
    #2-1 MUSCLE
    if {![info exists NoMUSCLE]} {
	set Muscle [file join $OutputDIR ${NameOfFile}.muscle]

	lappend L_MSA $Muscle

	if {![info exists NoMUSCLEFound] && ![file exists $Muscle]} {
	    set MUSCLE_cmd "$PATHMUSCLE -in $FastaFile -out $Muscle -quiet"
	    #puts $MUSCLE_cmd
	    if {[catch {eval exec $MUSCLE_cmd} Message]} {
		puts stderr "stderr $Message" 
		if {[file exists $Muscle]} {file delete -force $Muscle}
	    } else {
		#puts "Done MUSCLE: $FastaFile"
	    }
	}
    }

    #2-2 MAFFT 
    if {![info exists NoMAFFT]} {
	set Mafft [file join $OutputDIR ${NameOfFile}.mafft]
	
	lappend L_MSA $Mafft
	
	if {![info exists NoMAFFTFound] && ![file exists $Mafft]} {
	    set MAFFT_cmd "$PATHMAFFT $FastaFile"
	    #puts $MAFFT_cmd
	    set  CODE [catch {eval exec $MAFFT_cmd > $Mafft} Message]
	    if {$CODE != "1"} {
		puts stderr "stderr $Message" 
		if {[file exists $Mafft]} {file delete -force $Mafft}
	    } else {
		#puts "Done MAFFT: $FastaFile"
	    }
	}
    }
	
    #2-3 CLUSTALO 
    if {![info exists NoCLUSTALO]} {
	set Clustalo [file join $OutputDIR ${NameOfFile}.clustalo]

	lappend L_MSA $Clustalo

	if {![info exists NoCLUSTALOFound] && ![file exists $Clustalo]} {
	    set CLUSTALO_cmd "$PATHCLUSTALO -i $FastaFile -o $Clustalo --force"
	    #puts $CLUSTALO_cmd
	    if {[catch {eval exec $CLUSTALO_cmd} Message]} {
		puts stderr "stderr $Message" 
		if {[file exists $Clustalo]} {file delete -force $Clustalo}
	    } else {
		#puts "Done CLUSTALO: $FastaFile"
	    }
	}
    }

    #2-4 RASCAL
    foreach MSAFile $L_MSA {
	set RASCALOUT [file join ${MSAFile}.rascal]
	
	lappend L_MSAR $RASCALOUT

	if {![file exists $MSAFile] || [file exists $RASCALOUT]} {continue}
	
	set RASCAL_cmd "$PATHRASCAL $MSAFile $RASCALOUT"
	#puts $RASCAL_cmd
	if {[catch {eval exec $RASCAL_cmd} Message]} {
	    puts stderr "$Message" 
	    if {[file exists $RASCALOUT]} {file delete -force $RASCALOUT}
	} else {
	    #puts "Done RASCAL: $MSAFile"
	}
    }
    
    set L_MSA [concat $L_MSA $L_MSAR]

    #2-5 NORMD
    set L_NORMD {}

    foreach MSAFile $L_MSA {
	set NORMDOUT [file join ${MSAFile}.normd]
	
	lappend L_NORMD $NORMDOUT

	if {![file exists $MSAFile] || [file exists $NORMDOUT]}  {continue}

	set NORMD_cmd "$PATHNORMD $MSAFile"
	#puts $NORMD_cmd
	if {[catch {eval exec $NORMD_cmd > $NORMDOUT} Message]} {
	    puts stderr "$Message" 
	    if {[file exists $NORMDOUT]} {file delete -force $NORMDOUT}
	} else {
	    #puts "Done NORMD: $MSAFile"
	}
    }
    
    set L_NORMD_AVAILABLE {}
    foreach F $L_NORMD {
	if {![file exists $F]} {
            puts "One normd file does not exists yet ([file tail $F])"
        } else {
            lappend L_NORMD_AVAILABLE $F
        }
    }

    #3- Selection of the best MSA computed

    set  BEST [Get_Best_NORMD $L_NORMD_AVAILABLE]
    if {$BEST!= "-1"} {
		regsub ".normd$" $BEST "" BEST
	
		#Creating either a link or copy of the file
		if {![file exists $OutputFile]} {
			cd [file dirname $OutputFile]
	
			file link -symbolic [file tail $OutputFile] [file tail $BEST]
			#file link -symbolic $OutputFile $BEST
			#file copy -force $BEST $Output
		} else {
			puts "$OutputFile already done."
		}
    } else {
	puts "Problem with one of the normd files"
    }
    
    return "$OutputFile"
}

proc DisplayUsage {} {
    
    puts "" 
    puts "AQUA - v1.3"
    puts ""
    puts "Usage is:"
    puts ""
    puts "./AQUA.tcl <FastaFile> <OutputDIR> (optional)"
    puts ""
    puts "<FastaFile> is your sequence inputfile in fasta format."
    puts "<OutputDIR> is the directory where to store the results. Default is setup to the same directory as <FastaFile>."
    puts ""

    return ""
}

proc Get_Best_NORMD {L_File} {

    #Return the File with the Best NORMD Score

    set BestScore "-999999999999999"
    set BestFile ""

    set Counter 0

    foreach File $L_File {
	if {![file exists $File]} {continue}

	#puts "$File"
	incr Counter

	set  Score [Parser_NORMD $File]
	#if {$Score=="" || ![regexp {[0-9]+} $Score]} {return "-1"}
	if {$Score=="" || ![regexp {[0-9]+} $Score]} {continue}

	if {$Score > $BestScore} {
	    set BestScore $Score
	    set BestFile  $File
	}
    }
    if {$BestFile==""} {set BestFile [lindex $L_File 0]}

    if {$Counter < 1} {return "-1"}

    return $BestFile
}

proc Parser_NORMD {FileIn} {

    #Returns NORMD score

    if {! [file exists $FileIn]} {return ""}

    set NormdScore ""
    set F [open "$FileIn"]
    while {[gets $F Line]>=0} {
	if {! [regexp -nocase {[a-z0-9]+} $Line] || [regexp -nocase {^[\#]} $Line]} {continue}
	#scan $Line "%s %s" Name NormdScore
	scan $Line "%s" NormdScore
    }
    close $F

    return $NormdScore
}

proc CheckFastaFormat {FileIn} {

    #Quick basic check if the input file is a fasta file

    set nbLine 0
    set nbSuperiorTo 0

    set F [open "$FileIn"]
    while {[gets $F Line]>=0} {
	if {! [regexp -nocase {[a-z0-9]+} $Line]} {continue}

	incr nbLine

	if {[string range $Line 0 0] == ">"} {incr nbSuperiorTo}

	#At least one in the first 3 lines otherwise quit
	if {$nbSuperiorTo=="0" && $nbLine>3} {break}
	#If 2 then ok enough quit
	if {$nbSuperiorTo>2} {break}
    }
    close $F

    return $nbSuperiorTo
}

Main

