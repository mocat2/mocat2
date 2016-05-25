#! /bin/sh

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

if [ $# -lt 3 ]; then
        echo "sh $0 [scr folder] [correct result] [name] [output]"
        exit;
fi
#############################################################################
### echo $3."corrcted_sacftig" >> $4 ### Commented out to support previous file structures
perl $1/MOCATAssemblyRevision_contig_stat.pl $2 0:500:1000 >> $4
