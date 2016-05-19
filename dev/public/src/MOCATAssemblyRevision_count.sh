if [ $# -ne 2 ];then 
    echo "Usage $0 [scr folder] [outdir]"
    exit
fi

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

######################################################################################################
find $2/correct_Result/*log >$2/revised.log.list
find $2/break_contig/*num >$2/contig_break_result.log.list
######################################################################################################
perl $1/MOCATAssemblyRevision_break_result_stat.pl  $2/contig_break_result.log.list $2/contig_break_result.stat

######################################################################################################
perl $1/MOCATAssemblyRevision_correction_count.pl $2/revised.log.list $2/revision.count.stat

######################################################################################################
perl $1/MOCATAssemblyRevision_combine.pl $2/revision.count.stat $2/contig_break_result.stat $2/contig_correct.stat

######################################################################################################
rm $2/revised.log.list
rm $2/contig_break_result.log.list
rm $2/contig_break_result.stat
rm $2/revision.count.stat
