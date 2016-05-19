#!/bin/bash
# 
# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.
#
# PBS parameters
#PBS -S /bin/bash
#PBS -N MOCATJobArray
#PBS -v MY_PREFIX=MOCATJob,MY_SUFFIX=sh
#PBS -o ./logs/other/MOCATJob_PBSARRAY.SAMPLE.OUTPUT.log
#PBS -e ./logs/other/MOCATJob_PBSARRAY.SAMPLE.ERROR.log
#PBS -l select=ncpus=CPUS:mem=15gb
#PBS -J 1-MAXJOBS

# Job environment parameters
MY_HOST=`hostname`
MY_DATE=`date`

cd $PBS_O_WORKDIR
echo -e "STARTING JOB ID PBSJOBID : $MY_PREFIX.$PBS_ARRAY_INDEX.$MY_SUFFIX\ton $MY_HOST\t@ $MY_DATE" >> ./logs/PBSARRAY/startstop/MOCATJob_PBSARRAY.startstop.log
chmod 740 ./RJB.$PBS_ARRAY_INDEX.sh &&
./RJB.$PBS_ARRAY_INDEX.sh
ERR=$?
if [ $ERR -ne 0 ]
then
    MY_DATE=`date`
    echo -e "FINISHED JOB ID PBSJOBID : $MY_PREFIX.$PBS_ARRAY_INDEX.$MY_SUFFIX\ton $MY_HOST\t@ $MY_DATE & EXIT STATUS $ERR" >> ./logs/PBSARRAY/startstop/MOCATJob_PBSARRAY.startstop.log
    echo "Error $ERR. PBSERRORBYMOCAT! Something went wrong when executing script. Please check log files." >> ./logs/PBSARRAY/commands/MOCATJob_PBSARRAY.PBSDATE.command.log;
    exit $ERR
fi

MY_DATE=`date`
echo -e "FINISHED JOB ID PBSJOBID : $MY_PREFIX.$PBS_ARRAY_INDEX.$MY_SUFFIX\ton $MY_HOST\t@ $MY_DATE & EXIT STATUS $ERR" >> ./logs/PBSARRAY/startstop/MOCATJob_PBSARRAY.startstop.log
