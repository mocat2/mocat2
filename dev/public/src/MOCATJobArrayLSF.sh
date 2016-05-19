#!/bin/bash
# 
# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.
#
# LSF parameters
# --not working-- #BSUB -J "MOCATJobArray[1-MAXJOBS]"
# --not working-- #BSUB -n CPUS
# --not working-- #BSUB -K
# --not working-- #BSUB -o CWD/logs/other/MOCATJob_LSFARRAY.SAMPLE.OUTPUT.%I.log
# --not working-- #BSUB -e CWD/logs/other/MOCATJob_LSFARRAY.SAMPLE.ERROR.%I.log

# Job environment parameters
MY_HOST=`hostname`
MY_DATE=`date`
echo -e "STARTING JOB ID LSFJOBID : MOCATJob.$LSB_JOBINDEX.sh\ton $MY_HOST\t@ $MY_DATE" >> CWD/logs/LSFARRAY/startstop/MOCATJob_LSFARRAY.startstop.log
chmod 740 CWD/RJB.$LSB_JOBINDEX.sh &&
CWD/RJB.$LSB_JOBINDEX.sh
ERR=$?
if [ $ERR -ne 0 ]
then
    MY_DATE=`date`
    echo -e "FINISHED JOB ID LSFJOBID : MOCATJob.$LSB_JOBINDEX.sh\ton $MY_HOST\t@ $MY_DATE & EXIT STATUS $ERR" >> CWD/logs/LSFARRAY/startstop/MOCATJob_LSFARRAY.startstop.log
    echo "$LSB_JOBINDEX.failed" >> CWD/logs/LSF_specific/LSFARRAY.LSFDATE.status
    exit $ERR
fi
MY_DATE=`date`
echo -e "FINISHED JOB ID LSFJOBID : MOCATJob.$LSB_JOBINDEX.sh\ton $MY_HOST\t@ $MY_DATE & EXIT STATUS $ERR" >> CWD/logs/LSFARRAY/startstop/MOCATJob_LSFARRAY.startstop.log
echo "$LSB_JOBINDEX.completed" >> CWD/logs/LSF_specific/LSFARRAY.LSFDATE.status
