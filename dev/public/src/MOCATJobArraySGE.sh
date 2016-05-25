#!/bin/bash
# 
# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.
#
# SGE parameters
#$ -cwd
#$ -S /bin/bash
#$ -N MOCATJobArray
#$ -v MY_PREFIX=MOCATJob,MY_SUFFIX=sh
#$ -pe smp 8

# Job environment parameters
MY_HOST=`hostname`
MY_DATE=`date`
echo -e "STARTING JOB ID SGEJOBID : $MY_PREFIX.$SGE_TASK_ID.$MY_SUFFIX\ton $MY_HOST\t@ $MY_DATE" >> CWD/logs/SGEARRAY/startstop/MOCATJob_SGEARRAY.startstop.log

chmod 740 ./RJB.$SGE_TASK_ID.sh && ./RJB.$SGE_TASK_ID.sh

ERR=$?
MY_DATE=`date`
echo -e "FINISHED JOB ID SGEJOBID : $MY_PREFIX.$SGE_TASK_ID.$MY_SUFFIX\ton $MY_HOST\t@ $MY_DATE & EXIT STATUS $ERR" >> CWD/logs/SGEARRAY/startstop/MOCATJob_SGEARRAY.startstop.log
exit $ERR
