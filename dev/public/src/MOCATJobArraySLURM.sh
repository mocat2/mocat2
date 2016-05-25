#!/bin/bash
# 
# This code is part of the MOCAT analysis pipeline
# This code is released under GNU GPL v3.

#SBATCH --job-name=NAME
#SBATCH --export=ALL
#SBATCH --distribution=block
#SBATCH --mincpus=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-MAXJOBS

    MY_HOST=`hostname`
MY_DATE=`date`
echo -e "STARTING JOB ID SLURMJOBID : MOCATJob.$SLURM_ARRAY_TASK_ID.sh\ton $MY_HOST\t@ $MY_DATE" >> CWD/logs/SLURMARRAY/startstop/MOCATJob_SLURMARRAY.startstop.log
chmod 740 CWD/RJB.$SLURM_ARRAY_TASK_ID.sh &&
CWD/RJB.$SLURM_ARRAY_TASK_ID.sh
ERR=$?
if [ $ERR -ne 0 ]
then
    MY_DATE=`date`
    echo -e "FINISHED JOB ID SLURMJOBID : MOCATJob.$SLURM_ARRAY_TASK_ID.sh\ton $MY_HOST\t@ $MY_DATE & EXIT STATUS $ERR" >> CWD/logs/SLURMARRAY/startstop/MOCATJob_SLURMARRAY.startstop.log
    echo "$SLURM_ARRAY_TASK_ID.failed" >> CWD/logs/SLURM_specific/SLURMARRAY.SLURMDATE.status
    exit $ERR
fi
MY_DATE=`date`
echo -e "FINISHED JOB ID SLURMJOBID : MOCATJob.$SLURM_ARRAY_TASK_ID.sh\ton $MY_HOST\t@ $MY_DATE & EXIT STATUS $ERR" >> CWD/logs/SLURMARRAY/startstop/MOCATJob_SLURMARRAY.startstop.log
echo "$SLURM_ARRAY_TASK_ID.completed" >> CWD/logs/SLURM_specific/SLURMARRAY.SLURMDATE.status

