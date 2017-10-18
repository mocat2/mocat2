#!/bin/bash

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

# USER EMAIL (enter your email address here)
EMAIL=''

# Inital configuration
CWD=`pwd`
USER=`whoami`
SUBMITDATE=`date`
PROJECT=`echo $CWD | sed 's/.*\///'`

# Load arguments
LOAD_SF=0
LOAD_CFG=0
LOADED_ANY=1
SAMPLE_FILE=""
CFG_FILE=""
for argument in "$@"
do
	
	# Check
	if [ "$LOADED_ANY" == "1" ]; then
		LOADED_ANY=0
	fi
	
	# Load
	if [ "$LOAD_SF" == "1" ]; then
		SAMPLE_FILE="$argument"
		LOAD_SF=0
	fi
	if [ "$LOAD_CFG" == "1" ]; then
		CFG_FILE="$argument"
		LOAD_CFG=0
	fi

	# Parse
	if [ "$argument" == "-sf" ]; then
		LOAD_SF=1
		LOADED_ANY=2
	fi
	if [ "$argument" == "-cfg" ]; then
		LOAD_CFG=1
		LOADED_ANY=2
	fi
done

# Welcome message
echo ""
echo " ##############################################################################"
echo " #                     WELCOME TO THE MOCAT EXECUTER v1.3                     #"
echo " ##############################################################################"
echo ""
echo "   This shell script is used to execute a number of MOCAT commands in a row."
echo "   Typically this is used to process raw reads up to final taxonomic or mOTU"
echo "   profiles. Of course you can process each step individually using MOCAT.pl"
echo "   but we have created this software for your ease to execute these commands"
echo "   with ease without prior knowledge of how to run MOCAT.    Ps. Remember to"
echo "   setup your email in the script file: MOCAT/src/runMOCAT.sh         Enjoy!"

echo ""
echo "          Usage: runMOCAT.sh [-sf SAMPLE_FILE -cfg CONFIG_FILE]"
echo ""

# Exit if incorrect parameters
if [ "$LOADED_ANY" == "0" ]; then
	echo " If any other options than the ones above are entered,"
	echo " this software will not launch. Please specify only"
	echo " the correct options, or no options at all."
	echo " If no options are specified, this software will attempt"
	echo " to locate default sample files and config files."
	exit 1
fi

# Check for CFG and sample file
if [ "$CFG_FILE" == "" ]; then
	CFG_FILE="MOCAT.cfg"
fi

# Set full path
CFG_FILE=`readlink -f "$CFG_FILE"`

if [ -f "$CFG_FILE" ];
then
	echo -n ""
else
	echo " ERROR: Could not find $CFG_FILE"
	echo " Unless a config file is specified, MOCAT.cfg is expected to exist in the current folder."
	exit 1
fi

rm -f TMP.runMOCAT.sample_files TMP.runMOCAT.files TMP.runMOCAT.folders
if [ "$SAMPLE_FILE" == "" ]; then
	echo " SAMPLE_FILE not specified with option -sf SAMPLE_FILE"
	echo " Looking for valid sample files in the current folder:"
	echo "   Getting files..."
	find . -maxdepth 1 -xtype f -size -500k -size +0 | sed 's/\.\///' | grep -v '.summary$' | grep -v '~$' > TMP.runMOCAT.files
	echo "   Getting folders..."
	find . -maxdepth 1 -type d | sed 's/\.\///' > TMP.runMOCAT.folders
	echo -n "   Processing files"
	for file in `cat TMP.runMOCAT.files`; do
		echo -n "."
		MATCHING=`cat $file | fgrep -x -f - TMP.runMOCAT.folders | grep -c .`
		TOTAL=`grep -c . $file`
		if [ "$TOTAL" == "$MATCHING" ]; then
			echo "$file" >> TMP.runMOCAT.sample_files
		fi
	done
	echo ""
	
	# Select SAMPLE_FILE
	echo ""
	LIST=`cat TMP.runMOCAT.sample_files 2>/dev/null | sed 's/^/   - /g'`
	if [ "$LIST"  == "" ]; then
		echo " There are no valid sample files in this folder."
		echo ""
		echo " Have you created the following:"
		echo "  1. A folder for each sample containing .fq or .fq.gz files for each lane"
		echo "     The pair 1 and pair 2 files should be on the format: LANE.1.fq(.gz) and LANE.2.fq(.gz)"
		echo "  2. A sample file containing the name of the samples to analyze, one per line"
		echo "     This is the names of the folders, in which you have placed the fq files"
		exit 1
	fi
	echo " SELECT A SAMPLE FILE:"
	echo "$LIST"
	echo ""
	echo -n " ENTER SAMPLE FILE: "
	read SAMPLE_FILE < /dev/tty
fi

# Check sample file
SAMPLE_FILE=`readlink -f "$SAMPLE_FILE"`
if [ -f "$SAMPLE_FILE" ]
then
	echo -n ""
else
	echo ""
	echo " ERROR: Could not find $SAMPLE_FILE"
	echo " Please specify a correct sample file."
	exit 1
fi

# Load defaults
rm -f TMP.runMOCAT.sample_files TMP.runMOCAT.files TMP.runMOCAT.folders
DEFAULT_READS=''
MODE=`grep 'MOCAT_data_type' $CFG_FILE | sed s'/.*:\s*\(\S*\).*/\1/'`
SCRIPTS=`grep 'MOCAT_dir' $CFG_FILE | sed s'/.*:\s*\(\S*\).*/\1/'`
SCRIPTS="$SCRIPTS/scripts"

# Set parameters
SAMPLES="$SAMPLE_FILE"
NO_OF_SAMPLES=`wc -l $SAMPLES | sed s'/ .*//'`

# Prepare folders
mkdir -p $CWD/logs/status $CWD/logs/tmp

# Select numer of samples in each group
if [ "$GROUPING" = "" ]; then
	GROUPING=`echo " $NO_OF_SAMPLES / 5 + 1" | bc`
fi

echo ""
echo " SAMPLE FILE: $SAMPLE_FILE"
echo " CONFIG FILE: $CFG_FILE"


### BUILD EXECUTABLE ###

# Set variables
IDENTIFIER_DATE=`date "+%Y%m%d_%H%M%S"`
DIRECTORY="$CWD"
SAMPLE_FILE="$SAMPLES"
ALL_EXECUTE=""
CFG="$CFG_FILE"
STATUS="$CWD/logs/processing/"
TMPDIR="$CWD/logs/tmp/"

# Get input from user, and create tmp/TMP.job files
if [ "$EXECUTE" = "" ]; then
	AVAILABLE=`find $SCRIPTS/ -xtype f | sort | grep -v '~'`

	echo ""
	echo " AVAILABLE SCRIPTS:"
	COUNTER=1
	for f in `echo $AVAILABLE`
	do
		g=`echo $f | sed 's/.*scripts\///'`
		h=`echo $f | cut -f 7 -d'/'`
		TEXT=`head -1 $f | sed 's/^# //'`
		echo -e " $COUNTER: $g\n    $TEXT"
		COUNTER=`expr $COUNTER + 1`
	done
	echo -ne "\n STEP TO EXECUTE (enter number): "
	read EXECUTE < /dev/tty
	if echo "$EXECUTE" | egrep -q '^[0-9]+$'; then
		COUNTER=1
		for f in `echo $AVAILABLE`
		do 
			if [ "$COUNTER" = "$EXECUTE" ]; then
				g=`echo $f | sed 's/.*scripts\///'`
				EXECUTE="$g"
			fi
		COUNTER=`expr $COUNTER + 1`
		done
	fi
fi

# We only allow one line from now on
if [ -f "$SCRIPTS/$EXECUTE" ];
then
	echo -n ""
else
	echo "Script file $SCRIPTS/$EXECUTE does not exists."
	exit 1
fi
EXE=`basename $EXECUTE`
ALL_EXECUTE="$EXE"

# Override server, and if so, edit final file
SERVER=`hostname`

# Process the script in EXECUTE to build the commands
NAME=`head -1 $SCRIPTS/$EXECUTE | sed 's/^# //'`
echo "" >> $TMPDIR/$IDENTIFIER_DATE.TMP.jobs
echo "echo \`date\` \":: STARTING $NAME\" >> \$DIRECTORY/PROGRESS &&" >> $TMPDIR/$IDENTIFIER_DATE.TMP.jobs
echo "" >> $TMPDIR/$IDENTIFIER_DATE.TMP.jobs

IFS=$'\n'
for LINE in `grep '^ *\[command].*' $SCRIPTS/$EXECUTE`
do
	CURRENTSERVER=`echo $LINE | sed 's/.*\[\(.*\)\].*/\1/'`
	CMD1=`echo $LINE | sed 's/.*\[command\][ ]*-\([^ ]*\) *\([^ ^-]*\).*/\1 \2/' | sed 's/ /_/g' | sed 's/_$//'`
	COMMAND=`echo $LINE | sed 's/.*\][ ]*//'`
	echo "echo \`date\` \":: starting $CMD1\" >> \$DIRECTORY/PROGRESS &&" >> $TMPDIR/$IDENTIFIER_DATE.TMP.jobs
	echo "cd \$DIRECTORY && MOCAT.pl -sf $SAMPLE_FILE -cfg \$CFG/MOCAT.cfg $COMMAND &&" >> $TMPDIR/$IDENTIFIER_DATE.TMP.jobs
	echo "echo \`date\` \":: completed $CMD1\" >> \$DIRECTORY/PROGRESS &&" >> $TMPDIR/$IDENTIFIER_DATE.TMP.jobs
	echo "rm \$STATUS/\$IDENTIFIER/status/BASE.* &&" >> $TMPDIR/$IDENTIFIER_DATE.TMP.jobs
	echo "touch \$STATUS/\$IDENTIFIER/status/BASE.$CMD1$CMD2.completed &&" >> $TMPDIR/$IDENTIFIER_DATE.TMP.jobs
	echo ""  >> $TMPDIR/$IDENTIFIER_DATE.TMP.jobs
done

echo "" >> $TMPDIR/$IDENTIFIER_DATE.TMP.jobs
echo "rm \$STATUS/\$IDENTIFIER/status/BASE.* &&" >> $TMPDIR/$IDENTIFIER_DATE.TMP.jobs
echo "touch \$STATUS/\$IDENTIFIER/status/BASE.finished &&" >> $TMPDIR/$IDENTIFIER_DATE.TMP.jobs
echo "echo \`date\` \":: FINISHED $NAME\" >> \$DIRECTORY/PROGRESS &&" >> $TMPDIR/$IDENTIFIER_DATE.TMP.jobs
echo "" >> $TMPDIR/$IDENTIFIER_DATE.TMP.jobs
unset IFS

# Set the options
IFS=$'\n'
for LINE in `grep '^\{option}' $SCRIPTS/$EXECUTE`
do
	COUNTER=1
	CHANGE=`echo $LINE | cut -f 1 -d'=' | sed -e 's/{option}//' -e s'/ //g'`
	OPTIONS=`echo $LINE | cut -f 2 -d'=' | sed -e 's/\s*:\s*/\n/g' -e 's/^\s*//' -e 's/\s*$//'`
	echo ""
	echo " Please specify $CHANGE to use:"
	for OPTION in $OPTIONS
	do
echo "  $COUNTER: $OPTION"
COUNTER=`expr $COUNTER + 1`
	done
	echo -en "\n Select $CHANGE (enter number): "
	read SELECTED < /dev/tty
	if echo "$SELECTED" | egrep -q '^[0-9]+$'; then
COUNTER=1
for OPTION in $OPTIONS
do 
	if [ "$COUNTER" = "$SELECTED" ]; then
SELECTED="$OPTION"
	fi
	COUNTER=`expr $COUNTER + 1`
done
	fi
	sed -i "s/${CHANGE}/${SELECTED}/" $TMPDIR/$IDENTIFIER_DATE.TMP.jobs
done
unset IFS


# Set new identifier variable
SAMPLES_BASE=`basename $SAMPLES`
IDENTIFIER="$USER.$SAMPLES_BASE.$IDENTIFIER_DATE.$ALL_EXECUTE"
mkdir -p $TMPDIR/$IDENTIFIER

# Generate inital part of final file
echo "TMPDIR='$TMPDIR'"       >> $TMPDIR/$IDENTIFIER/TMP.init
echo "IDENTIFIER='$IDENTIFIER' &&"       >> $TMPDIR/$IDENTIFIER/TMP.init
echo "DIRECTORY='$DIRECTORY' &&"         >> $TMPDIR/$IDENTIFIER/TMP.init
echo "CFG='$TMPDIR/$IDENTIFIER' &&" >> $TMPDIR/$IDENTIFIER/TMP.init
echo "STATUS='$STATUS' &&"         >> $TMPDIR/$IDENTIFIER/TMP.init
echo "ORIGINAL_SAMPLE_FILE=$SAMPLE_FILE &&"       >> $TMPDIR/$IDENTIFIER/TMP.init
echo "USER=$USER &&"            >> $TMPDIR/$IDENTIFIER/TMP.init
echo ''                >> $TMPDIR/$IDENTIFIER/TMP.init

# Mail Part
echo "INFORMATION:"                   >> $TMPDIR/$IDENTIFIER/TMP.mail
echo "$PROJECT | $SUBMITDATE | $SAMPLES_BASE ($NO_OF_SAMPLES samples) | $ALL_EXECUTE | $SERVER" >> $TMPDIR/$IDENTIFIER/TMP.mail
echo "Complete log files: $DIRECTORY/logs/processing/$IDENTIFIER/log"       >> $TMPDIR/$IDENTIFIER/TMP.mail
echo "Quick overview:  $STATUS/$IDENTIFIER/progress/$SAMPLES_BASE.progress"           >> $TMPDIR/$IDENTIFIER/TMP.mail

# Generate last part of final file
echo ''                    >  $TMPDIR/$IDENTIFIER/TMP.exit
echo 'echo `date` ":: FINISHED PIPELINE " >> $DIRECTORY/PROGRESS'    >> $TMPDIR/$IDENTIFIER/TMP.exit
echo 'echo `date` ":: PIPELINE HAS NOW STOPPED " >> $DIRECTORY/PROGRESS'   >> $TMPDIR/$IDENTIFIER/TMP.exit
echo 'if [ -f $STATUS/$IDENTIFIER/status/BASE.finished ];'         >> $TMPDIR/$IDENTIFIER/TMP.exit
echo 'then'                   >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' echo OK > $STATUS/$IDENTIFIER/status/BASE.finished'         >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' echo "Part COUNTER / TOTALCOUNTER has finished successfully" >> $TMPDIR/$IDENTIFIER/TMP.finished' >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' FINISHED=`grep -c . $TMPDIR/$IDENTIFIER/TMP.finished`' >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' if [ "$FINISHED" == "TOTALCOUNTER" ]; '        >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' then'                  >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' echo -e "\nTHIS EMAIL INDICATES THAT ALL JOBS FINISHED SUCCESSFULLY FOR THIS SUBMISSION" >> $TMPDIR/$IDENTIFIER/TMP.mail' >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' TIME=`date`'                 >> $TMPDIR/$IDENTIFIER/TMP.exit
echo "  mail -s \"SUCCESS : $PROJECT [$SAMPLES] ($ALL_EXECUTE) on $SERVER @ \$TIME\" \"$EMAIL\" < $TMPDIR/$IDENTIFIER/TMP.mail" >> $TMPDIR/$IDENTIFIER/TMP.exit
echo '  echo `date` ":: Sent success email " >> $DIRECTORY/PROGRESS'    >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' fi'                 >> $TMPDIR/$IDENTIFIER/TMP.exit
echo 'else'                   >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' rm $STATUS/$IDENTIFIER/status/BASE.*'           >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' touch $STATUS/$IDENTIFIER/status/BASE.FAILED'         >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' cat $DIRECTORY/SAMPLES > $STATUS/$IDENTIFIER/status/BASE.FAILED'     >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' echo -e "\nJOB OUTPUT:\nFILE: PROGRESS" >> $TMPDIR/$IDENTIFIER/TMP.mail.COUNTER' >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' cat $DIRECTORY/logs/processing/$IDENTIFIER/progress/BASE* >> $TMPDIR/$IDENTIFIER/TMP.mail.COUNTER ' >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' echo -e "\nMOCAT OUTPUT:" >> $TMPDIR/$IDENTIFIER/TMP.mail.COUNTER' >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' cat $DIRECTORY/logs/processing/$IDENTIFIER/log/BASE* >> $TMPDIR/$IDENTIFIER/TMP.mail.COUNTER ' >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' cat $TMPDIR/$IDENTIFIER/TMP.mail $TMPDIR/$IDENTIFIER/TMP.mail.COUNTER > $TMPDIR/$IDENTIFIER/TMP.mail.COUNTER.final' >> $TMPDIR/$IDENTIFIER/TMP.exit
echo ' TIME=`date`'                 >> $TMPDIR/$IDENTIFIER/TMP.exit
echo " cat $TMPDIR/$IDENTIFIER/TMP.mail.COUNTER.final | mail -s \"FAILED : $PROJECT [$SAMPLES_BASE] ($ALL_EXECUTE) on $SERVER @ \$TIME\" \"$EMAIL\"" >> $TMPDIR/$IDENTIFIER/TMP.exit
echo 'fi'                  >> $TMPDIR/$IDENTIFIER/TMP.exit

# Paste files
cat $TMPDIR/$IDENTIFIER/TMP.init  $TMPDIR/$IDENTIFIER_DATE.TMP.jobs $TMPDIR/$IDENTIFIER/TMP.exit > $TMPDIR/$IDENTIFIER/TMP.exe

# Edit executable: set reads for the different projects
echo ""
echo " Preparing scripts..."

# Do the actual editing of the file
cp $CFG_FILE $TMPDIR/$IDENTIFIER/MOCAT.cfg
sed -i "s/READS/${DEFAULT_READS}/g" $TMPDIR/$IDENTIFIER/TMP.exe

### BUILD AND SUBMIT JOBS ###

# Prepare
mkdir -p $DIRECTORY/logs/processing/$IDENTIFIER/jobs $DIRECTORY/logs/processing/$IDENTIFIER/log $DIRECTORY/logs/processing/$IDENTIFIER/progress $STATUS/$IDENTIFIER/status/

echo " Submitting..."
cp $SAMPLES $DIRECTORY/logs/processing/$IDENTIFIER/jobs
base="$SAMPLES_BASE"
COUNTER=1
TOTALCOUNTER=1
sed "s|SAMPLES|logs/processing/${IDENTIFIER}/jobs/${base}|" $TMPDIR/$IDENTIFIER/TMP.exe > $DIRECTORY/logs/processing/$IDENTIFIER/jobs/$base.process
sed -i "s|PROGRESS|logs/processing/${IDENTIFIER}/progress/${base}.progress|" $DIRECTORY/logs/processing/$IDENTIFIER/jobs/$base.process
sed -i "s|BASE|${base}|" $DIRECTORY/logs/processing/$IDENTIFIER/jobs/$base.process
sed -i "s|TOTALCOUNTER|$TOTALCOUNTER|g" $DIRECTORY/logs/processing/$IDENTIFIER/jobs/$base.process
sed -i "s|COUNTER|$COUNTER|g" $DIRECTORY/logs/processing/$IDENTIFIER/jobs/$base.process
sed -i "s|ADDITIONAL_PARAMETERS|${ADDITIONAL_PARAMETERS}|" $DIRECTORY/logs/processing/$IDENTIFIER/jobs/$base.process
sh $DIRECTORY/logs/processing/$IDENTIFIER/jobs/$base.process >$DIRECTORY/logs/processing/$IDENTIFIER/log/$base.log 2>$DIRECTORY/logs/processing/$IDENTIFIER/log/$base.log.error &
touch $STATUS/$IDENTIFIER/status/$base.queued_or_started

# Print final message
echo ""
echo " ##############################################################################"
echo " #  YOUR JOBS HAVE BEEN SUBMITTED TO THE CURRENT MACHINE YOU ARE LOGGED INTO  #"
echo " ##############################################################################"
echo "   "
echo "$PROJECT | $SUBMITDATE | $SAMPLES_BASE ($NO_OF_SAMPLES samples) | $ALL_EXECUTE | $SERVER"
echo "Complete log files: $DIRECTORY/logs/processing/$IDENTIFIER/log/"
echo "Quick overview:	 $STATUS/$IDENTIFIER/progress/$SAMPLES_BASE.progress"
echo "   "

