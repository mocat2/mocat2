########################################################
### THIS SCRIPT PRODUCES mOTU AND TAXONOMIC PROFILES ###
########################################################



########################################################
# INTRODUCTION                                         #
########################################################

# BEFORE YOU RUN THIS YOU MUST DOWNLOAD THE DATASET
# FROM http://mocat.embl.de
# Link: Download HMP mock community metagenome (from MOCAT article; 500 MB)
# And extract it into the folder even_sample

# This script will run MOCAT to first generate high
# quality reads and then map these reads against the
# two databases to produce mOTU and NCBI species profiles.

# THIS IS A GENERIC SCRIPT THAT YOU CAN USE ON YOUR
# OWN METAGENOMICS DATASET TO GENERATE PROFILES.
# OR, you can execute runMOCAT.sh, which is a bash
# implementation to run generic MOCAT commands on a project

# Copy this script into your project folder (described below)
# and change the sample SAMPLE_FILE and OUTPUT_FOLDER as desired



########################################################
# DEFINE VARIABLES                                     #
########################################################

# Set variables used later in the script
SAMPLE_FILE='sample'
OUTPUT_FOLDER='RESULTS'



########################################################
# DESCRIPTION OF A PROJECT FOLDER                      #
########################################################

# The required setup looks like this:
# In this folder you need:
# 1. MOCAT.cfg
# 2. a sample file (sample) in our example below

# You also need one folder for each sample mentioned in the
# sample file. Each of these folders should contain files
# called .fq.gz or .fq. If the samples are processed using
# paired end sequencing method, the lanes should be called:
# XXX.1.fq.gz and XXX.2.fq.gz. It is OK to have many different
# lane files in each sample folder. 

# An example (with paried end) could look like this:
# ____________________________________
# This folder: /usr/me/project/
# ____________________________________
# Required files:
# /usr/me/project/MOCAT.cfg
# /usr/me/project/sample
# ____________________________________
# Content of sample file:
# SAMPLE_1
# SAMPLE_2
# ____________________________________
# Sample folders:
# /usr/me/project/SAMPLE_1/lane1.1.fq
# /usr/me/project/SAMPLE_1/lane1.2.fq
# /usr/me/project/SAMPLE_1/lane2.1.fq
# /usr/me/project/SAMPLE_1/lane2.2.fq
# /usr/me/project/SAMPLE_1/lane3.1.fq
# /usr/me/project/SAMPLE_1/lane4.2.fq
#
# /usr/me/project/SAMPLE_2/laneID.1.fq
# /usr/me/project/SAMPLE_2/laneID.2.fq

# Or (single end) like this:
# ____________________________________
# This folder: /usr/me/project/
# ____________________________________
# Required files:
# /usr/me/project/MOCAT.cfg
# /usr/me/project/sample
# ____________________________________
# Content of sample file:
# even_sample
# ____________________________________
# Sample folders:
# /usr/me/project/even_sample/SRR172902.fq.gz

# NOTE:
# Change the 'MOCAT_paired_end' variable in the
# config file to 'yes' or 'no'.



########################################################
# EXECUTE MOCAT TO GENERATE PROFILES                   #
########################################################

# Initial sample processing #
MOCAT.pl -sf $SAMPLE_FILE -rtf

# Generate mOTU profiles #
MOCAT.pl -sf $SAMPLE_FILE -s mOTU.v1.padded -r reads.processed -identity 97 -extracted_files
MOCAT.pl -sf $SAMPLE_FILE -f mOTU.v1.padded -r reads.processed -identity 97
MOCAT.pl -sf $SAMPLE_FILE -p mOTU.v1.padded -r reads.processed -identity 97 -mode mOTU -o $OUTPUT_FOLDER

# Generate taxonomic profiles #
MOCAT.pl -sf $SAMPLE_FILE -s RefMG.v1.padded -r mOTU.v1.padded -e -identity 97
MOCAT.pl -sf $SAMPLE_FILE -f RefMG.v1.padded -r mOTU.v1.padded -e -identity 97
MOCAT.pl -sf $SAMPLE_FILE -p RefMG.v1.padded -r mOTU.v1.padded -e -identity 97 -mode NCBI -previous_db_calc_tax_stats_file -o $OUTPUT_FOLDER

# Done
echo "-----------------------------------------------------------------"
echo "The generated mOTU and taxonomic profiles should be available in:"
echo "$OUTPUT_FOLDER"
echo "-----------------------------------------------------------------"


# END OF SCRIPT #

