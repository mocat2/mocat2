

###############################################
### This script runs a gene catalog example ###
###############################################




########################################################
# INTRODUCTION                                         #
########################################################

# This script assumes that high quality reads have
# already been generated, that these have been assembled
# and that gene have been predicted on these assemblies.

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
SAMPLE_FILE='2samples'



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
# OPTION A:
# /usr/me/project/SAMPLE_1/lane1.1.fq
# /usr/me/project/SAMPLE_1/lane1.2.fq
# /usr/me/project/SAMPLE_1/lane2.1.fq
# /usr/me/project/SAMPLE_1/lane2.2.fq
# /usr/me/project/SAMPLE_1/lane3.1.fq
# /usr/me/project/SAMPLE_1/lane3.2.fq
# OPTION B:
# /usr/me/project/SAMPLE_1/lane1.pair.1.fq
# /usr/me/project/SAMPLE_1/lane1.pair.2.fq
# /usr/me/project/SAMPLE_1/lane1.single.fq
# /usr/me/project/SAMPLE_1/lane2.pair.1.fq
# /usr/me/project/SAMPLE_1/lane2.pair.2.fq
# /usr/me/project/SAMPLE_1/lane2.single.fq
# ____________________________________
# This folder: /usr/me/project/
# ____________________________________
# Required files:
# /usr/me/project/MOCAT.cfg
# /usr/me/project/2samples
# ____________________________________
# Content of sample file:
# MH0002
# MH0004
# ____________________________________



########################################################
# PRE-EXECUTED MOCAT COMMANDS                          #
########################################################

# FOR THE PURPOSE OF THIS EXAMPLE,
# WE ASSUME THESE HAVE ALREADY BEEN EXECUTED

# Initial sample processing
# MOCAT.pl -sf $SAMPLE_FILE -rtf

# High quality reads
# MOCAT.pl -sf $SAMPLE_FILE -s hg19 -screened_files -identity 90

# Assemble
# MOCAT.pl -sf $SAMPLE_FILE -a hg19

# Predict genes
# MOCAT.pl -sf $SAMPLE_FILE -gp assembly



########################################################
# MAKE GENE CATALOG                                    #
########################################################

# Extract the genes, cluster them and generate the catalog and
# also a padded version of the catalog. The padded version of the catalog
# will automatically be copied to the MOCAT data folder.
# Here assembly type is required and should be either 'assembly' or 'assembly.revised'
# Because, in this example, we do not revise the assembly, we use 'assembly'
# -r is set to hg19 because we screened reads against the hg19 catalog
# -cpu is set to 10 as an example. We recommend using e.g. 40 here if possible
# the default number of CPUs used is 8.

MOCAT.pl -sf $SAMPLE_FILE -make_gene_catalog -r hg19 -assembly_type assembly -cpu 10 &&

# The created gene catalog can be found in ./GENE_CATALOGS/2samples/catalog
# The created gene catalog with padded sequences (up to 100 bp at each end of the genes)
# can be found in ./GENE_CATALOGS/2samples/padded_catalog
# The padded catalog can also be found in the MOCAT data folder



########################################################
# ANNOTATE GENE CATALOG                                #
########################################################

# This will annotate an existing gene catalog
# By default it is expected that the gene catalog to annotate has been
# constructed using MOCAT in the step above, this menas as input file
# [current directory]/GENE_CATALOGS/$SAMPLE_FILE/catalog/$SAMPLE_FILE.faa
# will be used. It is also possible to specify a file after '-annotate_gene_catalog'
# as such: -annotate_gene_catalog /path/to/my/file/genecat
# But this is not advised as general usage
# Note that the input gene catalog has to be in amino acid sequence.
# When profiles are created, and reads are mapped to to a gene catalog,
# the nucleotide sequences are used.
# The field '-r hg19' does not have to be specified
# Again, we use -cpu 10 as an example

MOCAT.pl -sf $SAMPLE_FILE -annotate_gene_catalog_or_file -blasttype blastp -cpu 10 &&

# The annotations for the  gene catalog can be found in ./GENE_CATALOG_ANNOTATIONS/2samples/
# and in the MOCAT data folder



########################################################
# CREATE FUNCITONAL PROFILES BASED ON CATALOG          #
########################################################

# Now we have clustered and annotated the genes into an annotated
# gene catalog, and the .padded version of the gene catalog, as well
# as the funcitonal profiles have been imported to the data folder.
# Now we can map the high quality reads against the padded catalog
# We do this at a 95% sequence identity threshold.

MOCAT.pl -sf $SAMPLE_FILE -s $SAMPLE_FILE.padded -r hg19 -identity 95 &&
MOCAT.pl -sf $SAMPLE_FILE -f $SAMPLE_FILE.padded -r hg19 -identity 95 &&
MOCAT.pl -sf $SAMPLE_FILE -p $SAMPLE_FILE.padded -r hg19 -identity 95 -mode functional &&

echo DONE


# END OF SCRIPT #

