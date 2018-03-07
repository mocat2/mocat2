######################################################
### THIS SHELL SCRIPT IS PART OF THE MOCAT PIPELINE ##
######################################################


################ IMPORTANT! #################
# TO RUN THIS PIPELINE YOU NEED TO DOWNLOAD #
# THE OPTIONAL COMPONENT USEARCH!           #
# PLEASE SEE THE README!                    #
#############################################


#########################################################
# BEFORE YOU RUN THIS YOU MUST DOWNLOAD THE DATASET     #
# FROM http://mocat.embl.de                             #
# Link: Download HMP mock community metagenome          #
# And extract it into the folder even_sample            #
#                                                       #
# AND ALSO THE HMP mock community reference genomes     #
# FROM http://mocat.embl.de                             #
# Link: Download HMP mock community reference genomes   #
# And extract it into the MOCAT/data folder             #
#########################################################


###########################################
# This script run the MOCAT commands on   #
# the MOCK COMMUNITY, that were run to    #
# produce the reults in the article.      #
# Please note that some additional steps, #
# such as running blastall on the genes   #
# against the reference database were run #
# outside of the pipeline. However, the   #
# objective of this script is to only     #
# reproduce the steps inside the pipeline #
# to illustrate how MOCAT can be run.     #
###########################################


### CHECK THAT USEARCH HAS BEEN INSTALLED ###############################
mocatpath=`grep 'MOCAT_dir' MOCAT.cfg | cut -f 2 -d":" | sed 's/ //'`   #
if [ -e $mocatpath/ext/usearch/usearch ]; then                          #
echo "*** Usearch INSTALLED ***"                                        #
chmod u+x $mocatpath/ext/usearch/usearch                                #
#########################################################################


### INITAL REQUIRED PROCESSING STEPS #########
# i) Quality trimming and filtering of reads #
MOCAT.pl -sf sample -rtf &&                  #
##############################################


### MAP HQ READS AGAINST REF DB AND CALCULATE COVERAGES ##############
# ii,a) Mapping of processed reads (high quality reads)              #
#       to the 22 references genomes of the MOCK COMMUNITY           #
#       Here, we could also specify the MOCAT command as:            #
#       MOCAT.pl -sf sample -s mock_community_ref -r reads.processed #
#       but we don't have to because this line in the config:        #
#       MOCAT_default_reads      : reads.processed                   #
MOCAT.pl -sf sample -s mock_community_ref &&                         #
                                                                     #
# ii,b) Mapped reads are passed through a filtering step             #
#       For the MOCK COMMUNITY this has no major effect,             #
#       because a) we already filtered on the specified              #
#       cut off in the previous step, b) we have single end          #
#       reads, and the 'filter_paired_end_filtering' is set          #
#       to 'no' in the config file. However, to be able to run       #
#       the next step, this step HAS TO BE RUN. It has been          #
#       implemented this way to be able to map at a lower            #
#       cut off and then filter at different cut offs.               #
MOCAT.pl -sf sample -f mock_community_ref &&                         #
                                                                     #
# ii,c) Calculate the base and read (insert) coverages               #
MOCAT.pl -sf sample -p mock_community_ref &&                         #
######################################################################


### FIRST REMOVE READS MATCHING ADAPTER SEQUENCES, THEN MAP THOSE #####
### READS AGAINST THE 22 REFERENCE GENOMES ############################
# ii,d) Align reads against the mock_adapters, please                 #
#       note that running this step requires you to install           #
#       the optional component Usearch, which requires a              #
#       seperate download and licence key. Info can be found          #
#       in the readme file how to do this.                            #
MOCAT.pl -sf sample -sff mock_adapters &&                             #
                                                                      #
# ii,e) Same as (ii,a) but using reads that have been screened        #
#       for adapter sequences (the adapter sequences have been        #
#       removed). Had we wanted to use the adapter sequences to       #
#       map against the database, we would've specifeid this:         #
#       MOCAT.pl -sf sample -s mock_community_ref -r mock_adapters -e #
MOCAT.pl -sf sample -s mock_community_ref -r mock_adapters &&         #
                                                                      #
# ii,f) Same as (ii,b)                                                #
MOCAT.pl -sf sample -f mock_community_ref -r mock_adapters &&         #
                                                                      #
# ii,g) Same as (ii,c)                                                #
MOCAT.pl -sf sample -p mock_community_ref -r mock_adapters &&         #
#######################################################################


### ASSEMBLE READS SCREENED FOR ADAPTERS AND THE PREDICT GENES ######
# iii) Assembly of the reads that have been screened for adapters   #
MOCAT.pl -sf sample -a -r mock_adapters &&                          #
                                                                    #
# iv) Had we wanted to do an assembly revision here, the command    #
#     would have been:                                              #
#     MOCAT.pl -sf sample -ar -r mock_adapters                      #
#     BUT WE CANNOT DO ASSEMBLY REVISION because the MOCK COMMUNITY #
#     only has single end reads and not paired end reads.           #
                                                                    #
# v) Gene prediction on scaftigs created in the assembly step       #
#    Here we could also predict genes on contigs. This is changed   #
#    in the config file.                                            #
#    NOTE! BY DEFAULT PRODIGAL IS USED TO PREDICT GENES! THIS CAN   #
#    BE CHANGED IN THE CONFIG FILE, BUT THEN METAGENEMARK HAS TO    #
#    BE MANUALLY DOWNLOADED! IN THE ARTICLE WE USED METAGENEMARK    #
#    BUT THE TWO PRODUCE VERY SIMILAR RESULTS                       #
MOCAT.pl -sf sample -gp assembly -r mock_adapters &&                #
#####################################################################

### GET STATISTICS ##########
MOCAT.pl -sf sample -ss  && #
#############################

echo DONE

### IF USEARCH WAS NOT INSTALLED CORRECTLY ############################################
else                                                                                  #
echo "Usearch could not be found at this location:"                                   #
echo "$mocatpath/ext/usearch/usearch"                                                 #
echo "Do the following:"                                                              #
echo "1. Download Usearch from http://www.drive5.com/usearch/nonprofit_form.html";    #
echo "2. Rename the downloaded file to 'usearch'";                                    #
echo "3. Copy it to $mocatpath/ext/usearch";                                          #
echo "4. Run chmod u+x $mocatpath/ext/usearch/usearch to make it executable";         #
echo "5. Re-run this script"                                                          #
fi                                                                                    #
#######################################################################################

