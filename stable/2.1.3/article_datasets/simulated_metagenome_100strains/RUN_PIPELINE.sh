######################################################
### THIS SHELL SCRIPT IS PART OF THE MOCAT PIPELINE ##
######################################################


################ IMPORTANT! #################
# TO RUN THIS PIPELINE YOU NEED TO DOWNLOAD #
# THE OPTIONAL COMPONENT USEARCH!           #
# PLEASE SEE THE README!                    #
#############################################


################################################################
# YOU MUST ALSO DOWNLOAD                                       #
#                                                              #
# 1. THE DATASET FROM http://mocat.embl.de                     #
#    Link: Download simulated metagenome (100 strains)         #
#    And extract it into the folder 100strains                 #
#                                                              #
# 2. THE REFERENCE STRAINS FROM http://mocat.embl.de           #
#    Link: Download 100 simulated reference strains            #
#    And extract it into the MOCAT/data folder                 #
################################################################


###########################################
# This script run the MOCAT commands on   #
# the SIMULATED METAGENOME, that were run #
# to produce the reults in the article.   #
# Please note that some additional steps, #
# such as running blastall on the genes   #
# against the reference database were run #
# outside of the pipeline. However, the   #
# objective of this script is to only     #
# reproduce the steps inside the pipeline #
# to illustrate how MOCAT can be run.     #
###########################################


### INITAL REQUIRED PROCESSING STEPS #########
# i) Quality trimming and filtering of reads #
MOCAT.pl -sf sample -rtf &&                  #
##############################################


### MAP HQ READS AGAINST REF DB AND CALCULATE COVERAGES ##################
# ii,a) Mapping of processed reads (high quality reads)                  #
#       to the 22 references genomes of the SIMULATED METAGENOME         #
#       Here, we could also specify the MOCAT command as:                #
#       MOCAT.pl -sf sample -s 100_strains_references -r reads.processed #
#       but we don't have to because this line in the config:            #
#       MOCAT_default_reads      : reads.processed                       #
MOCAT.pl -sf sample -s 100_strains_references &&                         #
                                                                         #
# ii,b) Mapped reads are passed through a filtering step                 #
#       For the SIMULATED METAGENOME this does have an effect,           #
#       because  we have paired-end reads, and the                       #
#       'filter_paired_end_filtering' is set to 'yes'                    #
#       in the config file. This paired-end filter removes a read        #
#       if its mate maps to a different reference with a higher % ID     #
MOCAT.pl -sf sample -f 100_strains_references &&                         #
                                                                         #
# ii,c) Calculate the base and read (insert) coverages                   #
MOCAT.pl -sf sample -p 100_strains_references &&                         #
##########################################################################


### ASSEMBLE READS SCREENED FOR ADAPTERS AND THE PREDICT GENES ######
# iii) Assembly of the HQ reads. Note that the assembly is done     #
#      using the quality trimmed and filtered reads. Ie, this       #
#      step is entirely decoupled form the previous steps, where    #
#      reads were mapped against the database of genomes.           #
MOCAT.pl -sf sample -a &&                                           #
                                                                    #
# iv) Had we wanted to do an assembly revision here, the command    #
#     would have been:                                              #
#     MOCAT.pl -sf sample -ar                                       #
#     BUT we have chosen not to do this to keep consistency with    #
#     the analysis of the single-end MOCK COMMUNITY                 #
                                                                    #
# v) Gene prediction on scaftigs created in the assembly step       #
#    Here we could also predict genes on contigs. This is changed   #
#    in the config file.                                            #
#    NOTE! BY DEFAULT PRODIGAL IS USED TO PREDICT GENES! THIS CAN   #
#    BE CHANGED IN THE CONFIG FILE, BUT THEN METAGENEMARK HAS TO    #
#    BE MANUALLY DOWNLOADED! IN THE ARTICLE WE USED METAGENEMARK    #
#    BUT THE TWO PRODUCE VERY SIMILAR RESULTS                       #
MOCAT.pl -sf sample -gp assembly &&                                 #
#####################################################################

### GET STATISTICS ##########
MOCAT.pl -sf sample -ss  && #
#############################

echo DONE

# END OF FILE #

