




for xGROUPING in "GROUPING=Phenotype GROUPS=Lean_to_normal|Overweight_to_obese" "GROUPING=Sex GROUPS=M|F" "GROUPING=Diet GROUPS=Hi_Protein_low_CHO|Low_Protein_Hi_CHO" "GROUPING=Breed GROUPS=BEA|LAB" "GROUPING=Cohort GROUPS=1|2"
do



# This script will run wilcoxon tests between the two groups to see which ones are
# significant and then generate heatmaps of coverage and horizontal coverages

# This is a MOCAT module extension script
# You can use this script as a template
# to cerate your own modules. The current
# 'rules' to creating a new module are>
# 1. create a new folder in the MOCAT/mod folder
# 2. create a script inside the folder with the
#    <same name>.sh
# There are a number of environmental PATHs that
# are set automatically within MOCAT that can be
# accessed inside thes emodule script, a short
# list of these are:
# MC_SF = the sample file
# MC_TMP = the temp folder
# MC_BIN = the bin folder
# MC_MOD = this module folder, that means this
#          folder can act as a bin folder, where
#          you can store executables or results
# MC_OUT = the final output folder whwre final
#          results should be saved
# MC_R   = the -r option
# MC_DB  = the database DB eg (MOCAT.pl -p <DB>)
#          in the shortened format
# MC_ID  = the identity cutoff
# MC_CPU = number of CPUs
# MC_OPT = Additional options passed to the
#          MOCAT -options <OPT> command
#          these should preferably be on the form
#          -options OPT1=a OPT2=b ...
# MC_M  = -mode option
# MC_WD = working directory
# MC_E  = yes or no whether the -e option is activated
# MC_MODE = solexaqa or fastx
# MC_MAPPING = allbest, unique or random
# MC_LEN = specfied length
# MC_FUNCT_MAP = if it exists, the MOCAT funcitonal map file
#                it is important that the funcitonal map is not
#                in gz format

##########################################################################################################################
# Define variables, these will later be set inside MOCAT
##########################################################################################################################
MC_FUNCT_MAP=/g/bork1/kultima/MOCAT/data/dogs65v1.RGC.padded.functional.map
MC_SF=dogs.65
MC_OUT=/g/bork5/mocat/DOGS/ANALYSES/wilcox
MC_TMP=/g/bork5/mocat/DOGS/tmp
MC_R=screened.adapter.on.cm3
MC_DB=dogs65v1.RGC.padded
MC_ID=95
MC_M=functional
MC_OPT="TYPE=base NORM=mm.dist.among.unique.norm CATEGORIES=UNKNOWN:ARDB_drug|UNKNOWN:ARDB_superclass|UNKNOWN:CARD_aro|UNKNOWN:CARD_drug|UNKNOWN:E4abe|UNKNOWN:ICEberg|UNKNOWN:K62_ko|UNKNOWN:K62_module|UNKNOWN:K62_pathway|UNKNOWN:K74ppve_ko|UNKNOWN:K74ppve_module|UNKNOWN:K74ppve_pathway|UNKNOWN:Pfam28_clan|UNKNOWN:Pfam28_pf|UNKNOWN:biocyc|UNKNOWN:dbcan|UNKNOWN:ep74_module|UNKNOWN:ep74_reaction|UNKNOWN:mvir|UNKNOWN:prophages|UNKNOWN:resfam_drug|UNKNOWN:resfam_id|UNKNOWN:seed_categoryA|UNKNOWN:seed_categoryB|UNKNOWN:seed_enzyme|UNKNOWN:seed_pathway|UNKNOWN:vFam METADATA=/g/bork/kultima/metadata/dogs/dogs ANNOTATION_FILE='' $xGROUPING"

MC_WD=/g/bork5/mocat/DOGS/
MC_MODE=solexaqa
MC_MAPPING=allbest
MC_LEN=45
MC_ID=95
MC_BIN=/g/bork8/kultima/GIT/MOCAT/master/public/bin/
MC_MOD=/g/bork8/kultima/GIT/MOCAT/master/public/mod/wilcox_heatmap/
filter_low_abundant_features=0.000001






if [ "$MC_E" == "yes" ]; then
    MC_R="extracted.$MC_R"
else
    MC_R="screened.$MC_R"
fi
categories=`echo $MC_OPT | sed 's/ /\n/g' | grep '^CATEGORIES=' | cut -f 2 -d'=' | sed 's/|/\n/g'` &&
type=`echo $MC_OPT | sed 's/ /\n/g' | grep '^TYPE=' | cut -f 2 -d'='` &&
norm=`echo $MC_OPT | sed 's/ /\n/g' | grep '^NORM=' | cut -f 2 -d'='` &&
metadata=`echo $MC_OPT | sed 's/ /\n/g' | grep '^METADATA=' | cut -f 2 -d'='` &&
grouping=`echo $MC_OPT | sed 's/ /\n/g' | grep '^GROUPING=' | cut -f 2 -d'='` &&
groups=`echo $MC_OPT | sed 's/ /\n/g' | grep '^GROUPS=' | cut -f 2 -d'='` &&
annotation_file=`echo $MC_OPT | sed 's/ /\n/g' | grep '^ANNOTATION_FILE=' | cut -f 2 -d'='` &&
DATE=`date +\%Y\%b\%d_\%H\%M\%S` &&
DATE1=$DATE &&

# HARD CODED OPTIONS FOR THIS SCRIPT
# '{if($4<=0.05){print}} :: gets the significant categories based on FDR because that's the fourth column of that file


##########################################################################################################################
# Extract categories including genes, including generating files for the first R session
##########################################################################################################################
echo "=== PART 1 ===" &&
mkdir -p $MC_OUT/$DATE/data $MC_OUT/$DATE/settings $MC_OUT/$DATE/tables $MC_OUT/$DATE/figures &&
echo -e "index\tfile\tinsert_or_base\tnorm_type\ttax_level\tsettings\tall\tzipped\tload" > $MC_OUT/$DATE/settings/R.files &&
index=1 &&
for category in $categories    #"gene:gene"
do
main=`echo $category | cut -f 1 -d':'` &&
main2=`echo $category | cut -f 1 -d':'` &&
sub=`echo $category | cut -f 2 -d':'` &&
load=yes &&
if [ "$main" == "gene" ]; then
    main2="" &&
    MC_M2=$MC_M &&
    MC_M='gene' &&
    load=no
fi
ZIP="$MC_WD/PROFILES/$MC_M.profiles/$MC_SF/$main2/$MC_SF.$main.profiles.$MC_R.on.$MC_DB.$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID.zip"
zip1="$MC_SF.$MC_M.profile.$MC_R.on.$MC_DB.$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID.$type.$norm.$sub"
zip2="$MC_SF.$MC_M.profile.$MC_R.on.$MC_DB.$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID.horizontal.$sub"
if [ -e $ZIP ]; then
    OK=1
else
    echo -e "FAILED : $ZIP does not exist"
    exit 1
fi
$MC_BIN/unzip -p $ZIP $zip1 > $MC_OUT/$DATE/data/$main.$sub
if [ "$?" == "0" ]; then
    echo "OK : $MC_BIN/unzip -p $ZIP $zip1 > $MC_OUT/$DATE/data/$main.$sub"
else 
    echo -e "FAILED : $MC_BIN/unzip -p $ZIP $zip1 > $MC_OUT/$DATE/data/$main.$sub"
    exit 1
fi
if [ "$main" == "gene" ]; then
    MC_M=$MC_M2
fi
$MC_BIN/unzip -p $ZIP $zip2 > $MC_OUT/$DATE/data/$main.$sub.horizontal
if [ "$?" == "0" ]; then
    echo "OK : $MC_BIN/unzip -p $ZIP $zip2 > $MC_OUT/$DATE/data/$main.$sub.horizontal"
else 
    echo -e "FAILED : $MC_BIN/unzip -p $ZIP $zip2 > $MC_OUT/$DATE/data/$main.$sub.horizontal"
    exit 1
fi
if [ "$main" == "gene" ]; then
    MC_M=$MC_M2
fi
echo -e "$index\t$MC_OUT/$DATE/data/$main.$sub\t$type\t$norm\t$sub\t$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID\t$type.$norm.$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID\tno\t$load" >> $MC_OUT/$DATE/settings/R.files &&
index=`expr $index + 1` 
done

##########################################################################################################################
# Setup the first R session and run it
##########################################################################################################################
echo "=== PART 2 ==="
echo -e "type\tsetting" > $MC_OUT/$DATE/settings/R.settings
echo -e "working_dir\t$MC_OUT/$DATE/settings" >> $MC_OUT/$DATE/settings/R.settings
echo -e "tables_dir\t$MC_OUT/$DATE/tables" >> $MC_OUT/$DATE/settings/R.settings
echo -e "data_dir\t$MC_OUT/$DATE/data" >> $MC_OUT/$DATE/settings/R.settings
echo -e "figures_dir\t$MC_OUT/$DATE/figures" >> $MC_OUT/$DATE/settings/R.settings
echo -e "downsample\tnone" >> $MC_OUT/$DATE/settings/R.settings
echo -e "load_functions\t$MC_MOD/load.functions.R" >> $MC_OUT/$DATE/settings/R.settings
echo -e "filter\t0" >> $MC_OUT/$DATE/settings/R.settings
echo -e "filter_fraction\t0" >> $MC_OUT/$DATE/settings/R.settings
echo -e "filter_low_abundant_features\t$filter_low_abundant_features" >> $MC_OUT/$DATE/settings/R.settings
echo -e "filter_prevalence_of_features\t0" >> $MC_OUT/$DATE/settings/R.settings
echo -e "normalize\tdecostand:total" >> $MC_OUT/$DATE/settings/R.settings
echo -e "load_groups\tyes" >> $MC_OUT/$DATE/settings/R.settings
echo -e "load_metadata\t$metadata" >> $MC_OUT/$DATE/settings/R.settings
echo -e "load_feature_annotations\tno" >> $MC_OUT/$DATE/settings/R.settings
echo -e "grouping\t$grouping" >> $MC_OUT/$DATE/settings/R.settings
echo -e "groups\t$groups" >> $MC_OUT/$DATE/settings/R.settings
echo -e "function\t$MC_MOD/wilcox.R" >> $MC_OUT/$DATE/settings/R.settings
Rscript $MC_MOD/run.R $MC_OUT/$DATE/settings/R.settings
if [ "$?" == "0" ]; then
    echo "OK : Finished wilcoxon tests"
else 
    echo "FAILED : Running wilcoxon tests"
    exit 1
fi

##########################################################################################################################
# Process the wilcox results and run the second R session, here we load the horizontal coverages
##########################################################################################################################
# Now that we have generated these results, we want to check the genes that are in common for all these levels
# and also extract these genes from the actual abundance files, and horizontal files
echo "=== PART 3 ==="
significant=`awk 'NR>1' $MC_OUT/$DATE/tables/*.p-values.comparison*.table | awk -F"," '{if($4<=0.05){print}}' | sed 's/,/\t/g' | grep -c .`
if [ "$significant" == 0 ]; then
    touch $MC_OUT/$DATE/data/categories.signif.genes
    echo "$xGROUPING"
    echo "================================================================="
    echo "Pipeline part 1-2 finished successfully. However:"
    echo "All categories have in total $significant significant entities"
    echo "Will exit now because pipeline may fail downstream."
    echo "Data can be found here:"
    echo "$MC_OUT/$DATE/tables/*.p-values.comparison*.table"
    echo "The following was executed to test:"
    echo "awk 'NR>1' $MC_OUT/$DATE/tables/*.p-values.comparison*.table | awk -F"," '{if(\$4<=0.05){print}}' | sed 's/,/\t/g' | grep -c ."
    echo "================================================================="
    #exit 1
else
    echo "All categories have in total $significant significant entities"
    awk 'NR>1' $MC_OUT/$DATE/tables/*.p-values.comparison*.table | awk -F"," '{if($4<=0.05){print}}' | sed 's/,/\t/g' > $MC_OUT/$DATE/data/significant.categories
    awk 'NR>1' $MC_OUT/$DATE/tables/*.p-values.comparison*.table | awk -F"," '{if($4<=0.05){print}}' | sed 's/,/\t/g' | cut -f 1 | fgrep -wf - $MC_FUNCT_MAP > $MC_OUT/$DATE/data/categories.signif.genes
    head -1 $MC_FUNCT_MAP > $MC_OUT/$DATE/data/categories.signif.genes.header
    cat $MC_OUT/$DATE/data/categories.signif.genes.header $MC_OUT/$DATE/data/categories.signif.genes > $MC_OUT/$DATE/data/categories.signif.genes2
    mv $MC_OUT/$DATE/data/categories.signif.genes2 $MC_OUT/$DATE/data/categories.signif.genes
fi
for FILE in $MC_OUT/$DATE/data/gene.gene $MC_OUT/$DATE/data/gene.gene.horizontal
do
echo "Extracting from $FILE"
cut -f 1 $MC_OUT/$DATE/data/categories.signif.genes | sed 's/$/\t/' | fgrep -f - $FILE > $FILE.significant
head -6 $FILE > $FILE.header
cat $FILE.header  $FILE.significant >  $FILE.significant2
mv  $FILE.significant2  $FILE.significant
done
echo "Extractions completed"










    echo "$xGROUPING"
    echo "================================================================="
    echo "Pipeline part 1-2 finished successfully."
    echo "All categories have in total $significant significant entities"
    echo "Data can be found here:"
    echo "$MC_OUT/$DATE/tables/*.p-values.comparison*.table"
    echo "The following was executed to test:"
    echo "awk 'NR>1' $MC_OUT/$DATE/tables/*.p-values.comparison*.table | awk -F"," '{if(\$4<=0.05){print}}' | sed 's/,/\t/g' | grep -c ."
    echo "================================================================="
    #exit 0




done # end xGROUPING

exit 0






















# Step two (deep_analysis script) of this analysis may require the fractions of the genes
$MC_BIN/MOCATFraction.pl -in $MC_OUT/$DATE/data/gene.gene -out $MC_OUT/$DATE/data/gene.gene.fraction

# We have to load the horizontal coverages, note that now we are starting to mix DATE1 and DATE2
DATE=`date +\%Y\%b\%d_\%H\%M\%S`
DATE2=$DATE
mkdir -p $MC_OUT/$DATE/data $MC_OUT/$DATE/settings $MC_OUT/$DATE/tables $MC_OUT/$DATE/figures
echo -e "index\tfile\tinsert_or_base\tnorm_type\ttax_level\tsettings\tall\tzipped\tload" > $MC_OUT/$DATE/settings/R.files
index=1
for category in $categories "gene:gene"
do
main=`echo $category | cut -f 1 -d':'`
main2=`echo $category | cut -f 1 -d':'`
sub=`echo $category | cut -f 2 -d':'`
load=yes
end=""
if [ "$main" == "gene" ]; then
    end=".significant"
    echo -e "$index\t$MC_OUT/$DATE1/data/$main.$sub.horizontal$end\t$type\t$norm\t$sub.horizontal\t$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID\t$type.$norm.$MC_MODE.$sub.horizontal.$MC_MAPPING.l$MC_LEN.p$MC_ID\tno\t$load" >> $MC_OUT/$DATE/settings/R.files
    index=`expr $index + 1`
    echo -e "$index\t$MC_OUT/$DATE1/data/$main.$sub$end\t$type\t$norm\t$sub\t$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID\t$type.$norm.$MC_MODE.$sub.$MC_MAPPING.l$MC_LEN.p$MC_ID\tno\t$load" >> $MC_OUT/$DATE/settings/R.files
else
    echo -e "$index\t$MC_OUT/$DATE1/data/$main.$sub.horizontal$end\t$type\t$norm\t$sub.horizontal\t$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID\t$type.$norm.$MC_MODE.$sub.horizontal.$MC_MAPPING.l$MC_LEN.p$MC_ID\tno\t$load" >> $MC_OUT/$DATE/settings/R.files
    index=`expr $index + 1`
fi
done

# Generate settings
echo "=== PART 4 ==="
echo -e "type\tsetting" > $MC_OUT/$DATE/settings/R.settings
echo -e "working_dir\t$MC_OUT/$DATE/settings" >> $MC_OUT/$DATE/settings/R.settings
echo -e "tables_dir\t$MC_OUT/$DATE/tables" >> $MC_OUT/$DATE/settings/R.settings
echo -e "data_dir\t$MC_OUT/$DATE/data" >> $MC_OUT/$DATE/settings/R.settings
echo -e "figures_dir\t$MC_OUT/$DATE/figures" >> $MC_OUT/$DATE/settings/R.settings
echo -e "downsample\tnone" >> $MC_OUT/$DATE/settings/R.settings
echo -e "load_functions\t$MC_MOD/load.functions.R" >> $MC_OUT/$DATE/settings/R.settings
echo -e "filter\t0" >> $MC_OUT/$DATE/settings/R.settings
echo -e "filter_fraction\t0" >> $MC_OUT/$DATE/settings/R.settings
echo -e "filter_low_abundant_features\t0" >> $MC_OUT/$DATE/settings/R.settings # this has changed from DATE1; this was set to 0 while above it was set to 0.00001. Why? Because we are looking at plotting the genes.
echo -e "filter_prevalence_of_features\t0" >> $MC_OUT/$DATE/settings/R.settings
echo -e "normalize\tnone" >> $MC_OUT/$DATE/settings/R.settings # this has changed from DATE1
echo -e "load_groups\tyes" >> $MC_OUT/$DATE/settings/R.settings
echo -e "load_metadata\t$metadata" >> $MC_OUT/$DATE/settings/R.settings
echo -e "load_feature_annotations\t$annotation_file" >> $MC_OUT/$DATE/settings/R.settings # this has changed from DATE1
echo -e "grouping\t$grouping" >> $MC_OUT/$DATE/settings/R.settings
echo -e "function\t$MC_MOD/plot.R" >> $MC_OUT/$DATE/settings/R.settings # this has changed from DATE1
echo -e "load_file\t$MC_OUT/$DATE1/data/session.RData" >> $MC_OUT/$DATE/settings/R.settings # this has changed from DATE1
echo -e "significant_gene\t$MC_OUT/$DATE1/data/categories.signif.genes" >> $MC_OUT/$DATE/settings/R.settings # this has changed from DATE1

# include the list of the significant categories
for FILE in $MC_OUT/$DATE1/tables/*.p-values.comparison*.table
do
part=`echo $FILE | sed 's/.p-values.comparison.*//' | sed 's|.*/||'`
echo -e "significant_$part\t$FILE" >> $MC_OUT/$DATE/settings/R.settings # this has changed from DATE1
done


# Run the R script
Rscript $MC_MOD/run.R $MC_OUT/$DATE/settings/R.settings
if [ "$?" == "0" ]; then
    echo "OK : Plotting"
else 
    echo "FAILED : Plotting"
    exit 1
fi

echo "$xGROUPING"
echo "================================================================="
echo "Pipeline finished successfully. Now you can launch R and run:"
echo ">load('$MC_OUT/$DATE/data/session.RData')"
echo "Please use the $MC_MOD/deep_analysis.R script as a starting point"
echo ""
echo "Figures for the significant categories are here:   $MC_OUT/$DATE/figures"
echo "Significant categories in tabular format are here: $MC_OUT/$DATE1/tables"
echo "================================================================="








