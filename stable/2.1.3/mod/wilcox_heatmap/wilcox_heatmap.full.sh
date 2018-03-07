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
# MC_DATE = the DATE identifier for the MOCAT job

##########################################################################################################################
# Define variables, these will later be set inside MOCAT
##########################################################################################################################

# Specify if extracted or screened reads are used
if [ "$MC_E" == "yes" ]; then
    MC_R="extracted.$MC_R"
else
    MC_R="screened.$MC_R"
fi

# check if these already defined, then save them to over write whatever was in cfg
export restore=0
export SELECT=false
if [ "$filter_var" != "" ]; then
    export restore=1
    export filter_var_safe="$filter_var"
    export allowed_set_safe="$allowed_set"
fi
# Load additional configs from cfg file, there will be some errors for some lines, but we dont care
source $MC_MODCFG 2>/dev/null

if [ "$restore" == 1 ]; then
    export filter_var="$filter_var_safe"
    export SELECT="true"
    export allowed_set="{$allowed_set_safe}"
fi

# We need to process some of the input variables to relate to what we select
if [ "$MC_TAX" == "mOTU" ]; then
    categories="mOTU:mOTU.cluster"
fi
if [ "$MC_TAX" == "NCBI" ]; then
    categories=NCBI:phylum|NCBI:class|NCBI:order|NCBI:family|NCBI:genus|NCBI:species|NCBI:specI_cluster
fi
if [ "$MC_TAX" == "NA" ]; then
if [ "$MC_M" == "functional" ]; then
    ZIP="$MC_WD/PROFILES/$MC_M.profiles/$MC_SF/$MC_SF.$MC_M.profiles.$MC_R.on.$MC_DB.$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID.zip"
    if [ -e "$ZIP" ]; then
	categories=`unzip -l $ZIP | sed 's/.*\.//' | grep -P -v '^[0-9]*\s+' | grep -v '^----' | grep -v '^zip$' | sort -u | sed 's/^/functional:/'`
	if [ "$?" == "0" ]; then
	    echo "OK : extracting categories"
	else 
	    echo -e "FAILED : unzip -l $ZIP | sed 's/.*\.//' | grep -P -v '^[0-9]*\s+' | grep -v '^----' | grep -v '^zip$' | sort -u | sed 's/^/functional:/'"
	    exit 1
	fi
    else
	echo "ERROR & EXIT: Missing $ZIP"
	exit 1
    fi
else
echo "ERROR & EXIT: Expected MC_M=functional but MC_M=$MC_M"
exit 1
fi
fi

# set some more variables
metadata=`echo $MC_OPT | sed 's/ /\n/g' | grep '^METADATA=' | cut -f 2 -d'='` &&
grouping=`echo $MC_OPT | sed 's/ /\n/g' | grep '^GROUPING=' | cut -f 2 -d'='` &&
groups=`echo $MC_OPT | sed 's/ /\n/g' | grep '^GROUPS=' | cut -f 2 -d'='` &&
DATE="$MC_DATE"
#DATE=`date +\%Y\%b\%d_\%H\%M\%S` &&
if [ -s "$metadata" ]; then
echo "OK : $metadata"
else
echo "ERROR & EXIT: Missing $metadata"
fi
column=` head -1 $metadata | sed 's/\t/\n/g' | grep -n . | grep $grouping | cut -f 1 -d':'` &&
if [ "$column" == "" ]; then
    echo "ERROR & EXIT: It seems like the metadata file does not contain the column name '$grouping' on the first line"
    exit 1
else
    echo "OK : grouping"
fi
gr1=`echo $groups | cut -f 1 -d","` &&
gr2=`echo $groups | cut -f 2 -d","` &&
DATE1=$DATE &&
main=`echo $categories | head -1 | cut -f 1 -d':' | tr '[A-Z]' '[a-z]'` && # main should be the same the whole time
if [ "$MC_NAME" == "" ]; then
    echo "OK : name not set"
else
    MC_NAME=".$MC_NAME"
    echo "OK : name set"
fi
export MAINDATADIR="$MC_OUT/$MC_SF.$main.$grouping$MC_NAME.$DATE" &&



##########################################################################################################################
# Extract data files
##########################################################################################################################
echo "=== PART 1 ===" &&
mkdir -p $MAINDATADIR/data $MAINDATADIR/settings $MAINDATADIR/tables $MAINDATADIR/figures &&
if [ "$?" == "0" ]; then
    echo "OK : creating folders"
else 
    echo -e "FAILED : mkdir -p $MAINDATADIR/data $MAINDATADIR/settings $MAINDATADIR/tables $MAINDATADIR/figures"
    exit 1
fi
echo -e "index\tfile\tinsert_or_base\tnorm_type\ttax_level\tsettings\tall\tzipped\tload" > $MAINDATADIR/settings/R.files &&
index=1 &&
for category in $categories
do
main=`echo $category | cut -f 1 -d':' | tr '[A-Z]' '[a-z]'` &&
main2=`echo $category | cut -f 1 -d':'` &&
sub=`echo $category | cut -f 2 -d':'` &&
load=yes &&
if [ "$main" == "gene" ]; then
    main2="" &&
    MC_M2=$MC_M &&
    MC_M='gene' &&
    load=no
fi
if [ "$MC_M" == "taxonomic" ]; then
    MC_M3=$main
fi
if [ "$MC_M" == "functional" ]; then
    main2="" &&
    MC_M3=$MC_M
fi
if [ "$main" == "motu" ]; then
ZIP="$MC_WD/PROFILES/$MC_M.profiles/$MC_SF/$main2/$MC_SF.$main.profile.$MC_R.on.$MC_DB.$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID.insert.mm.dist.among.unique.scaled.mOTU.clusters" &&
sub=mOTU.clusters &&
cp $ZIP $MAINDATADIR/data/$main.$sub
if [ "$?" == "0" ]; then
    echo "OK : cp $ZIP $MAINDATADIR/data/$main.$sub"
else 
    echo -e "FAILED : cp $ZIP $MAINDATADIR/data/$main.$sub"
    exit 1
fi
else
ZIP="$MC_WD/PROFILES/$MC_M.profiles/$MC_SF/$main2/$MC_SF.$main.profiles.$MC_R.on.$MC_DB.$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID.zip"
zip1="$MC_SF.$MC_M3.profile.$MC_R.on.$MC_DB.$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID.$type.$norm.$sub"
# SIAMCAT currently do not use the horizontal values, but could in theory be used
#zip2="$MC_SF.$MC_M3.profile.$MC_R.on.$MC_DB.$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID.horizontal.$sub"
if [ -e $ZIP ]; then
    OK=1
else
    echo -e "FAILED : $ZIP does not exist"
    exit 1
fi
$MC_BIN/unzip -p $ZIP $zip1 > $MAINDATADIR/data/$main.$sub
if [ "$?" == "0" ]; then
    echo "OK : $MC_BIN/unzip -p $ZIP $zip1 > $MAINDATADIR/data/$main.$sub"
else 
    echo -e "FAILED : $MC_BIN/unzip -p $ZIP $zip1 > $MAINDATADIR/data/$main.$sub"
    exit 1
fi
fi

# This is an important bit. Here we treat functional in a specific way, because to get fractions we divide by the "mapped"
# As of now, the other methods are not supported. But it should be just taking fractions when including -1 and then removing -1, 
# possibly all this can be done within SIAMCAT.
if [ "$MC_M" == "functional" ]; then
export normalize="none"
grep -v '^#' $MAINDATADIR/data/$main.$sub | awk -F"\t" '{if(NR==1){out=""; for(i=3;i<=NF;i++){out=out"\t"$i}; print $2out}else{print}}' | grep -Pv '^NA\t' | grep -Pv '^unassigned\t' | perl -F"\t" -lane 'chomp(@F); if(m/^mapped\t/){@D=@F} else{  if($D[0]){@G=();$G[0]=$F[0]; foreach $f (1..scalar @F-1){$G[$f]=$F[$f]/$D[$f]}; print join "\t", @G }else{print $_}}' > $MAINDATADIR/data/feat_$main.$sub.tsv
if [ "$?" == "0" ]; then
    echo "OK : awk & perl filter"
else 
    echo -e "FAILED : awk & perl filter"
    exit 1
fi
else
if [ "$MC_M" == "taxonomic" ]; then
if [ "$main" == "ncbi" ]; then
export normalize="decostand:total"
grep -v '^#' $MAINDATADIR/data/$main.$sub | awk -F"\t" '{if(NR==1){out=""; for(i=3;i<=NF;i++){out=out"\t"$i}; print $2out}else{print}}' > $MAINDATADIR/data/feat_$main.$sub.tsv
if [ "$?" == "0" ]; then
    echo "OK : awk"
else 
    echo -e "FAILED : awk"
    exit 1
fi
else
if [ "$main" == "motu" ]; then
export normalize="decostand:total"
grep -v '^#' $MAINDATADIR/data/$main.$sub | awk -F"\t" '{if(NR==1){out=""; for(i=3;i<=NF;i++){out=out"\t"$i}; print $2out}else{print}}' > $MAINDATADIR/data/feat_$main.$sub.tsv
if [ "$?" == "0" ]; then
    echo "OK : awk"
else 
    echo -e "FAILED : awk"
    exit 1
fi
else
echo "mode=$MC_M main=$main currently not supported."
exit 1
fi
fi
else
echo "mode=$MC_M main=$main currently not supported."
exit 1
fi
fi

# This could be added later if needed
if [ "$main" == "gene" ]; then
    MC_M=$MC_M2
    echo "ERROR & EXIT: SIAMCAT was not designed to run using genes"
    exit 1
fi

# Here we would extract the horizontal coverages
# $MC_BIN/unzip -p $ZIP $zip2 > $MAINDATADIR/data/$main.$sub.horizontal
# if [ "$?" == "0" ]; then
#     echo "OK : $MC_BIN/unzip -p $ZIP $zip2 > $MAINDATADIR/data/$main.$sub.horizontal"
# else 
#     echo -e "FAILED : $MC_BIN/unzip -p $ZIP $zip2 > $MAINDATADIR/data/$main.$sub.horizontal"
#     exit 1
# fi

if [ "$main" == "gene" ]; then
    MC_M=$MC_M2
fi
echo -e "$index\t$MAINDATADIR/data/feat_$main.$sub.tsv\t$type\t$norm\t$sub\t$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID\t$type.$norm.$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID\tno\t$load" >> $MAINDATADIR/settings/R.files &&
index=`expr $index + 1` &&
echo "$main.$sub" >> $MAINDATADIR/settings/tags.files
done
echo "Part 2 complete"




##########################################################################################################################
# Setup the first R session and run it
##########################################################################################################################
echo "=== PART 2 ==="
echo -e "type\tsetting" > $MAINDATADIR/settings/R.settings
echo -e "working_dir\t$MAINDATADIR/settings" >> $MAINDATADIR/settings/R.settings
echo -e "tables_dir\t$MAINDATADIR/tables" >> $MAINDATADIR/settings/R.settings
echo -e "data_dir\t$MAINDATADIR/data" >> $MAINDATADIR/settings/R.settings
echo -e "figures_dir\t$MAINDATADIR/figures" >> $MAINDATADIR/settings/R.settings
echo -e "downsample\tnone" >> $MAINDATADIR/settings/R.settings
echo -e "load_functions\t$MC_MOD/load.functions.R" >> $MAINDATADIR/settings/R.settings
echo -e "filter\t0" >> $MAINDATADIR/settings/R.settings
echo -e "filter_fraction\t0" >> $MAINDATADIR/settings/R.settings
echo -e "filter_low_abundant_features\t$filter_low_abundant_features" >> $MAINDATADIR/settings/R.settings
echo -e "filter_prevalence_of_features\t0" >> $MAINDATADIR/settings/R.settings
echo -e "normalize\t$normalize" >> $MAINDATADIR/settings/R.settings
echo -e "load_groups\tyes" >> $MAINDATADIR/settings/R.settings
echo -e "load_metadata\t$metadata" >> $MAINDATADIR/settings/R.settings
echo -e "load_feature_annotations\tno" >> $MAINDATADIR/settings/R.settings
echo -e "grouping\t$grouping" >> $MAINDATADIR/settings/R.settings
echo -e "groups\t$groups" >> $MAINDATADIR/settings/R.settings
echo -e "function\t$MC_MOD/wilcox.R" >> $MAINDATADIR/settings/R.settings
echo -e "select_samples\t$SELECT" >> $MAINDATADIR/settings/R.settings
echo -e "filter_var\t$filter_var" >> $MAINDATADIR/settings/R.settings
echo -e "allowed_set\t$allowed_set" >> $MAINDATADIR/settings/R.settings
$pathtoR/Rscript $MC_MOD/run.R $MAINDATADIR/settings/R.settings
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
significant=`awk 'NR>1' $MAINDATADIR/tables/*.p-values.comparison*.table | awk -F"," '{if($4<=0.05){print}}' | sed 's/,/\t/g' | grep -c .`
if [ "$significant" == 0 ]; then
    touch $MAINDATADIR/data/categories.signif.genes
    echo "$xGROUPING"
    echo "================================================================="
    echo "Pipeline part 1-2 finished successfully. However:"
    echo "All categories have in total $significant significant entities"
    echo "Will exit now because pipeline may fail downstream."
    echo "Data can be found here:"
    echo "$MAINDATADIR/tables/*.p-values.comparison*.table"
    echo "The following was executed to test:"
    echo "awk 'NR>1' $MAINDATADIR/tables/*.p-values.comparison*.table | awk -F"," '{if(\$4<=0.05){print}}' | sed 's/,/\t/g' | grep -c ."
    echo "================================================================="
    exit 1
else
    echo "All categories have in total $significant significant entities"
    awk 'NR>1' $MAINDATADIR/tables/*.p-values.comparison*.table | awk -F"," '{if($4<=0.05){print}}' | sed 's/,/\t/g' > $MAINDATADIR/data/significant.categories
    awk 'NR>1' $MAINDATADIR/tables/*.p-values.comparison*.table | awk -F"," '{if($4<=0.05){print}}' | sed 's/,/\t/g' | cut -f 1 | fgrep -wf - $MC_FUNCT_MAP > $MAINDATADIR/data/categories.signif.genes
    head -1 $MC_FUNCT_MAP > $MAINDATADIR/data/categories.signif.genes.header
    cat $MAINDATADIR/data/categories.signif.genes.header $MAINDATADIR/data/categories.signif.genes > $MAINDATADIR/data/categories.signif.genes2
    mv $MAINDATADIR/data/categories.signif.genes2 $MAINDATADIR/data/categories.signif.genes
fi
for FILE in $MAINDATADIR/data/gene.gene $MAINDATADIR/data/gene.gene.horizontal
do
echo "Extracting from $FILE"
cut -f 1 $MAINDATADIR/data/categories.signif.genes | sed 's/$/\t/' | fgrep -f - $FILE > $FILE.significant
head -6 $FILE > $FILE.header
cat $FILE.header  $FILE.significant >  $FILE.significant2
mv  $FILE.significant2  $FILE.significant
done
echo "Extractions completed"
echo "================================================================="
echo "Pipeline part 1-2 finished successfully."
echo "All categories have in total $significant significant entities"
echo "Data can be found here:"
echo "$MAINDATADIR/tables/*.p-values.comparison*.table"
echo "The following was executed to test:"
echo "awk 'NR>1' $MAINDATADIR/tables/*.p-values.comparison*.table | awk -F"," '{if(\$4<=0.05){print}}' | sed 's/,/\t/g' | grep -c ."
echo "================================================================="
exit 0






















# Step two (deep_analysis script) of this analysis may require the fractions of the genes
$MC_BIN/MOCATFraction.pl -in $MAINDATADIR/data/gene.gene -out $MAINDATADIR/data/gene.gene.fraction

# We have to load the horizontal coverages, note that now we are starting to mix DATE1 and DATE2
DATE=`date +\%Y\%b\%d_\%H\%M\%S`
DATE2=$DATE
mkdir -p $MAINDATADIR/data $MAINDATADIR/settings $MAINDATADIR/tables $MAINDATADIR/figures
echo -e "index\tfile\tinsert_or_base\tnorm_type\ttax_level\tsettings\tall\tzipped\tload" > $MAINDATADIR/settings/R.files
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
    echo -e "$index\t$MAINDATADIR1/data/$main.$sub.horizontal$end\t$type\t$norm\t$sub.horizontal\t$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID\t$type.$norm.$MC_MODE.$sub.horizontal.$MC_MAPPING.l$MC_LEN.p$MC_ID\tno\t$load" >> $MAINDATADIR/settings/R.files
    index=`expr $index + 1`
    echo -e "$index\t$MAINDATADIR1/data/$main.$sub$end\t$type\t$norm\t$sub\t$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID\t$type.$norm.$MC_MODE.$sub.$MC_MAPPING.l$MC_LEN.p$MC_ID\tno\t$load" >> $MAINDATADIR/settings/R.files
else
    echo -e "$index\t$MAINDATADIR1/data/$main.$sub.horizontal$end\t$type\t$norm\t$sub.horizontal\t$MC_MODE.$MC_MAPPING.l$MC_LEN.p$MC_ID\t$type.$norm.$MC_MODE.$sub.horizontal.$MC_MAPPING.l$MC_LEN.p$MC_ID\tno\t$load" >> $MAINDATADIR/settings/R.files
    index=`expr $index + 1`
fi
done

# Generate settings
echo "=== PART 4 ==="
echo -e "type\tsetting" > $MAINDATADIR/settings/R.settings
echo -e "working_dir\t$MAINDATADIR/settings" >> $MAINDATADIR/settings/R.settings
echo -e "tables_dir\t$MAINDATADIR/tables" >> $MAINDATADIR/settings/R.settings
echo -e "data_dir\t$MAINDATADIR/data" >> $MAINDATADIR/settings/R.settings
echo -e "figures_dir\t$MAINDATADIR/figures" >> $MAINDATADIR/settings/R.settings
echo -e "downsample\tnone" >> $MAINDATADIR/settings/R.settings
echo -e "load_functions\t$MC_MOD/load.functions.R" >> $MAINDATADIR/settings/R.settings
echo -e "filter\t0" >> $MAINDATADIR/settings/R.settings
echo -e "filter_fraction\t0" >> $MAINDATADIR/settings/R.settings
echo -e "filter_low_abundant_features\t0" >> $MAINDATADIR/settings/R.settings # this has changed from DATE1; this was set to 0 while above it was set to 0.00001. Why? Because we are looking at plotting the genes.
echo -e "filter_prevalence_of_features\t0" >> $MAINDATADIR/settings/R.settings
echo -e "normalize\tnone" >> $MAINDATADIR/settings/R.settings # this has changed from DATE1
echo -e "load_groups\tyes" >> $MAINDATADIR/settings/R.settings
echo -e "load_metadata\t$metadata" >> $MAINDATADIR/settings/R.settings
echo -e "load_feature_annotations\t$annotation_file" >> $MAINDATADIR/settings/R.settings # this has changed from DATE1
echo -e "grouping\t$grouping" >> $MAINDATADIR/settings/R.settings
echo -e "function\t$MC_MOD/plot.R" >> $MAINDATADIR/settings/R.settings # this has changed from DATE1
echo -e "load_file\t$MAINDATADIR1/data/session.RData" >> $MAINDATADIR/settings/R.settings # this has changed from DATE1
echo -e "significant_gene\t$MAINDATADIR1/data/categories.signif.genes" >> $MAINDATADIR/settings/R.settings # this has changed from DATE1

# include the list of the significant categories
for FILE in $MAINDATADIR1/tables/*.p-values.comparison*.table
do
part=`echo $FILE | sed 's/.p-values.comparison.*//' | sed 's|.*/||'`
echo -e "significant_$part\t$FILE" >> $MAINDATADIR/settings/R.settings # this has changed from DATE1
done


# Run the R script
$pathtoR/Rscript $MC_MOD/run.R $MAINDATADIR/settings/R.settings
if [ "$?" == "0" ]; then
    echo "OK : Plotting"
else 
    echo "FAILED : Plotting"
    exit 1
fi

echo "$xGROUPING"
echo "================================================================="
echo "Pipeline finished successfully. Now you can launch R and run:"
echo ">load('$MAINDATADIR/data/session.RData')"
echo "Please use the $MC_MOD/deep_analysis.R script as a starting point"
echo ""
echo "Figures for the significant categories are here:   $MAINDATADIR/figures"
echo "Significant categories in tabular format are here: $MAINDATADIR/tables"
echo "================================================================="








