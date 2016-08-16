# MOCAT annotate gene catalog module

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

# Load additional configs from cfg file, there will be some errors for some lines, but we dont care
grep '^export' $MC_MODCFG | source /dev/stdin 2>/dev/null &&

# Define
export LC_ALL=C
FOLDER="$MC_WD/GENE_CATALOG_ANNOTATIONS/$MC_SF" &&
LOCAL_TMP="$MC_WD/GENE_CATALOG_ANNOTATIONS/$MC_SF/temp" &&
PROJECT_FOLDER="$MC_WD" &&
DATA="$MC_MOD/data" &&
BLASTTYPE=`echo $MC_OPT | sed 's/ /\n/g' | grep '^BLASTTYPE=' | cut -f 2 -d'='` &&


# Prepare
echo 'Preparing...' &&
rm -fr $FOLDER &&
mkdir -p $FOLDER &&
mkdir -p $LOCAL_TMP &&
cd $FOLDER &&

if [ "$MC_DB" == "" ]; then
if [ -e "$MC_WD/GENE_CATALOGS/$MC_SF/catalog/$MC_SF.faa" ]; then
echo "EXISTS FILE=$MC_WD/GENE_CATALOGS/$MC_SF/catalog/$MC_SF.faa"
FILE="$MC_WD/GENE_CATALOGS/$MC_SF/catalog/$MC_SF.faa"
else
echo "MISSING FILE=$MC_WD/GENE_CATALOGS/$MC_SF/catalog/$MC_SF.faa"
exit 1
fi
else
if [ -e "$MC_DB" ]; then
echo "EXISTS FILE=$MC_DB"
FILE="$MC_DB"
else
echo "MISSING FILE=$MC_DB"
exit 1
fi
fi

echo "Running DIAMOND..." &&
NAME=`echo $FILE | sed 's%.*/%%'` &&
NAME2=`echo $NAME | sed 's/.faa$//'` &&
FILEFOLDER=`dirname $FILE` &&
ln -s $FILE $NAME &&
CWD=`pwd` &&
{  time $MC_BIN/diamond $BLASTTYPE -v --sensitive -p $MC_CPU -t $LOCAL_TMP  -d $DATA/eggNOG -k 10 -q $CWD/$NAME -a $CWD/$NAME2-eggNOG >$CWD/$NAME2-eggNOG.log 2>$CWD/$NAME2-eggNOG.log && $MC_BIN/diamond view -a $CWD/$NAME2-eggNOG -o $CWD/$NAME2-eggNOG.blastp && $MC_SRC/filterBlastReport.pl --flavor NCBI --topbits=0 --bits=60 $CWD/$NAME2-eggNOG.blastp  > $CWD/$NAME2-eggNOG.filtered ; } 2>$CWD/$NAME2-eggNOG.time

echo "OG mapping..." &&
BRH=$MC_SRC/find_best_hit.pl &&
OG=$MC_SRC/og_mapping.py &&
LEN=$MC_MOD/data/eggNOG.len &&
OGFILE=$MC_MOD/data/eggNOG.OG &&
$BRH -i $CWD/$NAME2-eggNOG.filtered -f NCBI > $CWD/$NAME2-eggNOG.besthit &&
$OG -b $CWD/$NAME2-eggNOG.filtered -t $CWD/$NAME2-eggNOG.besthit -g $OGFILE -p $LEN -f NCBI 2> $CWD/$NAME2-eggNOG.OG.err > $CWD/$NAME2-eggNOG.OG.pre &&
egrep -v ' NONE ' $CWD/$NAME2-eggNOG.OG.pre > $CWD/$NAME2-eggNOG.OG &&
rm $CWD/$NAME2-eggNOG.OG.pre &&
FILE=$CWD/$NAME2-eggNOG.OG &&
FILE2=$CWD/$NAME2-eggNOG.annotation &&
echo -e "#gene\teggNOG_OG" > $FILE2
cut -f 1,3 $FILE | grep -v '^#' | perl -ane '
chomp(@F);
foreach $i (1 .. scalar @F-1){
$h{$F[0]}{$i}{$F[$i]}=1
};
END{
foreach $g (keys %h) {
print "$g";
foreach $i (1 .. scalar @F-1) {
print "\t";
foreach $j (keys %{$h{$g}{$i}}){
print "$j|";
}
}
print "\n"
}
}' | sed 's/|$//' | sed 's/\t\t/\tNA\t/' | sed 's/\t\t/\tNA\t/' | sed 's/\t\t/\tNA\t/' | sed 's/|\t/\t/' | sort -k 1,1 >> $FILE2 &&




for DB in `ls $MC_MOD/data/*.HMM $MC_MOD/data/*.DIAMOND 2>/dev/null | sed 's/.*eggNOG.\(.*\.[DIAMOND|HMM]\)$/\1/'`
do
DB2=`echo $DB | cut -f 1 -d'.'` &&
HMM=`echo $DB | grep -c HMM` &&
sleep 2
if [ "$HMM" == "0" ]; then
echo "MAKING $NAME2 $DB2" &&
grep -v '#' $CWD/$NAME2-eggNOG.filtered | cut -f 1,2,9,10,12 | sort -k2,2 > $NAME2.$DB2.tmp1 &&
join -1 2 -2 1 -o 1.1 1.2 1.3 1.4 1.5 2.1 2.2 2.3 2.4 2.5 -e "xNAx" $NAME2.$DB2.tmp1 $MC_MOD/data/eggNOG-$DB | sed 's/ /\t/g' > $NAME2.$DB2.tmp3 &&
echo -en ""
fi # end HMM=0
if [ "$HMM" == "1" ]; then
echo "MAKING $NAME2 $DB2 HMM" &&
grep -v '#' $CWD/$NAME2-eggNOG.filtered | cut -f 1,2,9,10,12 | sort -k2,2 > $NAME2.$DB2.tmp1 &&
join -1 2 -2 2 -o 1.1 1.2 1.3 1.4 1.5 2.1 -e "xNAx" $NAME2.$DB2.tmp1 $MC_MOD/data/eggNOG-$DB | sed 's/ /\t/g' > $NAME2.$DB2.tmp3 &&
echo -en ""
fi # end HMM=1
if [ "$DB" == "ARDB.DIAMOND" ]; then scores="60:0.9"; fi
if [ "$DB" == "CARD.DIAMOND" ]; then scores="60:0.9"; fi
if [ "$DB" == "DBETH.DIAMOND" ]; then scores="100:0.9"; fi
if [ "$DB" == "KEGG.DIAMOND" ]; then scores="60:0.9"; fi
if [ "$DB" == "PATRIC.DIAMOND" ]; then scores="100:0.9"; fi
if [ "$DB" == "Pfam.HMM" ]; then scores="60:0.9"; fi
if [ "$DB" == "VFDB.DIAMOND" ]; then scores="100:0.9"; fi
if [ "$DB" == "Victors.DIAMOND" ]; then scores="100:0.9"; fi
if [ "$DB" == "MetaCyc.DIAMOND" ]; then scores="100:0.9"; fi
if [ "$DB" == "dbCAN.HMM" ]; then scores="60:0.9"; fi
if [ "$DB" == "DrugBank.DIAMOND" ]; then scores="100:0.9"; fi
if [ "$DB" == "MvirDB.DIAMOND" ]; then scores="100:0.9"; fi
if [ "$DB" == "Prophages.DIAMOND" ]; then scores="100:0.9"; fi
if [ "$DB" == "Resfams.HMM" ]; then scores="60:0.9"; fi
if [ "$DB" == "SEED.DIAMOND" ]; then scores="60:0.9"; fi
if [ "$DB" == "Superfamily.HMM" ]; then scores="60:0.9"; fi
if [ "$DB" == "vFam.DIAMOND" ]; then scores="100:0.9"; fi
export score=`echo $scores | cut -f 1 -d':'` &&
export overlap=`echo $scores | cut -f 2 -d':'` &&
if [ "$HMM" == "1" ]; then
perl -wlane '
BEGIN{use List::Util qw(min max)}
print if ($F[4] >= $ENV{score} );
' $NAME2.$DB2.tmp3 | cut -f 1,6 | perl -F"\t" -ane '
chomp(@F);
foreach $i (1 .. scalar @F-1){
@l = split /%/, $F[$i];
for $j (0 .. scalar @l-1){
$k=$j+1;
@l2 = split /\|/, $l[$j];
foreach $l2 (@l2) {
$h{$F[0]}{$k}{$l2}=1;
}
}
};
END{
foreach $g (sort keys %h) {
print "$g";
foreach $i (1 .. scalar keys %{$h{$g}}) {
print "\t";
$counter=0;
foreach $j (sort keys %{$h{$g}{$i}}){
$counter++;
if ( ($counter > 1 && $j ne "NA") || $counter==1 ) {
print "$j|";
}
}
}
print "\n"
}
}' | sed 's/|$//' | sed 's/\t\t/\tNA\t/' | sed 's/\t\t/\tNA\t/' | sed 's/\t\t/\tNA\t/' | sed 's/|\t/\t/' | sed 's/|\t/\t/g' | sed 's/|\t/\t/g' |  sort -k1,1 > $NAME2.$DB2.$scores.tmp5 &&
echo -en ""
fi # end HMM=1
if [ "$HMM" == "0" ]; then
perl -wlane '
BEGIN{use List::Util qw(min max)}
$overlap=max(0, min($F[3], $F[8]) - max($F[2], $F[7])) / min( $F[3]-$F[2],$F[8]-$F[7] );
print if (($F[4] >= $ENV{score} && $F[9] >= $ENV{score} && $overlap >= $ENV{overlap}));
' $NAME2.$DB2.tmp3 | sort -k7,7 > $NAME2.$DB2.$scores.tmp4 &&
join -1 7 -2 1 -o 1.1 2.2 -e "xNAx" $NAME2.$DB2.$scores.tmp4 $MC_MOD/data/$DB2.MAP | sed 's/ /\t/g' | sort -u | perl -F"\t" -ane '
chomp(@F);
foreach $i (1 .. scalar @F-1){
@l = split /%/, $F[$i];
for $j (0 .. scalar @l-1){
$k=$j+1;
@l2 = split /\|/, $l[$j];
foreach $l2 (@l2) {
$h{$F[0]}{$k}{$l2}=1;
}
}
};
END{
foreach $g (sort keys %h) {
print "$g";
foreach $i (1 .. scalar keys %{$h{$g}}) {
print "\t";
$counter=0;
foreach $j (sort keys %{$h{$g}{$i}}){
$counter++;
if ( ($counter > 1 && $j ne "NA") || $counter==1 ) {
print "$j|";
}
}
}
print "\n"
}
}' | sed 's/|$//' | sed 's/\t\t/\tNA\t/' | sed 's/\t\t/\tNA\t/' | sed 's/\t\t/\tNA\t/' | sed 's/|\t/\t/' | sed 's/|\t/\t/g' > $NAME2.$DB2.$scores.tmp5 &&
echo -en ""
fi # end HMM=0
echo -e "\t$DB2\t" | fgrep -f - $MC_MOD/data/MAP | sort -k3,3 | cut -f 1 | sed ':a;N;$!ba;s/\n/ /g' | sed 's/^/#gene\t/' | sed 's/ /\t/g' > $NAME2.$DB2.header &&
cat $NAME2.$DB2.header $NAME2.$DB2.$scores.tmp5 | sort > $NAME2-$DB2.annotation &&
echo -en ""
done # end DB




echo "Removing..." &&
#rm $CWD/*tmp[0-9] $CWD/*header &&

echo "Joining..." &&
join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 2.2 2.3 $NAME2-eggNOG.annotation $NAME2-ARDB.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 2.2 2.3 - $NAME2-CARD.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 2.2 2.3 - $NAME2-DBETH.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 2.2 2.3 2.4 2.5 2.6 2.7 - $NAME2-DrugBank.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 2.2 - $NAME2-ICEberg.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 2.2 2.3 2.4 - $NAME2-KEGG.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 2.2 - $NAME2-MetaCyc.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 2.2 - $NAME2-MvirDB.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.20 2.2 2.3 2.4 - $NAME2-PATRIC.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.20 1.21 1.22 1.23 2.2 2.3 - $NAME2-Pfam.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.20 1.21 1.22 1.23 1.24 1.25 2.2 - $NAME2-Prophages.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.20 1.21 1.22 1.23 1.24 1.25 1.26 2.2 - $NAME2-Resfams.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.20 1.21 1.22 1.23 1.24 1.25 1.26 1.27 2.2 2.3 2.4 2.5 - $NAME2-SEED.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.20 1.21 1.22 1.23 1.24 1.25 1.26 1.27 1.28 1.29 1.30 1.31 2.2 2.3 - $NAME2-Superfamily.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.20 1.21 1.22 1.23 1.24 1.25 1.26 1.27 1.28 1.29 1.30 1.31 1.32 1.33 2.2 2.3 2.4 - $NAME2-VFDB.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.20 1.21 1.22 1.23 1.24 1.25 1.26 1.27 1.28 1.29 1.30 1.31 1.32 1.33 1.34 1.35 1.36 2.2 2.3 2.4 - $NAME2-Victors.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.20 1.21 1.22 1.23 1.24 1.25 1.26 1.27 1.28 1.29 1.30 1.31 1.32 1.33 1.34 1.35 1.36 1.37 1.38 1.39 2.2 2.3 - $NAME2-dbCAN.annotation | sed 's/ /\t/g' | join -1 1 -2 1 -a 1 -a 2 -e "NA" -o 0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.20 1.21 1.22 1.23 1.24 1.25 1.26 1.27 1.28 1.29 1.30 1.31 1.32 1.33 1.34 1.35 1.36 1.37 1.38 1.39 1.40 1.41 2.2 - $NAME2-vFam.annotation | sed 's/ /\t/g' | perl -F"\t" -lane 'foreach (@F) {s/^NA$//}; print join "\t", @F' > $NAME2.functional.map &&

echo "Cleaning up..." &&
mkdir -p temp &&
mv $NAME *tmp[1-5] *.header *.annotation *OG.err *OG *blastp *daa *log *time *filtered *besthit temp &&
echo "You can remove all data in the temp folder." &&
echo "Copying..." &&
cp $NAME2.functional.map $MC_DATA &&
ln -fs $MC_DATA/$NAME2.functional.map $MC_DATA/$NAME2.padded.functional.map &&
NUM=`ls $MC_WD/GENE_CATALOGS/$MC_SF/padded_catalog/*coord | grep -c .`
if [ "$NUM" == 1 ]; then
echo -en ""
else
ln -fs $MC_DATA/$NAME2.functional.map $MC_DATA/$NAME2.padded.1-$NUM.functional.map &&
echo -en ""
fi

echo "DONE"