# MOCAT make gene catalog module

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

# Load additional configs from cfg file, there will be some errors for some lines, but we dont care
grep '^export' $MC_MODCFG | source /dev/stdin 2>/dev/null &&

# Define
SAMPLE_FILE="$MC_SF" &&
FOLDER="$MC_WD/GENE_CATALOGS/$MC_SF" &&
CLUSTERS="$MC_SF.clustered.0.95.fna" &&
PROJECT_FOLDER="$MC_WD" &&
FINAL_CATALOG_NAME="$MC_SF" &&
SAMPLE_FILE_WITH_SLASH_AT_END="$MC_SF.slash" &&
ASSEMBLY_FORMAT="$MC_ASS.$MC_R.$MC_MODE.K*" &&
ASSEMBLY_FORMAT2="$ASSEMBLY_FORMAT*scaftig.gz" &&

# Prepare
echo 'Preparing...' &&
rm -fr $FOLDER && 
mkdir -p $FOLDER && 
cd $FOLDER &&
sed 's/$/\//' $PROJECT_FOLDER/$SAMPLE_FILE > $SAMPLE_FILE_WITH_SLASH_AT_END &&
echo 'Checking...' &&
for f in `cat $PROJECT_FOLDER/$SAMPLE_FILE`
do
ls $PROJECT_FOLDER/$f/gene.prediction.$ASSEMBLY_FORMAT*/*fna 2>/dev/null
done > files.genes &&
n=`grep -c . $PROJECT_FOLDER/$SAMPLE_FILE` &&
n2=`grep -c . files.genes` &&
if [ "$n" == "$n2" ]; then
echo -en ""
else 
echo "Some predicted genes are missing. Check it." &&
exit 1
fi &&

## Paste
echo 'Pasting...' &&
cat `cat files.genes` > $CLUSTERS.IN &&

## Cluster
echo 'Clustering...' &&
$MC_BIN/cd-hit-est -i $CLUSTERS.IN -c 0.95 -T $MC_CPU -M 0 -G 0 -aS 0.9 -g 1 -r 1 -d 0 -o $CLUSTERS.OUT > $CLUSTERS.log &&
CLUSTERS="$CLUSTERS.OUT" &&

# Filter shorter than 100bp
echo 'Filtering reads <100bp...' &&
$MC_SRC/fasta2tabdel.pl $CLUSTERS &&
awk -F"\t" '{if(length($2)>=100){print}}' $CLUSTERS.tabdel | sed 's/\t/\n/' > step.2.filtered.seeds &&
grep -c '>' step.2.filtered.seeds > stats.2.seeds.after.filtering &&
grep '>' step.2.filtered.seeds | sed 's/.*length:\([0-9]*\) .*/\1/' | sort -nu > stats.2.seed.lengths &&

# Check min length
echo 'Checking...' &&
MINLEN=`head -1 stats.2.seed.lengths` &&
if [ "$MINLEN" \< "100" ]; then
    echo "Min length in step.2.filtered.seeds is less than 100bp." &&
    exit 1
fi &&

# Get scaftigs
echo 'Getting scaftigs...' &&
find $PROJECT_FOLDER -name "*$ASSEMBLY_FORMAT2" > step.3.scaftig.files &&
fgrep -f $SAMPLE_FILE_WITH_SLASH_AT_END step.3.scaftig.files > step.3.scaftig.files.filtered &&
echo 'Pasting scaftigs...' &&
zcat `cat step.3.scaftig.files.filtered` > step.3.scaftigs &&
echo 'Indexing scaftigs...' &&
$MC_BIN/cdbfasta step.3.scaftigs 2> /dev/null &&

# Get genes to grep
echo 'Getting genes...' &&
grep '>' step.2.filtered.seeds | perl -lane 'm/>(\S+)_gene(\S+) strand:([+-]) start:(\d+) stop:(\d+) length:(\d+) start_codon:([a-z]+) stop_codon:([a-z]+) gene_type:([a-z]*)/; $name=$1; $geneID=$2; $gene="$1_gene$2"; $length=$6; $start=$4; $stop=$5; $strand=$3;  $newstart = $start - 100; if ( $newstart < 1 ) { $newstart = 1; }; $newstop = $stop + 100; $pos_start = $start - $newstart + 1; $pos_stop = $pos_start + $length - 1; print "$gene\t$name\t$newstart\t$newstop\t$geneID\t$pos_start\t$pos_stop"' > step.4.to.grep &&
		
# we need to use my modified version of cdbyank to get the correct IDs
echo 'Extracting genes...' &&
cut -f 2,3,4,5,6,7 step.4.to.grep | $MC_BIN/cdbyank step.3.scaftigs.cidx -R -Q -x > step.4.seeds.padded &&

echo 'Getting stats...' &&
grep -c '>' step.4.seeds.padded > stats.4.no.of.seeds &&

c1=`md5sum stats.2.seeds.after.filtering | cut -f 1 -d' '` &&
c2=`md5sum stats.4.no.of.seeds | cut -f 1 -d' '` &&
if [ "$c1" == "$c2" ]; then

# Split the padded db
echo 'Splitting DB...' &&
$MC_SRC/falen.pl -in step.4.seeds.padded -out step.4.seeds.padded.len &&
perl -lane 'BEGIN{$p=1}chomp; $sum+=$F[1]; if($sum>=3500000000){$sum=0; print $F[0]}' step.4.seeds.padded.len > step.4.seeds.padded.split &&
perl -lane 'BEGIN{open IN, "<step.4.seeds.padded.split"; while(<IN>){chomp; $h{$_}=1};$p=1; close IN; open OUT, ">step.6.padded.1"}; foreach $k (keys %h){if($F[0]=~m/^>$k/){$p++; close OUT; open OUT, ">step.6.padded.$p"}}; print OUT $_;' step.4.seeds.padded &&
grep -c '>' step.6.padded.* > stats.6.no.of.seeds &&

# Generate links to unpadded catalog
echo 'Copying...' &&
mkdir -p catalog &&
cp step.2.filtered.seeds catalog/$FINAL_CATALOG_NAME.fna &&
cd catalog &&
ln -s $FINAL_CATALOG_NAME.fna $FINAL_CATALOG_NAME &&

# Generate AA file
echo 'Generating AA file...' &&
$MC_BIN/nt2aa < $FINAL_CATALOG_NAME.fna > $FINAL_CATALOG_NAME.faa &&
cd .. &&

# Links to padded catalog
echo 'Copying...' &&
mkdir -p padded_catalog &&
for f in step.6.padded.*
do
h=`echo $f | sed 's/step.6.padded.//'` &&
cp $f padded_catalog/$FINAL_CATALOG_NAME.padded.$h &&
echo -en ""
done

echo 'Indexing...' &&
for f in padded_catalog/*padded*
do
$MC_BIN/2bwt-builder $f &&
echo -en ""
done

echo "Making lengths file..." &&
for f in padded_catalog/*padded.[0-9]
do
$MC_SRC/falen.pl -in $f -out $f.len &&
echo -en ""
done



# Make coord file
echo 'Making coord files...' &&
grep '>' step.4.seeds.padded | sed 's/>//' > step.7.padded.coord &&
for f in padded_catalog/*padded.[0-9]
do
g=`echo $f | sed 's/\.padded\.\([0-9]*\)/.padded.\1.coord/'` &&
grep '>' $f | sed 's/>//' | sed 's/$/\t/' | cut -f 1 | fgrep -wf - step.7.padded.coord > $g &&
echo -en ""
done

# Renaming if only 1 part
echo "Renaming and moving..."
p=`ls padded_catalog/*padded.[0-9] | grep -c .` &&
if [ "$p" == "1" ]; then
cd padded_catalog &&
echo -en ""
for f in *padded.1 *padded.1.index*
do
g=`echo $f | sed 's/padded\.1/padded/'` &&
mv $f $g &&
echo -en ""
done
for f in *padded.1.coord
do
g=`echo $f | sed 's/padded.1.coord/padded.coord/'` &&
mv $f $g &&
echo -en ""
done
for f in *padded.1.len
do
g=`echo $f | sed 's/padded.1.len/padded.len/'` &&
mv $f $g &&
echo -en ""
done
fi
echo "Copying..." &&
D=`pwd` &&
cp $D/* $MC_DATA/ &&

# Moving temp files
cd .. &&
mkdir temp &&
mv step.* temp &&
mv $MC_SF.clust* temp &&
mv *slash temp &&
mv files.genes temp &&
mv stats.* temp &&

echo "Done."

else
echo "ERROR & EXIT: Number of padded sequences do not match number of initial sequences in catalog" &&
exit 1
fi
