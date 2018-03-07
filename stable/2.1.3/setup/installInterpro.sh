# THIS SCRIPT INSTALLS AND CONFIGURES INTERPRO ON A LOCAL SYSTEM #


# Configure
INST='/g/bork5/kultima/MOCAT/external/InterPro'
mkdir -p $INST
cd $INST


# Download and extract InterProScan
FILE='interproscan-5.2-45.0-64-bit.tar.gz'
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.2-45.0/$FILE
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.2-45.0/$FILE.md5
CHECK=`md5sum -c $FILE.md5 | cut -f 2 -d':' | sed 's/ //g'`
if [ "$CHECK" != "OK" ]; then
  echo "ERROR: $FILE md5sum was incorrect. Please restart installation."
  exit 0
fi
tar -pxvzf $FILE


# Panther download and setup
cd $INST/$FILE/data
FILE='panther-data-8.1.tar.gz'
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/$FILE
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/$FILE.md5
CHECK=`md5sum -c $FILE.md5 | cut -f 2 -d':' | sed 's/ //g'`
if [ "$CHECK" != "OK" ]; then
  echo "ERROR: $FILE md5sum was incorrect. Please restart installation."
  exit 0
fi
tar -pxvzf $FILE
