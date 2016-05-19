soap2sam=../../public/src/MOCATFilter_soap2sam.py
DEPENDENCIES=${soap2sam}
INPUTS=(input2k.soap)
STDOUT=output.sam
function execute () {
    python MOCATFilter_soap2sam.py -v MIN_LEN=45 -v MIN_AS=95 < input2k.soap
}

