soap2sam=../../public/src/MOCATFilter_soap2sam.py
DEPENDENCIES=${soap2sam}
INPUTS=(../soap2sam.basic/input2k.soap)
STDOUT=output97.sam
function execute () {
    python $(basename $soap2sam) -v MIN_LEN=45 -v MIN_AS=97 < input2k.soap
}

