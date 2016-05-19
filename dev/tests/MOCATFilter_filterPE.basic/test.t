filterPE=../../public/src/MOCATFilter_filterPE.pl
DEPENDENCIES=${filterPE}
INPUTS=(2k.sam)
STDOUT=expected.txt
function execute () {
    perl $(basename ${filterPE}) < ${INPUTS[0]}
}
