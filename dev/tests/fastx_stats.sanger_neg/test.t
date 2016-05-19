INPUTS=(negSanger.fq)
OUTPUTS=(output.txt)
DEPENDENCIES=(../../public/bin/fastx_quality_stats)
function execute () {
    ./fastx_quality_stats -o output.txt < negSanger.fq
}
