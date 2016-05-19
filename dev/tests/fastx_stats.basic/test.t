INPUTS=(gutSimulation3.4k.1.fq.gz)
OUTPUTS=(output.txt)
DEPENDENCIES=(../../public/bin/fastx_quality_stats)
function execute () {
    gunzip -c ${INPUTS[1]} | ./fastx_quality_stats -o output.txt
}
