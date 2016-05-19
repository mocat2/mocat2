filter_stats=../../public/src/MOCATFilter_stats.pl
DEPENDENCIES=(${filter_stats})
INPUTS=(istats.txt test.sam)
OUTPUTS=(istats.txt.filtered.l40.p95.stats)
STDOUT=stdout.txt
function execute () {
    perl $(basename ${filter_stats}) --format SAM --identity 95 --length 40 --stats ${INPUTS[0]} < test.sam
}
