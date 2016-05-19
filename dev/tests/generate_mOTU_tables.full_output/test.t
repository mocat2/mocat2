mOTU_tables=../../public/src/MOCATPasteTaxonomyCoverageFiles_generate_mOTU_tables.pl
DEPENDENCIES=${mOTU_tables}
INPUTS=(../generate_mOTU_tables.basic/mOTU.nr.padded.motu.linkage.map ../generate_mOTU_tables.basic/mOTU_tables.table.gz)
function execute () {
    perl $(basename ${mOTU_tables})               --wd ${PWD} --map $(basename ${INPUTS[1]}) --table $(basename ${INPUTS[2]}) --prefix motus_test
    perl $(basename ${mOTU_tables}) --full-output --wd ${PWD} --map $(basename ${INPUTS[1]}) --table $(basename ${INPUTS[2]}) --prefix motus_test_full
}

function check () {
    (diff -u motus_test.annotated.mOTU.clusters.tab motus_test_full.annotated.mOTU.clusters.tab |grep '^+' | grep -v 0$) && return 1
    (diff -u motus_test.annotated.mOTU.clusters.fraction.tab motus_test_full.annotated.mOTU.clusters.fraction.tab |grep '^+' | grep -v 0$) && return 1
    (diff -u motus_test.mOTU.clusters.tab motus_test.mOTU.clusters.tab | grep '^+' | grep -v 0$) && return 1
    (diff -u motus_test.mOTU.clusters.fraction.tab motus_test.mOTU.clusters.fraction.tab | grep '^+' | grep -v 0$) && return 1
    return 0

}
