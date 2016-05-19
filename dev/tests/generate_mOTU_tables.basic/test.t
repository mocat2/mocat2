mOTU_tables=../../public/src/MOCATPasteTaxonomyCoverageFiles_generate_mOTU_tables.pl
DEPENDENCIES=${mOTU_tables}
INPUTS=(mOTU.nr.padded.motu.linkage.map mOTU_tables.table.gz)
OUTPUTS=(motus_test.annotated.mOTU.clusters.fraction.tab motus_test.annotated.mOTU.clusters.tab motus_test.mOTU.clusters.fraction.tab motus_test.mOTU.clusters.tab)
function execute () {
    perl $(basename ${mOTU_tables}) --wd ${PWD} --map ${INPUTS[1]} --table ${INPUTS[2]} --prefix motus_test
}
