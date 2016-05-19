INPUTS=(gutSimulation3.4k.1.fq.bz2 gutSimulation3.4k.2.fq.bz2)
STDOUT=stdout.txt
DEPENDENCIES=(../../public/bin/fastq_trim_filter_v5_EMBL expected.single.fq.gz expected.pair.1.fq.gz expected.pair.2.fq.gz)
function execute () {
    ./fastq_trim_filter_v5_EMBL -m solexaqa -a gutSimulation3.4k.1.fq.bz2 -b gutSimulation3.4k.2.fq.bz2 -o output -Q 64 -f 5 -2 7
}

function check() {
     zdiff expected.single.fq.gz output.single.fq.gz || return 1
     zdiff expected.pair.1.fq.gz output.pair.1.fq.gz || return 1
     zdiff expected.pair.2.fq.gz output.pair.2.fq.gz || return 1
}

