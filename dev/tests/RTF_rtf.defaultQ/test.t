INPUTS=(gutSimulation3.4k.1.fq.gz)
STDOUT=stdout.txt
DEPENDENCIES=(../../public/bin/fastq_trim_filter_v5_EMBL expected.single.fq.gz expected.pair.1.fq.gz expected.pair.2.fq.gz)
function execute () {
    ./fastq_trim_filter_v5_EMBL -m solexaqa -a gutSimulation3.4k.1.fq.gz -o output -f 5
}

function check() {
     zdiff expected.single.fq.gz output.single.fq.gz || return 1
     zdiff expected.pair.1.fq.gz output.pair.1.fq.gz || return 1
     zdiff expected.pair.2.fq.gz output.pair.2.fq.gz || return 1
}

