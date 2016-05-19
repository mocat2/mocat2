INPUTS=(t1.fq)
STDOUT=stdout.txt
DEPENDENCIES=(../../public/bin/fastq_trim_filter_v5_EMBL expected.single.fq.gz expected.pair.1.fq.gz expected.pair.2.fq.gz)
function execute () {
    ./fastq_trim_filter_v5_EMBL -m fastx -a t1.fq -o output -Q 33 -f 3
}

function check() {
     zdiff expected.single.fq.gz output.single.fq.gz || return 1
     zdiff expected.pair.1.fq.gz output.pair.1.fq.gz || return 1
     zdiff expected.pair.2.fq.gz output.pair.2.fq.gz || return 1
}

