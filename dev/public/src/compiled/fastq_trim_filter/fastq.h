#include <inttypes.h>
#include <iostream>
#include <string>
#include <vector>
#include <cstdio>

struct fastq_read {
    std::string header;
    std::string seq;
    std::string quals;
    unsigned size() const { return quals.size(); }
    void clear() {
        seq.clear();
        quals.clear();
    }
};

const int nr_quality_values = 128;
enum fastq_type { sanger, solexa };

bool read_fastq(std::istream& input, fastq_read& read, int& lineno);
