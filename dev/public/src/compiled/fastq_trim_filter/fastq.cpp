#include "fastq.h"
#include <sstream>
inline
bool skip_eol(std::istream& in) {
    return in.ignore(4096, '\n');
}
bool read_fastq(std::istream& input, fastq_read& read, int& lineno) {
    if (!std::getline(input, read.header) ||
        !std::getline(input, read.seq) ||
        !skip_eol(input) ||
        !std::getline(input, read.quals)) {
            return false;
    }
    if (read.seq.size() != read.quals.size()) {
        std::ostringstream err;
        err << "line " << (lineno + 1) << ": FASTQ Format error: quality string is not of the same size as the sequence";
        throw err.str();
    }
    lineno += 4;
    return true;
}
