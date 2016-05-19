#include "popenstream.h"
#include "fastq.h"
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <assert.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/replace.hpp>

#if __GXX_EXPERIMENTAL_CXX0X__
    #include <unordered_map>
    typedef std::unique_ptr<std::istream> unique_istream_ptr;
    typedef std::unordered_map<std::string, fastq_read> queue_type;
#else
    #include <tr1/unordered_map>
    typedef std::auto_ptr<std::istream> unique_istream_ptr;
    typedef std::tr1::unordered_map<std::string, fastq_read> queue_type;
#endif

const char* version = "FASTQ trim/filter. EMBL v5 (2015-03-26)";


unique_istream_ptr open_stream(std::string fname) {
    if (fname.rfind(".gz") == fname.size() - 3) {
        return unique_istream_ptr(new ipopengzstream(fname.c_str()));
    }
    if (fname.rfind(".bz2") == fname.size() - 4) {
        return unique_istream_ptr(new ipopenbzstream(fname.c_str()));
    }
    return unique_istream_ptr(new std::ifstream(fname.c_str()));
}


struct trim_filter_stats {
    trim_filter_stats()
        :reads(0)
        ,sum(0)
        ,max(-1)
        ,inserts(0)
        { }

    void maybe_output_read(std::ostream& output, const fastq_read& r) {
        const int rs = r.size();
        if (!rs) return;
        if (rs > max) { max = rs; }
        sum += rs;
        ++reads;
        ++inserts;
        output
            << r.header << '\n'
            << r.seq << "\n+\n"
            << r.quals << '\n';
    }
    void maybe_output_readpair(std::ostream& out1, const fastq_read& r1, std::ostream& out2, const fastq_read& r2) {
        this->maybe_output_read(out1, r1);
        this->maybe_output_read(out2, r2);
        --inserts;
    }

    long long int reads;
    long long int sum;
    long long int max;
    long long int inserts;
};

enum trim_mode_t { solexaqa, fastx };
std::ostream& operator << (std::ostream& out, const trim_mode_t m) {
    return out <<
        ((m == solexaqa) ? "solexaqa" : "fastx");
}


struct trim_options {
    trim_options(trim_mode_t m, int f, int mQ, int mL, int min_fastx_p)
        :mode(m)
        ,start_pos(f)
        ,min_qual(mQ)
        ,min_len(mL)
        ,min_fastx_p(min_fastx_p)
        { }
    const enum trim_mode_t mode;
    const int start_pos;
    const int min_qual;
    const int min_len;
    const int min_fastx_p;
};

std::vector<std::string> header_string_split(const std::string str, const std::string sep) {
    std::vector<std::string> tokens;
    tokens.reserve(13);
    size_t start = 0;
    size_t pos = 0;
    while ((pos = str.find_first_of(sep, start)) != std::string::npos) {
        tokens.push_back(str.substr(start, pos-start));
        start = pos + 1;
    }
    if (start < str.length()) {
        tokens.push_back(str.substr(start));
    } else if (start == str.length()) {
        tokens.push_back("0");
    }
    return tokens;
}
std::string join(const std::vector<std::string>& tokens, const std::string sep) {
    std::string res;
    bool first_element = true;
    for (std::vector<std::string>::const_iterator first = tokens.begin(), past = tokens.end();
                    first != past;
                    ++first) {
        if (!first_element) res += sep;
        res += *first;
        first_element = false;
    }
    return res;
}


std::string header_stem(fastq_read& seq, int lineno=-1) {
    using boost::replace_all_copy;
    assert(seq.size());
    static bool warned_before = false;
    const size_t slash = seq.header.find('/'); // FIND /1 | /2 at end
    if (slash == seq.header.size() - 2) {
        return replace_all_copy(seq.header.substr(0, slash), " ", "_");
    }
    std::vector<std::string> tokens = header_string_split(seq.header, ": ");
    if (tokens.size() == 11) {
        if (tokens[8] == "Y") seq.clear();
        if (tokens[7] == "1" || tokens[7] == "2") {
            tokens[7] = "0";
            return join(tokens, ":");
        } else {
            std::cerr << "Cannot handle this header type: " << seq.header << "\n";
            seq.clear();
            return "";
        }
    } else if (tokens.size() == 6 || tokens.size() == 7) {
        return tokens[0]+":"+tokens[2]+":"+tokens[3] +":"+tokens[4]+":"+tokens[5]+"#0";
    } else {
        if (!warned_before && lineno >= 0) {
            std::cerr <<
                "WARNING: FastQ header format unknown on line " << lineno << ": Skipping header normalization step.\n";
            warned_before = true;
        }
        return replace_all_copy(seq.header, " ", "_");
    }
}

void fastx_trim_discard(fastq_read& read, const trim_options& opts) {
    int end = read.seq.size() - 1;
    while (end >= opts.start_pos && read.quals[end] < opts.min_qual) --end;

    const int len = end - opts.start_pos + 1;
    assert(len >= 0);
    if (end < opts.min_len) {
        read.clear();
        return;
    }
    int above = 0;
    for (int i = opts.start_pos; i <= end; ++i) {
        if (read.quals[i] >= opts.min_qual) ++above;
    }
    above *= 100;
    above /= len;
    if (above < opts.min_fastx_p) {
        read.clear();
    } else {
        read.seq = read.seq.substr(opts.start_pos, len);
        read.quals = read.quals.substr(opts.start_pos, len);
    }
}
void solexaqa_trim_discard(fastq_read& read, const trim_options& opts) {
    int s = opts.start_pos;
    int max_before = -1;
    int best_s = -1;
    const int size = read.seq.size();
    for (int i = opts.start_pos; i < (size+1); ++i) {
        const int q = (i < size ? read.quals[i] : -1);
        if (q <= opts.min_qual) {
            if ((i - s) > max_before) {
                max_before = (i - s);
                best_s = s;
            }
            s = (i+1);
        }
    }
    if (max_before < opts.min_len) {
        read.clear();
    } else {
        read.seq = read.seq.substr(best_s, max_before);
        read.quals = read.quals.substr(best_s, max_before);
    }
}

void trim_discard(fastq_read& read, const trim_options& opts) {
    if (opts.mode == solexaqa) solexaqa_trim_discard(read, opts);
    else fastx_trim_discard(read, opts);
}

trim_filter_stats rtf_se(std::istream& a,
            std::ostream& single,
            const bool header_transform,
            const trim_options& opts) {
    fastq_read ra;
    trim_filter_stats st;
    int lineno = 0;

    while (read_fastq(a, ra, lineno)) {
        trim_discard(ra, opts);
        if (ra.size()) {
            if (header_transform) ra.header = header_stem(ra, lineno);
            ra.header += "/1";
        }
        st.maybe_output_read(single, ra);
    }

    return st;
}


trim_filter_stats rtf_pe_simple(std::istream& a,
            std::istream& b,
            std::istream* c,
            std::ostream& single,
            std::ostream& pair1,
            std::ostream& pair2,
            const trim_options& opts1,
            const trim_options& opts2,
            const trim_options& opts3
            ) {
    trim_filter_stats st;
    fastq_read reads[2];
    int lines[2];
    lines[0] = lines[1] = 0;

    while (bool(a) && bool(b)) {
        const bool ok0 = read_fastq(a, reads[0], lines[0]);
        const bool ok1 = read_fastq(b, reads[1], lines[1]);
        if (ok0 != ok1) {
            std::cerr << "Unexpected end of file (paired-end files are not of the same size?)\n";
            return st;
        }
        if (!ok0) return st;
        trim_discard(reads[0], opts1);
        trim_discard(reads[1], opts1);
        reads[0].header += "/1";
        reads[1].header += "/2";

        if (reads[0].size() && reads[1].size()) {
            st.maybe_output_readpair(pair1, reads[0], pair2, reads[1]);
        } else if (reads[0].size() || reads[1].size()) {
            st.maybe_output_read(single, (reads[0].size() ? reads[0] : reads[1]));
        }
    }
    lines[0] = 0;
    while (c && read_fastq(*c, reads[0], lines[0])) {
        trim_discard(reads[0], opts3);
        st.maybe_output_read(single, reads[0]);
    }
    return st;
}
trim_filter_stats rtf_pe(std::istream& a,
            std::istream& b,
            std::istream* c,
            std::ostream& single,
            std::ostream& pair1,
            std::ostream& pair2,
            const bool header_transform,
            const trim_options& opts1,
            const trim_options& opts2,
            const trim_options& opts3
            ) {
    if (!header_transform) return rtf_pe_simple(a, b, c, single, pair1, pair2, opts1, opts2, opts3);
    fastq_read reads[2];
    queue_type rqueues[2];
    std::istream* inputs[2] = { &a, &b };
    trim_filter_stats st;
    int lineno[2];
    lineno[0] = lineno[1] = 0;

    while (bool(a) || bool(b)) {
        const int ri = ((rqueues[0].size() <= rqueues[1].size()) && bool(a)) ? 0 : 1;

        if (!read_fastq(*inputs[ri], reads[ri], lineno[ri])) {
            if (ri == 1) break;
            else continue;
        }
        std::string mirrorstr = header_stem(reads[ri]);
        reads[ri].header = mirrorstr + (ri == 0 ? "/1" : "/2");
        queue_type::iterator mirror = rqueues[1-ri].find(mirrorstr);
        if (mirror != rqueues[1-ri].end()) {
            reads[1-ri] = mirror->second;
            rqueues[1-ri].erase(mirror);
        } else {
            rqueues[ri][mirrorstr] = reads[ri];
            continue;
        }
        trim_discard(reads[0], opts1);
        trim_discard(reads[1], opts2);
        if (reads[0].size() && reads[1].size()) {
            st.maybe_output_readpair(pair1, reads[0], pair2, reads[1]);
        } else if (reads[0].size() || reads[1].size()) {
            st.maybe_output_read(single, (reads[0].size() ? reads[0] : reads[1]));
        }
    }
    for (int mi = 0; mi < 2; ++mi) {
        for (queue_type::iterator it = rqueues[mi].begin(), past = rqueues[mi].end();
                it != past;
                ++it) {
            fastq_read r = it->second;
            trim_discard(r, (mi == 0 ? opts1: opts2));
            st.maybe_output_read(single, r);
        }
    }
    while (bool(a) && read_fastq(a, reads[0], lineno[0])) {
        reads[0].header = header_stem(reads[0]) + "/1";
        trim_discard(reads[0], opts1);
        st.maybe_output_read(single, reads[0]);
    }
    if (c) {
        lineno[0] = 0;
        while (read_fastq(*c, reads[0], lineno[0])) {
                trim_discard(reads[0], opts3);
                st.maybe_output_read(single, reads[0]);
        }
    }
    return st;
}
int main(int argc, char* argv[]) {
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "produce version string")
        ("fast", "use fast gzip compression (faster computation, but larger output files)")
        ("input1,a", po::value<std::string>(), "input file 1 [required argument]")
        ("input2,b", po::value<std::string>(), "input file 2")
        ("input3,c", po::value<std::string>(), "input file 3: singleton file which will be added to single")
        ("prefix,o", po::value<std::string>(), "output file prefix [required argument]")
        ("method,m", po::value<std::string>()->default_value("fastx"), "algorithm (solexaqa or fastx)")
        ("min-length,l", po::value<int>()->default_value(30), "min length")
        ("min-quality,q", po::value<int>()->default_value(20), "min quality")
        ("min-fastx-p,p", po::value<int>()->default_value(50), "Min percentage of above quality bases for fastx filtering")
        ("quality-offset,Q", po::value<int>()->default_value(64), "quality offset")
        ("trim1,f", po::value<int>()->default_value(1), "First base of mate 1")
        ("trim2,2", po::value<int>()->default_value(1), "First base of mate 2")
        ("trim3,3", po::value<int>()->default_value(1), "First base of singleton")
        ("header-transform", po::value<bool>()->default_value(true), "Normalize FASTQ header")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("version")) {
        std::cout << version << std::endl;
        return 0;
    }
    if (!vm.count("input1") || !vm.count("prefix")) {
        std::cout << version << std::endl;
        std::cout << desc << std::endl;
        // 64 is EX_USAGE in sysexits.h (inlined to not rely on a non-portable include)
        return 64;
    }

    try {
        const std::string aname = vm["input1"].as<std::string>();
        const std::string prefix = vm["prefix"].as<std::string>();

        const int min_len = vm["min-length"].as<int>();
        const int min_qual = vm["min-quality"].as<int>();
        const int q_offset = vm["quality-offset"].as<int>();

        const bool header_transform = vm["header-transform"].as<bool>();

        // Subtract one to adjust for 1-based indexing:
        const int f1 = vm["trim1"].as<int>() - 1;
        const int f2 = vm["trim2"].as<int>() - 1;
        const int f3 = vm["trim3"].as<int>() - 1;
        const int min_fastx_p = vm["min-fastx-p"].as<int>();
        const char* maybe_fast = "";
        const trim_mode_t mode = (vm["method"].as<std::string>() == "fastx" ? fastx : solexaqa);

        if (vm.count("fast")) { maybe_fast = "--fast"; }
        unique_istream_ptr a = open_stream(aname);
        opopengzstream single(prefix + ".single.fq.gz", maybe_fast);
        opopengzstream pair1(prefix + ".pair.1.fq.gz", maybe_fast);
        opopengzstream pair2(prefix + ".pair.2.fq.gz", maybe_fast);
        trim_filter_stats st;
        if (vm.count("input2")) {
            unique_istream_ptr b = open_stream(vm["input2"].as<std::string>());
            unique_istream_ptr c;
            if (vm.count("input3")) {
                c = open_stream(vm["input3"].as<std::string>());
            }
            st = rtf_pe(*a, *b, c.get(), single, pair1, pair2, header_transform,
                                trim_options(mode, f1, min_qual + q_offset, min_len, min_fastx_p),
                                trim_options(mode, f2, min_qual + q_offset, min_len, min_fastx_p),
                                trim_options(mode, f3, min_qual + q_offset, min_len, min_fastx_p));
        } else {
            st = rtf_se(*a, single, header_transform, trim_options(mode, f1, min_qual + q_offset, min_len, min_fastx_p));
        }
        std::cout << version
                  << ": using \"" << mode << "\" method to trim and filter reads, ";
        if (q_offset == 64) std::cout << "quality format: \"solexa\"" << std::endl;
        else if (q_offset == 33) std::cout << "quality format: \"sanger\"" << std::endl;

        if (!st.reads) {
            std::cout << "SAMPLE ERROR & EXIT: (from fastq_trim_filter) "
                        "Sample seem to not have any reads passing read trim filter. "
                        "Something wrong?\n";
            return 1;
        }
        std::cout << st.reads << std::endl;
        std::cout << st.sum << std::endl;
        std::cout << st.max << std::endl;
        std::cout << st.inserts << std::endl;

    } catch (std::bad_alloc&) {
        std::cerr << "Out of memory error.\n";
        return 2;
    } catch (const char* error) {
        std::cerr << error << '\n';
        return 2;
    } catch (std::string error) {
        std::cerr << error << '\n';
        return 2;
    }
    return 0;
}

