#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdio>

const char* version = "Quality statistics [EMBL Version 2015-01-26]";

inline
bool skip_eol(std::istream& in) {
    return in.ignore(4096, '\n');
}

const int min_quality_value = -10;

struct position_stats {
    position_stats()
        :nuc_counts(6)
        ,qual_counts(256)
        { }
    std::vector<int> nuc_counts;
    std::vector<int> qual_counts;
};

int nucleotide_index(char c) {
    switch(c) {
        case 'a':
        case 'A': return 0;
        case 'c':
        case 'C': return 1;
        case 'g':
        case 'G': return 2;
        case 't':
        case 'T': return 3;
        case 'n':
        case 'N': return 4;
        default: return 5;
    }
}


void stats(std::istream& input, std::ostream& output, const int quality_offset) {
    unsigned size = 0;
    std::vector<position_stats> statistics;
    int nucleotide_table[256];
    for (int i = 0; i < 256; ++i) nucleotide_table[i] = nucleotide_index(char(i));

    std::string seq;
    std::string qualstr;
    int line = 0;
    while (skip_eol(input) && std::getline(input, seq) && skip_eol(input) && std::getline(input, qualstr)) {
        const unsigned s = seq.size();
        if (qualstr.size() != s) {
            std::ostringstream err;
            err << "FASTQ Format error in lines " << line << "-" << (line + 3)
                    << ": sequence has length " << s << ", and quality string has length "
                    << qualstr.size() << "." << std::endl;
            throw err.str();
        }
        line += 4;
        while (s > size) {
            statistics.push_back(position_stats());
            ++size;
        }
        for (unsigned i = 0; i != s; ++i) {
            const int nuc_index = nucleotide_table[ int(seq[i]) ];
            ++statistics[i].nuc_counts[nuc_index];
            ++statistics[i].qual_counts[int(qualstr[i])];
        }
    }
    for (unsigned i = 0; i != statistics.size(); ++i) {
        std::vector<int> rescaled(256);
        for (int qv = 0; qv < quality_offset + min_quality_value; ++qv) {
            if (statistics[i].qual_counts[qv]) {
                std::cerr << "Found quality value below expected offset!\n";
                std::cerr << "Please check -Q option!\n";
            }
        }
        for (int qv = 0; qv < (256 - quality_offset - min_quality_value); ++qv) {
              rescaled[qv] = statistics[i].qual_counts[qv + quality_offset + min_quality_value];
        }
        statistics[i].qual_counts.swap(rescaled);
    }
    output << "column\tcount\tmin\tmax\tsum\tmean\tQ1\tmed\tQ3\tIQR\tlW\trW\tA_Count\tC_Count\tG_Count\tT_Count\tN_Count\tMax_count\n";
    int count0 = -1;
    const int negative_quality = min_quality_value-1;
    for (unsigned i = 0; i != size; ++i) {
        int count = 0;
        int minQ = 256;
        int maxQ = negative_quality;
        long long int sumQ = 0;
        int medQ = negative_quality;
        int q1Q  = negative_quality;
        int q3Q  = negative_quality;
        for (int qvi = 0; qvi != int(statistics[i].qual_counts.size()); ++qvi) {
            const int curc = statistics[i].qual_counts[qvi];
            const int qv = qvi + min_quality_value;
            if (curc) {
                if (qv < minQ) minQ = qv;
                if (qv > maxQ) maxQ = qv;
                count += curc;
                sumQ += static_cast<long long>(curc) * static_cast<long long>(qv);
            }
        }
        if (i == 0) count0 = count;
        int accum = 0;
        for (int qvi = 0; qvi != int(statistics[i].qual_counts.size()); ++qvi) {
            accum += statistics[i].qual_counts[qvi];

            const int qv = qvi + min_quality_value;
            if ((q1Q == negative_quality) && accum > (count / 4)) q1Q = qv;
            if ((medQ == negative_quality) && accum > (count / 2)) medQ = qv;
            if ((q3Q == negative_quality) && accum > (3 * count / 4)) q3Q = qv;
        }
        int lWQ = std::max<int>(q1Q - (q3Q - q1Q)*3/2, minQ);
        int rWQ = std::min<int>(q3Q + (q3Q - q1Q)*3/2, maxQ);

        char meanbuffer[128];
        std::snprintf(meanbuffer, sizeof(meanbuffer), "%3.2f", double(sumQ)/count);

        output
            << (i+1) << '\t' // column  = column number (1 to 36 for a 36-cycles read solexa file)
            << count << '\t' // count   = number of bases found in this column.
            << minQ << '\t' // min     = Lowest quality score value found in this column.
            << maxQ << '\t' // max     = Highest quality score value found in this column.
            << sumQ << '\t' // sum     = Sum of quality score values for this column.
            << meanbuffer << '\t' // mean    = Mean quality score value for this column.
            << q1Q << '\t' // Q1  = 1st quartile quality score.
            << medQ << '\t' // med = Median quality score.
            << q3Q << '\t' // Q3  = 3rd quartile quality score.
            << (q3Q - q1Q) << '\t' // IQR = Inter-Quartile range (Q3-Q1).
            << lWQ << '\t' // lW  = 'Left-Whisker' value (for boxplotting).
            << rWQ << '\t' // rW  = 'Right-Whisker' value (for boxplotting).
            << statistics[i].nuc_counts[0] << '\t' // A_Count = Count of 'A' nucleotides found in this column.
            << statistics[i].nuc_counts[1] << '\t' // C_Count = Count of 'C' nucleotides found in this column.
            << statistics[i].nuc_counts[2] << '\t' // G_Count = Count of 'G' nucleotides found in this column.
            << statistics[i].nuc_counts[3] << '\t' // T_Count = Count of 'T' nucleotides found in this column.
            << statistics[i].nuc_counts[4] << '\t' // N_Count = Count of 'N' nucleotides found in this column.
            << count0 << std::endl; // max-count = max. number of bases (in all cycles)
    }
}

int parse_args(int argc, char* argv[], std::string& iname, std::string& oname, int& qoffset) {
    const char* usage = " [-i INPUT] [-o OUTPUT]\n";
    iname = "-";
    oname = "-";
    int argi = 1;

    while (argi < argc) {
        if (argv[argi] == std::string("-v") || argv[argi] == std::string("--version")) {
            std::cout << version << std::endl;
            // 64 is EX_USAGE in sysexits.h (inlined to not rely on a non-portable include)
            return 64;
        } else if (argv[argi] == std::string("-h") || argv[argi] == std::string("--help")) {
            std::cout << version << std::endl;
            std::cout << "\nUsage:\n\t"
                    << argv[0] << usage << std::endl;
            // 64 is EX_USAGE in sysexits.h (inlined to not rely on a non-portable include)
            return 64;
        } else if (argv[argi][0] == '-') {
            ++argi;
            if (argi == argc) {
                std::cerr << argv[0] << usage;
                std::cerr << "\n[Wrong number of arguments]\n";
                return 1;
            }
            std::string value = argv[argi];
            if (argv[argi-1] == std::string("-i")) iname = value;
            else if (argv[argi-1] == std::string("-o")) oname = value;
            else if (argv[argi-1] == std::string("-Q")) {
                std::string argvalue(argv[argi]);
                std::istringstream in(argvalue);
                if (!(in >> qoffset)) {
                    std::cerr << "Could not parse -Q option value as integer.\n";
                    return 1;
                }
                if (qoffset < 0) {
                    std::cerr << "Qoffset must be a positive number.\n";
                    return 1;
                }
            } else {
                std::cerr << usage;
                std::cerr << "\n[Only -i and -o are accepted as arguments]\n";
                return 1;
            }
            ++argi;
        } else {
            std::cerr << usage;
            return 1;
        }
    }
    return 0;
}

int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::string iname, oname;
    int qoffset = 64;
    int r = parse_args(argc, argv, iname, oname, qoffset);
    if (r) return r;
    try {
        std::istream* input;
        std::ostream* output;
        if (iname == "-") {
            input = &std::cin;
        } else {
            input = new std::ifstream(iname.c_str());
        }
        if (oname == "-") {
            output = &std::cout;
        } else {
            output = new std::ofstream(oname.c_str());
        }
        stats(*input, *output, qoffset);
        if (iname != "-") delete input;
        if (oname != "-") delete output;
    } catch (std::string error) {
        std::cerr << error << '\n';
        return 2;
    } catch (const char* error) {
        std::cerr << error << '\n';
        return 2;
    }
    return 0;
}

