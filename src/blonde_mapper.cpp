#include <iostream>
#include <string>
#include <getopt.h>
#include <vector>
#include <unordered_map>
#include <time.h>
#include <cmath>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "blonde_alignment.h"
#include "blonde_minimizers.h"

#define VERSION "v0.1.3"

namespace blonde {

constexpr int LENGTH_LIMIT = 5000;

static int help_flag = 0;         /* Flag set by �--help�.    */
static int version_flag = 0;      /* Flag set by �--version�. */
static int algorithm = 0;         /* Flag set by �-a�. */
static int match_cost = 1;        /* Flag set by �-m�. */
static int mismatch_cost = -1;    /* Flag set by �-n�. */
static int gap_cost = -1;         /* Flag set by �-g�. */
static double frequency = 0.001;  /* Flag set by �-f�. */
static int kmer_len = 15;         /* Flag set by �-k�. */
static int window_len = 5;        /* Flag set by �-w�. */

class Sequence {
public:
    Sequence(
        const char* name, std::uint32_t name_len,
        const char* data, std::uint32_t data_len) : name_(name, name_len), data_(data, data_len) {
    }

public:
    std::string name_, data_;
};

class SequenceFastq : public Sequence {
public:
    SequenceFastq(
        const char* name,    std::uint32_t name_len,
        const char* data,    std::uint32_t data_len,
        const char* quality, std::uint32_t quality_len) : Sequence(name, name_len, data, data_len), quality_(quality, quality_len) {
    }

public:
    std::string quality_;
};

void printReferenceGenomesInfo(const std::vector<std::unique_ptr<Sequence>>& genomes) {
    std::cerr << "Reference genome sequences:\n";
    for (int i = 0; i < int(genomes.size()); i++) {
        std::cerr << "\t" << genomes[i]->name_ << " , length = " << genomes[i]->data_.size() << '\n';
    }
}

void printFragmentsInfo(const std::vector<std::unique_ptr<Sequence>>& fragments, const bool fastq = false) {
    uint64_t length_sum = 0;
    std::vector<size_t> lengths(fragments.size());
    for (int i = 0; i < int(fragments.size()); i++) {
        lengths[i] = fragments[i]->data_.size();
        length_sum += lengths[i];
    }
    sort(lengths.begin(), lengths.end(), std::greater<size_t>());

    uint64_t N50 = -1, tmp_sum = 0;
    for (int i = 0; i < int(fragments.size()); i++) {
        tmp_sum += lengths[i];
        if (tmp_sum * 2 >= length_sum) {
            N50 = lengths[i];
            break;
        }
    }
    if (fastq) {
        std::cout << "FASTQ fragments:\n";
    } else {
        std::cout << "FASTA fragments:\n";
    }
    std::cerr << "Number of fragments: " << fragments.size() << '\n';
    std::cerr << "Average length: " << length_sum * 1.0 / fragments.size() << '\n';
    std::cerr << "N50 length: " << N50 << '\n';
    std::cerr << "Minimal length: " << lengths.back() << '\n';
    std::cerr << "Maximal length: " << lengths.front() << '\n';
}

void printMinimizerInfo(const std::vector<std::unique_ptr<Sequence>>& sequences) {
    std::unordered_map<unsigned int, unsigned int> minimizers_map;
    int num_of_singletons = 0;
    for (int i = 0; i < int(sequences.size()); i++) {
        std::vector<minimizers::Kmer> minimizers = minimizers::Minimize(sequences[i]->data_.c_str(), 
                                                                        sequences[i]->data_.size(), 
                                                                        kmer_len, window_len);
        for (int j = 0; j < int(minimizers.size()); j++) {
            minimizers_map[std::get<0>(minimizers[j])]++;
        }
    }
    
    std::vector<unsigned int> minimizers_vec;
    minimizers_vec.reserve(minimizers_map.size());
    for(auto& entry : minimizers_map) {
        if (entry.second == 1) num_of_singletons++;
        minimizers_vec.push_back(entry.second);
    }
    sort(minimizers_vec.begin(), minimizers_vec.end(), std::greater<>());
    size_t minimizers_to_skip = std::ceil(minimizers_map.size() * frequency);
    if (minimizers_to_skip >= minimizers_map.size()) minimizers_to_skip = minimizers_map.size() - 1;

    std::cout << "Number of distinct minimizers: " << minimizers_map.size() << std::endl;
    std::cout << "Fraction of singletons: " << ((double) num_of_singletons / minimizers_map.size()) << std::endl;
    std::cout << "Number of occurrences of the most frequent minimizer with top " << frequency * 100
              << "% most frequent ignored: " << minimizers_vec[minimizers_to_skip] << std::endl;
}

void processGenomes(
    const std::vector<std::unique_ptr<Sequence>>& genomes,
    const std::vector<std::unique_ptr<Sequence>>& fragments,
    const bool fastq = false) {

    printReferenceGenomesInfo(genomes);
    printFragmentsInfo(fragments, fastq);

    std::vector<Sequence> short_fragments;
    for (const auto& p : fragments) {
        if (p->data_.size() < LENGTH_LIMIT) {
            short_fragments.push_back(*p);
        }
    }

    const Sequence first_fragment = short_fragments[rand() % short_fragments.size()];
    const Sequence second_fragment = short_fragments[rand() % short_fragments.size()];

    unsigned int first_len = first_fragment.data_.size();
    unsigned int second_len = second_fragment.data_.size();

    std::cout << "\nAligning two random sequences: \n";
    std::cout << "Target sequence:\n";
    std::cout << second_fragment.name_ << "\n"; 
    std::cout << "Query Sequence:\n";
    std::cout << first_fragment.name_ << "\n";

    std::string cigar;
    unsigned int target_begin;

    int64_t align_score = alignment::Align(first_fragment.data_.c_str(),
                                           first_len,
                                           second_fragment.data_.c_str(),
                                           second_len,
                                           (alignment::AlignmentType)algorithm,
                                           match_cost,
                                           mismatch_cost,
                                           gap_cost,
                                           &cigar,
                                           &target_begin);
    std::cout << "Align result: " << align_score << std::endl;
    std::cout << "Cigar str: " << cigar << std::endl;
    std::cout << "Target begin: " << target_begin << std::endl;

    //Minimizers from reference genome and fragments
    std::cout << "\nMINIMIZER INFO:\n";
    std::cout << "Reference genome:\n";
    printMinimizerInfo(genomes);
    std::cout << "\nFragments:\n";
    printMinimizerInfo(fragments);
}

/* Modificiran primjer https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html */
/* Pojasnjenje primjera https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Options.html */

const std::string HELP_MESSAGE = "blonde_mapper usage: \n\n"
                                 "flags: \n"
                                 "-h or --help     prints help message \n"
                                 "-v or --version  prints version      \n"
                                 "-a <type>        sets the alignment type. 0 for local (default), 1 for global, 2 for semi-global \n"
                                 "-m <X>           sets the match cost to X. default is 1 \n"
                                 "-n <X>           sets the mismatch cost to X. default is 1 \n"
                                 "-g <X>           sets the gap cost to X. default is 1 \n"
                                 "-k <X>           sets the kmer length to to X. default is 15 \n"
                                 "-w <X>           sets the minimizer window length to X. default is 5 \n"
                                 "-f <X>           sets the fraction of top minimizers to ignore to X. default is 0.001 \n"
                                 "\nblonde_mapper takes two filenames as command line arguments: \n"
                                 "The first file will contain a reference genome in FASTA format,\n"
                                 "while the second file will contain a set of fragments in either\n"
                                 "FASTA or FASTQ format.\n\n";
}

int main (int argc, char **argv) {
    using namespace blonde;
    srand (time(NULL)); /* initialize random seed: */
    int c;              /* result variable for getopt_long function */

    while (1) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                {"help",    no_argument, &help_flag,    1},
                {"version", no_argument, &version_flag, 1},
                {"algorithm", required_argument, &algorithm, 1},
                {"match", required_argument, &match_cost, 1},
                {"mismatch", required_argument, &mismatch_cost, 1},
                {"gap", required_argument, &gap_cost, 1},
                {"kmerlen", required_argument, &kmer_len, 1},
                {"windowlen", required_argument, &window_len, 1},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "hva:m:n:g:k:w:f:",
                        long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
        case 0:
            break;

        case 'h':
            help_flag = 1;
            break;

        case 'v':
            version_flag = 1;
            break;

        case 'a':
            algorithm = std::stoi(optarg);
            break;

        case 'm':
            match_cost = std::stoi(optarg);
            break;

        case 'n':
            mismatch_cost = -std::stoi(optarg);
            break;

        case 'g':
            gap_cost = -std::stoi(optarg);
            break;

        case 'k':
            kmer_len = std::stoi(optarg);
            break;

        case 'w':
            window_len = std::stoi(optarg);
            break;

        case 'f':
            frequency = std::stod(optarg);
            break;

        case '?':
            /* getopt_long already printed an error message. */
            break;

        default:
            abort ();
        }
    }

    if (help_flag) {
        std::cout << HELP_MESSAGE;
    } else if (version_flag) {
        std::cout << VERSION << std::endl;
    } if (optind < argc) {
        auto genome_parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[optind]);
        auto genomes = genome_parser->Parse(-1);
        try {
            auto fragment_parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[optind + 1]);
            auto fragments = fragment_parser->Parse(-1);
            processGenomes(genomes, fragments);
        } catch (std::invalid_argument e) {
            auto fragment_parser = bioparser::Parser<SequenceFastq>::Create<bioparser::FastqParser>(argv[optind + 1]);

            // parse in chunks
            std::vector<std::unique_ptr<Sequence>> fragments;
            std::uint32_t chunk_size = 500 * 1024 * 1024;  // 500 MB
            for (auto t = fragment_parser->Parse(chunk_size); !t.empty(); t = fragment_parser->Parse(chunk_size)) {
                fragments.insert(
                    fragments.end(),
                    std::make_move_iterator(t.begin()),
                    std::make_move_iterator(t.end()));
            }
            processGenomes(genomes, fragments, true);
        }
    }

    return 0;
}

