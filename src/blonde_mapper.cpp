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

#define VERSION "v0.1.5"

namespace blonde {

constexpr int LENGTH_LIMIT = 5000;
constexpr int EPSILON_BAND = 500;
constexpr int MAPPING_QUALITY = 255;

static int help_flag = 0;         /* Flag set by �--help�.    */
static int version_flag = 0;      /* Flag set by �--version�. */
static int algorithm = 0;         /* Flag set by �-a�. */
static int match_cost = 1;        /* Flag set by �-m�. */
static int mismatch_cost = -1;    /* Flag set by �-n�. */
static int gap_cost = -1;         /* Flag set by �-g�. */
static double frequency = 0.001;  /* Flag set by �-f�. */
static int kmer_len = 15;         /* Flag set by �-k�. */
static int window_len = 5;        /* Flag set by �-w�. */
static int cigar_flag = 0;        /* Flag set by �-c�. */
static int thread_num = 1;        /* Flag set by �-t�. */

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

// q_begin, q_end, t_begin, t_end
using CandidateRegion = std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>;
// relative strand 0 = same, 1 = different, i - i' / i + i', i', 
using Match = std::tuple<bool, unsigned int, unsigned int>;

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

void makeIndex(
    const std::unique_ptr<Sequence> sequence, 
    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>>& index) {
    
    std::vector<minimizers::Kmer> minimizers = minimizers::Minimize(sequence->data_.c_str(), 
                                                                    sequence->data_.size(), 
                                                                    kmer_len, window_len);
    for (auto minimizer : minimizers) {
        index[std::get<0>(minimizer)].emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
    }
}

void makeCleanedReferenceIndex(
    const std::unique_ptr<Sequence> sequence, 
    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>>& index) {
    
    std::cout << "Making reference index..." << std::endl;
    makeIndex(sequence, index);
    std::vector<std::pair<unsigned int, unsigned int>> minimizers_occurance;
    minimizers_occurance.reserve(index.size());
    int num_of_singletons = 0;
    for (auto& entry : index) {
        if (entry.second.size() == 1) num_of_singletons++;
        minimizers_occurance.emplace_back(std::make_pair(entry.second.size(), entry.first));
    }
    std::cout << "Sorting reference minimizers by occurance..." << std::endl;
    sort(minimizers_occurance.begin(), minimizers_occurance.end(), std::greater<std::pair<unsigned int, unsigned int>>());
    size_t minimizers_to_skip = std::ceil(index.size() * frequency);
    if (minimizers_to_skip >= index.size()) minimizers_to_skip = index.size() - 1;
    
    std::cout << "Number of distinct minimizers: " << index.size() << std::endl;
    std::cout << "Fraction of singletons: " << ((double) num_of_singletons / index.size()) << std::endl;
    std::cout << "Number of occurrences of the most frequent minimizer with top " << frequency * 100
              << "% most frequent ignored: " << minimizers_occurance[minimizers_to_skip] << std::endl;
    
    std::cout << "Removing too frequent minimizers from reference index..." << std::endl;
    for (int i = 0; i < minimizers_to_skip; i++) {
        index.erase(minimizers_occurance[i].second);
    }
    std::cout << "Done." << std::endl;
}

// INSPIRED BY https://www.geeksforgeeks.org/construction-of-longest-monotonically-increasing-subsequence-n-log-n/?ref=rp
int GetCeilIndex(std::vector<Match>& curr_cluster, std::vector<unsigned int>& T, int l, int r, 
                 unsigned int key) 
{ 
    while (r - l > 1) { 
        int m = l + (r - l) / 2; 
        if (std::get<2>(curr_cluster[T[m]]) >= key) 
            r = m; 
        else
            l = m; 
    } 
  
    return r; 
} 

// INSPIRED BY https://www.geeksforgeeks.org/construction-of-longest-monotonically-increasing-subsequence-n-log-n/?ref=rp
void longestIncreasingSubsequence(
    std::vector<Match>& curr_cluster,
    std::vector<Match>& result) {
    
    // Add boundary case, when array n is zero 
    // Depend on smart pointers 
    // curr_cluster = arr in comments
    size_t n = curr_cluster.size();
    std::vector<unsigned int> tailIndices(n, 0); // Initialized with 0 
    std::vector<unsigned int> prevIndices(n, -1); // initialized with -1 
  
    unsigned int len = 1; // it will always point to empty location 
    for (int i = 1; i < int(n); i++) { 
        if (std::get<2>(curr_cluster[i]) < std::get<2>(curr_cluster[tailIndices[0]])) { 
            // new smallest value 
            tailIndices[0] = i; 
        } 
        else if (std::get<2>(curr_cluster[i]) > std::get<2>(curr_cluster[tailIndices[len - 1]])) { 
            // arr[i] wants to extend largest subsequence 
            prevIndices[i] = tailIndices[len - 1]; 
            tailIndices[len++] = i; 
        } 
        else { 
            // arr[i] wants to be a potential condidate of 
            // future subsequence 
            // It will replace ceil value in tailIndices 
            int pos = GetCeilIndex(curr_cluster, tailIndices, -1, 
                                   len - 1, std::get<2>(curr_cluster[i])); 
  
            prevIndices[i] = tailIndices[pos - 1]; 
            tailIndices[pos] = i; 
        } 
    } 
  
    std::vector<unsigned int> lisIndexes;
    lisIndexes.reserve(len);
    for (int i = tailIndices[len - 1]; i >= 0; i = prevIndices[i]) 
        lisIndexes.emplace_back(i);
    std::reverse(lisIndexes.begin(), lisIndexes.end());
    for(auto i : lisIndexes) {
        result.emplace_back(curr_cluster[i]);
    }
    return;
}

unsigned int getQueryLoc(Match& x) {
    if (std::get<0>(x)) {
        return std::get<1>(x) - std::get<2>(x);
    } else {
        return std::get<1>(x) + std::get<2>(x);
    }
}

unsigned int getRefLoc(Match& x) {
    return std::get<2>(x);
}

unsigned int getRelStrandString(Match& x) {
    return std::get<0>(x) ? '+' : '-';
}

int kmerOverlap(unsigned int prev, unsigned int curr) {
    if (curr - prev < kmer_len) {
        return curr - prev;
    } else {
        return kmer_len;
    }
}

int col10Approx(std::vector<Match>& lis_of_cluster) {
    int approx_val_q = kmer_len;
    int approx_val_r = kmer_len;
    unsigned int prev_q = getQueryLoc(lis_of_cluster[0]);
    unsigned int prev_r = getRefLoc(lis_of_cluster[0]);
    bool first = true;
    for (auto match : lis_of_cluster) {
        if (first) {
            first = false;
            continue;
        }
        unsigned int curr_q = getQueryLoc(match);
        unsigned int curr_r = getRefLoc(match);
        approx_val_q += kmerOverlap(prev_q, curr_q);
        approx_val_r += kmerOverlap(prev_r, curr_r);
        prev_q = curr_q;
        prev_r = curr_r;
    }
    return std::min(approx_val_r, approx_val_q);
}

void splitIntoClustersAndFindBest(
    std::vector<Match>& matches,
    std::vector<Match>& match_cluster) {

    std::vector<Match> curr_cluster;
    int max_col10 = 0;
    curr_cluster.emplace_back(matches[0]);
    for (int i = 1; i <= int(matches.size()); i++) {
        if (i == int(matches.size()) ||
            std::get<0>(matches[i]) != std::get<0>(matches[i - 1]) ||
            std::get<1>(matches[i]) - std::get<1>(matches[i - 1]) >= EPSILON_BAND) {
            
            std::vector<Match> lis_of_cluster;
            longestIncreasingSubsequence(curr_cluster, lis_of_cluster);
            int curr_col10 = col10Approx(lis_of_cluster);
            if (curr_col10 > max_col10) {
                max_col10 = curr_col10;
                match_cluster.clear();
                match_cluster = lis_of_cluster;
            }
            curr_cluster.clear()
        }
        curr_cluster.emplace_back(matches[i]);
    }
}

bool matchComparator(const Match& a, const Match& b) {
    if (std::get<0>(a) != std::get<0>(b)) {
        return !std::get<0>(a);
    } else {
        if(std::get<1>(a) != std::get<1>(b)) {
            return std::get<1>(a) < std::get<1>(b);
        } else {
            return std::get<2>(a) < std::get<2>(b);
        }
    }
}

void bestMatchCluster(
    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>>& fragment_index,
    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>>& reference_index,
    std::vector<Match>& match_cluster) {
    
    std::vector<Match> matches;
    for (auto f_entry : fragment_index) {
        if (reference_index.count(f_entry.first) != 0) {
            for (auto f_location : f_entry.second) {
                for (auto r_location : reference_index[f_entry.first]) {
                    bool diff_strand = f_location.second ^ r_location.second;
                    unsigned int relative_position;
                    if (diff_strand) {
                        relative_position = f_location.first + r_location.first;
                    } else {
                        relative_position = f_location.first - r_location.first;
                    }
                    matches.emplace_back(std::make_tuple(diff_strand, relative_position, r_location.first));
                }
            }
        }
    }

    if (!matches.empty()) {
        sort(matches.begin(), matches.end(), matchComparator);
        splitIntoClustersAndFindBest(matches, match_cluster);
    }

}

std::string getPaf(
    const std::unique_ptr<Sequence> fragment,
    const std::unique_ptr<Sequence> reference,
    std::vector<Match>& match_cluster) {

    unsigned int q_begin = getQueryLoc(match_cluster.front());
    unsigned int q_end = getQueryLoc(match_cluster.back()) + kmer_len;
    unsigned int t_begin = getRefLoc(match_cluster.front());
    unsigned int t_end = getRefLoc(match_cluster.back()) + kmer_len;

    std::string result = fragment->name_ +                            '\t' +
                         fragment->data_.size() +                     '\t' +
                         std::to_string(q_begin) +                    '\t' +
                         std::to_string(q_end) +                      '\t' +
                         getRelStrandString(match_cluster.front()) +  '\t' +
                         reference->name_ +                           '\t' +
                         reference->data_.size() +                    '\t' +
                         std::to_string(t_begin) +                    '\t' +
                         std::to_string(t_end) +                      '\t' +
                         std::to_string(col10Approx(match_cluster)) + '\t';
    
    std::string cigar;
    if (cigar_flag) {
        unsigned int target_begin;

        alignment::Align(fragment.data_.c_str() + q_begin,
                         q_end - q_begin,
                         fragment.data_.c_str() + t_begin,
                         t_end - t_begin,
                         (alignment::AlignmentType) algorithm,
                         match_cost,
                         mismatch_cost,
                         gap_cost,
                         &cigar,
                         &target_begin);
        int col10 = 0;
        int col11 = 0;
        int num = 0;
        for (auto x : cigar) {
            if (std::isdigit(x)) {
                num *= 10;
                num += std::atoi(&x);
            } else {
                if (x == '=') col10 += num;
                if (x != 'S') col11 += num;
                num = 0;
            }
        }
        result += std::to_string(col10) + '\t' +
                  std::to_string(col11) + '\t';
    } else {
        int col11_approx = std::max(q_end - q_begin, t_end - t_begin);
        resut += std::to_string(col10Approx(match_cluster)) + '\t' +
                 std::to_string(col11_approx) +               '\t';
    }

    result += std::to_string(MAPPING_QUALITY) + '\t';
    if (cigar_flag) result += "cg:Z:" + cigar;
    return result;
}

std::string mapFragmentsToReference(
    const std::vector<std::unique_ptr<Sequence>>& fragments,
    const std::unique_ptr<Sequence> reference,
    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>>& reference_index,
    int fragments_begin, int fragments_end) {  // vidi jel zadnja dva parametra moraju biti reference il nesto 

    std::string result = "";
    for (int i = fragments_begin; i < fragments_end; i++) {
        std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> fragment_index;
        makeIndex(fragments[i], fragment_index);
        std::vector<Match> match_cluster;
        bestMatchCluster(fragment_index, reference_index, match_cluster);
        if(match_cluster.empty()) {
            // hmmm idk smisli nes
            continue;
        }
        std::string paf = getPaf(fragments[i], reference, match_cluster);
        result += paf + "\n";
    }
    return result;
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
                                 "-c <X>           if present the cigar string will be calculated \n"
                                 "-t <X>           sets the number of threads, default is 1 \n"
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

        case 'c':
            cigar_flag = 1;
            break;
        
        case 't':
            thread_num = std::stoi(optarg);
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

