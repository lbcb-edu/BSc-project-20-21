#include <iostream>
#include <getopt.h>
#include "include/bioparser/include/bioparser/parser.hpp"
#include "include/bioparser/include/bioparser/fastq_parser.hpp"
#include "include/bioparser/include/bioparser/fasta_parser.hpp"
#include "include/thread_pool/include/thread_pool/thread_pool.hpp"
#include "include/thread_pool/include/thread_pool/semaphore.hpp"
#include "include/brown_alignment.hpp"
#include "include/brown_minimizer.hpp"
#include <stdlib.h>
#include <time.h>
#include <map>
#include <set>
#include <ctype.h>

#define VERSION "v1.0"
#define DEFAULT_KMER_LENGTH 15
#define DEFAULT_WINDOW_LENGTH 5
#define DEFAULT_MINIMIZER_FREQUENCY 0.001

struct Sequence {
  public:
    std::string sequenceName;
    std::string sequenceSequence;   
    std::string sequenceQuality;

    Sequence( 
        const char* name, std::uint32_t nameLength,
        const char* sequence, std::uint32_t sequenceLength,
        const char* quality = nullptr, std::uint32_t qualityLength = 0) :
            sequenceName(name, nameLength),
            sequenceSequence(sequence, sequenceLength) {
                if (quality != nullptr) {
                    sequenceQuality = std::string(quality, qualityLength);
                }
            }
};

std::vector<std::unique_ptr<Sequence>> referenceGenom;
std::vector<std::unique_ptr<Sequence>> fragments;
std::vector<std::tuple<unsigned int, unsigned int, bool>> reference_minimizers;
brown::AlignmentType type;
int match, gap, mismatch;
unsigned int kmer_length = DEFAULT_KMER_LENGTH, window_length = DEFAULT_WINDOW_LENGTH;
int thread_number;
double frequency = DEFAULT_MINIMIZER_FREQUENCY;
bool cigar_flag;


static option long_options[] = {
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'v'},
    {0, 0, 0, 0}
};

static std::string help = "brown_mapper \n"
                            "-h or --help for help\n"
                            "-v or --version for programs version\n\n"
                            "Alignment arguments:\n"
                            "-m   value for matching\n"
                            "-n   value for mismatching\n"
                            "-g   value for gap\n"
                            "-a   alignment type (GLOBAL, SEMIGLOBAL or LOCAL)\n"
                            "-c   cigar string enabled\n\n"
                            "Minimizer arguments:\n"
                            "-k   k-mer length (default: " + std::to_string(DEFAULT_KMER_LENGTH) + ")\n"
                            "-w   window length (default: " + std::to_string(DEFAULT_WINDOW_LENGTH) + ")\n"
                            "-f   top f frequent minimizers that will not be taken in account (default: " 
                            + std::to_string(DEFAULT_MINIMIZER_FREQUENCY) + ")\n"
                            "Program accepts two files as floating arguments.\n"
                            "The first file needs to contain a reference genome in FASTA format.\n"
                            "The second file needs to contain a set of fragments\n"
                            "in either FASTA or FASTQ format.";

static std::string genomLine = "\n-------------------------\n"
                                "      REFERENCE GENOM     \n"
                                "-------------------------\n\n";

static std::string fragmentLine = "\n-------------------------\n"
                                "       FRAGMENTS     \n"
                                "-------------------------\n\n";


bool compareFragmentLengths(std::unique_ptr<Sequence>& s1, std::unique_ptr<Sequence>& s2 ) {
    return s1->sequenceSequence.length() > s2->sequenceSequence.length();
}

void printsFragmentsStats(std::vector<std::unique_ptr<Sequence>>& fragments) {

    int sumLength = 0;
    sort(fragments.begin(), fragments.end(), compareFragmentLengths);
    for (int i = 0; i < fragments.size(); i++ ) 
        sumLength += fragments[i]->sequenceSequence.length();
    
    int N50sum = 0;
    int N50;
    for (int i = 0; i < fragments.size(); i++ ) {
        N50sum += fragments[i]->sequenceSequence.length();
        if (N50sum > 0.5 * sumLength) {
            N50 = fragments[i]->sequenceSequence.length();
            break;
        }
    }
        
    std::cerr << "Number of sequences: " <<  fragments.size() << std::endl;
    std::cerr << "Total length of all fragments: " << sumLength << std::endl;
    std::cerr << "Largest fragment: " << fragments[0]->sequenceName << std::endl;
    std::cerr << "      length: " << fragments[0]->sequenceSequence.length() << std::endl;
    std::cerr << "Smallest fragment: " << fragments[fragments.size() - 1]->sequenceName << std::endl;
    std::cerr << "      length: " << fragments[fragments.size() - 1]->sequenceSequence.length() << std::endl;
    std::cerr << "Average length: " << ((double) sumLength) / fragments.size() << std::endl;
    std::cerr << "N50 length: " << N50 << std::endl << std::endl;
}

bool compareMinimizerOccurence(std::pair<unsigned int, unsigned int>& a, std::pair<unsigned int, unsigned int>& b) { 
    return a.second < b.second; 
}

bool comparePairOfMatches(std::pair<unsigned int, unsigned int>& a, std::pair<unsigned int, unsigned int>& b) {
    if (a.second != b.second)
        return a.second < b.second;

    return a.first < b.first;
}

void cutOffTooFrequentMinimizers(std::vector<std::tuple<unsigned int, unsigned int, bool>>& minimizers,
                                std::string name, bool ref_genom) {
    std::map<unsigned int, unsigned int> minimizers_map;
    unsigned int singletons = 0;
    for(unsigned int i = 0; i < minimizers.size(); i++) {
        unsigned int currentMinimizer = std::get<0>(minimizers[i]);
        if(minimizers_map.count(currentMinimizer) == 0) {
            minimizers_map.insert(std::pair<unsigned int,unsigned int>(currentMinimizer, 1));
            singletons++;
        } else { 
            minimizers_map[currentMinimizer]++;
            if (minimizers_map[currentMinimizer] == 2)
                singletons--;
        }
    }
    std::cerr << "Finished counting occurences of minimizers.." << std::endl << std::endl;
    unsigned int fth_minimizer_index = minimizers_map.size() - 1 - frequency * minimizers_map.size();
    
    std::vector<std::pair<unsigned, unsigned int>> sorted_map_values(minimizers_map.begin(), minimizers_map.end());
    sort(sorted_map_values.begin(), sorted_map_values.end(), compareMinimizerOccurence);

    std::string seq = ref_genom == true? "reference genom" : "fragment"; 
    std::cerr << "Total number of minimizers in " << seq << " " << name << ": " << minimizers.size() << "." << std::endl;
    std::cerr << "The fraction of singletone minimizers in " << seq << ":  " << ((double) singletons) / minimizers.size()
                << "." << std::endl;
    std::cerr << "The number of occurrences of the most frequent minimizer when the top " << frequency 
                << " frequent minimizers are not taken in account: " << sorted_map_values[fth_minimizer_index].second << std::endl << std::endl;
    
    std::vector<unsigned int> unwanted_minimizers;
    for (unsigned int i = fth_minimizer_index; i < sorted_map_values.size(); i++) {
        unwanted_minimizers.push_back(sorted_map_values[i].first);
    }
    std::cerr << "Found unwanted minimizers...\n"
                << "Deleting them..." << std::endl;


    unsigned int original_size = minimizers.size();
    for (unsigned int i = 0; i < minimizers.size(); i++) {
        if (std::find(unwanted_minimizers.begin(), unwanted_minimizers.end(), std::get<0>(minimizers[i])) != unwanted_minimizers.end()) {
            minimizers.erase(minimizers.begin() + i);
        }
        /*for (unsigned int j = 0; j < unwanted_minimizers.size(); j++) {
            if (std::get<0>(minimizers[i]) == unwanted_minimizers[j]) {
                minimizers.erase(minimizers.begin() + i);
                cuurent_size = minimizers.size();
                break;
            }
        }*/
    }

    /*for (unsigned int i = 0; i < unwanted_minimizers.size(); i++) {
        std::tuple<unsigned int, unsigned int, bool> fake_tuple = std::make_tuple(unwanted_minimizers[i], 0, false);
        std::remove_if(minimizers.begin(), minimizers.end(), [fake_tuple] (std::tuple<unsigned int, unsigned int, bool> n) {
                                                                return std::get<0>(fake_tuple) == std::get<0>(n);
                                                                });
    }*/

}

unsigned int GetCeilIndex(std::vector<std::pair<unsigned int, unsigned int>>& matches,
                            std::vector<unsigned int>& T, unsigned int l, unsigned int r,
                            std::pair<unsigned int, unsigned int>key) { 
    while (r - l > 1) { 
        int m = l + (r - l) / 2; 
        if (matches[T[m]].second > key.second) {
            r = m;
            continue;
        } 
        else if (matches[T[m]].second < key.second) {
            l = m;
            continue;
        }

        if (matches[T[m]].first >= key.first) 
            r = m; 
        else
            l = m;  
    } 
  
    return r; 
} 
  
unsigned int LongestIncreasingSubsequence(std::vector<std::pair<unsigned int, unsigned int>>& matches,
                                    std::vector<unsigned int>& tail_indices,
                                    std::vector<unsigned int>& prev_indices) { 

    unsigned int len = 1; 

    std::cerr << "Founding LIS..." << std::endl;
    for (int i = 1; i < matches.size(); i++) { 
        if (matches[i].first < matches[tail_indices[0]].first) { 
            tail_indices[0] = i; 
        } 
        else if (matches[i].first > matches[tail_indices[len - 1]].first
                && matches[i].second != matches[prev_indices[tail_indices[len - 1]]].second) { 
            prev_indices[i] = tail_indices[len - 1]; 
            tail_indices[len++] = i; 
        } 
        else { 
            unsigned int pos = GetCeilIndex(matches, tail_indices, -1, len - 1, matches[i]); 
            prev_indices[i] = tail_indices[pos - 1]; 
            tail_indices[pos] = i; 
        } 
    } 
  
    std::cerr << "LIS compleated" << std::endl;
    return len; 
} 


void mapping(Sequence fragment) {
    std::vector<std::tuple<unsigned int, unsigned int, bool>> fragment_minimizers = brown::Minimize(fragment.sequenceSequence.c_str(), 
                                                                                                    fragment.sequenceSequence.length(),
                                                                                                    kmer_length,
                                                                                                    window_length);
    
    cutOffTooFrequentMinimizers(fragment_minimizers, fragment.sequenceName, false);
    
    std::vector<std::pair<unsigned int, unsigned int>> original_matches;
    std::vector<std::pair<unsigned int, unsigned int>> revcompl_matches;

    std::cerr << "Matching reference and fragment minimizers..." << std::endl;
    
    for (unsigned int i = 0; i < reference_minimizers.size(); i++) {
        
        for (unsigned int j = 0; j < fragment_minimizers.size(); j++) {
            //std::cerr << "Uso u petlju " << i * fragment_minimizers.size() + j << " put\n";
            if (std::get<0>(fragment_minimizers[j]) == std::get<0>(reference_minimizers[i])) {
                original_matches.push_back(std::make_pair(std::get<1>(reference_minimizers[i]), std::get<1>(fragment_minimizers[j])));
            }
            
            unsigned int reverse_complement = brown::getReversedComplKmerValue(std::get<0>(fragment_minimizers[j]), kmer_length);
            if (reverse_complement == std::get<0>(reference_minimizers[i])) {
                revcompl_matches.push_back(std::make_pair(std::get<1>(reference_minimizers[i]), 
                                                        fragment.sequenceSequence.length() - std::get<1>(fragment_minimizers[j]) - kmer_length));
            }
        }

    }

    std::cerr << "Matching compleated" << std::endl;

    std::sort(original_matches.begin(), original_matches.end(), comparePairOfMatches);
    std::sort(revcompl_matches.begin(), revcompl_matches.end(), comparePairOfMatches);

    std::cerr << "Sorting matches compleated" << std::endl;
    std::vector<unsigned int> tail_indices_original(original_matches.size(), 0); 
    std::vector<unsigned int> prev_indices_original(original_matches.size(), UINT32_MAX);

    std::vector<unsigned int> tail_indices_revcompl(original_matches.size(), 0); // Initialized with 0 
    std::vector<unsigned int> prev_indices_revcompl(original_matches.size(), UINT32_MAX); // initialized with -1 

    unsigned int lis_length_original = LongestIncreasingSubsequence(original_matches, tail_indices_original, prev_indices_original);
    unsigned int lis_length_revcompl = LongestIncreasingSubsequence(revcompl_matches, tail_indices_revcompl, prev_indices_revcompl);

    std::vector<std::pair<unsigned int, unsigned int>> lis_vector;
    char strand;
    if (lis_length_original > lis_length_revcompl) {
        strand = '+';
        for (int i = tail_indices_original[lis_length_original - 1]; i >= 0; i = prev_indices_original[i]) {
            lis_vector.push_back(original_matches[i]);
        }
    } else { 
        strand = '-';
        for (int i = tail_indices_revcompl[lis_length_revcompl - 1]; i >= 0; i = prev_indices_revcompl[i]) {
            lis_vector.push_back(revcompl_matches[i]);
        }
    }

    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    int align_result = brown::Align(referenceGenom[0]->sequenceSequence.c_str() + lis_vector.back().first,
                                    lis_vector.front().first - lis_vector.back().first,
                                    fragment.sequenceSequence.c_str() + lis_vector.back().second,
                                    lis_vector.front().second - lis_vector.back().second,
                                    type, match, mismatch, gap, cigar, target_begin);

    unsigned int block_length = 0, num_of_matches = 0;
    for (unsigned int i = 0; i < (*cigar).length(); i++) {
        if (isdigit((*cigar)[i])) {
            std::string numbers;
            
            while (isdigit((*cigar)[i]))
                numbers.append((std::to_string((*cigar)[i++])));

            block_length += std::stoi(numbers);
            if ((*cigar)[i] == 'M')
                num_of_matches += std::stoi(numbers);
        }
    }

    std::cout << referenceGenom[0]->sequenceName
                << "\t" << referenceGenom[0]->sequenceSequence.length()
                << "\t" << lis_vector.back().first
                << "\t" << lis_vector.front().first
                << "\t" << strand
                << "\t" << fragment.sequenceName
                << "\t" << fragment.sequenceSequence.length()
                << "\t" << lis_vector.back().second
                << "\t" << lis_vector.front().second
                << "\t" << num_of_matches
                << "\t" << block_length
                << "\t" << 255 - (((double) num_of_matches) / block_length) * 255;

    if (cigar_flag)
        std::cout << "\t" << *cigar;
    
    std::cout << std::endl;

    delete cigar;
    delete target_begin;

}



int main(int argc, char* argv[]) { 

    bool match_flag;
    bool mismatch_flag;
    bool gap_flag;
    bool type_flag;

    int c;
    while ((c = getopt_long(argc, argv, "m:g:n:a:k:w:f:t:hvc", long_options, 0)) != -1) {
        switch (c){
            case 'h' :
                std::cerr << help << std::endl;
                return 0;
            case 'v' :
                std::cerr << VERSION << std::endl;
                return 0;
            case 'm':
                match = atoi(optarg);
                std::cerr << "Match is : " << optarg << std::endl;
                match_flag = true;
                break;
            case 'n' :
                mismatch = atoi(optarg);
                std::cerr << "Mismatch is : " << optarg << std::endl;
                mismatch_flag = true;
                break;
            case 'g' :
                gap = atoi(optarg);
                std::cerr << "Gap is : " << optarg << std::endl;
                gap_flag = true;
                break;
            case 'a' :
                type = static_cast<brown::AlignmentType>(atoi(optarg));
                std::cerr << "Allignment type is : " << optarg << std::endl;
                type_flag = true;
                break;
            case 'k' :
                kmer_length = atoi(optarg);
                std::cerr << "K-mer length is: " << kmer_length << std::endl; //ovo steka
                break;
            case 'w':
                window_length = atoi(optarg);
                std::cerr << "Window length is: " << window_length << std::endl; //ovo steka
                break;
            case 'f' :
                frequency = atof(optarg);   //ovo steka
                std::cerr << "Top f frequent minimizers that will not be taken in account: " << frequency << std::endl;
                break;
            case 't':
                thread_number = atoi(optarg);
                std::cerr << "Number of threads: " << thread_number << std::endl;
                break;
            case 'c':
                cigar_flag = true;
                std::cerr << "Cigar flag raised" << std::endl;
                break;
            case '?' :
                if (optopt == 'm' || optopt == 'n' || optopt == 'g')
                    std::cerr << "Option -" << optopt << " requires an argument." << std::endl;
                else 
                    std::cerr << "Unknown option." << std::endl;
                exit(EXIT_FAILURE); 
            default:
                exit(EXIT_FAILURE);
        }
    }


    if (optind < argc) {
        std::string file1 = argv[optind++];
        std::string file2 = argv[optind];

        if (file1.compare(file1.length() - 6, 6, ".fasta") == 0 ||
            file1.compare(file1.length() - 3, 3, ".fa") == 0 ||
            file1.compare(file1.length() - 4, 4, ".fna") == 0 ||
            file1.compare(file1.length() - 4, 4, ".ffn") == 0 ||
            file1.compare(file1.length() - 4, 4, ".faa") == 0 ||
            file1.compare(file1.length() - 4, 4, ".frn") == 0) {
                auto p = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(file1);
                referenceGenom = p->Parse(-1);    
        } else {
            std::cerr << "First file needs to be in FASTA format!" << std::endl;
            return 1;
        }
        
        if (file2.compare(file2.length() - 6, 6, ".fastq") == 0 ||
            file2.compare(file2.length() - 3, 3, ".fq") == 0) {
                auto p = bioparser::Parser<Sequence>::Create<bioparser::FastqParser>(file2);
                std::uint32_t chunk_size = 500 * 1024 * 1024;  // 500 MB
                //std::cout << "ude parsirat fragmente" << std::endl;
                for (auto t = p->Parse(chunk_size); !t.empty(); t = p->Parse(chunk_size)) {
                    fragments.insert(
                        fragments.end(),
                        std::make_move_iterator(t.begin()),
                        std::make_move_iterator(t.end()));
                }
        } else if (file2.compare(file2.length() - 6, 6, ".fasta") == 0 ||
            file2.compare(file2.length() - 3, 3, ".fa") == 0 ||
            file2.compare(file2.length() - 4, 4, ".fna") == 0 ||
            file2.compare(file2.length() - 4, 4, ".ffn") == 0 ||
            file2.compare(file2.length() - 4, 4, ".faa") == 0 ||
            file2.compare(file2.length() - 4, 4, ".frn") == 0) {
                auto p = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[optind]);
                fragments = p->Parse(-1);
        } else {
            std::cerr << "Second file needs to be in FASTA or FASTQ format!" << std::endl;
            return 1;
        }

        std::cerr << genomLine;
        std::cerr << "Reference genom name: " << referenceGenom.front()->sequenceName << std::endl;
        std::cerr << "Reference genom length: "<< referenceGenom.front()->sequenceSequence.length() << std::endl;
        
        std::cerr << fragmentLine;
        printsFragmentsStats(fragments);
        
        if (match_flag == false || mismatch_flag == false || gap_flag == false || type_flag == false) {
            std::cerr << "Some arguments for alignment are missing!" << std::endl;
            exit(EXIT_FAILURE);
        }

        std::string* cigar = new std::string();
        unsigned int* target_begin = new unsigned int();
        srand(time(0));
        int fragment_postion1, fragment_postion2;
        do {
            fragment_postion1 = rand() % fragments.size();
        } while (fragments[fragment_postion1]->sequenceSequence.length() < 5000);

        do {
            fragment_postion2 = rand() % fragments.size();
        } while (fragments[fragment_postion2]->sequenceSequence.length() < 5000);

        int result = brown::Align(fragments[fragment_postion1]->sequenceSequence.c_str(), 
                                    fragments[fragment_postion1]->sequenceSequence.length(),
                                    fragments[fragment_postion2]->sequenceSequence.c_str(), 
                                    fragments[fragment_postion2]->sequenceSequence.length(),
                                    type, match, mismatch, gap, cigar, target_begin);        
        std::cerr << std::endl << "Alignment results for two randomly chossen genom fragments in second file: " << std::endl
                    << "    Alignment score: " << result << std::endl
                    << "    Begining of alignemnt is on position " << *target_begin << " of target genom." << std::endl;

        if (cigar_flag)
            std::cerr << "    Cigar string of alignment: " << *cigar << std::endl;
        
        delete cigar;
        delete target_begin;

        
        reference_minimizers = brown::Minimize(referenceGenom[0]->sequenceSequence.c_str(), 
                                                referenceGenom[0]->sequenceSequence.length() , 
                                                kmer_length, window_length);
        
        cutOffTooFrequentMinimizers(reference_minimizers, referenceGenom[0]->sequenceName, true);

        thread_pool::ThreadPool pool(thread_number);

        std::vector<std::future<void>> futures;
        for (int i = 0; i < fragments.size(); i++) {
            futures.emplace_back(pool.Submit(mapping, std::cref(*(fragments[i]))));
        }
        for (auto& it : futures) {
          it.wait();
        }
        
        
    } else {
        throw "Files missing\n";
    }

    return 0;
}
