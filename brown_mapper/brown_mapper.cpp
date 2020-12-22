#include <iostream>
#include <getopt.h>
#include "bioparser/include/bioparser/parser.hpp"
#include "bioparser/include/bioparser/fastq_parser.hpp"
#include "bioparser/include/bioparser/fasta_parser.hpp"
#include "brown_alignment.hpp"
#include "brown_minimizer.hpp"
#include <stdlib.h>
#include <time.h>
#include <map>

#define VERSION "v0.5"
#define DEFAULT_KMER_LENGTH 15
#define DEFAULT_WINDOW_LENGTH 5
#define DEFAULT_MINIMIZER_FREQUENCY 0.001

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
                            "-g   value for gap\n\n"
                            "Minimizizer arguments:\n"
                            "-k   k-mer length (default: " + std::to_string(DEFAULT_KMER_LENGTH) + ")\n"
                            "-w   window length (default: " + std::to_string(DEFAULT_WINDOW_LENGTH) + ")\n"
                            "-f   top f frequent minimizers that will not be taken in account (default:" 
                            + std::to_string(DEFAULT_MINIMIZER_FREQUENCY) + ")\n\n"
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



int main(int argc, char* argv[]) {

    brown::AlignmentType type;
    int match, gap, mismatch;
    int kmer_length = DEFAULT_KMER_LENGTH, window_length = DEFAULT_WINDOW_LENGTH;
    double frequency = DEFAULT_MINIMIZER_FREQUENCY; 

    bool match_flag;
    bool mismatch_flag;
    bool gap_flag;
    bool type_flag;

    int c;
    while ((c = getopt_long(argc, argv, "m:g:n:a:k::w::f::hv", long_options, 0)) != -1) {
        switch (c){
            case 'h' :
                std::cerr << help << std::endl;
                return 0;
            case 'v' :
                std::cerr << VERSION << std::endl;
                return 0;
            case 'm':
                match = atoi(optarg);
                std::cout << "Match is : " << optarg << std::endl;
                match_flag = true;
                break;
            case 'n' :
                mismatch = atoi(optarg);
                std::cout << "Mismatch is : " << optarg << std::endl;
                mismatch_flag = true;
                break;
            case 'g' :
                gap = atoi(optarg);
                std::cout << "Gap is : " << optarg << std::endl;
                gap_flag = true;
                break;
            case 'a' :
                type = static_cast<brown::AlignmentType>(atoi(optarg));
                std::cout << "Allignment type is : " << optarg << std::endl;
                type_flag = true;
                break;
            case 'k' :
                kmer_length = atoi(optarg);
                std::cout << "K-mer length is: " << kmer_length << std::endl; //ovo steka
                break;
            case 'w':
                window_length = atoi(optarg);
                std::cout << "Window length is: " << window_length << std::endl; //ovo steka
                break;
            case 'f' :
                frequency = atof(optarg);   //ovo steka
                std::cout << "Top f frequent minimizers that will not be taken in account: " << frequency << std::endl;
                break;
            case '?' :
                if (optopt == 'm' || optopt == 'n' || optopt == 'g')
                    std::cout << "Option -" << optopt << " requires an argument." << std::endl;
                else 
                    std::cout << "Unknown option." << std::endl;
                exit(EXIT_FAILURE); 
            default:
                exit(EXIT_FAILURE);
        }
    }


    //std:: cout << "evo me 2" << std::endl;
    if (optind < argc) {
        //std::cout << "ide pravit stringove";
        std::string file1 = argv[optind++];
        std::string file2 = argv[optind];

        std::vector<std::unique_ptr<Sequence>> referenceGenom;
        if (file1.compare(file1.length() - 6, 6, ".fasta") == 0 ||
            file1.compare(file1.length() - 3, 3, ".fa") == 0 ||
            file1.compare(file1.length() - 4, 4, ".fna") == 0 ||
            file1.compare(file1.length() - 4, 4, ".ffn") == 0 ||
            file1.compare(file1.length() - 4, 4, ".faa") == 0 ||
            file1.compare(file1.length() - 4, 4, ".frn") == 0) {
                auto p = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[1]);
                referenceGenom = p->Parse(-1);                
        } else {
            std::cerr << "First file needs to be in FASTA format!" << std::endl;
            return 1;
        }
        
        std::vector<std::unique_ptr<Sequence>> fragments;
        if ((file2.compare(file2.length() - 6, 6, ".fastq") == 0 ||
            file2.compare(file2.length() - 3, 3, ".fq") == 0)) {
                auto p = bioparser::Parser<Sequence>::Create<bioparser::FastqParser>(argv[2]);
                std::uint32_t chunk_size = 500 * 1024 * 1024;  // 500 MB
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
                auto p = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[2]);
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
            std::cout << "Some arguments for alignment are missing!" << std::endl;
            exit(EXIT_FAILURE);
        }

        std::string cigar;
        unsigned int target_begin;
        srand(time(0));
        int fragment_postion1, fragment_postion2;
        do {
            fragment_postion1 = rand() % fragments.size();
        } while (referenceGenom[fragment_postion1]->sequenceSequence.length() < 5000);

        do {
            fragment_postion2 = rand() % fragments.size();
        } while (referenceGenom[fragment_postion2]->sequenceSequence.length() < 5000);

        /*int result = brown::Align(fragments[fragment_postion1]->sequenceSequence.c_str(), 
                                    fragments[fragment_postion1]->sequenceSequence.length(),
                                    fragments[fragment_postion2]->sequenceSequence.c_str(), 
                                    fragments[fragment_postion2]->sequenceSequence.length(),
                                    type, match, mismatch, gap, &cigar, &target_begin);*/

        char sequence[referenceGenom[0]->sequenceSequence.length()+1];
        strcpy(sequence, referenceGenom[0]->sequenceSequence.c_str());
        std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers=
        brown::Minimize(sequence, referenceGenom[0]->sequenceSequence.length()+1 , kmer_length, window_length);
        int id=1; // na kraju ce biti broj razlicitih minimizera
        int size=minimizers.size();
        std::map<int, int> mapaMinimizera; //prvo minimizer, pa koliko se pojavljuje
        for(int i=0;i<size;i++) {
            int currentMinimizer=std::get<0>(minimizers[i]);
            if(mapaMinimizera.count(currentMinimizer) == 0) {
                mapaMinimizera.insert(std::pair<int,int>(currentMinimizer, 1));
            }
            else {
                mapaMinimizera[currentMinimizer]++;
            }
        }
        int numOfMinimizers=mapaMinimizera.size();
        int numOfSingletons=0;
        int fthMinimizer;
        std::map<int, int>::iterator itr;
        std::multimap<int, int> mapByFrequency;
        for(itr = mapaMinimizera.begin(); itr != mapaMinimizera.end(); itr++) {
            if(itr->second == 1) numOfSingletons++;
            mapByFrequency.insert(std::pair<int, int>(itr->second, itr->first));
        }
        double fractionOfSingletons=(double) numOfSingletons/size;
        size=mapByFrequency.size();
        int counter=0;
        for(itr = mapByFrequency.begin(); itr != mapByFrequency.end(); itr++) {
            if (counter == size - frequency -1) {
                fthMinimizer = itr->second;
                break;
            }
        }
        std::cout << numOfMinimizers << " " << fractionOfSingletons << " " << fthMinimizer << std::endl;
        
        
        
    }
    
    


    return 0;
}
