#include <iostream>
#include <getopt.h>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "orange_alignment/orange_alignment.h"
#include "orange_minimizer/orange_minimizer.h"

using namespace std;
#define VERSION "v0.1.0"

static int help_flag = 0;
static int version_flag = 0;
static int alg_flag = 0;
static int match_flag = 1;
static int mismatch_flag = -1;
static int gap_flag = -1;
static int kmer_len = 15;
static int window_len = 5;
static double frequency = 0.001;
//std::unordered_map<Kmer, unsigned int> minimizer_count;

void printHelp(){
    std::cout << "\norange_mapper usage:\n"
    "./orange_mapper [OPTION]\n"
    "./orange_mapper [FILE1][FILE2]\n\n"
    "Options :\n"
    "-h or --help\t Prints this help message\n"
    "-v or --version\t Prints mapper version\n"
    "-a or --algorithm\t select alignment type: 0 for global(default), 1 for local, 2 for semi-global alignment\n"
    "-m or --match\t Define match cost\n"
    "-n or --mismatch\t Define mismatch cost\n"
    "-g or --gap\t Define gap cost\n\n"
    "Mapper accepts two files arguments.\n"
    "The first file contains a reference genome in FASTA format.\n"
    "The second file contains a set of fragments in either FASTA or FASTQ format.\n\n"
    "Mapper parses two files and displays statistics for them which are:\n"
    "names of sequences in the reference file and their lengths,\n"
    "number of sequences in the fragments file, their average length,\n"
    "N50 length, minimal and maximal length, alignment score.\n\n";
}

void printVersion(){
    std::cout << VERSION << std::endl;
}

class Sequence {
public:
    Sequence(  // required arguments
        const char* name, std::uint32_t name_len,
        const char* data, std::uint32_t data_len):
        names(name, name_len), 
        datas(data, data_len)
    {
    }
    string names, datas;
};

void get_minimizers(const vector<unique_ptr<Sequence>> &fragments, unsigned int kmer_len, unsigned int window_len, double f){

	unsigned int sequence_len;
	vector <tuple<unsigned int, unsigned int, bool>> distinct;
	vector <tuple<unsigned int, unsigned int, bool>> all;
	char *sequence;
        for (int i = 0; i < fragments.size(); i++) {
            strcpy(sequence, (fragments[i]->datas).c_str());
        }
	orange::Minimizer *m = new orange::Minimizer(sequence, sequence_len, kmer_len, window_len);
	distinct = m->Minimize(sequence, sequence_len, kmer_len, window_len);
	all = m->all_minimizers(sequence, sequence_len, kmer_len, window_len);
	cout << "Number of distinct minimizers: " << distinct.size();
	cout << "Fraction of singletons: "<< distinct.size()/all.size();
	return;
}

int calculate_score(int algorithm, int gap, int match, int mismatch, const vector<unique_ptr<Sequence>> &fragments){
	int index1, index2;
	string cigar;
    	unsigned int target_begin;

	srand((unsigned)time(NULL));
	while(1){

		index1 = rand()%(fragments.size()) + 0;
		if(fragments[index1]->datas.length() < 5000){

			break;
			
		}
	}
	while(1){

		index2 = rand()%(fragments.size()) + 0;
		if(fragments[index2]->datas.length() < 5000){

			break;
		}
	}
	int query_len = fragments[index1]->datas.length();
        char query[query_len]; 
    	strcpy(query, fragments[index1]->datas.c_str());
    	int target_len = fragments[index2]->datas.length();
    	char target[target_len];
    	strcpy(query, fragments[index2]->datas.c_str());
	orange::AlignmentType type;
	if(algorithm == 0){
		type = orange::AlignmentType::global;
	}else if(algorithm == 1){
		type = orange::AlignmentType::local;
	}
	else if(algorithm == 2){
		type = orange::AlignmentType::semiGlobal;
	}
    orange::Alignment* a = new orange::Alignment(query, sizeof(query), target, sizeof(target), type, match, mismatch, gap, &cigar, &target_begin);
    return a->Align(query, query_len, target, target_len, type, match, mismatch, gap, &cigar, &target_begin);
}

int main(int argc, char *argv[]){
    int i;
    
   static struct option long_options[] = {
            {"help", no_argument, &help_flag, 1},
            {"version", no_argument, &version_flag, 1},
            {"algorithm", required_argument, &alg_flag, 1},
            {"match", required_argument, &match_flag, 1},
            {"mismatch", required_argument, &mismatch_flag, 1},
            {"gap", required_argument, &gap_flag, 1},
            {"kmer", required_argument, &kmer_len, 1},
            {"window", required_argument, &window_len, 1},
            {0, 0, 0, 0}    
    };
    
       while((i = getopt_long(argc, argv, "hva:m:n:g", long_options, nullptr)) != -1){
        
        switch (i)
        {
        case 0:
            break;

        case 'h':
            help_flag = 1;
            break;

        case 'v':
            version_flag = 1;
            break;
        case 'a':
            alg_flag = stoi(optarg);
            break;
        case 'm':
            match_flag = stoi(optarg);
            break;
        case 'n':
            mismatch_flag = stoi(optarg);
            break;
        case 'g':
            gap_flag = stoi(optarg);
            break;
        case 'k':
            kmer_len = stoi(optarg);
            break;
        case 'w':
            window_len = stoi(optarg);
            break;
       	case 'f':
            frequency = stoi(optarg);
            break;    
        case '?':
            break;    
        default:
            abort();
        }
    }
    if(help_flag){
        printHelp();
    } else if(version_flag){
        printVersion();
    }

    if(optind < argc){
        auto genome = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[1]);
        auto genome_parsed = genome->Parse(-1);

        auto fragments = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[2]);
        auto fragments_parsed = fragments->Parse(-1);

        //print names of genome from reference file and their lengths
        cerr<<"Sequnces in the reference genome file: \n";
        for(int i = 0; i < genome_parsed.size(); i++){
            cerr<< genome_parsed[i]-> names << " , length= "<<genome_parsed[i]->datas.length()<<"\n";
        }

        //number of sequences in the fragments file, average length
        cerr<<"Number of sequences in fragments: "<<fragments_parsed.size()<<"\n";
        uint64_t sumFragmentsLen = 0;
        for(int i = 0; i < fragments_parsed.size();i++){
            sumFragmentsLen += fragments_parsed[i]->datas.length();
        }
        uint64_t avgLen = sumFragmentsLen/fragments_parsed.size();
        cerr<<"Average length: "<< avgLen<<"\n";
        //N50 length
        vector<size_t> len(fragments_parsed.size());
        for(int i = 0 ; i < fragments_parsed.size(); i++){
            len[i] = fragments_parsed[i]->datas.length();
        }
        //sorting in descending order(biggest to smallest)
        sort(len.begin(), len.end(), greater<size_t>());
        uint64_t n50Len = 0;
        int i = 0;
        do{
            n50Len += len[i];
            i++;
        }while(n50Len < (sumFragmentsLen/2));
        cerr<<"N50 length = "<<n50Len<<"\n";

        //minimal length - last element
        cerr<<"Minimal length = "<<len.back()<<"\n";

        //maximal length - first element
        cerr<<"Maximal length = "<<len.front()<<"\n";
         std::vector<std::unique_ptr<Sequence>> fdok;
            std::uint32_t chunk_size = 500 * 1024 * 1024;  // 500 MB
            for (auto t = fragments->Parse(chunk_size); !t.empty(); t = fragments->Parse(chunk_size)) {
                fdok.insert(
                    fdok.end(),
                    std::make_move_iterator(t.begin()),
                    std::make_move_iterator(t.end()));
            }
        //alignment score 
        int score = calculate_score(alg_flag, gap_flag, match_flag, mismatch_flag, fragments_parsed);
       // cerr<<"Alignment score = "<<score<<"\n";

    	get_minimizers(fdok, kmer_len, window_len, frequency);
	
    }
    return 0;
}
