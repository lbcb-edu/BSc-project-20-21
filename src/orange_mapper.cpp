#include <iostream>
#include <getopt.h>
#include <vector>

#include "bioparser/bioparser/include/bioparser/fasta_parser.hpp"
#include "bioparser/bioparser/include/bioparser/fastq_parser.hpp"
#include "bioparser/bioparser/include/bioparser/parser.hpp"
using namespace std;
#define VERSION "v0.1.0"

static int help_flag;
static int version_flag;
void printHelp(){
    std::cout << "\norange_mapper usage:\n"
    "./orange_mapper [OPTION]\n"
    "./orange_mapper [FILE1][FILE2]\n\n"
    "Options :\n"
    "-h or --help\t Prints this help message\n"
    "-v or --version\t Prints mapper version\n\n"
    "Mapper accepts two files arguments.\n"
    "The first file contains a reference genome in FASTA format.\n"
    "The second file contains a set of fragments in either FASTA or FASTQ format.\n\n"
    "Mapper parses two files and displays statistics for them which are:\n"
    "names of sequences in the reference file and their lengths,\n"
    "number of sequences in the fragments file, their average length,\n"
    "N50 length, minimal and maximal length.\n\n";
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

int main(int argc, char *argv[]){
    int i;
    
    static struct option long_options[] = {
            {"help", no_argument, &help_flag, 1},
            {"version", no_argument, &version_flag, 1},
            {0, 0, 0, 0}    
    };
    
    while((i = getopt_long(argc, argv, "hv", long_options, nullptr)) != -1){
        
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
    }
    return 0;
}
