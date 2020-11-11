#include <iostream>
#include <getopt.h>
#include "bioparser/include/bioparser/parser.hpp"
#include "bioparser/include/bioparser/fastq_parser.hpp"
#include "bioparser/include/bioparser/fasta_parser.hpp"

#define VERSION "v0.1"

static option long_options[] = {
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'v'},
    {0, 0, 0, 0}
};

static std::string help = "brown_mapper \n"
                            "-h or --help for help\n"
                            "-v or --version for programs version\n"
                            "Program accepts two files as floating arguments.\n"
                            "The first file needs to contain a reference genome in FASTA format.\n"
                            "The second file needs to contain a set of fragments\n"
                            "in either FASTA or FASTQ format.";

/*fastq
static struct Sequence {  // or any other name
 public:
  Sequence(  // required arguments
      const char*, std::uint32_t,
      const char*, std::uint32_t,
      const char*, std::uint32_t) {
    // implementation
  }
}
auto p = bioparser::Parser<Sequence>::Create<bioparser::FastqParser>(path);

// parse in chunks
std::vector<std::unique_ptr<Sequence>> s;
std::uint32_t chunk_size = 500 * 1024 * 1024;  // 500 MB
for (auto t = p->parse(chunk_size); !t.empty(); t = p->parse(chunk_size)) {
  s.insert(
      s.end(),
      std::make_move_iterator(t.begin()),
      std::make_move_iterator(t.end()));
}
*/

//fasta
struct Sequence {  // or any other name
  public:
    std::string sequenceName;
    std::string sequenceSequence;   
    std::string sequenceQuality;

    Sequence(  // required arguments
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
//auto p = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(path);

// parse whole file
//auto s = p->Parse(-1);




int main(int argc, char* argv[]) {

    //std::cout << "evo me1" << std::endl;

    int c = getopt_long(argc, argv, "hv", long_options, 0);
    if (c != -1) {
        switch (c){
            case 'h' :
                std::cout << help << std::endl;
                return 0;;
            case 'v' :
                std::cout << VERSION << std::endl;
                return 0;
            default:
                return 0;
        }
    }
    //std:: cout << "evo me 2" << std::endl;
    if (argc == 2) { //promjeni u 3
        //std::cout << "ide pravit stringove";
        std::string file1 = argv[1];
        //std::string file2 = argv[2];
        if (file1.compare(file1.length() - 6, 6, ".fasta") == 0 ||
            file1.compare(file1.length() - 3, 3, ".fa") == 0) {
                auto p = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[1]);
                std::vector<std::unique_ptr<Sequence>> s = p->Parse(-1);
                std::cout << s.front()->sequenceName << std::endl;
                std::cout << s.front()->sequenceSequence.length() << std::endl;
                
        } else {
            std::cerr << "First file needs to be in FASTA format!" << std::endl;
        }

        /*if ((file2.compare(file2.length() - 6, 6, ".fastq") == 0 ||
            file2.compare(file2.length() - 3, 3, ".fq") == 0)) {
                auto p = bioparser::Parser<Sequence>::Create<bioparser::FastqParser>(argv[2]);
                std::vector<std::unique_ptr<Sequence>> s;
                std::uint32_t chunk_size = 500 * 1024 * 1024;  // 500 MB
                for (auto t = p->Parse(chunk_size); !t.empty(); t = p->Parse(chunk_size)) {
                     s.insert(
                        s.end(),
                        std::make_move_iterator(t.begin()),
                        std::make_move_iterator(t.end()));
                }
        }*/

    }

    return 0;
}
