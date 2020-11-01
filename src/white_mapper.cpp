#include <iostream>
#include <getopt.h>
#include "projectControl.h"
#include "../bioparser/include/bioparser/fasta_parser.hpp"
#include "../bioparser/include/bioparser/fastq_parser.hpp"

static int version_req;
static int help_req;

void version_print(){
    std::cout << PROJECT_VER << std::endl;
}

//Checks if passed arguments are in fasta and fastq formats
bool checkArgs() {
	std::string extension_fasta[5] = [".fasta", ".fna", ".ffn", ".faa", ".frn"];
	char* p = strrchr(argv[0], '.');
	int found = 0;
	if (p) {
		for (int i = 0; i < 5; i++) {
			if (strcmp(p, extension[i]) {
				found++;
					break;
			}
		}
	}

	if (found != 1) {
		return false;
	}

	std::string extension_fastq[7] = [".fasta", ".fna", ".ffn", ".faa", ".frn", ".fastq", ".fq"];
	p = strrchr(argv[1], '.');
	if (p) {
		for (int i = 0; i < 7; i++) {
			if (strcmp(p, extension[i]) {
				return true;
			}
		}
	}
	return false;
}
	

void help_print(){
    std::cout << "\n" PROJECT_NAME " usage:\n"
    "Two options :\n"
    "-h or --help\t Prints help message\n"
    "-v or --version\t Prints version\n\n"
    "Mapper accepts two files arguments.\n"
    "The first file  FASTA format.\n"
    "The second file  FASTA or FASTQ format.\n\n"
    "Mapper shows these statistics :\n"
    "names of sequences in the reference file and their lengths,\n"
    "number of sequences in the fragments file,\n"
    "their average length,\n"
    "N50 length,\n"
    "minimal and maximal length.\n\n";
}

int main(int argc, char *argv[]){

    int i;
    
    static struct option long_options[] = {
            {"help", no_argument, &help_req, 1},
            {"version", no_argument, &version_req, 1},
            {0, 0, 0, 0}    
    };
    
    while ((i = getopt_long(argc, argv, "hv", long_options, nullptr)) != -1){
        
        switch (i)
        {
        case 0:
            break;

        case 'h':
            help_req = true;
            break;

        case 'v':
            version_req = true;
            break;

        case '?':
            break;    
        
        default:
            abort();
        }
    }

    if (version_req){
        version_print();
    } else if(help_req){
        help_print();
    }

	if (argc != 2) {
		std::cerr << "error: invalid input, please include exactly two files";
		return 1;
	}

	if (!checkArgs()) {
		std::cerr << "error: invalid file format, please pass two files in FASTA and FASTQ formats in that order"
			return 1;
	}

	struct Sequence {  // or any other name
	public:
		Sequence(  // required arguments
			const char*, std::uint32_t,
			const char*, std::uint32_t) {
			// implementation
		}
	}
	auto p = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(path);

	// parse whole file
	auto s = p->Parse(-1);
	return 0;
}