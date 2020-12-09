#include <iostream>
#include <getopt.h>
#include "projectControl.h"
#include "../bioparser/include/bioparser/fasta_parser.hpp"
#include "../bioparser/include/bioparser/fastq_parser.hpp"
#include "white_alignment.h"

static int version_req;
static int help_req;
static int align_algorithm = 0;
static int match_cost = 1;
static int mismatch_cost = -1;
static int gap_cost = -1;

//gets the project version from projectControl.h
void version_print(){
    std::cout << PROJECT_VER << std::endl;
}

//basic definition of a Sequence structure
struct Sequence { 
	public:
		Sequence(
			const char* name, std::uint32_t nameLength,
			const char* data, std::uint32_t dataLength) {
			this->name = name;
			this->data = data;
			this->dataLength = dataLength;
			this->nameLength = nameLength;
		}

		std::string getName() {
			return this->name;
		}

		std::string getData() {
			return this->data;
		}

		std::uint32_t getNameLength() {
			return this->nameLength;
		}

		std::uint32_t getDataLength() {
			return this->dataLength;
		}

	private:
		const char* name;
		const char* data;
		std::uint32_t nameLength;
		std::uint32_t dataLength;
	};

//Checks if passed arguments are in fasta and fastq formats
bool checkArgs(char *argv[]) {
	const char* extension_fasta[5] = { ".fasta", ".fna", ".ffn", ".faa", ".frn" };
	char* p = strrchr(argv[1], '.');
	int found = 0;
	if (p) {
		for (int i = 0; i < 5; i++) {
			if (strcmp(p, extension_fasta[i])) {
				found++; //if the first file checks out
				break; //break this loop and check the second file
			}
		}
	}

	if (found != 1) { //in case the for loop has come to an end and the file doesn't pass the check
		return false;
	}

	const char* extension_fastq[7] = {".fasta", ".fna", ".ffn", ".faa", ".frn", ".fastq", ".fq"};
	p = strrchr(argv[2], '.');
	if (p) {
		for (int i = 0; i < 7; i++) {
			if (strcmp(p, extension_fastq[i])) {
				return true; //if the second file checks out, then both files are good and the function returns true
			}
		}
	}
	return false; //if the second file doesn't check out, the function returns false
}

//calculates the index of N50-th member
int calculateN50(std::vector<size_t> fragmentVector, int sum) {
	int temp_sum = 0;
	for (int i = fragmentVector.size(); i > 0; i--) {
		temp_sum += fragmentVector[i];
		if (temp_sum > sum/2) {
			return i;
		}
	}
	return -1; //the program should never get to this, but it was necessary to implement because of compilation issues
}
	
//prints the standard help message
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

int calcAlignment(int size, const std::vector<std::unique_ptr<Sequence>> &fragment_list, std::string* cigar, unsigned int* target_begin) {
	std::cout << "check0\n";
	srand (time(NULL));
	int query_index;
	int target_index;
	white::AlignmentType align_type;
	std::cout << "check0.25\n";
	do {
		query_index = rand() % (size);
		target_index = rand() % (size);
	} while (fragment_list[query_index] -> getDataLength() > 5000 &&
		fragment_list[target_index] -> getDataLength() > 5000);
	std::cout << "check0.5\n";
	switch (align_algorithm) {
	case 0:
		align_type = white::AlignmentType::GLOBAL;
		break;
	case 1:
		align_type = white::AlignmentType::LOCAL;
		break;
	case 2:
		align_type = white::AlignmentType::SEMIGLOBAL;
		break;
	default:
		break;
	}
	std::cout << "check1\n";
	white::Aligner aligner = white::Aligner(fragment_list[query_index] -> getData().c_str(),
		fragment_list[query_index] -> getDataLength(),
		fragment_list[target_index] -> getData().c_str(),
		fragment_list[target_index] -> getDataLength(),
		match_cost, mismatch_cost, gap_cost, cigar, target_begin);
	
		int align_score = aligner.Align (align_type);
	return align_score;
};

int main(int argc, char *argv[]){

    int opt;

	int align_algorithm = 0;
	int match_cost = 1;
	int mismatch_cost = -1;
	int gap_cost = -1;
    
    static struct option long_options[] = {
            {"help", no_argument, &help_req, 1},
            {"version", no_argument, &version_req, 1},
			{"algorithm", required_argument, nullptr, 'a'},
			{"match_cost", required_argument, nullptr, 'm'},
			{"mismatch_cost", required_argument, nullptr, 'n'},
			{"gap_cost", required_argument, nullptr, 'g'},
            {0, 0, 0, 0}    
    };

	
    
    while ((opt = getopt_long(argc, argv, "a:m:n:g:hv", long_options, nullptr)) != -1){
        
        switch (opt)
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

		case 'a':
			align_algorithm = atoi(optarg);
			break;

		case 'm':
			match_cost = atoi(optarg);
			break;

		case 'n':
			mismatch_cost = atoi(optarg);
			break;

		case 'g':   
			gap_cost = atoi(optarg);
			break;

        default:
            abort();
        }
    }

	//prints help or version and ends the program, depending on which option flag was input as an argument
    if (version_req){
        version_print();
		return 0;
    } else if(help_req){
        help_print();
		return 0;
	}

	//if the input wasn't help or version flag, then it checks if we have given two arguments
	if (argc != 3) {
		std::cerr << "error: invalid input, please include exactly two files\n";
		return 1;
	}

	//if the given two arguments aren't in valid formats, the program recieves an error
	if (!checkArgs(argv)) {
		std::cerr << "error: invalid file format, please pass two files in FASTA and FASTQ formats in that order\n";
		return 1;
	}

	

	//parsing the first file
	auto genome = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[1]);
	auto g = genome->Parse(-1);
	int g_size = (int)g.size();

	//parsing the second file
	auto fragments = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[2]);
	auto f = fragments->Parse(-1);
	int f_size = (int)f.size();
	
	int sum = 0;

	//names of sequences in the reference file and their lengths
	/* for (int i = 0; i < g_size; i++) {
		std::cerr << "Name of sequence: " << g[i]->getName() << "\n" << "Length of sequence: " << g[i]->getDataLength() << "\n\n";
	} */

	//number of sequences in the fragments file
	for (int i = 0; i < f_size; i++) {
		sum += f[i]->getDataLength();
	}
	std::cerr << "Number of sequences: " << sum << "\n\n";

	float avg_size = sum / f_size;
	std::cerr << "Average length of fragments: " << avg_size << "\n\n";


	//filling a vector with necessary data
	std::vector<size_t> fragmentVector;
	for (int i = 0; i < f_size; i++) {
		fragmentVector.push_back(f[i]->getDataLength());
	}

	std::sort(fragmentVector.begin(), fragmentVector.end()); //sorting a vector

	//calculates at which place is N50 located, and outputs it
	int N50 = calculateN50(fragmentVector, sum);
	std::cerr << "N50 length: " << fragmentVector[N50] << "\n\n";

	//prints the minimum and the maximum value from the fragment vector
	std::cerr << "Minimum: " << fragmentVector.front() << "\n\n";
	std::cerr << "Maximum: " << fragmentVector.back() << "\n\n";

	std::string *cigar;
	unsigned int *target_begin;
	int alignment_score = calcAlignment(f_size, f, cigar, target_begin);

	 std::cout << "Alignment score: " << alignment_score << '\n'
            << "Target begin index: " << target_begin << "\n\n"
            << "CIGAR: " << cigar << std::endl;

	return 0;
}