#include <iostream>
#include <getopt.h>
#include <unordered_map>
#include "projectControl.h"
#include "../bioparser/include/bioparser/fasta_parser.hpp"
#include "../bioparser/include/bioparser/fastq_parser.hpp"
#include "white_minimizers.hpp"

#include "white_alignment.hpp"

static int version_req;
static int help_req;
static int align_algorithm = 0;
static int match_cost = 1;
static int mismatch_cost = -1;
static int gap_cost = -1;
static double ignored_fraction = 0.001;
static bool calc_cigar = false;
static unsigned int kmer_len = 15;
static unsigned int window_len = 5;
static int thread_count = 1;

//basic definition of a Sequence structure
struct Sequence
{
public:
	Sequence(
		const char *name, std::uint32_t nameLength,
		const char *data, std::uint32_t dataLength,
		const char *quality = nullptr, std::uint32_t qualityLength = 0) : name(name, nameLength),
																		  data(data, dataLength) {}

	std::string getName()
	{
		return this->name;
	}

	std::string getData()
	{
		return this->data;
	}

private:
	std::string name;
	std::string data;
	std::string quality;
};

//Checks if passed arguments are in fasta and fastq formats
bool checkArgs(char *argv[])
{
	const char *extension_fasta[5] = {".fasta", ".fna", ".ffn", ".faa", ".frn"};
	char *p = strrchr(argv[1], '.');
	int found = 0;
	if (p)
	{
		for (int i = 0; i < 5; i++)
		{
			if (strcmp(p, extension_fasta[i]))
			{
				found++; //if the first file checks out
				break;	 //break this loop and check the second file
			}
		}
	}

	if (found != 1)
	{ //in case the for loop has come to an end and the file doesn't pass the check
		return false;
	}

	const char *extension_fastq[7] = {".fasta", ".fna", ".ffn", ".faa", ".frn", ".fastq", ".fq"};
	p = strrchr(argv[2], '.');
	if (p)
	{
		for (int i = 0; i < 7; i++)
		{
			if (strcmp(p, extension_fastq[i]))
			{
				return true; //if the second file checks out, then both files are good and the function returns true
			}
		}
	}
	return false; //if the second file doesn't check out, the function returns false
}

//calculates the index of N50-th member
int calculateN50(std::vector<size_t> fragmentVector, int sum)
{
	int temp_sum = 0;
	for (int i = fragmentVector.size(); i > 0; i--)
	{
		temp_sum += fragmentVector[i];
		if (temp_sum > sum / 2)
		{
			return i;
		}
	}
	return -1; //the program should never get to this, but it was necessary to implement because of compilation issues
}

//prints the standard help message
void help_print()
{
	std::cout << "\n" PROJECT_NAME " usage: white_mapper [options] <reference-genome> <fragments>\n\n"
				 "\t<reference-genome>\n\t reference genome in FASTA format\n"
				 "\t<fragments>\n\t set of fragments in FASTA or FASTQ format\n\n"
				 "	options:\n"
				 "		-h, --help\t Prints help message\n"
				 "		-v, --version\t Prints version\n"
				 "	---------------------------------------"
				 "		-a, --algorithm <int>\n"
				 "			alignment algorithm:\n"
				 "			 0 - LOCAL/Smith-Waterman (default)\n"
				 "			 1 - GLOBAL/Needleman-Wunsch\n"
				 "			 2 - SEMI-GLOBAL\n"
				 "		-m, --match_cost <int>\n"
				 "			sets match cost during alignment\n"
				 "			default: 1\n"
				 "		-n, --mismatch_cost <int>\n"
				 "			sets mismatch cost during alignment\n"
				 "			default: -1\n"
				 "		-g, --gap_cost <int>\n"
				 "			sets linear gap cost during alignment\n"
				 " 			default: -1\n"
				 "		-k, --kmer_len <int>\n"
				 "			sets kmer length for minimizers\n"
				 "			default: 15\n"
				 "		-w, --window_len <int>\n"
				 "			sets window length for minimizers\n"
				 "			default: 5\n"
				 "		-f, --ingored_fraction <double> in range [0, 1)\n"
				 "			sets the percentage of top minimizers to ignore\n"
				 "			default: 0.001 (0.1%)\n"
				 "		-c, --cigar\n"
				 "			mapper calculates and prints CIGAR string\n"
				 "		-t, --thread <int>\n"
				 "			sets the number of threads\n"
				 "			default: 1"
				 "Mapper shows following statistics :\n"
				 "names of sequences in the reference file and their lengths,\n"
				 "number of sequences in the fragments file,\n"
				 "their average length,\n"
				 "N50 length,\n"
				 "minimal and maximal length.\n\n";
}

//gets the project version from projectControl.h
void version_print()
{
	std::cout << PROJECT_VER << std::endl;
}
int calcAlignment(int size, const std::vector<std::unique_ptr<Sequence>> &fragment_list, std::string *cigar, unsigned int *target_begin)
{
	srand(time(NULL));
	int query_index;
	int target_index;
	white::AlignmentType align_type;
	;
	do
	{
		query_index = rand() % (size);
		target_index = rand() % (size);
	} while (fragment_list[query_index]->getData().size() > 5000 &&
			 fragment_list[target_index]->getData().size() > 5000);
	switch (align_algorithm)
	{
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

	white::Aligner *aligner = new white::Aligner(fragment_list[query_index]->getData().c_str(),
												 fragment_list[query_index]->getData().size(),
												 fragment_list[target_index]->getData().c_str(),
												 fragment_list[target_index]->getData().size(),
												 match_cost, mismatch_cost, gap_cost, cigar, target_begin);

	int align_score = aligner->Align(align_type);
	return align_score;
};

int main(int argc, char *argv[])
{

	int opt;

	static struct option long_options[] = {
		{"help", no_argument, &help_req, 1},
		{"version", no_argument, &version_req, 1},
		{"algorithm", required_argument, nullptr, 'a'},
		{"match_cost", required_argument, nullptr, 'm'},
		{"mismatch_cost", required_argument, nullptr, 'n'},
		{"gap_cost", required_argument, nullptr, 'g'},
		{"ignored_fraction", required_argument, nullptr, 'f'},
		{"kmer_len", required_argument, nullptr, 'k'},
		{"window_len", required_argument, nullptr, 'w'},
		{"cigar", required_argument, nullptr, 'c'},
		{"thread", required_argument, nullptr, 't'},
		{0, 0, 0, 0}};

	while ((opt = getopt_long(argc, argv, "a:m:n:g:k:f:w:hv:", long_options, nullptr)) != -1)
	{

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

		case 'f':
			ignored_fraction = atof(optarg);
			break;

		case 'k':
			kmer_len = atoi(optarg);
			break;

		case 'w':
			window_len = atoi(optarg);
			break;

		case 'c':
			calc_cigar = true;
			break;

		case 't':
			thread_count = atoi(optarg);
			break;

		default:
			abort();
		}
	}

	//prints help or version and ends the program, depending on which option flag was input as an argument
	if (version_req)
	{
		version_print();
		return 0;
	}
	else if (help_req)
	{
		help_print();
		return 0;
	}

	//if the input wasn't help or version flag, then it checks if we have given two arguments
	if (argc != 3)
	{
		std::cerr << "error: invalid input, please include exactly two files\n";
		return 1;
	}

	//if the given two arguments aren't in valid formats, the program recieves an error
	if (!checkArgs(argv))
	{
		std::cerr << "error: invalid file format, please pass two files in FASTA and FASTQ formats in that order\n";
		return 1;
	}
	
	//parsing the first file
	auto ref = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[argc - 2]);
	auto reference = ref->Parse(-1);
	int reference_size = (int)reference.size();

	//parsing the second file
	auto frag = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[argc - 1]);
	auto fragments = frag->Parse(-1);
	int fragments_size = (int)fragments.size();

	//names of sequences in the reference file and their lengths
	for (int i = 0; i < reference_size; i++)
	{
		std::cerr << "Name of sequence: " << reference[i]->getName() << "\n"
				  << "Length of sequence: " << reference[i]->getData().size() << "\n\n";
	}

	//number of sequences in the fragments file
	int sum = 0;
	for (int i = 0; i < fragments_size; i++)
	{
		sum += fragments[i]->getData().size();
	}
	float avg_size = sum / fragments_size;

	std::vector<size_t> fragmentVector;
	for (int i = 0; i < fragments_size; i++)
	{
		fragmentVector.push_back(fragments[i]->getData().size());
	}
	std::sort(fragmentVector.begin(), fragmentVector.end());

	int N50 = calculateN50(fragmentVector, sum);

	std::cerr << "Number of sequences: " << fragments_size << "\n";
	std::cerr << "Total length of all fragments: " << sum << "\n";
	std::cerr << "Minimum length: " << fragmentVector.front() << "\n";
	std::cerr << "Maximum length: " << fragmentVector.back() << "\n";
	std::cerr << "Average length of fragments: " << avg_size << "\n";
	std::cerr << "N50 length: " << fragmentVector[N50] << "\n\n";

	std::string cigar;
	unsigned int target_begin;
	int alignment_score = calcAlignment(fragments_size, fragments, &cigar, &target_begin);
	std::cout << "Alignment score: " << alignment_score << std::endl;
	std::cout << "Target begin index: " << target_begin << std::endl;
	if (calc_cigar)
		std::cout << "Cigar: " << cigar << std::endl;

	std::unordered_map<unsigned int, unsigned int> minimizers_map;
	int singletons_count = 0;
	for (auto &seq : reference)
	{
		std::cout << "test"
				  << "\n\n";
		auto minimizer_vector =
			white::Minimize(seq->getData().c_str(), seq->getData().size(), kmer_len, window_len);

		for (auto &minimizer : minimizer_vector)
		{
			minimizers_map[std::get<0>(minimizer)]++;
		}
	}

	std::vector<unsigned int> minimizer_sorted(minimizers_map.size());
	for (auto &minimizer : minimizers_map)
	{
		if (minimizer.second == 1)
		{
			singletons_count++;
		}
		minimizer_sorted.push_back(minimizer.second);
	}

	std::sort(minimizer_sorted.begin(), minimizer_sorted.end(), std::greater<>());

	std::cout << "Number of distinct minimizers: " << minimizers_map.size() << std::endl;
	std::cout << "Fraction of singletons: " << (double)singletons_count / minimizers_map.size() << std::endl;
	std::cout << "Number of occurrences of the most frequent minimizer with top " 
		<< ignored_fraction * 100 << "% most frequent ignored: " 
		<< minimizer_sorted[std::ceil(minimizers_map.size() * ignored_fraction)] 
		<< std::endl;

	return 0;
}