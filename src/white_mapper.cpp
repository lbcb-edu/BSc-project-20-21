#include <iostream>
#include <getopt.h>
#include <unordered_map>

#include "projectControl.h"

#include "../bioparser/include/bioparser/fasta_parser.hpp"
#include "../bioparser/include/bioparser/fastq_parser.hpp"
#include "../thread_pool/include/thread_pool/thread_pool.hpp"

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
static std::string separator = "----------------------------------------\n";

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
bool checkArgs(char *argv[], int argc)
{
	const char *extension_fasta[5] = {".fasta", ".fna", ".ffn", ".faa", ".frn"};
	char *p = strrchr(argv[argc - 2], '.');
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
	p = strrchr(argv[argc - 1], '.');
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
				 "		-v, --version\t Prints version\n		"
			  << separator << "		-a, --algorithm <int>\n"
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
							  "			default: 1\n\n"
							  "Mapper shows following statistics :\n"
							  "names of sequences in the reference file and their lengths,\n"
							  "number of sequences in the fragments file,\n"
							  "their average length,\n"
							  "N50 length,\n"
							  "minimal and maximal length.\n\n"
			  << separator;
}

//gets the project version from projectControl.h
void version_print()
{
	std::cout << PROJECT_VER << std::endl;
}

int calcAlignment(int size, const std::vector<std::unique_ptr<Sequence>> &fragment_list,
				  std::string *cigar, unsigned int *target_begin)
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

void createMinimizerIndex(const std::unique_ptr<Sequence> &seq,
						  std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> &index)
{
	//const char* seq_test = "TACGTACCGTA";
	auto minimizers_vector = //white::Minimize(seq_test, 11, kmer_len, window_len);
		white::Minimize(seq->getData().c_str(), seq->getData().size(), kmer_len, window_len);
	for (auto& vec : minimizers_vector) {
		std::cout << std::get<0>(vec) << ", " << std::get<1>(vec) << ", " << (std::get<2>(vec) ? "1" : "0") << "\n";
	}
	
	for (auto min : minimizers_vector)
	{
		index[std::get<0>(min)]
			.emplace_back(std::make_pair(std::get<1>(min), std::get<2>(min)));
	}
}

void createReferenceIndex(const std::unique_ptr<Sequence> &seq,
						  std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> &index)
{
	std::cerr << "Creating minimizer index for reference genome... ";
	createMinimizerIndex(seq, index);

	std::vector<std::pair<unsigned int, unsigned int>> distinct_minimizer_occurances;
	distinct_minimizer_occurances.reserve(index.size());

	int singleton_count = 0;
	for (auto &min : index)
	{
		if (min.second.size() == 1)
			singleton_count++;
		distinct_minimizer_occurances.emplace_back(std::make_pair(min.second.size(), min.first));
	}

	//sort minimizers by occurance descending
	sort(distinct_minimizer_occurances.begin(), distinct_minimizer_occurances.end(),
		 std::greater<std::pair<unsigned int, unsigned int>>());

	std::cerr << "Done!\n\n";

	int ignore_top_minimizers = std::ceil(index.size() * ignored_fraction);
	if (ignore_top_minimizers >= index.size())
		ignore_top_minimizers = index.size() - 1;
	std::cerr << "Reference index information:\n";
	std::cerr << "	Number of distinct minimizers: " << index.size() << std::endl;
	std::cerr << "	Fraction of singletons: " << (double)singleton_count / index.size() << std::endl;
	std::cerr << "	Number of occurrences of the most frequent minimizer with top "
			  << ignored_fraction * 100 << "% most frequent ignored: "
			  << distinct_minimizer_occurances[ignore_top_minimizers].first
			  << std::endl;

	//remove most frequently occuring minimizers from index
	for (int i = 0; i < ignore_top_minimizers; i++)
	{
		index.erase(distinct_minimizer_occurances[i].second);
	}
}

std::string mapToReference(const std::vector<std::unique_ptr<Sequence>> &fragments,
						   const std::unique_ptr<Sequence> &reference,
						   std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> &reference_index,
						   int fragments_begin, int fragments_end)
{
}

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

	while ((opt = getopt_long(argc, argv, "a:m:n:g:k:f:w:c:t:hv:", long_options, nullptr)) != -1)
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
	if (argc < 3)
	{
		std::cerr << "error: invalid input, please include exactly two files\n";
		return 1;
	}

	//if the given two arguments aren't in valid formats, the program recieves an error
	if (!checkArgs(argv, argc))
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
		std::cerr << "\nReference genome information:\n"
				  << "	Name of sequence: " << reference[i]->getName() << "\n"
				  << "	Length of sequence: " << reference[i]->getData().size() << "\n\n";
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

	std::cerr << "Fragment information:\n";
	std::cerr << "	Number of sequences: " << fragments_size << "\n";
	std::cerr << "	Total length of all fragments: " << sum << "\n";
	std::cerr << "	Minimum length: " << fragmentVector.front() << "\n";
	std::cerr << "	Maximum length: " << fragmentVector.back() << "\n";
	std::cerr << "	Average length of fragments: " << avg_size << "\n";
	std::cerr << "	N50 length: " << fragmentVector[N50] << "\n"
			  << separator;

	//ALIGNMENT
	std::string cigar;
	unsigned int target_begin;
	std::cout << "Aligning two random sequences...\n\n";
	int alignment_score = calcAlignment(fragments_size, fragments, &cigar, &target_begin);
	std::cout << "Alignment score: " << alignment_score << std::endl;
	std::cout << "Target begin index: " << target_begin << std::endl;
	if (calc_cigar)
		std::cout << "Cigar: " << cigar << std::endl;

	//MINIMIZERS
	std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> reference_index;
	createReferenceIndex(reference[0], reference_index);
	std::ios_base::sync_with_stdio(false);

	int fragments_in_thread = ceil(fragments_size / thread_count);
	int frag_begin = 0;
	thread_pool::ThreadPool thread_pool{};
	std::vector<std::future<std::string>> futures;

	for (int i = 0; i < (thread_count - 1); i++)
	{
		futures.emplace_back(
			thread_pool.Submit(mapToReference,
							   std::ref(fragments),
							   std::ref(reference[0]),
							   std::ref(reference_index),
							   frag_begin, frag_begin + fragments_in_thread));
		frag_begin += fragments_in_thread;
	}
	futures.emplace_back(
		thread_pool.Submit(mapToReference,
						   std::ref(fragments),
						   std::ref(reference[0]),
						   std::ref(reference_index),
						   frag_begin, frag_begin + fragments_in_thread));

	for (auto &fut : futures)
	{
		std::string s = fut.get();
		std::cout << s;
	}
	return 0;
}