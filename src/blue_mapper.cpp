#include "blue_mapper.hpp"

#include <getopt.h>

#include <iostream>
#include <random>

#include "alignment/blue_alignment.hpp"
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"

static constexpr option options[] = {
    {"algorithm", required_argument, nullptr, 'a'},
    {"match", required_argument, nullptr, 'm'},
    {"mismatch", required_argument, nullptr, 'n'},
    {"gap", required_argument, nullptr, 'g'},
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}};

void Version() {
  std::cout << "v" << blue_mapper_VERSION_MAJOR << "."
            << blue_mapper_VERSION_MINOR << "." << blue_mapper_VERSION_PATCH
            << std::endl;
}

void Help() {
  std::cout << "usage: blue_mapper [options] <reference-genome> <fragments>\n"
               "\n"
               "  <reference-genome>\treference genome in FASTA format"
               " (can be compressed with gzip)\n"
               "  <fragments>       \tset of fragments in FASTA/FASTQ format"
               "(can be compressed with gzip)\n"
               "\n"
               "  options:\n"
               "    -a, --algorithm <int>\n"
               "      default: 0\n"
               "      alignment algorithm:\n"
               "        0 - local (Smith-Waterman)\n"
               "        1 - global (Needleman-Wunsch)\n"
               "        2 - semi-global\n"
               "    -m, --match <int>\n"
               "      default: 5\n"
               "      score for matching bases\n"
               "    -n, --mismatch <int>\n"
               "      default: -4\n"
               "      score for mismatching bases\n"
               "    -g, --gap <int>\n"
               "      default: -8\n"
               "      linear gap penalty\n"
               "    -v, --version\n"
               "      Print version number\n"
               "    -h, --help\n"
               "      Print usage information\n";
}

void InvalidExtension(const std::string& file) {
  std::cerr
      << "Error: file " << file
      << " has unsupported format extension\n"
         "  Valid extensions for FASTA: .fasta, .fa, .fna, .ffn, .ffa, .frn\n"
         "  Valid extensions for FASTQ: .fastq, .fq\n"
         "\n"
         "The files can also be compressed with gzip\n";
}

class Sequence {
 public:
  std::string name_;
  std::string data_;
  std::string quality_;  // optional

  Sequence(const char* name, std::uint32_t name_len, const char* data,
           std::uint32_t data_len)
      : name_(name, name_len), data_(data, data_len) {}

  Sequence(const char* name, std::uint32_t name_len, const char* data,
           std::uint32_t data_len, const char* quality,
           std::uint32_t quality_len)
      : name_(name, name_len),
        data_(data, data_len),
        quality_(quality, quality_len) {}
};

bool HasExtension(const std::string& file,
                  const std::vector<std::string>& extensions) {
  auto has_suffix = [](const std::string& str, const std::string& suffix) {
    return str.size() >= suffix.size() ? str.compare(str.size() - suffix.size(),
                                                     suffix.size(), suffix) == 0
                                       : false;
  };

  for (auto& ext : extensions) {
    if (has_suffix(file, ext) || has_suffix(file, ext + ".gz")) return true;
  }
  return false;
}

bool IsFasta(const std::string& file) {
  return HasExtension(file, {".fasta", ".fa", ".fna", ".ffn", ".faa", ".frn"});
}

bool IsFastq(const std::string& file) {
  return HasExtension(file, {".fastq", ".fq"});
}

void PrintStats(const std::vector<std::unique_ptr<Sequence>>& fragments) {
  std::vector<size_t> lengths(fragments.size());
  uint64_t length_sum = 0;
  for (int i = 0; i < fragments.size(); i++) {
    length_sum += fragments[i]->data_.size();
    lengths[i] = fragments[i]->data_.size();
  }

  sort(lengths.rbegin(), lengths.rend());  // descending

  // calculate N50
  uint64_t N50 = 0;
  uint64_t partial_sum = 0;
  for (int i = 0; i < fragments.size(); i++) {
    partial_sum += lengths[i];
    if (partial_sum * 2 >= length_sum) {
      N50 = lengths[i];
      break;
    }
  }

  std::cerr << "Fragments\n"
            << "  Number of sequences: " << fragments.size() << '\n'
            << "  Average length:      "
            << length_sum / (double)fragments.size() << '\n'
            << "  N50 length:          " << N50 << '\n'
            << "  Minimal Length:      " << lengths.back() << '\n'
            << "  Maximal Length:      " << lengths[0] << "\n\n";
}

template <template <class> class T>
std::vector<std::unique_ptr<Sequence>> Parse(const std::string& file) {
  return bioparser::Parser<Sequence>::Create<T>(file)->Parse(-1);
}

int main(int argc, char* argv[]) {
  int opt;
  int8_t match_cost = 5;
  int8_t mismatch_cost = -4;
  int8_t gap_cost = -8;
  int8_t algorithm = 0;

  const char* opt_string = "a:m:n:g:hv";
  while ((opt = getopt_long(argc, argv, opt_string, options, nullptr)) != -1) {
    switch (opt) {
      case 'a': algorithm = atoi(optarg); break;
      case 'm': match_cost = atoi(optarg); break;
      case 'n': mismatch_cost = atoi(optarg); break;
      case 'g': gap_cost = atoi(optarg); break;
      case 'v': Version(); return 0;
      case 'h': Help(); return 0;
      default: return 1;
    }
  }

  if (argc != optind + 2) {
    std::cerr << "Error: wrong number of arguments." << std::endl;
    Help();
    return 1;
  }

  std::string reference_file = argv[optind];
  std::string fragments_file = argv[optind + 1];

  std::vector<std::unique_ptr<Sequence>> fragments;
  std::unique_ptr<Sequence> reference;

  std::cout << "Parsing..." << std::endl;
  try {
    if (IsFasta(reference_file)) {
      reference =
          std::move(Parse<bioparser::FastaParser>(reference_file).front());
    } else {
      InvalidExtension(reference_file);
      return 1;
    }

    if (IsFasta(fragments_file)) {
      fragments = Parse<bioparser::FastaParser>(fragments_file);
    } else if (IsFastq(fragments_file)) {
      fragments = Parse<bioparser::FastqParser>(fragments_file);
    } else {
      InvalidExtension(fragments_file);
      return 1;
    }
  } catch (const std::invalid_argument& exception) {
    // unable to open or corrupted
    std::cerr << exception.what() << std::endl;
    return 1;
  }

  std::cerr << "Reference genome\n"
            << "  Name:   " << reference->name_ << '\n'
            << "  Length: " << reference->data_.size() << "\n\n";

  PrintStats(fragments);

  // align two random sequences
  std::cout << "Aligning two random sequences..." << std::endl;
  std::vector<Sequence*> short_fragments;  // fragments with length < 5000
  for (auto& ptr : fragments)
    if (ptr->data_.size() < 5000) short_fragments.push_back(ptr.get());

  std::vector<Sequence*> out;
  std::sample(short_fragments.begin(), short_fragments.end(),
              std::back_inserter(out), 2, std::mt19937{std::random_device{}()});

  std::cout << "  Target\n"
            << "    Name:   " << out[0]->name_ << '\n'
            << "    Length: " << out[0]->data_.size() << '\n'
            << "  Query\n"
            << "    Name:   " << out[1]->name_ << '\n'
            << "    Length: " << out[1]->data_.size() << "\n\n";

  std::string& target = out[0]->data_;
  std::string& query = out[1]->data_;
  std::string cigar;
  unsigned int target_begin;
  int alignment_score =
      blue::Align(query.c_str(), query.size(), target.c_str(), target.size(),
                  static_cast<blue::AlignmentType>(algorithm), match_cost,
                  mismatch_cost, gap_cost, &cigar, &target_begin);

  std::cout << "  Alignment score:    " << alignment_score << '\n'
            << "  Target begin index: " << target_begin << "\n\n"
            << "  CIGAR: " << cigar << std::endl;

  return 0;
}
