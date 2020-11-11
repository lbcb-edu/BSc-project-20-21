#include "blue_mapper.hpp"

#include <getopt.h>

#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"

static constexpr option options[] = {{"version", no_argument, nullptr, 'v'},
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
               "  options\n"
               "    -v, --version\tPrint version number\n"
               "    -h, --help   \tPrint usage information\n";
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
            << "  Maximal Length:      " << lengths[0] << '\n';
}

template <template <class> class T>
std::vector<std::unique_ptr<Sequence>> Parse(const std::string& file) {
  return bioparser::Parser<Sequence>::Create<T>(file)->Parse(-1);
}

int main(int argc, char* argv[]) {
  int opt;
  while ((opt = getopt_long(argc, argv, "hv", options, nullptr)) != -1) {
    switch (opt) {
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

  return 0;
}
