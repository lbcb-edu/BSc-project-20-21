#include "blue_mapper.hpp"

#include <getopt.h>

#include <algorithm>
#include <ios>
#include <iostream>
#include <map>
#include <random>
#include <unordered_map>

#include "thread_pool/thread_pool.hpp"
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "blue_alignment.hpp"
#include "blue_minimizers.hpp"
#include "matcher.hpp"


static constexpr option options[] = {
    {"algorithm", required_argument, nullptr, 'a'},
    {"match", required_argument, nullptr, 'm'},
    {"mismatch", required_argument, nullptr, 'n'},
    {"gap", required_argument, nullptr, 'g'},
    {"kmer_len", required_argument, nullptr, 'k'},
    {"window_len", required_argument, nullptr, 'w'},
    {"ignored_fraction", required_argument, nullptr, 'f'},
    {"cigar", no_argument, nullptr, 'c'},
    {"threads", required_argument, nullptr, 't'},
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
               " (can be compressed with gzip)\n"
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
               "    -k, --kmer_len <int>\n"
               "      default: 15\n"
               "      kmer length\n"
               "    -w, --window_len <int>\n"
               "      default: 5\n"
               "      window length\n"
               "    -f, --ignored_fraction <double> in range [0,1)\n"
               "      default: 0.001\n"
               "      fraction of most frequent minimizers to be ignored\n"
               "    -c, --cigar\n"
               "      print CIGAR string\n"
               "    -t, --threads\n"
               "      default: 1\n"
               "      number of threads\n"
               "    -v, --version\n"
               "      print version number\n"
               "    -h, --help\n"
               "      print usage information\n";
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

void MinimizerStats(std::vector<std::unique_ptr<Sequence>>& sequences,
                    int8_t kmer_len, int8_t window_len,
                    double ignored_fraction) {
  std::unordered_map<unsigned int, unsigned int> minimizer_count;

  for (auto& seq : sequences) {
    auto minimizers = blue::Minimize(seq->data_.c_str(), seq->data_.size(),
                                     kmer_len, window_len);

    for (auto kmer : minimizers) minimizer_count[std::get<0>(kmer)]++;
  }

  std::vector<unsigned int> counts;
  counts.reserve(minimizer_count.size());

  unsigned int singletons = 0;
  for (auto pair : minimizer_count) {
    if (pair.second == 1) singletons++;
    counts.push_back(pair.second);
  }
  std::sort(counts.begin(), counts.end());

  unsigned int ignore_size = ignored_fraction * counts.size();
  std::cout << "  Number of distinct minimizers: " << minimizer_count.size()
            << '\n'
            << "  Fraction of singletons:        "
            << singletons / (double)minimizer_count.size() << '\n'
            << "  Number of occurrences of the most frequent minimizer: "
            << counts[counts.size() - ignore_size - 1] << "\n\n";
}

template <template <class> class T>
std::vector<std::unique_ptr<Sequence>> Parse(const std::string& file) {
  return bioparser::Parser<Sequence>::Create<T>(file)->Parse(-1);
}

// creates index of distinct minimizers, ignoring too frequent sequences
MinimizerIndex CreateMinimizerIndex(std::unique_ptr<Sequence>& sequence,
                                    int8_t kmer_len, int8_t window_len,
                                    double ignored_fraction) {
  auto minimizers = blue::Minimize(
      sequence->data_.c_str(), sequence->data_.size(), kmer_len, window_len);

  MinimizerIndex index;
  for (auto& [kmer, pos, origin] : minimizers)
    index[kmer].emplace_back(pos, origin);

  unsigned ignore_size = ignored_fraction * index.size();
  std::multimap<unsigned, unsigned> ignored;  // count -> kmer

  for (auto& [kmer, vec] : index) {
    ignored.emplace(vec.size(), kmer);
    if (ignored.size() > ignore_size)
      ignored.erase(ignored.begin());  // first elem is smallest
  }

  for (auto& [size, kmer] : ignored) index.erase(kmer);

  return index;
}

// Longest Increasing Subsequence | O(n*logn) | based on patience sort
Subsequence LIS(std::vector<Match>& matches, bool type) {
  if (matches.size() == 0) return {0, 0, 0, -1};

  std::vector<unsigned> tail;  // stores indexes of current end-elements
  std::vector<int> prev(matches.size(), -1);
  // stores 'backpointers' to reconstruct LIS

  auto comp = [&matches, &type](auto const& i, auto const& j) {
    return type ? matches[i].ref_pos > matches[j].ref_pos   // different strand
                : matches[i].ref_pos < matches[j].ref_pos;  // same strand
  };

  for (unsigned i = 0; i < matches.size(); i++) {
    auto it = std::lower_bound(tail.begin(), tail.end(), i, comp);

    if (it != tail.begin()) prev[i] = *(it - 1);
    if (it == tail.end()) {
      tail.push_back(i);
    } else {
      *it = i;
    }
  }

  unsigned beg_pos;
  for (int i = tail.back(); i >= 0; i = prev[i]) beg_pos = i;

  return {beg_pos, tail.back(), tail.size(), type};
}

// maps fragment to reference, prints mapping in PAF
void Map(MinimizerIndex reference_index, std::unique_ptr<Sequence>& reference,
         std::unique_ptr<Sequence>& fragment, int8_t kmer_len,
         int8_t window_len, blue::AlignmentType alignment_type, int match,
         int mismatch, int gap, bool print_cigar) {
  std::vector<blue::Kmer> frag_minimizers = blue::Minimize(
      fragment->data_.c_str(), fragment->data_.size(), kmer_len, window_len);

  // matches[0] - same strand, matches[1] - different strand
  std::vector<std::vector<Match>> matches(2);

  for (auto [frag_val, frag_pos, frag_origin] : frag_minimizers) {
    if (reference_index.count(frag_val)) {
      for (auto [ref_pos, ref_origin] : reference_index[frag_val])
        matches[frag_origin ^ ref_origin].push_back({frag_pos, ref_pos});
    }
  }

  auto comp = [](Match const& a, Match const& b) {
    return a.frag_pos < b.frag_pos;
  };

  Subsequence best = {0, 0, 0, -1};
  for (unsigned i = 0; i < matches.size(); i++) {
    if (matches[i].size() == 0) continue;

    sort(matches[i].begin(), matches[i].end(), comp);
    Subsequence s = LIS(matches[i], i);
    if (s.size > best.size) best = s;
  }

  int type = best.type;  // 0 - same strand, 1 - different
  if (type == -1) return;

  // raw indexes of the region to be aligned
  auto query_start = matches[type][best.beg].frag_pos;
  auto target_start = matches[type][best.beg].ref_pos;
  auto query_end = matches[type][best.end].frag_pos;
  auto target_end = matches[type][best.end].ref_pos;

  // lengths used for alignment
  auto query_align_len = query_end + kmer_len - query_start;
  auto target_align_len =
      (target_end + kmer_len - target_start) * (type ? -1 : 1);

  // apply orientation to indexes
  target_end = (type ? target_start : target_end + kmer_len);
  target_start = (type ? target_start - target_align_len : target_start);
  query_end += kmer_len;
  // TODO: granice +- kmer_len ?

  // clang-format off
  std::string paf = fragment->name_ + '\t'
    + std::to_string(fragment->data_.size()) + '\t'
    + std::to_string(query_start) + '\t'
    + std::to_string(query_end) + '\t'
    + (type ? '-' : '+') + '\t'
    + reference->name_ + '\t'
    + std::to_string(reference->data_.size()) + '\t'
    + std::to_string(target_start) + '\t'
    + std::to_string(target_end) + '\t';
  // clang-format on

  auto complement = [](char base) {
    switch (base) {
      case 'A': return 'T';
      case 'T': return 'A';
      case 'C': return 'G';
      case 'G': return 'C';
    }
    throw("[mapper:complement] Invalid base given as argument");
  };

  std::string cigar;
  if (print_cigar) {
    std::string query_rc;
    if (type) {  // compute query reverse complement
      for (int i = query_end; i >= query_start; i--) {
        query_rc.push_back(complement(fragment->data_[i]));
      }
    }

    std::cout << query_rc.size() << " " << query_align_len << " "
              << target_align_len << std::endl;
    unsigned _;  // not used
    blue::Align(
        (type ? query_rc.c_str() : fragment->data_.c_str() + query_start),
        query_align_len, reference->data_.c_str() + target_start,
        target_align_len, alignment_type, match, mismatch, gap, &cigar, &_);

    unsigned match_count = 0;
    unsigned total = 0;  // all operations count
    unsigned buff = 0;
    for (char c : cigar) {
      if (std::isdigit(c)) {
        buff *= 10;
        buff += c - '0';
      } else {
        if (c == '=') match_count += buff;
        total += buff;
        buff = 0;
      }
    }

    paf += std::to_string(match_count) + '\t' + std::to_string(total) + '\t';
  } else {
    // unsigned aprox10 = 0, aprox11 = 0;  // TODO aproximate columns 10 & 11
    // paf += std::to_string(aprox10) + '\t' + std::to_string(aprox11) + '\t';
  }

  unsigned quality = 255;
  paf += std::to_string(quality) + '\t';

  if (print_cigar) paf += "cg:Z:" + cigar;

  std::cout << paf << std::endl;
}

int main(int argc, char* argv[]) {
  int8_t match_cost = 5;
  int8_t mismatch_cost = -4;
  int8_t gap_cost = -8;
  int8_t algorithm = 0;
  int8_t kmer_len = 15;
  int8_t window_len = 5;
  int8_t num_of_threads = 1;
  double ignored_fraction = 0.001;
  bool print_cigar = false;

  const char* opt_string = "a:m:n:g:k:w:f:t:vhc";
  int opt;
  while ((opt = getopt_long(argc, argv, opt_string, options, nullptr)) != -1) {
    switch (opt) {
      case 'a': algorithm = atoi(optarg); break;
      case 'm': match_cost = atoi(optarg); break;
      case 'n': mismatch_cost = atoi(optarg); break;
      case 'g': gap_cost = atoi(optarg); break;
      case 'k': kmer_len = atoi(optarg); break;
      case 'w': window_len = atoi(optarg); break;
      case 'f': ignored_fraction = atof(optarg); break;
      case 't': num_of_threads = atof(optarg); break;
      case 'c': print_cigar = true; break;
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
  if (ignored_fraction < 0 || ignored_fraction >= 1) {
    std::cerr << "Error: ignored_fraction out of range [0,1). " << std::endl;
    return 1;
  }

  std::string reference_file = argv[optind];
  std::string fragments_file = argv[optind + 1];

  std::vector<std::unique_ptr<Sequence>> fragments;
  std::vector<std::unique_ptr<Sequence>> reference;

  std::cout << "Parsing..." << std::endl;
  try {
    if (IsFasta(reference_file)) {
      reference = Parse<bioparser::FastaParser>(reference_file);
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

  std::cerr << "Reference genome\n";
  for (auto& ref : reference) {
    std::cerr << "  Name:   " << ref->name_ << '\n'
              << "  Length: " << ref->data_.size() << "\n\n";
  }

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
            << "  CIGAR: " << cigar << "\n\n";

  // MINIMIZER STATS
  std::cout << "Reference minimizers" << std::endl;
  MinimizerStats(reference, kmer_len, window_len, ignored_fraction);
  std::cout << "Fragments minimizers" << std::endl;
  MinimizerStats(fragments, kmer_len, window_len, ignored_fraction);

  // MAPPING
  std::ios_base::sync_with_stdio(false);


  thread_pool::ThreadPool thread_pool(num_of_threads);
  std::vector<std::future<void>> futures;

  blue::AlignmentType alignment_type = static_cast<blue::AlignmentType>(algorithm);

  for (auto& ref : reference) {
    MinimizerIndex reference_index =
        CreateMinimizerIndex(ref, kmer_len, window_len, ignored_fraction);
    
    for (auto& frag : fragments) {
      futures.emplace_back(thread_pool.
                    Submit(Map,std::ref(reference_index),std::ref(ref),std::ref(frag),
                    std::ref(kmer_len),std::ref(window_len), std::ref(alignment_type),
                    std::ref(match_cost),std::ref(mismatch_cost),std::ref(gap_cost),std::ref(print_cigar)));
    }

    for (const auto& it : futures) {
      it.wait();
    }
  }
  
  return 0;
}
