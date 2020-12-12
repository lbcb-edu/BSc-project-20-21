#pragma once
#include <tuple>
#include <vector>

namespace blue {

using Kmer = std::tuple<unsigned int, unsigned int, bool>;

struct KmerHashFunction {
  size_t operator()(const blue::Kmer& kmer) const {
    return std::hash<unsigned int>()(std::get<0>(kmer)) ^
           std::hash<int>()(std::get<1>(kmer)) ^
           std::hash<bool>()(std::get<2>(kmer));
  }
};

// returns minimizers for given sequence
std::vector<Kmer> Minimize(const char* sequence, unsigned int sequence_len,
                           unsigned int kmer_len, unsigned int window_len);

}  // namespace blue
