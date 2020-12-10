#pragma once
#include <tuple>
#include <vector>

namespace blue {

using Kmer = std::tuple<unsigned int, unsigned int, bool>;

// returns minimizers for given sequence
std::vector<Kmer> Minimize(const char* sequence, unsigned int sequence_len,
                           unsigned int kmer_len, unsigned int window_len);

}  // namespace blue
