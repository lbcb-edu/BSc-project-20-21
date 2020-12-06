#ifndef BLONDE_MINIMIZERS_H_
#define BLONDE_MINIMIZERS_H_

#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <deque>

namespace blonde {
namespace minimizers {

using Kmer = std::tuple<unsigned int, unsigned int, bool>;

std::vector<Kmer> Minimize(
    const char* sequence, unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len);

}
}

#endif
