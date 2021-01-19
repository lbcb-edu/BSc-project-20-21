#pragma once
#include <unordered_map>
#include <vector>

using MinimizerIndex =  // kmer -> vector<{position, origin}>
    std::unordered_map<unsigned, std::vector<std::pair<unsigned, bool>>>;

struct Match {
  unsigned int frag_pos;
  unsigned int ref_pos;
};

struct Subsequence {  // holds positions
  unsigned beg;
  unsigned end;
  size_t size;
  int type;
};
