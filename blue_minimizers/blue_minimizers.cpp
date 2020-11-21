#include "blue_minimizers.hpp"

#include <bitset>
#include <deque>
#include <ios>
#include <iostream>
#include <string>
#include <unordered_set>

// encodes given char as 2 bits
static unsigned int Encode(char c, bool even) {
  switch (c) {
    case 'C': return even ? 3 : 0;
    case 'A': return even ? 2 : 1;
    case 'T': return even ? 1 : 2;
    case 'G': return even ? 0 : 3;
  }
  throw "[blue::Minimize] Error: Invalid base";
}

struct KmerHashFunction {
  size_t operator()(const blue::Kmer& kmer) const {
    return std::hash<unsigned int>()(std::get<0>(kmer)) ^
           std::hash<int>()(std::get<1>(kmer)) ^
           std::hash<bool>()(std::get<2>(kmer));
  }
};


//minimizer anchored to the start of the sequence
void StartMinimizer(const char* sequence,  
                    unsigned int sequence_len,
                    unsigned int kmer_len,
                    unsigned int window_len,
                    std::unordered_set<blue::Kmer, KmerHashFunction> &minimizers) {

  unsigned int pop_front_mask = 
      (1 << (2 * kmer_len - 2)) - 1;  // 0 for the first kmer char, 1 for others
  
  unsigned int first_kmer = 0;
  
  int i;
  for (i = 0; i < kmer_len; i++)
      first_kmer = (first_kmer << 2) | Encode(sequence[i], i % 2 == 0);
  
  int u = window_len-1;

  blue::Kmer current;
  std::deque<blue::Kmer> window;

  current = {first_kmer, 0, true};
  window.push_back(current);
  auto& [kmer, pos, original] = current;

  minimizers.insert(window.front()); //add the first kmer to the list of minimizers

  for (; i < u + kmer_len - 1; i++) {
    // create kmer
    kmer &= pop_front_mask;                                // delete first char
    kmer = (kmer << 2) | Encode(sequence[i], i % 2 == 0);  // add new char
    pos++;

    while (!window.empty() && current < window.back()) window.pop_back();
    window.push_back(current);
    minimizers.insert(window.front());
  }
}


//minimizer anchored to the end of the sequence
void EndMinimizer(const char* sequence,
                  unsigned int sequence_len,
                  unsigned int kmer_len,
                  unsigned int window_len,
                  std::unordered_set<blue::Kmer, KmerHashFunction> &minimizers) {

  unsigned int first_kmer =  0;
  
  int i;
  for (i = sequence_len-1; i > sequence_len-1 - kmer_len; i--)
      first_kmer = (first_kmer >> 2) | (Encode(sequence[i], i % 2 == 0) << 4);

  int u = window_len-1;
  
  blue::Kmer current;
  std::deque<blue::Kmer> window;

  current = {first_kmer, sequence_len - kmer_len, true};
  window.push_back(current);
  auto& [kmer, pos, original] = current;

  minimizers.insert(window.front()); //add the last kmer to the list of minimizers

  for (; i > sequence_len - kmer_len - u; i--) {
    // create kmer
    kmer = (kmer >> 2) | (Encode(sequence[i], i % 2 == 0) << 4);  // add new char
    pos--;

    while (!window.empty() && current < window.back()) window.pop_back();
    window.push_back(current);
    minimizers.insert(window.front());
  }
}


// TODO: reverse complement
std::vector<blue::Kmer> blue::Minimize(const char* sequence,
                                       unsigned int sequence_len,
                                       unsigned int kmer_len,
                                       unsigned int window_len) {

  std::unordered_set<blue::Kmer, KmerHashFunction> minimizers;

  StartMinimizer(sequence, sequence_len, kmer_len, window_len, minimizers);
  EndMinimizer(sequence, sequence_len, kmer_len, window_len, minimizers);

  unsigned int pop_front_mask = 
      (1 << (2 * kmer_len - 2)) - 1;  // 0 for the first kmer char, 1 for others

  blue::Kmer current;
  std::deque<blue::Kmer> window;
  unsigned int first_kmer = 0;

  int i;
  for (i = 0; i < kmer_len; i++)
    first_kmer = (first_kmer << 2) | Encode(sequence[i], i % 2 == 0);

  current = {first_kmer, 0, true};
  window.push_back(current);
  auto& [kmer, pos, original] = current;

  // add first w kmers to window
  for (; i < window_len + kmer_len - 1; i++) {
    // create kmer
    kmer &= pop_front_mask;                                // delete first char
    kmer = (kmer << 2) | Encode(sequence[i], i % 2 == 0);  // add new char
    pos++;

    while (!window.empty() && current < window.back()) window.pop_back();
    window.push_back(current);
  }

  // loop over the rest
  for (; i < sequence_len; i++) {
    minimizers.insert(window.front());  // save old window minimizer

    // remove elements not in current window
    while (!window.empty() &&
           std::get<1>(window.front()) + (window_len + kmer_len - 1) <= i)
      window.pop_front();

    // create kmer
    kmer &= pop_front_mask;                                // delete first char
    kmer = (kmer << 2) | Encode(sequence[i], i % 2 == 0);  // add new char
    pos++;

    while (!window.empty() && current < window.back()) window.pop_back();
    window.push_back(current);
  }

  minimizers.insert(window.front());

  return std::vector<blue::Kmer>(minimizers.begin(), minimizers.end());
}

// Minimizer = lexicographically smallest substring of size kmer_len in window
// of size window_len + kmer_len - 1
