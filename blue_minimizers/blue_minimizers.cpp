#include "blue_minimizers.hpp"

#include <bitset>
#include <deque>
#include <ios>
#include <iostream>
#include <string>
#include <unordered_set>

class Minimizer {
 private:
  const char* sequence;
  unsigned int sequence_len;
  unsigned int kmer_len;
  unsigned int window_len;
  unsigned int pop_front_mask;  // 0 for the first kmer char, 1 for others
  std::deque<blue::Kmer> window;
  std::unordered_set<blue::Kmer, blue::KmerHashFunction> minimizers;
  blue::Kmer current;    // current kmer
  blue::Kmer current_c;  // current kmer reverse complement
  int u;                 // max window size for end minimizers

  unsigned int Complement(unsigned int base) { return 3 - base; }

  // encodes given char as 2 bits
  unsigned int Encode(char c) {
    switch (c) {
      case 'C': return 0;
      case 'A': return 1;
      case 'T': return 2;
      case 'G': return 3;
    }
    throw "[blue::Minimize] Error: Invalid base";
  }

  // adds to window kmer whose last char is at given position
  void AddKmerToWindow(int position) {
    auto& [kmer, pos, original] = current;
    auto& [kmer_c, pos_c, original_c] = current_c;
    // create kmer
    kmer &= pop_front_mask;                           // delete first char
    kmer = (kmer << 2) | Encode(sequence[position]);  // add new char
    pos++;
    // create reverse complement kmer
    kmer_c &= pop_front_mask;  // delete first char
    kmer_c <<= 2;              // make room for new char
    kmer_c |= Complement(
        Encode(sequence[sequence_len - position - 1]));  // add new char
    pos_c++;

    while (!window.empty() && std::min(current, current_c) < window.back())
      window.pop_back();
    window.push_back(std::min(current, current_c));
  }

  // minimizers anchored at the end of the sequence
  // current and current_c should be first kmer when calling this method
  void StartMinimizer() {
    blue::Kmer old_current = current;
    blue::Kmer old_current_c = current_c;

    window.push_back(std::min(current, current_c));
    minimizers.insert(std::min(current, current_c));

    for (int i = kmer_len; i < u + kmer_len - 1; i++) {
      AddKmerToWindow(i);
      minimizers.insert(window.front());
    }
    // reset global variables
    current = old_current;
    current_c = old_current_c;
    window.clear();
  }

  // minimizers anchored to the end of the sequence
  // current and current_c should be last kmer when calling this method
  void EndMinimizer() {
    window.clear();
    window.push_back(std::min(current, current_c));
    minimizers.insert(std::min(current, current_c));

    auto& [kmer, pos, original] = current;
    auto& [kmer_c, pos_c, original_c] = current_c;

    for (int i = sequence_len - 1 - kmer_len; i > sequence_len - kmer_len - u;
         i--) {
      // create kmer
      kmer = (kmer >> 2) | (Encode(sequence[i]) << 4);  // add new char
      pos--;

      // create reverse complement kmer
      kmer_c = (kmer_c >> 2) |  // add new char
               (Complement(Encode(sequence[sequence_len - i - 1])) << 4);
      pos_c--;

      while (!window.empty() && std::min(current, current_c) < window.back())
        window.pop_back();
      window.push_back(std::min(current, current_c));
      minimizers.insert(window.front());
    }
  }

 public:
  Minimizer(const char* sequence, unsigned int sequence_len,
            unsigned int kmer_len, unsigned int window_len)
      : sequence(sequence),
        sequence_len(sequence_len),
        kmer_len(kmer_len),
        window_len(window_len) {
    pop_front_mask = (1 << (2 * kmer_len - 2)) - 1;
    u = window_len - 1;
  }

  std::vector<blue::Kmer> Minimize() {
    unsigned int first_kmer = 0;
    unsigned int first_kmer_c = 0;  // reverse complement first kmer

    // compute first kmer
    int i;
    for (i = 0; i < kmer_len; i++) {
      first_kmer = (first_kmer << 2) | Encode(sequence[i]);
      first_kmer_c = (first_kmer_c << 2) |
                     Complement(Encode(sequence[sequence_len - i - 1]));
    }

    current = {first_kmer, 0, true};
    current_c = {first_kmer_c, 0, false};  // reverse complement

    StartMinimizer();

    window.push_back(std::min(current, current_c));

    // add first window_len kmers to window
    for (; i < window_len + kmer_len - 1; i++) AddKmerToWindow(i);

    // loop over the rest
    for (; i < sequence_len; i++) {
      minimizers.insert(window.front());  // save old window minimizer

      // remove elements not in current window
      while (!window.empty() &&
             std::get<1>(window.front()) + (window_len + kmer_len - 1) <= i)
        window.pop_front();

      AddKmerToWindow(i);
    }
    minimizers.insert(window.front());

    // current and current_c hold last kmer
    EndMinimizer();

    return std::vector<blue::Kmer>(minimizers.begin(), minimizers.end());
  }
};

std::vector<blue::Kmer> blue::Minimize(const char* sequence,
                                       unsigned int sequence_len,
                                       unsigned int kmer_len,
                                       unsigned int window_len) {
  return Minimizer(sequence, sequence_len, kmer_len, window_len).Minimize();
}
