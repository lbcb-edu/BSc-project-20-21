
#include "blonde_minimizers.h"

namespace std {
using blonde::minimizers::Kmer;

template <>
struct hash<Kmer> {
    inline size_t operator()(const Kmer& kmer) const { // maybe some other hash function would be better ?
        return hash<int>()(std::get<0>(kmer) ^ std::get<1>(kmer) ^ std::get<2>(kmer));
    }
};
}

namespace blonde {
namespace minimizers {  

int get_base_mask(const char base) { // base encoding that matches up with the lexicographical ordering of the base pairs
    switch(base) {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    case 'T':
        return 3;
    }
    throw "invalid base pair";
}

int complement_mask(int mask) {
    return 3 - mask;
}

void load_next_base(const char base, unsigned int kmer_len, unsigned int lengthwise_mask, unsigned int& kmer, unsigned int& comp_kmer) {
    unsigned int base_mask = get_base_mask(base);
    unsigned int comp_base_mask = complement_mask(base_mask);

    kmer <<= 2;
    kmer |= base_mask;
    kmer &= lengthwise_mask;

    comp_base_mask <<= kmer_len * 2;
    comp_kmer |= comp_base_mask;
    comp_kmer >>= 2;
}

std::vector<Kmer> get_all_kmers(const char* sequence, unsigned int sequence_len, unsigned int kmer_len) {
    if (kmer_len > 16) throw "only offers support for 1-16 mers"; // jer za svaki nukleotid koristimo 2 bita a int imma 32 bita
    if (kmer_len > sequence_len) throw "kmer length longer than sequence!";
    
    std::vector<Kmer> result;
    result.reserve(sequence_len - kmer_len + 1);
    unsigned int lengthwise_mask = ~0U;
    lengthwise_mask >>= 32 - (kmer_len * 2);
    unsigned int kmer = 0;
    unsigned int comp_kmer = 0;

    //load first kmer
    int i = 0;
    for (; i < kmer_len; i++) load_next_base(sequence[i], kmer_len, lengthwise_mask, kmer, comp_kmer);
    
    // store first kmer and load and store all other kmers
    for(; i <= sequence_len; i++) {
        if(comp_kmer < kmer) {
            result.emplace_back(comp_kmer, i - kmer_len, false);
        } else {
            result.emplace_back(kmer, i - kmer_len, true);
        }
        if(i < sequence_len) load_next_base(sequence[i], kmer_len, lengthwise_mask, kmer, comp_kmer);
    }
    
    // std::cout << "KMERS:\n";
    // for(Kmer mer : result) {
    //     std::cout << std::hex << std::get<0>(mer);
    //     std::cout << " ";
    //     std::cout << std::hex << std::get<1>(mer);
    //     std::cout << " ";
    //     std::cout << std::hex << std::get<2>(mer);
    //     std::cout << " ";
    //     std::cout << std::endl;
    // }

    return result;
}

bool Kmer_compare(const Kmer& first, const Kmer& second) {
    return std::get<0>(first) < std::get<0>(second);
}

std::vector<Kmer> Minimize(
    const char* sequence, unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len) {

    std::unordered_set<Kmer> minimizers;
    std::vector<Kmer> kmers = get_all_kmers(sequence, sequence_len, kmer_len);

    // Add beginning and end minimizers (not sure if we need to do this or not)
    Kmer minimal_kmer_begin = kmers[0];
    Kmer minimal_kmer_end = kmers[sequence_len - kmer_len];
    minimizers.insert(minimal_kmer_begin);
    minimizers.insert(minimal_kmer_end);
    for (int i = 1; i < window_len - 1; i++) {
        if (Kmer_compare(kmers[i], minimal_kmer_begin)) 
            minimal_kmer_begin = kmers[i];
        if (Kmer_compare(kmers[sequence_len - kmer_len - i], minimal_kmer_end)) 
            minimal_kmer_end = kmers[sequence_len - kmer_len - i];
        minimizers.insert(minimal_kmer_begin);
        minimizers.insert(minimal_kmer_end);
    }

    //Add normal minimizers
    std::deque<Kmer> kmers_in_window;
    for(int i = 0; i < window_len; i++) {
        kmers_in_window.push_back(kmers[i]);
    }
    minimizers.insert(*std::min_element(kmers_in_window.begin(), kmers_in_window.end(), Kmer_compare));
    for(int i = window_len; i < kmers.size(); i++) {
        kmers_in_window.pop_front();
        kmers_in_window.push_back(kmers[i]);
        minimizers.insert(*std::min_element(kmers_in_window.begin(), kmers_in_window.end(), Kmer_compare));
    }

    std::vector<Kmer> result(minimizers.begin(), minimizers.end());
    return result;
}

}
}