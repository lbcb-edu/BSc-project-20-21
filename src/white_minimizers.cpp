#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <deque>
#include <unordered_set>
#include <algorithm>
#include "white_minimizers.hpp"

class Minimizer
{
private:
    const char *sequence;
    unsigned int sequence_len;
    unsigned int kmer_len;
    unsigned int window_len;
    std::tuple<unsigned int, unsigned int, bool> current_kmer;
    std::tuple<unsigned int, unsigned int, bool> current_kmer_c;
    unsigned int delete_first_char;

    struct kmer_hash
    {
        size_t operator()(const std::tuple<unsigned int, unsigned int, bool> &kmer) const
        {
            return std::hash<unsigned int>()(std::get<0>(kmer)) ^
                   std::hash<int>()(std::get<1>(kmer)) ^
                   std::hash<bool>()(std::get<2>(kmer));
        }
    };

    unsigned int mapLetter(char letter)
    {
        switch (letter)
        {
        case 'C':
            return 0;
        case 'A':
            return 1;
        case 'T':
            return 2;
        case 'G':
            return 3;
        default:
            return -1;
        }
    }

    unsigned int Complement(unsigned int base)
    {
        return 3 - base;
    }

    static bool compareKmers(std::tuple<unsigned int, unsigned int, bool> first, std::tuple<unsigned int, unsigned int, bool> second)
    {
        return std::get<0>(first) < std::get<0>(second);
    }

    void nextKmer(unsigned int seq_position)
    {
        unsigned int new_kmer_value = std::get<0>(current_kmer) & delete_first_char;
        new_kmer_value = (new_kmer_value << 2) | mapLetter(sequence[seq_position]);
        unsigned int new_kmer_pos = seq_position - kmer_len + 1;

        unsigned int new_kmer_c_value = std::get<0>(current_kmer) & delete_first_char;
        new_kmer_c_value = (new_kmer_c_value << 2) | mapLetter(sequence[seq_position]);

        current_kmer = {new_kmer_value, new_kmer_pos, true};
        current_kmer_c = {new_kmer_c_value, new_kmer_pos, false};
    }

public:
    Minimizer(const char *sequence, uint32_t sequence_len,
              unsigned int kmer_len, unsigned int window_len) : sequence(sequence),
                                                                sequence_len(sequence_len),
                                                                kmer_len(kmer_len),
                                                                window_len(window_len)
    {
        delete_first_char = (1 << (2 * kmer_len - 2)) - 1;
    }

    std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize()
    {
        unsigned int kmer_value = 0;
        unsigned int kmer_c_value = 0;
        std::deque<std::tuple<unsigned int, unsigned int, bool>> window;
        std::unordered_set<std::tuple<unsigned int, unsigned int, bool>, kmer_hash> minimizers;

        std::cout << "test1"
				  << "\n\n";

        //create first kmer
        for (int i = 0; i < kmer_len; i++)
        {
            kmer_value = (kmer_value << 2) | mapLetter(sequence[i]);
            kmer_c_value = (kmer_c_value << 2) | Complement(mapLetter(sequence[i]));
        }
        current_kmer = {kmer_value, 0, true};
        current_kmer_c = {kmer_c_value, 0, false};

        std::cout << "first kmer created"
				  << "\n\n";

        //save first kmer and kmer_c
        auto first_kmer = current_kmer;
        auto first_kmer_c = current_kmer_c;

        //beginning minimizers
        auto min_begin_kmer = std::min(current_kmer, current_kmer_c, compareKmers);
        for (int i = 1; i < window_len - 1; i++)
        {
            nextKmer(i);
            min_begin_kmer = std::min(std::min(current_kmer, current_kmer_c, compareKmers), min_begin_kmer, compareKmers);
            minimizers.insert(min_begin_kmer);
        }

        std::cout << "begin minimizers done"
				  << "\n\n";

        //beggining minimizers inserted, reset current_kmer and current_kmer_c to first_kmer
        current_kmer = first_kmer;
        current_kmer_c = first_kmer_c;

        //interior minimizers

        //add first kmer to window
        window.push_back(std::min(current_kmer, current_kmer_c, compareKmers));

        //fill window with first window_len kmers
        for (int i = kmer_len; i < window_len + kmer_len - 1; i++)
        {
            nextKmer(i);
            window.push_back(std::min(current_kmer, current_kmer_c, compareKmers));
        }
        minimizers.insert(*std::min_element(window.begin(), window.end(), compareKmers));

        std::cout << "window filled with first window_len kmers"
				  << "\n\n";
        //iterate through the rest of the sequence
        for (int i = window_len + kmer_len - 1; i < sequence_len; i++)
        {
            window.pop_front();
            nextKmer(i);
            window.push_back(std::min(current_kmer, current_kmer_c, compareKmers));
            minimizers.insert(*std::min_element(window.begin(), window.end(), compareKmers));
        }

        std::cout << "interior minimizers done"
				  << "\n\n";
        //current_kmer is last kmer in sequence - start end minimizer
        auto min_end_kmer = std::min(current_kmer, current_kmer_c, compareKmers);
        minimizers.insert(min_end_kmer);
        for (int i = sequence_len - kmer_len - 1; i > sequence_len - kmer_len - window_len + 1; i--)
        {
            unsigned int new_kmer_value = (std::get<0>(current_kmer) >> 2) | (mapLetter(sequence[i]) << 4);
            unsigned int new_kmer_c_value = (std::get<0>(current_kmer_c) >> 2) | (Complement(mapLetter(sequence[i]) << 4));
            unsigned int new_kmer_pos = i;

            current_kmer = {new_kmer_value, new_kmer_pos, true};
            current_kmer_c = {new_kmer_c_value, new_kmer_pos, false};

            min_end_kmer = std::min(std::min(current_kmer, current_kmer_c, compareKmers), min_end_kmer, compareKmers);
            minimizers.insert(min_end_kmer);
        }

        std::cout << "end minimizers done"
				  << "\n\n";
        return std::vector<std::tuple<unsigned int, unsigned int, bool>>(minimizers.begin(), minimizers.end());
    }
};

std::vector<std::tuple<unsigned int, unsigned int, bool>> white::Minimize(
    const char *sequence, uint32_t sequence_len,
    unsigned int kmer_len, unsigned int window_len)
{
    return Minimizer(sequence, sequence_len, kmer_len, window_len).Minimize();
}