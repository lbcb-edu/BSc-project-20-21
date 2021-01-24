#ifndef BROWN_MINIM_H
#include <iostream>
#include <vector>

#define BROWN_MINIM_H
    namespace brown {

            std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
                const char* sequence, unsigned int sequence_len,
                unsigned int kmer_len,
                unsigned int window_len);

        unsigned int getReversedComplKmerValue(unsigned int value, unsigned int length);
    }
#endif