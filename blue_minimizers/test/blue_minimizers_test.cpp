#include "blue_minimizers.hpp"
#include <gtest/gtest.h>

bool compareKmers(blue::Kmer kmer1, blue::Kmer kmer2) {
    return (std::get<1>(kmer1) < std::get<1>(kmer2));
}

TEST(BlueMinimizersTest, First) {
    const char* sequence = "ACGTGAC";
    unsigned int sequence_len = 7;
    unsigned int kmer_len = 3;
    unsigned int window_len = 3;
    std::vector<blue::Kmer> minimizers = blue::Minimize(sequence, sequence_len, kmer_len, window_len);
    int expected_values[4][3] = {{19,0,1},{14,1,1},{4,2,0},{14,4,0}};

    std::sort(minimizers.begin(),minimizers.end(),compareKmers);

    for (auto [i,j]= std::tuple{minimizers.begin(),0}; i != minimizers.end(); i++,j++ ) {
        EXPECT_EQ(expected_values[j][0], std::get<0>(*i));
        EXPECT_EQ(expected_values[j][1], std::get<1>(*i));
        EXPECT_EQ(expected_values[j][2], std::get<2>(*i));
    }
}