#include "blue_minimizers.hpp"
#include <gtest/gtest.h>

TEST(BlueMinimizersTest, First) {
    const char* sequence = "ACGTGAC";
    unsigned int sequence_len = 7;
    unsigned int kmer_len = 3;
    unsigned int window_len = 3;
    std::vector<blue::Kmer> minimizers = blue::Minimize(sequence, sequence_len, kmer_len, window_len);
    int expected_values[4][3] = {{14,4,0},{4,2,0},{14,1,1},{19,0,1}};
    int k;
    for (auto [i,j] = std::tuple{minimizers.begin(),0}; i != minimizers.end(); i++, j++) {
        k = 0;
        EXPECT_EQ(expected_values[j][k++], std::get<0>(*i));
        EXPECT_EQ(expected_values[j][k++], std::get<1>(*i));
        EXPECT_EQ(expected_values[j][k], std::get<2>(*i));
    }
}