#include <gtest/gtest.h>
#include "../include/brown_minimizer.hpp"

TEST(MinimizerTests, MinimizerTest1) {
    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = brown::Minimize("CACGATATCG", 10, 3, 4);
    EXPECT_EQ(minimizers.size(), 5);
    EXPECT_EQ(std::get<0>(minimizers[0]), 17);
    EXPECT_EQ(std::get<1>(minimizers[0]), 0);
    EXPECT_EQ(std::get<2>(minimizers[0]), true);

    EXPECT_EQ(std::get<0>(minimizers[1]), 6);
    EXPECT_EQ(std::get<1>(minimizers[1]), 1);
    EXPECT_EQ(std::get<2>(minimizers[1]), true);

    EXPECT_EQ(std::get<0>(minimizers[2]), 12);
    EXPECT_EQ(std::get<1>(minimizers[2]), 4);
    EXPECT_EQ(std::get<2>(minimizers[2]), true);

    EXPECT_EQ(std::get<0>(minimizers[3]), 13);
    EXPECT_EQ(std::get<1>(minimizers[3]), 6);
    EXPECT_EQ(std::get<2>(minimizers[3]), true);

    EXPECT_EQ(std::get<0>(minimizers[4]), 54);
    EXPECT_EQ(std::get<1>(minimizers[4]), 7);
    EXPECT_EQ(std::get<2>(minimizers[4]), true);

    /*for (int i = 0; i < minimizers.size(); i++) {
        std::cerr << "Value of " << i << ": " << std::get<0>(minimizers[i]) << " on position " << std::get<1>(minimizers[i]) << std::endl;
    }*/

}

TEST(MinimizerTests, MinimizerTest2) {
    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = brown::Minimize("TCGTATAGCTCAACT", 15, 3, 5);
    EXPECT_EQ(minimizers.size(), 6);
    EXPECT_EQ(std::get<0>(minimizers[0]), 54);
    EXPECT_EQ(std::get<1>(minimizers[0]), 0);
    EXPECT_EQ(std::get<2>(minimizers[0]), true);

    EXPECT_EQ(std::get<0>(minimizers[1]), 27);
    EXPECT_EQ(std::get<1>(minimizers[1]), 1);
    EXPECT_EQ(std::get<2>(minimizers[1]), true);

    EXPECT_EQ(std::get<0>(minimizers[2]), 12);
    EXPECT_EQ(std::get<1>(minimizers[2]), 4);
    EXPECT_EQ(std::get<2>(minimizers[2]), true);

    EXPECT_EQ(std::get<0>(minimizers[3]), 9);
    EXPECT_EQ(std::get<1>(minimizers[3]), 6);
    EXPECT_EQ(std::get<2>(minimizers[3]), true);

    EXPECT_EQ(std::get<0>(minimizers[4]), 1);
    EXPECT_EQ(std::get<1>(minimizers[4]), 11);
    EXPECT_EQ(std::get<2>(minimizers[4]), true);

    EXPECT_EQ(std::get<0>(minimizers[5]), 7);
    EXPECT_EQ(std::get<1>(minimizers[5]), 12);
    EXPECT_EQ(std::get<2>(minimizers[5]), true);

    /*for (int i = 0; i < minimizers.size(); i++) {
        std::cerr << "Value of " << i << ": " << std::get<0>(minimizers[i]) << " on position " << std::get<1>(minimizers[i]) << std::endl;
    }*/

}