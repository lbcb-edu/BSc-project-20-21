#include <gtest/gtest.h>
#include "../brown_alignment.hpp"

TEST(AlignTests, LocalAlignTest1) {
    int result = brown::Align("ACCTAAGG", 8, "GGCTCAATCA", 10, brown::LOCAL, 2, -1, -2, nullptr, nullptr);
    EXPECT_EQ(result, 6);
}

TEST(AlignTests, GlobalTest1) {
    int result = brown::Align("ATGGCCTC", 8, "ACGGCTC", 7, brown::GLOBAL, 1, -3, -4, nullptr, nullptr);
    EXPECT_EQ(result, -1);
}