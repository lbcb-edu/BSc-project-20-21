#include "../src/blonde_alignment.h"
#include <gtest/gtest.h>

TEST(BlondeAlignmentTests, Global) {
    const char* query = "GATTACA";
    const char* target = "GCATGCU";

    std::string cigar;
    unsigned int target_begin;

    int value = blonde::alignment::Align(query, strlen(query), target, strlen(target), blonde::alignment::AlignmentType::kGlobal, 1, -1, -1, &cigar, &target_begin);
    EXPECT_EQ(value, 0);
    EXPECT_EQ(cigar, "1=1D1=1I1=1X1=1X");
    EXPECT_EQ(target_begin, 0);
}
