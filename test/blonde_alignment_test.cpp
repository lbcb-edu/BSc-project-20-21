#include <gtest/gtest.h>
#include "../src/blonde_alignment.h"

TEST(BlondeAlignmentTests, Global_1) {
    const char* target = "GCATGCU";
    const char* query = "GATTACA";

    std::string cigar;
    unsigned int target_begin;

    int value = blonde::alignment::Align(query, strlen(query), target, strlen(target), blonde::alignment::AlignmentType::kGlobal, 1, -1, -1, &cigar, &target_begin);
    EXPECT_EQ(value, 0);
    EXPECT_EQ(cigar, "1=1D1=1I1=1X1=1X");
    EXPECT_EQ(target_begin, 0);
}

TEST(BlondeAlignmentTests, Global_2) {
    const char* target = "AGTCTTTGAT";
    const char* query = "GATTAGA";

    std::string cigar;
    unsigned int target_begin;

    int value = blonde::alignment::Align(query, strlen(query), target, strlen(target), blonde::alignment::AlignmentType::kGlobal, 1, -1, -1, &cigar, &target_begin);
    EXPECT_EQ(value, 0);
    EXPECT_EQ(cigar, "1D1=1D1X2=1X2=1D");
    EXPECT_EQ(target_begin, 0);
}