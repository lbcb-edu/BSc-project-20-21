#include "blue_alignment.hpp"
#include <gtest/gtest.h>

TEST(BlueAlignmentTests, Global) {
    const char* query = "GATTACA";
    unsigned int query_len = 7;
    const char* target = "GCATGCU";
    unsigned int target_len = 7;
    int match = 1;
    int mismatch = -1;
    int gap = -1;
    std::string cigar;
    unsigned int target_begin;

    int expected_value = 0;
    std::string expected_cigar = "1=1I1=1D1=1X1=1X";
    unsigned int expected_target_begin = 1;

    int value = blue::Align(query, query_len, target, target_len, blue::AlignmentType::kGlobal, match, mismatch, gap, &cigar, &target_begin);
    EXPECT_EQ(value, expected_value);
    EXPECT_EQ(cigar, expected_cigar);
    EXPECT_EQ(target_begin, expected_target_begin);

}

TEST(BlueAlignmentTests, Local) {
    const char* query = "ACCTAAGG";
    unsigned int query_len = 8;
    const char* target = "GGCTCAATCA";
    unsigned int target_len = 10;
    int match = 2;
    int mismatch = -1;
    int gap = -2;
    std::string cigar;
    unsigned int target_begin;

    int expected_value = 6;
    std::string expected_cigar = "2S2=1I2=2S";
    unsigned int expected_target_begin = 3;

    int value = blue::Align(query, query_len, target, target_len, blue::AlignmentType::kLocal, match, mismatch, gap, &cigar, &target_begin);
    EXPECT_EQ(value, expected_value);
    EXPECT_EQ(cigar, expected_cigar);
    EXPECT_EQ(target_begin, expected_target_begin);
}

TEST(BlueAlignmentTests, SemiGlobal) {
    const char* query = "CCTCGGTTA";
    unsigned int query_len = 9;
    const char* target = "GGTTAGAAAT";
    unsigned int target_len = 10;
    int match = 2;
    int mismatch = -1;
    int gap = -2;
    std::string cigar;
    unsigned int target_begin;

    int expected_value = 10;
    std::string expected_cigar = "4S5=";
    unsigned int expected_target_begin = 1;

    int value = blue::Align(query, query_len, target, target_len, blue::AlignmentType::kSemiGlobal, match, mismatch, gap, &cigar, &target_begin);
    EXPECT_EQ(value, expected_value);
    EXPECT_EQ(cigar, expected_cigar);
    EXPECT_EQ(target_begin, expected_target_begin);
}