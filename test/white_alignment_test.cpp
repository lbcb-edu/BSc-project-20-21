#include "white_alignment.h"
#include <gtest/gtest.h>

TEST(WhiteAlignmentTests, Global) {
    const char* query = "GATTACA";
    unsigned int query_len = 7;
    const char* target = "GCATGCU";
    unsigned int target_len = 7;
    int match_cost = 1;
    int mismatch_cost = -1;
    int gap_cost = -1;
    std::string cigar;
    unsigned int target_begin;

    int expected_value = 0;
    std::string expected_cigar = "1=1I1=1D1=1X1=1X";
    unsigned int expected_target_begin = 1;

    white::Aligner *aligner = new white::Aligner(query, query_len, target, target_len,
		match_cost, mismatch_cost, gap_cost, &cigar, &target_begin);
    
    int result = aligner -> Align(white::AlignmentType::GLOBAL);
    
    EXPECT_EQ(result, expected_value);
    EXPECT_EQ(cigar, expected_cigar);
    EXPECT_EQ(target_begin, expected_target_begin);
}

TEST(WhiteAlignmentTests, Local) {
    const char* query = "ACCTAAGG";
    unsigned int query_len = 8;
    const char* target = "GGCTCAATCA";
    unsigned int target_len = 10;
    int match_cost = 2;
    int mismatch_cost = -1;
    int gap_cost = -2;
    std::string cigar;
    unsigned int target_begin;

    int expected_value = 6;
    std::string expected_cigar = "2S2=1I2=2S";
    unsigned int expected_target_begin = 3;

    white::Aligner *aligner = new white::Aligner(query, query_len, target, target_len,
		match_cost, mismatch_cost, gap_cost, &cigar, &target_begin);
    
    int result = aligner -> Align(white::AlignmentType::LOCAL);
    
    EXPECT_EQ(result, expected_value);
    EXPECT_EQ(cigar, expected_cigar);
    EXPECT_EQ(target_begin, expected_target_begin);
}

TEST(WhiteAlignmentTests, SemiGlobal) {
    const char* query = "CCTCGGTTA";
    unsigned int query_len = 9;
    const char* target = "GGTTAGAAAT";
    unsigned int target_len = 10;
    int match_cost = 2;
    int mismatch_cost = -1;
    int gap_cost = -2;
    std::string cigar;
    unsigned int target_begin;
    
    int expected_value = 10;
    std::string expected_cigar = "4S5=";
    unsigned int expected_target_begin = 1;

    white::Aligner *aligner = new white::Aligner(query, query_len, target, target_len,
		match_cost, mismatch_cost, gap_cost, &cigar, &target_begin);
    
    int result = aligner -> Align(white::AlignmentType::SEMIGLOBAL);
    
    EXPECT_EQ(result, expected_value);
    EXPECT_EQ(cigar, expected_cigar);
    EXPECT_EQ(target_begin, expected_target_begin);
}
    