#include "../orange_mapper.cpp"
#include "gtest/gtest.h"
#include <string>
#include <zlib.h>

using namespace std;

//testing alignment
TEST(OrangeAlignmentTests, Global) {
    char* query = "GATTACA";
    unsigned int query_len = 7;
    char* target = "GCATGCU";
    unsigned int target_len = 7;
    int match = 1;
    int mismatch = -1;
    int gap = -1;
    std::string cigar;
    unsigned int target_begin;

    int expected_value = 0;
    std::string expected_cigar = "1=1I1=1D1=1X1=1X";
    unsigned int expected_target_begin = 1;
orange::Alignment classAlign (query, query_len, target, target_len,  orange::AlignmentType::global, match, mismatch, gap, &cigar, &target_begin);
    int value = classAlign.Align(query, query_len, target, target_len, orange::AlignmentType::global, match, mismatch, gap, &cigar, &target_begin);
    EXPECT_EQ(value, expected_value);
    //EXPECT_EQ(cigar, expected_cigar);
    //EXPECT_EQ(target_begin, expected_target_begin);

}

TEST(OrangeAlignmentTests, Local) {
    char* query = "ACCTAAGG";
    unsigned int query_len = 8;
    char* target = "GGCTCAATCA";
    unsigned int target_len = 10;
    int match = 2;
    int mismatch = -1;
    int gap = -2;
    std::string cigar;
    unsigned int target_begin;

    int expected_value = 6;
    std::string expected_cigar = "2S2=1I2=2S";
    unsigned int expected_target_begin = 3;
	orange::Alignment classAlign (query, query_len, target, target_len,  orange::AlignmentType::local, match, mismatch, gap, &cigar, 	&target_begin);	
    int value = classAlign.Align(query, query_len, target, target_len, orange::AlignmentType::local, match, mismatch, gap, &cigar, &target_begin);
    EXPECT_EQ(value, expected_value);
   // EXPECT_EQ(cigar, expected_cigar);
    //EXPECT_EQ(target_begin, expected_target_begin);
}

TEST(OrangeAlignmentTests, SemiGlobal) {
     char* query = "CCTCGGTTA";
    unsigned int query_len = 9;
     char* target = "GGTTAGAAAT";
    unsigned int target_len = 10;
    int match = 2;
    int mismatch = -1;
    int gap = -2;
    std::string cigar;
    unsigned int target_begin;

    int expected_value = 10;
    std::string expected_cigar = "4S5=";
    unsigned int expected_target_begin = 1;
orange::Alignment classAlign (query, query_len, target, target_len,  orange::AlignmentType::semiGlobal, match, mismatch, gap, &cigar, &target_begin);
    int value = classAlign.Align(query, query_len, target, target_len, orange::AlignmentType::semiGlobal, match, mismatch, gap, &cigar, &target_begin);
    EXPECT_EQ(value, expected_value);
   // EXPECT_EQ(cigar, expected_cigar);
   // EXPECT_EQ(target_begin, expected_target_begin);
}
