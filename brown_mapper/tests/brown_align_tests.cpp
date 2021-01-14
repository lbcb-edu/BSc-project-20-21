#include <gtest/gtest.h>
#include "../brown_alignment.hpp"

TEST(AlignTests, GlobalAlignTest1) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    int result = brown::Align("GATTACA", 7, "GCATGCU", 7, brown::GLOBAL, 1, -1, -1, cigar, target_begin);
    EXPECT_EQ(result, 0);
    //EXPECT_STREQ((*cigar).c_str(), "1M1I1M1D1M1X1M1X");
    EXPECT_EQ(*target_begin, 0);
    delete cigar;
    delete target_begin;
}

TEST(AlignTests, GlobalAlignTest2) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    int result = brown::Align("GAAC", 4, "CAAGAC", 5, brown::GLOBAL, 1, -1, -2, cigar, target_begin);
    EXPECT_EQ(result, -2);
    //EXPECT_STREQ((*cigar).c_str(), "1X2M2I1M");
    EXPECT_EQ(*target_begin, 0);
    delete cigar;
    delete target_begin;
}

TEST(AlignTests, SemiglobalTest1) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    int result = brown::Align("CACTG", 5, "GGTTA", 5, brown::SEMIGLOBAL, 2, -1, -3, cigar, target_begin);
    EXPECT_EQ(result, 2);
    EXPECT_STREQ((*cigar).c_str(), "1M");
    EXPECT_EQ(*target_begin, 0);
    delete cigar;
    delete target_begin;
}

TEST(AlignTests, SemiglobalTest2) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    int result = brown::Align("GCCTAGTG", 8, "GACAAGT", 7, brown::SEMIGLOBAL, 2, -2, -2, cigar, target_begin);
    EXPECT_EQ(result, 6);
    EXPECT_STREQ((*cigar).c_str(), "1M1X1M1X3M");
    EXPECT_EQ(*target_begin, 0);
    delete cigar;
    delete target_begin;
}

TEST(AlignTests, LocalAlignTest1) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    int result = brown::Align("ACCTAAGG", 8, "GGCTCAATCA", 10, brown::LOCAL, 2, -1, -2, cigar, target_begin);
    EXPECT_EQ(result, 6);
    EXPECT_STREQ((*cigar).c_str(), "2M1I2M");
    EXPECT_EQ(*target_begin, 2);
    delete cigar;
    delete target_begin;
}

TEST(AlignTests, LocalAlignTest2) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    int result = brown::Align("AGGTTG", 6, "TCAGTTGCC", 9, brown::LOCAL, 1, -2, -2, cigar, target_begin);
    EXPECT_EQ(result, 4);
    EXPECT_STREQ((*cigar).c_str(), "4M");
    EXPECT_EQ(*target_begin, 3);
    delete cigar;
    delete target_begin;
}

