#include <gtest/gtest.h>
#include "../src/blonde_alignment.h"

void testAlignmentScenario(
    const char* target, 
    const char* query, 
    blonde::alignment::AlignmentType type, 
    int value_expect, 
    std::string cigar_expect, 
    int target_begin_expect, 
    int match = 1, 
    int mismatch = -1, 
    int gap = -1) {
    std::string cigar;
    unsigned int target_begin;

    int value = blonde::alignment::Align(query, 
        strlen(query), target, strlen(target), type, match, mismatch, gap, &cigar, &target_begin);

    EXPECT_EQ(value, value_expect);
    EXPECT_EQ(cigar, cigar_expect);
    EXPECT_EQ(target_begin, target_begin_expect);
}

TEST(BlondeAlignmentTests, Global) {
    testAlignmentScenario("GCATGCU", "GATTACA", blonde::alignment::AlignmentType::kGlobal, 0, "1=1D1=1I1=1X1=1X", 0);
    testAlignmentScenario("AGTCTTTGAT", "GATTAGA", blonde::alignment::AlignmentType::kGlobal, 0, "1D1=1D1X2=1X2=1D", 0);   
}

TEST(BlondeAlignmentTests, Local) {
    testAlignmentScenario("AGTCTTTGAT", "GATTAGA", blonde::alignment::AlignmentType::kLocal, 3, "2S2=1X2=", 4);
    testAlignmentScenario("AGTAGT", "TAG", blonde::alignment::AlignmentType::kLocal, 3, "3=", 2);
}