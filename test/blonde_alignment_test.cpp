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
    blonde::alignment::AlignmentType global = blonde::alignment::AlignmentType::kGlobal;
    testAlignmentScenario("GCATGCU", "GATTACA", global, 0, "1=1D1=1I1=1X1=1X", 0);
    testAlignmentScenario("AGTCTTTGAT", "GATTAGA", global, 0, "1D1=1D1X2=1X2=1D", 0);   
}

TEST(BlondeAlignmentTests, Local) {
    blonde::alignment::AlignmentType local = blonde::alignment::AlignmentType::kLocal;
    testAlignmentScenario("AGTCTTTGAT", "GATTAGA", local, 3, "2S2=1X2=", 4);
    testAlignmentScenario("AGTAGT", "TAGA", local, 3, "3=1S", 2);
}

TEST(BlondeAlignmentTests, SemiGlobal) {
    blonde::alignment::AlignmentType semi = blonde::alignment::AlignmentType::kSemiGlobal;
    testAlignmentScenario("AGTCTTTGAT", "GATTAGA", semi, 3, "3=4S", 7);
    testAlignmentScenario("AGTAGT", "TAGA", semi, 2, "3=1X", 2);
}