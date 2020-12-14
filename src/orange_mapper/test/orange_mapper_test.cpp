#include "orange_mapper/orange_mapper.cpp"
#include "gtest/gtest.h"
#include <string>

using namespace std;

//testing alignment
TEST(AlignmentTest, Global){
	const char* query = {'G', 'T', 'A', 'C', 'C'};
	unsigned int query_len = sizeof(query);
        const char* target = {'G', 'A', 'T', 'A', 'C', 'G', 'T', 'T', 'A'};
	unsigned int target_len = sizeof(target);
        orange::Alignment::AlignmentType type = orange::Alignment::AlignmentType::global;
        int match = 0;
        int mismatch = 1;
        int gap = 1;
	string* cigar = nullptr;
	unsigned int* target_begin = nullptr;
	int alignmentScore = orange::Alignment::Align(query, query_len, target, target_len, type, match, mismatch, gap, cigar, target_begin);
	EXPECT_EQ(alignmentScore, 5);
}
TEST(AlignmentTest, Local){
	const char* query = {'G', 'A', 'T', 'C', 'A', 'T', 'A', 'T', 'T'};
	unsigned int query_len = sizeof(query);
        const char* target = {'T', 'C', 'G', 'T', 'A', 'G', 'C', 'G'};
	unsigned int target_len = sizeof(target);
        orange::Alignment::AlignmentType type = orange::Alignment::AlignmentType::local;
        int match = 2;
        int mismatch = -1;
        int gap = -2;
	string* cigar = nullptr;
	unsigned int* target_begin = nullptr;
	int alignmentScore = orange::Alignment::Align(query, query_len, target, target_len, type, match, mismatch, gap, cigar, target_begin);
	EXPECT_EQ(alignmentScore, 7);
}
TEST(AlignmentTest, SemiGlobal){
	const char* query = {'G', 'T', 'A', 'C', 'C'};
	unsigned int query_len = sizeof(query);
        const char* target = {'T', 'T', 'C', 'A', 'C', 'G', 'T', 'T', 'A'};
	unsigned int target_len = sizeof(target);
        orange::Alignment::AlignmentType type = orange::Alignment::AlignmentType::semiGlobal;
        int match = 1;
        int mismatch = -1;
        int gap = -1;
	string* cigar = nullptr;
	unsigned int* target_begin = nullptr;
	int alignmentScore = orange::Alignment::Align(query, query_len, target, target_len, type, match, mismatch, gap, cigar, target_begin);
	EXPECT_EQ(alignmentScore, 2);
}
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
