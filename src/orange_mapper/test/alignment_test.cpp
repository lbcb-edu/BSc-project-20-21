#include "gtest/gtest.h"
#include "orange_alignment/orange_alignment.cpp"
#include <vector>
#include<string>
#include<iostream>
#include<algorithm>

using namespace std;
using namespace orange;
TEST(AlignmentTest, Global) {
	char query[] = { 'G', 'T', 'A', 'C', 'C' };
	unsigned int query_len = sizeof(query);
	char target[] = { 'G', 'A', 'T', 'A', 'C', 'G', 'T', 'T', 'A' };
	unsigned int target_len = sizeof(target);
	orange::AlignmentType type = orange::AlignmentType::global;
	int match = 0;
	int mismatch = 1;
	int gap = 1;
	string* cigar = nullptr;
	unsigned int* target_begin = nullptr;
	orange::Alignment* a = new orange::Alignment(query, sizeof(query), target, sizeof(target), type, match, mismatch, gap);
	int alignmentScore = a->Align(query, query_len, target, target_len, type, match, mismatch, gap, cigar, target_begin);
	EXPECT_EQ(alignmentScore, 5);
}
TEST(AlignmentTest, Local) {
	char query[] = { 'G', 'A', 'T','C', 'A', 'T', 'A', 'T', 'T' };
	char target[] = { 'T', 'C', 'G', 'T', 'A', 'G', 'C', 'G' };
	unsigned int query_len = sizeof(query);
	unsigned int target_len = sizeof(target);
	orange::AlignmentType type = orange::AlignmentType::local;
	int match = 2;
	int mismatch = -1;
	int gap = -2;
	string* cigar = nullptr;
	unsigned int* target_begin = nullptr;
	orange::Alignment* a = new orange::Alignment(query, sizeof(query), target, sizeof(target), type, match, mismatch, gap, cigar, target_begin);
	int alignmentScore = a->Align(query, query_len, target, target_len, type, match, mismatch, gap, cigar, target_begin);
	EXPECT_EQ(alignmentScore, 7);
}
TEST(AlignmentTest, SemiGlobal) {
	char query[] = { 'G', 'T', 'A', 'C', 'C' };
	char target[] = { 'T', 'T', 'C', 'A', 'C', 'G', 'T', 'T', 'A' };
	unsigned int query_len = sizeof(query);
	unsigned int target_len = sizeof(target);
	orange::AlignmentType type = orange::AlignmentType::semiGlobal;
	int match = 1;
	int mismatch = -1;
	int gap = -1;
	string* cigar = nullptr;
	unsigned int* target_begin = nullptr;
	orange::Alignment* a = new orange::Alignment(query, sizeof(query), target, sizeof(target), type, match, mismatch, gap, cigar, target_begin);
	int alignmentScore = a->Align(query, query_len, target, target_len, type, match, mismatch, gap, cigar, target_begin);
	EXPECT_EQ(alignmentScore, 2);
}
int main(int argc, char** argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

