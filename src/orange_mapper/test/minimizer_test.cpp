#include "orange_minimizer/orange_minimizer.h"
#include "gtest/gtest.h"
#include <vector>
#include<string>
#include<iostream>
#include<algorithm>

using namespace std;
using namespace orange;
TEST(MinimizerTest, testOne) {
	char sequence[] = { 'T', 'G', 'A', 'C', 'G', 'T', 'A', 'C', 'A', 'T', 'G', 'G', 'A', 'C', 'A' };
	std::vector<tuple<unsigned int, unsigned int, bool>> mins;
	orange::Minimizer* m = new orange::Minimizer(sequence, sizeof(sequence), 5, 3);
	mins = m->Minimize(sequence, sizeof(sequence), 3, 3);
	unsigned int expectedResult[7][3] = { {2,10,0},{6,4,0},{6,7,1},{11,1,0},{11,11,0},{17,12,1},{18,0,0} };
	for (int j = 0; j <= sizeof(mins); j++) {
		EXPECT_EQ(std::get<0>(mins.at(j)), expectedResult[j][0]);
		EXPECT_EQ(std::get<1>(mins.at(j)), expectedResult[j][1]);
		EXPECT_EQ(std::get<2>(mins.at(j)), expectedResult[j][2]);
	}
}
int main(int argc, char** argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}