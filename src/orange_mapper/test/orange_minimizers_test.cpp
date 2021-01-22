#include "../orange_mapper.cpp"
#include "../orange_minimizer/orange_minimizer.h"
#include "gtest/gtest.h"
#include <string>

using Kmer = std::tuple<unsigned int, unsigned int, bool>;
using namespace std;

TEST(MinimizersTest, First){
    	const char* sequence = "ACGTGAC";
    	unsigned int s = 7;
    	unsigned int k = 3;
    	unsigned int w = 3;
	int i = 0;
	orange::Minimizer *m = new orange::Minimizer(sequence, s, k, w);
    	vector<Kmer> minimizers = m->Minimize(sequence, s, k, w);  
    	int expected_values[4][3] = {{19,0,1},{14,1,1},{4,2,0},{14,4,0}};

	sort(minimizers.begin(),minimizers.end());
	
	
	for (std::vector<Kmer>::iterator it = minimizers.begin() ; it != minimizers.end(); ++it, i++) {
        EXPECT_EQ(expected_values[i][0], std::get<0>(*it));
        EXPECT_EQ(expected_values[i][1], std::get<1>(*it));
        EXPECT_EQ(expected_values[i][2], std::get<2>(*it));
    }	      
}

