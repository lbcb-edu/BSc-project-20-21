#include "../src/white_minimizers.hpp"
#include <gtest/gtest.h>


TEST(WhiteMinimizersTests, First) {
    const char* seq = "CGTTAC";

    std::vector<std::tuple<unsigned int, unsigned int, bool>> result = white::Minimize(seq, 6, 3, 3);
    
    std::vector<std::tuple<unsigned int, unsigned int, bool>> expected = {
    std::make_tuple(19,4,0),
    std::make_tuple(14,0,1),
    std::make_tuple(36,3,1)
    };

    EXPECT_EQ(result,expected);

}

TEST(WhiteMinimizersTests, Second) {

    const char* seq = "TACGTACCGTA";

    std::vector<std::tuple<unsigned int, unsigned int, bool>> result = white::Minimize(seq, 11, 3, 3);
    
    std::vector<std::tuple<unsigned int, unsigned int, bool>> expected = {
	std::make_tuple(14,7,1),
    std::make_tuple(17,1,1),
    std::make_tuple(14,2,1),
    std::make_tuple(36,0,1),
    std::make_tuple(16,5,1),
    std::make_tuple(57,8,1),
    std::make_tuple(3,6,1)
    };

    EXPECT_EQ(result,expected);

}