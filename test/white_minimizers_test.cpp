#include "../src/white_minimizers.hpp"
#include <gtest/gtest.h>


TEST(WhiteMinimizersTests, First) {
    const char* seq = "CGTTAC";

    std::vector<std::tuple<unsigned int, unsigned int, bool>> result = white::Minimize(seq, 6, 3, 3);
    
    std::vector<std::tuple<unsigned int, unsigned int, bool>> expected = {
    std::make_tuple(5,1,0),
    std::make_tuple(14,0,1),
    std::make_tuple(27,3,0)
    };

    EXPECT_EQ(result,expected);

}

TEST(WhiteMinimizersTests, Second) {

    const char* seq = "TACGTACCGTA";

    std::vector<std::tuple<unsigned int, unsigned int, bool>> result = white::Minimize(seq, 11, 3, 3);
    
    std::vector<std::tuple<unsigned int, unsigned int, bool>> expected = {
    std::make_tuple(6,3,0),
    std::make_tuple(6,8,0),
    std::make_tuple(16,5,1),
    std::make_tuple(27,0,0)
    };

    EXPECT_EQ(result,expected);

}