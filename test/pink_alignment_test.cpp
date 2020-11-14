#include <iostream>
#include <gtest/gtest.h>
#include "pink_alignment.hpp"



TEST(ALignment_test, Global_test){
    //Arrange
    char query[] = {'G', 'T', 'A','C', 'C'};
    char target[] = {'G', 'A', 'T', 'A', 'C', 'G', 'T', 'T', 'A'};
    int match = 0, mismatch = 1, gap = 1;
    pink::AlignmentType type = pink::AlignmentType::global;
    //Act
    int result = pink::Align(query, sizeof(query), target, sizeof(target), type, match, mismatch, gap);
    //Assert
    EXPECT_EQ(result,5);
}

TEST(ALignment_test, Local_test){
    //Arrange
    char query[] = {'G', 'A', 'T','C', 'A', 'T', 'A', 'T', 'T'};
    char target[] = {'T', 'C', 'G', 'T', 'A', 'G', 'C', 'G'};
    int match = 2, mismatch = -1, gap = -2;
    pink::AlignmentType type = pink::AlignmentType::local;
    //Act
    int result = pink::Align(query, sizeof(query), target, sizeof(target), type, match, mismatch, gap);
    //Assert
    EXPECT_EQ(result,7);
}

TEST(ALignment_test, Semiglobal_test){
    //Arrange
    char query[] = {'G', 'T', 'A','C', 'C'};
    char target[] = {'T', 'T', 'C', 'A', 'C', 'G', 'T', 'T', 'A'};
    int match = 1, mismatch = -1, gap = -1;
    pink::AlignmentType type = pink::AlignmentType::semiglobal;
    //Act
    int result = pink::Align(query, sizeof(query), target, sizeof(target), type, match, mismatch, gap);
    //Assert
    EXPECT_EQ(result,2);
}

int main(int argc, char **argv){

    testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
