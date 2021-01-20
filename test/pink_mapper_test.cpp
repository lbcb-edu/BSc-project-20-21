#include <iostream>
#include <gtest/gtest.h>
#include "pink_alignment.hpp"
#include "pink_minimizers.hpp"
#include <string>
#include <vector>
#include <tuple>

TEST(Alignment_test, Global_test)
{
    //Arrange
    char query[] = {'G', 'T', 'A', 'C', 'C'};
    char target[] = {'G', 'A', 'T', 'A', 'C', 'G', 'T', 'T', 'A'};
    int match = 0, mismatch = 1, gap = 1;
    string cigar_cmp = "1=1I3=1X3I";
    string cigar;
    pink::AlignmentType type = pink::AlignmentType::global;
    //Act
    int result = pink::Align(query, sizeof(query), target, sizeof(target), type, match, mismatch, gap, &cigar);
    //Assert
    EXPECT_EQ(result, 5);
    //EXPECT_EQ(cigar, cigar_cmp);
}

TEST(Alignment_test, Local_test)
{
    //Arrange
    char query[] = {'G', 'A', 'T', 'C', 'A', 'T', 'A', 'T', 'T'};
    char target[] = {'T', 'C', 'G', 'T', 'A', 'G', 'C', 'G'};
    int match = 2, mismatch = -1, gap = -2;
    pink::AlignmentType type = pink::AlignmentType::local;
    //Act
    int result = pink::Align(query, sizeof(query), target, sizeof(target), type, match, mismatch, gap);
    //Assert
    EXPECT_EQ(result, 7);
}

TEST(Alignment_test, Semiglobal_test)
{
    //Arrange
    char query[] = {'G', 'T', 'A', 'C', 'C'};
    char target[] = {'T', 'T', 'C', 'A', 'C', 'G', 'T', 'T', 'A'};
    int match = 1, mismatch = -1, gap = -1;
    pink::AlignmentType type = pink::AlignmentType::semiglobal;
    //Act
    int result = pink::Align(query, sizeof(query), target, sizeof(target), type, match, mismatch, gap);
    //Assert
    EXPECT_EQ(result, 2);
}

TEST(Minimizers_test, First_minimizers_test)
{
    //Arrange
    char sequence[] = {'T', 'G', 'A', 'C', 'G', 'T', 'A', 'C', 'A', 'T', 'G', 'G', 'A', 'C', 'A'};
    vector<tuple<unsigned int, unsigned int, bool>>
        min_cmp{make_tuple(1, 3, 0), make_tuple(6, 4, 0), make_tuple(6, 7, 1), make_tuple(14, 10, 0), make_tuple(17, 12, 1), make_tuple(19, 2, 1), make_tuple(45, 0, 1)};
    //Act
    auto mins = pink::Minimize(sequence, sizeof(sequence), 3, 3);
    //Assert
    EXPECT_EQ(mins, min_cmp);
}

TEST(Minimizers_test, Second_minimizers_test)
{
    //Arrange
    char sequence[] = {'T', 'A', 'G', 'G', 'A', 'A', 'C', 'T', 'G', 'A', 'T', 'T', 'G', 'C', 'C', 'T', 'T', 'A'};
    vector<tuple<unsigned int, unsigned int, bool>>
        min_cmp{make_tuple(9, 14, 0), make_tuple(10, 13, 1), make_tuple(22, 5, 0), make_tuple(30, 9, 0), make_tuple(45, 6, 1), make_tuple(82, 4, 1), make_tuple(95, 1, 0), make_tuple(122, 10, 0), make_tuple(151, 0, 0)};
    //Act
    auto mins = pink::Minimize(sequence, sizeof(sequence), 4, 3);
    //Assert
    EXPECT_EQ(mins, min_cmp);
}

int main(int argc, char **argv)
{

    testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
