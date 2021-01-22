// orange_alignment.cpp : library which implements all three alignment algorithms:
//the Needleman-Wunsch algorithm for global alignment, the Smith-Waterman algorithm for local alignment 
//and semi-global algorithms used for suffix-prefix and prefix-suffix alignments.
//
#include "orange_alignment.h"

#include<iostream>
#include<algorithm>
#include<vector>
#include<string>
using namespace std;
namespace orange
{
    const char* query;
    unsigned int query_len;
    const char* target;
    unsigned int target_len;
    AlignmentType type;
    int match;
    int mismatch;
    int gap;
    std::string* cigar = nullptr;
    unsigned int* target_begin = nullptr;
    std::vector<std::vector<Cell>> mainMatrix;
    //generates cigar string
    string generateCigar(int maxRow, int maxColumn, string initialVal, vector<std::vector<Cell>>& mainMatrix) {
        for (int i = maxRow; i >= 0; i--) {
            for (int j = maxColumn; j >= 0; j--) {
                switch (mainMatrix.at(i).at(j).op) {
                case sNone:
                    break;
                case sMatch:
                    initialVal = '=' + initialVal;
                    break;
                case sMismatch:
                    initialVal = 'X' + initialVal;
                    break;
                case sDelete:
                    initialVal = 'D' + initialVal;
                    j++;
                    break;
                case sInsert:
                    initialVal = 'I' + initialVal;
                    i++;
                    break;
                default:
                    break;
                }
            }
        }
        int counter = 1;
        char letter = initialVal[counter - 1];
        string result = "";
        for (int i = 0; i < initialVal.size(); i++) {
            if (letter == initialVal[i]) {
                counter++;
            }
            else {
                result += (char)counter + letter;
                letter = initialVal[i];
                counter = 1;
            }
        }
        result += (char)counter + letter;
        return result;
    }
    void calculateMax(int i, int j, const char* query,
        const char* target, orange::AlignmentType type,
        int match,
        int mismatch,
        int gap, std::vector<std::vector<Cell>>& mainMatrix)
    {
        Cell minCell = { numeric_limits<int>::min(), sNone };
        int maxColumnIndex;
        int maxRowIndex;
        int maxSmithWaterman;
        if (type == local) {
            minCell.value = 0;
        }
        Cell left = { mainMatrix.at(i).at(j - 1).value + gap, sInsert };
        Cell up = { mainMatrix.at(i - 1).at(j).value + gap, sDelete };
        Cell diagonal;
        if (query[i - 1] == target[j - 1]) {
            diagonal = { mainMatrix.at(i - 1).at(j - 1).value + match, sMatch };
        }
        else {
            diagonal = { mainMatrix.at(i - 1).at(j - 1).value + mismatch, sMismatch };
        }
        if (type == local) {
            if (left.value < 0) {
                left.value = 0;
            }
            if (up.value < 0) {
                up.value = 0;
            }
            if (diagonal.value < 0) {
                diagonal.value = 0;
            }
            mainMatrix.at(i).at(j) = std::max({ minCell, diagonal, up, left },
                [](Cell a, Cell b) { return a.value < b.value; });
        }
        else if (type == global) {
            mainMatrix.at(i).at(j) = std::min({ diagonal, up, left },
                [](Cell a, Cell b) { return a.value < b.value; });
        }
        else if (type == semiGlobal) {
            mainMatrix.at(i).at(j) = std::max({ minCell, diagonal, up, left },
                [](Cell a, Cell b) { return a.value < b.value; });
        }

    }
    int Alignment::Align(
        const char* query, unsigned int query_len,
        const char* target, unsigned int target_len,
        orange::AlignmentType type,
        int match,
        int mismatch,
        int gap,
        std::string* cigar,
        unsigned int* target_begin)
    {
        int alignmentScore;
        //intialization
        if (type == global) {
            for (int i = 1; i <= target_len; i++) {
                mainMatrix[0][i] = { (mainMatrix.at(0).at(i - 1).value + gap), sInsert };
            }
            for (int j = 1; j <= query_len; j++) {
                mainMatrix[j][0] = { (mainMatrix.at(j - 1).at(0).value + gap), sDelete };
            }
        }
        else {
            //Smith-Waterman algorithm or Needleman-wunsch
            for (int i = 1; i < query_len; i++) {
                int val = (mainMatrix.at(i).at(0).value + gap);
                if (val < 0) {
                    val = 0;
                }
                mainMatrix.at(i).at(0) = { val, sDelete };
            }
            for (int i = 1; i < target_len; i++) {
                int val = mainMatrix.at(0).at(i).value + gap;
                if (val < 0) {
                    val = 0;
                }
                mainMatrix.at(0).at(i) = { val, sInsert };
            }
        }
        //matrix filling
        for (int i = 1; i <= query_len; i++) {
            for (int j = 1; j <= target_len; j++) {
                calculateMax(i, j, query, target, type, match, mismatch, gap, mainMatrix);
            }
        }
        int maxValue = mainMatrix.at(query_len).at(target_len).value;
        int maxRow = query_len;
        int maxColumn = target_len;
        if (type == local) {
            maxValue = (mainMatrix.at(0).at(0).value);
            for (int i = 0; i <= query_len; i++) {
                for (int j = 0; j <= target_len; j++) {
                    if (mainMatrix.at(i).at(j).value > maxValue) {
                        maxValue = mainMatrix.at(i).at(j).value;
                        maxRow = i;
                        maxColumn = j;
                    };
                }
            }
        }
        else if (type == semiGlobal) {
            maxValue = numeric_limits<int>::min();
            for (int i = 0; i <= query_len; i++) {
                if (mainMatrix.at(i).at(target_len).value > maxValue) {
                    maxValue = mainMatrix.at(i).at(target_len).value;
                    maxRow = i;
                }
            }
            for (int j = 0; j <= target_len; j++) {
                if (mainMatrix.at(query_len).at(j).value > maxValue) {
                    maxValue = mainMatrix.at(query_len).at(j).value;
                    maxColumn = j;
                }
            }
        }
        if (cigar) {
            string initialVal = "";
            string result = "";
            switch (type) {
            case global:
                *cigar = generateCigar(maxRow, maxColumn, initialVal, mainMatrix);
                *target_begin = 1;
                break;
            case (local || semiGlobal):
                for (int a = maxRow + 1; a <= query_len; a++) {
                    initialVal = "S" + initialVal;
                }
                result = generateCigar(maxRow, maxColumn, initialVal, mainMatrix);
                if (maxRow != 0) {
                    result = to_string(maxRow) + "S" + result;
                }
                *target_begin = maxColumn + 1;
                *cigar = result;
                break;
            default:
                break;
            }
        }
        return maxValue;
    }
}
