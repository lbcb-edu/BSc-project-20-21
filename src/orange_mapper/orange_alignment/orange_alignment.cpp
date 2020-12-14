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
        string Alignment::generateCigar(int maxRow, int maxColumn, string initialVal, std::vector<std::vector<Cell>> mainMatrix) {
            for (int i = maxRow; i >= 0; i--) {
                for (int j = maxColumn; j >= 0; j--) {
                    switch (mainMatrix[i][j].op) {
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
        void Alignment::calculateMax(int i, int j, const char* query,
            const char* target, orange::AlignmentType type,
            int match,
            int mismatch,
            int gap, std::vector<std::vector<Cell>> mainMatrix)
        {
            Cell minCell = { numeric_limits<int>::min(), sNone };
            int maxColumnIndex;
            int maxRowIndex;
            int maxSmithWaterman;
            if (type == local) {
                minCell.value = 0;
            }
            Cell left = { mainMatrix[i][j - 1].value + gap, sInsert };
            Cell up = { mainMatrix[i - 1][j].value + gap, sDelete };
            Cell diagonal;
            if (query[i] == target[j]) {
                diagonal = { mainMatrix[i - 1][j - 1].value + match, sMatch };
            }
            else {
                diagonal = { mainMatrix[i - 1][j - 1].value + mismatch, sMismatch };
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
            }
            mainMatrix[i][j] = std::max({ minCell, diagonal, up, left },
                [](Cell a, Cell b) { return a.value < b.value; });
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
                //Neeedleman-Wunsch algorithm
                for (int i = 1; i <= query_len; i++) {
                    mainMatrix[i][0] = { (mainMatrix[i][0].value + gap), sDelete };
                }
                for (int j = 1; j <= target_len; j++) {
                    mainMatrix[0][j] = { (mainMatrix[0][j].value + gap), sInsert };
                }
            }
            else if (type == local) {
                //Smith-Waterman algorithm
                for (int i = 1; i < query_len; i++) {
                    int val = mainMatrix[i][0].value + gap;
                    if (val < 0) {
                        val = 0;
                    }
                    mainMatrix[i][0] = { val, sDelete };
                }
                for (int i = 1; i < target_len; i++) {
                    int val = mainMatrix[0][i].value + gap;
                    if (val < 0) {
                        val = 0;
                    }
                    mainMatrix[0][i] = { val, sInsert };
                }
            }
            //matrix filling
            for (int i = 1; i <= query_len; i++) {
                for (int j = 1; j <= target_len; j++) {
                    calculateMax(i, j, query, target, type, match, mismatch, gap, mainMatrix);
                }
            }
            int maxValue = mainMatrix[query_len][target_len].value;
            int maxRow = query_len;
            int maxColumn = target_len;
            if (type == local) {
                int maxValue = 0;
                for (int i = 0; i <= query_len; i++) {
                    for (int j = 0; j <= target_len; j++) {
                        if (mainMatrix[i][j].value > maxValue) {
                            maxValue = mainMatrix[i][j].value;
                            maxRow = i;
                            maxColumn = j;
                        };
                    }
                }
            }
            else if (type == semiGlobal) {
                int maxValue = numeric_limits<int>::min();
                for (int i = 0; i <= query_len; i++) {
                    if (mainMatrix[i][target_len].value > maxValue) {
                        maxValue = mainMatrix[i][target_len].value;
                        maxRow = i;
                    }
                }
                for (int j = 0; j <= target_len; j++) {
                    if (mainMatrix[query_len][j].value > maxValue) {
                        maxValue = mainMatrix[query_len][j].value;
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
//};
}
