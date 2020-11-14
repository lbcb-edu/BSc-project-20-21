#include <iostream>
#include <math.h>
namespace white {

    enum AlignmentType {GLOBAL, LOCAL, SEMIGLOBAL}; //type of alignment
    enum direction {UP, LEFT, UPLEFT, NONE};

    struct cell {
        direction dir = NONE;
        int value;
    };

    int Align(const char* query, unsigned int query_len, const char* target, unsigned int target_len,
    AlignmentType type, int match, int mismatch, int gap, std::string* cigar = nullptr, unsigned int* target_begin = nullptr) {

        cell Mat[query_len + 1][target_len + 1]; //comparison matrix
        Mat[0][0].value = 0;
        Mat[0][0].dir = NONE;

        switch(type) {
            case GLOBAL:
                for (int i = 1; i < query_len; i++) {
                    Mat[0][i].value = i;
                    Mat[0][i].dir = LEFT;
                }
                for (int i = 1; i < target_len; i++) {
                    Mat[i][0].value = i;
                    Mat[i][0].dir = UP;
                }
                for (int i = 1; i < target_len; i++) {
                    for (int j = 1; j < query_len; j++) {
                        cell upleft = Mat[i-1][j-1];
                        cell up = Mat[i-1][j];
                        cell left = Mat[i][j-1];
                        int upleftVal = upleft.value;
                        if (query[j] == target[i]) {
                            upleftVal -= match;
                        }
                        else {
                            upleftVal += mismatch;
                        }
                        int min = upleftVal;
                        direction d = UPLEFT; 
                        int upVal = up.value + gap;
                        int leftVal = left.value + gap;
                        if (upVal < min) {
                            min = upVal;
                            d = UP;
                        }
                        if (leftVal < min) {
                            min = leftVal;
                            d = LEFT;
                        }
                        Mat[i][j].value = min;
                        Mat[i][j].dir = d;
                    }
                }
                return Mat[query_len][target_len].value;
            case LOCAL:
                for (int i = 1; i < query_len; i++) {
                    Mat[0][i].value = 0;
                    Mat[0][i].dir = NONE;
                }
                for (int i = 1; i < target_len; i++) {
                    Mat[i][0].value = 0;
                    Mat[i][0].dir = NONE;
                }
                std::pair <int, int> maximalCell(0, 0);
                int maxVal = 0;
                for (int i = 1; i < target_len; i++) {
                    for (int j = 1; j < query_len; j++) {
                        cell upleft = Mat[i-1][j-1];
                        cell up = Mat[i-1][j];
                        cell left = Mat[i][j-1];
                        int upleftVal = upleft.value;
                        if (query[j] == target[i]) {
                            upleftVal += match;
                        }
                        else {
                            upleftVal -= mismatch;
                        }
                        int max = upleftVal;
                        direction d = UPLEFT; 
                        int upVal = up.value - gap;
                        int leftVal = left.value - gap;
                        if (upVal > max) {
                            max = upVal;
                            d = UP;
                        }
                        if (leftVal > max) {
                            max = leftVal;
                            d = LEFT;
                        }
                        if (max < 0) {
                            max = 0;
                            d = NONE;
                        }
                        Mat[i][j].value = max;
                        Mat[i][j].dir = d;
                        if (max > maxVal) {
                            maxVal = max;
                            maximalCell.first = i;
                            maximalCell.second = j;
                        }
                    }
                }
                return maxVal;
            case SEMIGLOBAL:
                for (int i = 1; i < query_len; i++) {
                    Mat[0][i].value = 0;
                    Mat[0][i].dir = LEFT;
                }
                for (int i = 1; i < target_len; i++) {
                    Mat[i][0].value = 0;
                    Mat[i][0].dir = UP;
                }
                std::pair <int, int> maximalCell(0, 0);
                int maxVal = 0;
                for (int i = 1; i < target_len; i++) {
                    for (int j = 1; j < query_len; j++) {
                        cell upleft = Mat[i-1][j-1];
                        cell up = Mat[i-1][j];
                        cell left = Mat[i][j-1];
                        int upleftVal = upleft.value;
                        if (query[j] == target[i]) {
                            upleftVal += match;
                        }
                        else {
                            upleftVal -= mismatch;
                        }
                        int max = upleftVal;
                        direction d = UPLEFT; 
                        int upVal = up.value - gap;
                        int leftVal = left.value - gap;
                        if (upVal > max) {
                            max = upVal;
                            d = UP;
                        }
                        if (leftVal > max) {
                            max = leftVal;
                            d = LEFT;
                        }
                        Mat[i][j].value = max;
                        Mat[i][j].dir = d;
                    }
                } 
                for (int i = 0; i < target_len; i++) {
                    if (Mat[query_len][i].value > maxVal) {
                        maximalCell.first = query_len;
                        maximalCell.second = i;
                        maxVal = Mat[query_len][i].value;
                    }
                }
                for (int i = 0; i < query_len; i++) {
                    if (Mat[i][target_len].value > maxVal) {
                        maximalCell.first = target_len;
                        maximalCell.second = i;
                        maxVal = Mat[target_len][i].value;
                    }
                }
                return maxVal;
        }

    }
}

