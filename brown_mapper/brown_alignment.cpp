#include <iostream>
#include "brown_alignment.hpp"
#include <algorithm>
#include <math.h>
namespace brown {

    //enum AlignmentType {GLOBAL, LOCAL, SEMIGLOBAL};

    //int bzvz() {
    //    return 1;
    //}

    int Align(
        const char* query, unsigned int query_len,
        const char* target, unsigned int target_len,
        AlignmentType type,
        int match,
        int mismatch,
        int gap,
        std::string* cigar = nullptr,
        unsigned int* target_begin = nullptr) {
            
            int resultRow = 0;
            int resultColumn = 0;
            int m[query_len + 1][target_len + 1];

            if(type == GLOBAL) {
                m[0][0]=0;
                for (int i = 1; i < query_len + 1; i++) {
                    m[i][0]=i * gap;
                }
                for (int j = 1; j < target_len + 1; j++) {
                    m[0][j]=j * gap;
                }
                for (int i = 1; i < query_len+1; i++)
                    for (int j = 1; j < target_len+1; j++) {
                        int matchCost;
                        if (query[i] == target[j]) matchCost=m[i-1][j-1] + match;
                        else matchCost=m[i-1][j-1] + mismatch;
                        m[i][j]=std::max(std::max(matchCost, m[i][j-1] + gap), m[i-1][j] + gap);
                    }
                resultRow = query_len;
                resultColumn = target_len;

            }
            else if (type == LOCAL) {
                int maxCell = 0;
                m[0][0] = 0;
                resultColumn = 0;
                resultRow = 0;
                for (int i = 1; i < query_len + 1; i++) {
                    m[i][0] = 0;
                }
                for (int j = 1; j < target_len + 1; j++) {
                    m[0][j] = 0;
                }
                for (int i = 1; i < query_len + 1; i++)
                    for (int j = 1; j < target_len + 1; j++) {
                        int matchCost;
                        if (query[i]==target[j]) matchCost = m[i-1][j-1]+match;
                        else matchCost=m[i-1][j-1] + mismatch;
                        m[i][j]=std::max(std::max(0, matchCost), std::max(m[i][j-1] + gap, m[i-1][j] + gap));
                        //maxCell=std::max(maxCell, m[i][j]);
                        if (m[i][j] > maxCell) {
                            maxCell = m[i][j];
                            resultRow = i;
                            resultColumn = j;
                        }
                    }
            }
            else if (type == SEMIGLOBAL){
                //TODO
            }
            else {
                std::cerr << "Invalid AlignmentType" << std::endl;
                exit(EXIT_FAILURE);
            }

            if (target_begin != nullptr) {
                //TODO
            }

            if (cigar != nullptr) {
                //TODO
            }

            return m[resultRow][resultColumn];
    }
}
