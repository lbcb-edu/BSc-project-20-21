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
                m[0][0] = 0;
                int maxCell = 0;
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
                        if (m[i][j] > maxCell && (i==query_len || j==target_len)) {
                            maxCell = m[i][j];
                            resultRow = i;
                            resultColumn = j;
                        }
                    }
            }
            else {
                std::cerr << "Invalid AlignmentType" << std::endl;
                exit(EXIT_FAILURE);
            }

            if (target_begin != nullptr) {
                //TODO
            }

            if (cigar != nullptr) {
                std::string cigarBeta="";
                    while(type==GLOBAL && (resultRow != 0 && resultColumn != 0) || (type==LOCAL || type==SEMIGLOBAL) && (resultRow != 0 || resultColumn != 0)) {
                        if(!(resultRow==0 || resultColumn==0)) {
                            if (m[resultRow-1][resultColumn-1]+match==m[resultRow][resultColumn]) {
                                cigarBeta+="M";
                                resultColumn--;
                                resultRow--;
                            }
                            else if (m[resultRow-1][resultColumn-1]+mismatch==m[resultRow][resultColumn]) {
                                cigarBeta+="X";
                                resultColumn--;
                                resultRow--;
                            }
                        }
                        else if (resultRow!=0 && m[resultRow-1][resultColumn]+gap==m[resultRow][resultColumn]) {
                            cigarBeta+="I";
                            resultRow--;
                        }
                        else if(resultColumn!=0) {
                            cigarBeta="D";
                            resultColumn--;
                        }

                    }

                *cigar="";
                //cigar->append(cigarBeta.substr(0,1));
                char current=cigarBeta.at(0);
                int counter=0;
                while(!cigarBeta.empty()) {
                    if(cigarBeta.at(0) == current) {
                        counter++;
                        cigarBeta=cigarBeta.substr(1);

                    }
                    else {
                        cigar->append(std::string (1, current));
                        cigar->append(std::to_string(counter));
                        counter=0;
                        current=cigarBeta.at(0);
                    }

                }
            }

            return m[resultRow][resultColumn];
    }
}
