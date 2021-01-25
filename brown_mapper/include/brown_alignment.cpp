#include <iostream>
#include "brown_alignment.hpp"
#include <algorithm>
#include <math.h>
namespace brown {

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

            int **m = new int*[query_len + 1];
            for (unsigned int i = 0; i < query_len + 1; i++)
                m[i] = new int[target_len + 1];

            AlignmentDirection** direction = new AlignmentDirection*[query_len + 1];
            for (unsigned int i = 0; i < query_len + 1; i++)
                direction[i] = new AlignmentDirection[target_len + 1];
            
            if(type == GLOBAL) {
                m[0][0]=0;
                direction[0][0] = NONE;

                for (int i = 1; i < query_len + 1; i++) {
                    m[i][0] = i * gap;
                    direction[i][0] = DELETION;
                }
                for (int j = 1; j < target_len + 1; j++) {
                    m[0][j] = j * gap;
                    direction[0][j] = INSERTION;
                }
                
                for (int i = 1; i < query_len+1; i++)
                    for (int j = 1; j < target_len+1; j++) {
                        int matchCost;
                        if (query[i-1] == target[j-1]) matchCost=m[i-1][j-1] + match;
                        else matchCost = m[i-1][j-1] + mismatch;
                        m[i][j]=std::max(std::max(matchCost, m[i][j-1] + gap), m[i-1][j] + gap);
                        
                        if (m[i][j] == matchCost) {
                            if (matchCost == m[i-1][j-1] + match)
                                direction[i][j] = MATCH;
                            else 
                                direction[i][j] = MISMATCH;
                        } 
                        else if (m[i][j] == m[i][j-1] + gap){
                            direction[i][j] = INSERTION;
                        } 
                        else if (m[i][j] == m[i-1][j] + gap){
                            direction[i][j] = DELETION;
                        }
                        
                    }
                resultRow = query_len;
                resultColumn = target_len;
                //std::cout << "odradio je global align\n";
            }
            else if (type == LOCAL) {
                int maxCell = INT32_MIN;
                m[0][0] = 0;
                direction[0][0] = NONE;
                
                for (int i = 1; i < query_len + 1; i++) {
                    m[i][0] = 0;
                    direction[i][0] = NONE;
                }
                for (int j = 1; j < target_len + 1; j++) {
                    m[0][j] = 0;
                    direction[0][j] = NONE;
                }

                for (int i = 1; i < query_len + 1; i++)
                    for (int j = 1; j < target_len + 1; j++) {
                        int matchCost;
                        if (query[i - 1] == target[j - 1]) matchCost = m[i-1][j-1]+match;
                        else matchCost = m[i-1][j-1] + mismatch;
                        m[i][j] = std::max(std::max(0, matchCost), std::max(m[i][j-1] + gap, m[i-1][j] + gap));
                        if (m[i][j] > maxCell) {
                            maxCell = m[i][j];
                            resultRow = i;
                            resultColumn = j;
                        }

                        if (m[i][j] == 0) {
                            direction[i][j] = NONE;
                        }
                        else if (m[i][j] == matchCost) {
                            if (matchCost == m[i-1][j-1] + match)
                                direction[i][j] = MATCH;
                            else 
                                direction[i][j] = MISMATCH;
                        } 
                        else if (m[i][j] == m[i][j-1] + gap){
                            direction[i][j] = INSERTION;
                        } 
                        else if (m[i][j] == m[i-1][j] + gap){
                            direction[i][j] = DELETION;
                        }
                    }
            }
            else if (type == SEMIGLOBAL){
                m[0][0] = 0;
                direction[0][0] = NONE;
                int maxCell = INT32_MIN;

                for (int i = 1; i < query_len + 1; i++) {
                    m[i][0] = 0;
                    direction[i][0] = NONE;
                }
                for (int j = 1; j < target_len + 1; j++) {
                    m[0][j] = 0;
                    direction[0][j] = NONE;
                }

                for (int i = 1; i < query_len + 1; i++)
                    for (int j = 1; j < target_len + 1; j++) {
                        int matchCost;
                        if (query[i - 1] == target[j - 1]) matchCost = m[i-1][j-1]+match;
                        else matchCost = m[i-1][j-1] + mismatch;
                        m[i][j]=std::max(matchCost, std::max(m[i][j-1] + gap, m[i-1][j] + gap));

                        if (m[i][j] > maxCell && (i == query_len || j == target_len)) {
                            maxCell = m[i][j];
                            resultRow = i;
                            resultColumn = j;
                        }

                        if (m[i][j] == matchCost) {
                            if (matchCost == m[i-1][j-1] + match)
                                direction[i][j] = MATCH;
                            else 
                                direction[i][j] = MISMATCH;
                        } 
                        else if (m[i][j] == m[i][j-1] + gap){
                            direction[i][j] = INSERTION;
                        } 
                        else if (m[i][j] == m[i-1][j] + gap){
                            direction[i][j] = DELETION;
                        }
                    }
            }
            else {
                std::cerr << "Invalid AlignmentType" << std::endl;
                exit(EXIT_FAILURE);
            }
            int returnRow = resultRow;
            int returnColumn = resultColumn;
            
            std::string cigarBeta = "";
            //int targetLocal;
            if (cigar != nullptr || target_begin !=nullptr) {
                
                    while(direction[resultRow][resultColumn] != NONE) {
                    
                        switch (direction[resultRow][resultColumn]) {
                            case MATCH:
                                cigarBeta += "M";
                                resultColumn--;
                                resultRow--;
                                break;
                            case MISMATCH :
                                cigarBeta += "X";
                                resultRow--;
                                resultColumn--;
                                break;
                            case INSERTION :
                                cigarBeta += "I";
                                resultColumn--;
                                break;
                            case DELETION :
                                cigarBeta += "D";
                                resultRow--;
                                break;
                            default:
                                std::cerr << "Direction not allowed" << std::endl;
                                exit(EXIT_FAILURE);
                        }
                        
                    }
            }

            if(target_begin != nullptr) { 
                    if (type == GLOBAL) *target_begin = 0;
                    else if (type == SEMIGLOBAL) *target_begin = resultColumn;
                    else *target_begin = resultColumn;
            }

            if(cigar != nullptr) {
                *cigar = "";
                char current = cigarBeta.at(0);
                int counter = 0;
                while(!cigarBeta.empty()) {
                    if(cigarBeta.at(0) == current) {
                        counter++;
                        cigarBeta = cigarBeta.substr(1);

                    }
                    else {
                        cigar->append(std::string (1, current));
                        cigar->append(std::to_string(counter));
                        counter = 0;
                        current = cigarBeta.at(0);
                    }
                }
                cigar->append(std::string (1, current));
                cigar->append(std::to_string(counter));
                std::reverse((*cigar).begin(), (*cigar).end());
            }

            int result = m[returnRow][returnColumn];
            for (int i = query_len; i >= 0; i--) {
                delete[] m[i];
                delete direction[i];
            }

            delete[] direction;
            delete[] m;

            return result;
    }
    
}
