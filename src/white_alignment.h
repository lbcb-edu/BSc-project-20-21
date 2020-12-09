#include <iostream>
#include <math.h>
#include <vector>
#include <limits>

namespace white {

    enum AlignmentType {GLOBAL, LOCAL, SEMIGLOBAL}; //type of alignment
    enum Operation {kMatch = 0, kMismatch, kDelete, kInsert, kNone};

    struct Cell {
        Operation operation;
        int value;
    };

    class Aligner {

        public:
            const char* query;
            unsigned int query_len;
            const char* target;
            unsigned int target_len;
            int match;
            int mismatch;
            int gap;
            std::string* cigar;
            unsigned int* target_begin;
            std::vector<std::vector<Cell>> Mat;

        
        Aligner(const char* query, unsigned int query_len, const char* target,
          unsigned int target_len, int match, int mismatch, int gap,
          std::string* cigar, unsigned int* target_begin)
            :   query(query),
                query_len(query_len),
                target(target),
                target_len(target_len),
                match(match),
                mismatch(mismatch),
                gap(gap),
                cigar(cigar),
                target_begin(target_begin),
                Mat(query_len + 1,
                    std::vector<Cell>(target_len + 1, {Operation::kNone, 0})) {}

        int NeedlemanWunsch() {
            //std::cout << "checkNeedle1\n";
            for (int i = 1; i < query_len; i++) {
                Mat[i][0].value = i;
                Mat[i][0].operation = kDelete;
            }
            //std::cout << "checkNeedle2\n";
            for (int i = 1; i < target_len; i++) {
                Mat[0][i].value = i;
                Mat[0][i].operation = kInsert;
            }
            //std::cout << "checkNeedle3\n";
            for (int i = 1; i < query_len; i++) {
                for (int j = 1; j < target_len; j++) {
                    Cell upleft = Mat[i-1][j-1];
                    Cell up = Mat[i-1][j];
                    Cell left = Mat[i][j-1];
                    Operation op;
                    int upleftVal = upleft.value;
                    if (query[i] == target[j]) {
                        upleftVal += match;
                        op = kMatch;
                    }
                    else {
                        upleftVal += mismatch;
                        op = kMismatch;
                    }
                    int min = upleftVal;                       
                    int upVal = up.value + gap;
                    int leftVal = left.value + gap;
                    if (upVal < min) {
                        min = upVal;
                        op = kDelete;
                    }
                    if (leftVal < min) {
                        min = leftVal;
                        op = kInsert;
                    }
                    Mat[i][j].value = min;
                    Mat[i][j].operation = op;
                }
                //std::cout << "checkNeedle5\n";
            }
            //std::cout << "check3\n";
            if (cigar) {
                //std::cout << "checkexist\n";
                std::string starting_cigar = "";
                int row = query_len;
                int col = target_len;
                *cigar = CigarBuilder(row, col, starting_cigar);
                //std::cout << "checkcigarassign\n";
                *target_begin = 1;
                //std::cout << "check4\n";
            }

            return Mat[query_len][target_len].value;
        }

        int SmithWaterman() {
            for (int i = 1; i < query_len; i++) {
                    Mat[0][i].value = 0;
                    Mat[0][i].operation = kNone;
                }
                for (int i = 1; i < target_len; i++) {
                    Mat[i][0].value = 0;
                    Mat[i][0].operation = kNone;
                }

                std::pair <int, int> maximalCell(0, 0);
                int maxVal = 0;

                for (int i = 1; i < target_len; i++) {
                    for (int j = 1; j < query_len; j++) {
                        Cell upleft = Mat[i-1][j-1];
                        Cell up = Mat[i-1][j];
                        Cell left = Mat[i][j-1];
                        Operation op; 

                        int upleftVal = upleft.value;
                        if (query[j] == target[i]) {
                            upleftVal += match;
                            op = kMatch;
                        }
                        else {
                            upleftVal += mismatch;
                            op = kMismatch;
                        }
                        int max = upleftVal;
                        int upVal = up.value + gap;
                        int leftVal = left.value + gap;
                        if (upVal > max) {
                            max = upVal;
                            op = kDelete;
                        }
                        if (leftVal > max) {
                            max = leftVal;
                            op = kInsert;
                        }
                        if (max < 0) {
                            max = 0;
                            op = kNone;
                        }
                        Mat[i][j].value = max;
                        Mat[i][j].operation = op;
                        if (max > maxVal) {
                            maxVal = max;
                            maximalCell.first = i;
                            maximalCell.second = j;
                        }
                    }
                }

                if (cigar) {
                    ClippedCigarBuilder(maximalCell.first, maximalCell.second);
                }
                return maxVal;
            }  

        int SemiGlobal() {
            for (int i = 1; i < query_len; i++) {
                Mat[0][i].value = std::numeric_limits<int>::min();
                Mat[0][i].operation = kDelete;
            }
            for (int i = 1; i < target_len; i++) {
                Mat[i][0].value = std::numeric_limits<int>::min();
                Mat[i][0].operation = kInsert;
            }
            std::pair <int, int> maximalCell(0, 0);
            int maxVal = std::numeric_limits<int>::min();;
            for (int i = 1; i < target_len; i++) {
                for (int j = 1; j < query_len; j++) {
                    Cell upleft = Mat[i-1][j-1];
                    Cell up = Mat[i-1][j];
                    Cell left = Mat[i][j-1];
                    Operation op;
                    int upleftVal = upleft.value;
                    if (query[j] == target[i]) {
                        upleftVal += match;
                        op = kMatch;
                    }
                    else {
                        upleftVal -= mismatch;
                        op = kMismatch;
                    }
                    int max = upleftVal;
                    
                    int upVal = up.value - gap;
                    int leftVal = left.value - gap;
                    if (upVal > max) {
                        max = upVal;
                        op = kDelete;
                    }
                    if (leftVal > max) {
                        max = leftVal;
                        op = kInsert;
                    }
                    Mat[i][j].value = max;
                    Mat[i][j].operation = op;
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

            if (cigar) {
               ClippedCigarBuilder(maximalCell.first, maximalCell.second);
            }
            return maxVal;
        }

        std::string CigarBuilder(int& row, int& col, std::string starting_cigar) {
            //std::cout << "checkCigarBuilder1\n";
            while (Mat[row][col].operation != kNone) {
                switch (Mat[row][col].operation) {
                    case kDelete:
                        starting_cigar = "D" + starting_cigar;
                        row--;
                        break;
                    
                    case kInsert:
                        starting_cigar = "I" + starting_cigar;
                        col--;
                        break;
                    
                    case kMatch:
                        starting_cigar = "=" + starting_cigar;
                        row--;
                        col--;
                        break;

                    case kMismatch:
                        starting_cigar = "X" + starting_cigar;
                        row--;
                        col--;
                        break;
                    
                    default: 
                        break;
                }
            }
            //std::cout << "checkcigar2\n";

            std::string final_cigar = "";
            char symbol = starting_cigar[0];
            int count = 1;
            bool flag;
            for (int i = 1; i < starting_cigar.length(); i++) {
                bool flag = (symbol == starting_cigar[i]);
                if (flag) {
                    count++;
                } else {
                    final_cigar += std::to_string(count) + symbol;
                    symbol = starting_cigar[i];
                    count = 1;
                }
            }
            final_cigar += std::to_string(count) + symbol; 
            //std::cout << "checkfinalcigar\n"; 
            return final_cigar;

        }

        void ClippedCigarBuilder(int row, int col) {
            std::string starting_cigar = "";
            for (int i = row + 1; i < query_len; i++) {
                starting_cigar = "S" + starting_cigar;
            }
            std::string final_cigar = CigarBuilder(row, col, starting_cigar);
            if (row != 0) {
               final_cigar = std::to_string(row) + "S" + final_cigar; 
            }
            *target_begin = col + 1;
            *cigar = final_cigar;
        }

        int Align(AlignmentType type) {

            switch(type) {
                case GLOBAL:
                    return NeedlemanWunsch();
                case LOCAL:
                    return SmithWaterman();
                case SEMIGLOBAL:
                    return SemiGlobal();
                default:
                    return -1;    
            }
        }
    };   
}


