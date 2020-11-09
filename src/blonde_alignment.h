
#ifndef BLONDE_ALIGNMENT_H_
#define BLONDE_ALIGNMENT_H_

#include <string>

namespace blonde {
namespace alignment {

enum AlignmentType {
    kLocal,
    kGlobal,
    kSemiGlobal
};

enum SrcDirection {
    kLeft,
    kUp,
    kDiagonal,
    kNone
};

struct Cell {
public:
    std::int32_t score_;
    SrcDirection direction_;

    Cell() : score_(0), direction_(kNone) {}
};

class CellComputer{
private:
    const char* query;
    const char* target;
    unsigned int query_len;
    unsigned int target_len;
    std::vector<std::vector<Cell>>& table;
    int match, mismatch, gap;

public:
    CellComputer(
        const char* query, unsigned int query_len,
        const char* target, unsigned int target_len,
        std::vector<std::vector<Cell>>& table,
        int match, int mismatch, int gap
    ) : query(query), query_len(query_len), target(target), target_len(target_len), table(table), match(match), mismatch(mismatch), gap(gap) {

    }

    void computeCell(int table_row, int table_col) {
        if (table_row == 0 || table_col == 0)
            throw "table_row and table_col must be greater than 0";
        if(table_row > query_len)
            throw "table_row out of bounds";
        if (table_col > target_len)
            throw "table_col out of bounds";

        // Diagonal
        int diagonal_score = table[table_row - 1][table_col - 1].score_;
        diagonal_score += (query[table_row - 1] == target[table_col - 1]) ? match : mismatch;

        //Top
        int top_score = table[table_row - 1][table_col].score_ + gap;

        //Left
        int left_score = table[table_row][table_col - 1].score_ + gap;

        int max_score;
        SrcDirection dir;
        if (diagonal_score >= left_score && diagonal_score >= top_score) {
            dir = kDiagonal;
            max_score = diagonal_score;
        } else if (left_score >= diagonal_score && left_score >= top_score) {
            dir = kLeft;
            max_score = left_score;
        } else if(top_score >= left_score && top_score >= diagonal_score) {
            dir = kUp;
            max_score = top_score;
        }
        table[table_row][table_col].score_ = max_score;
        table[table_row][table_col].direction_ = dir;
    }

};


int semiGlobal(
    const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    int match,
    int mismatch,
    int gap,
    std::string* cigar = nullptr,
    unsigned int* target_begin = nullptr) {

    std::vector<std::vector<Cell>> table = std::vector<std::vector<Cell>> (query_len + 1, std::vector<Cell>(target_len + 1));
    
    for (int i = 1; i < query_len + 1; i++) {
        table[i][0].score_ = 0;
        table[i][0].direction_ = kUp;
    }
    for (int i = 1; i < target_len + 1; i++) {
        table[0][i].score_ = 0;
        table[0][i].direction_ = kLeft;
    }

    CellComputer computer = CellComputer(query, query_len, target, target_len, table, match, mismatch, gap);

    for (int i = 1; i < query_len + 1; i++) {
        for (int j = 1; j < target_len + 1; j++) {
            computer.computeCell(i, j);
        }
    }

    //Find Maximum in last row or column
    int maximum = table[0][target_len].score_;
    int max_indx_row = 0;
    int max_indx_col = target_len;
    for (int i = 0; i < query_len + 1; i++) {
        if (table[i][target_len].score_ > maximum) {
            maximum = table[i][target_len].score_;
            max_indx_row = i;
            max_indx_col = target_len;
        }
    }
    for (int i = 0; i < target_len + 1; i++) {
        if (table[query_len][i].score_ > maximum) {
            maximum = table[query_len][i].score_;
            max_indx_row = query_len;
            max_indx_col = i;
        }
    }

    if(cigar || target_begin) {
        int i = max_indx_row;
        int j = max_indx_col;
        int target_begin_result = j;
        std::string cigar_tmp = "";
        for(int k = i + 1; k < query_len + 1; k++) {
            cigar_tmp = "I" + cigar_tmp;
        }
        while (i != 0) {
            switch (table[i][j].direction_) {
            case kDiagonal:
                if(table[i-1][j-1].score_ + mismatch == table[i][j].score_) {
                    cigar_tmp = "X" + cigar_tmp;
                } else {
                    cigar_tmp = "=" + cigar_tmp;
                }
                i--;
                j--;
                break;
            
            case kUp:
                cigar_tmp = "I" + cigar_tmp;
                i--;
                break;

            case kLeft:
                cigar_tmp = "D" + cigar_tmp;
                j--;
                break;

            case kNone:
                throw "cell doesn't have a computed direction";
                break;
            }
            target_begin_result = j;
        }

        std::string cigar_result = "";
        char letter = cigar_tmp[0];
        int cnt = 1;
        for (int i = 1; i < cigar_tmp.size(); i++) {
            if (cigar_tmp[i] != letter) {
                cigar_result += std::to_string(cnt) + letter;
                letter = cigar_tmp[i];
                cnt = 1;
            } else {
                cnt++;
            }
        }
        if (cigar) {
            cigar_result += std::to_string(cnt) + letter;
            *cigar = cigar_result;
        }

        if (target_begin) {
            *target_begin = target_begin_result;
        }
    }
    return maximum;
}

int needlemanWunsch(
    const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    int match,
    int mismatch,
    int gap,
    std::string* cigar = nullptr,
    unsigned int* target_begin = nullptr) {

    std::vector<std::vector<Cell>> table = std::vector<std::vector<Cell>> (query_len + 1, std::vector<Cell>(target_len + 1));
    
    for (int i = 1; i < query_len + 1; i++) {
        table[i][0].score_ = i * gap;
        table[i][0].direction_ = kUp;
    }
    for (int i = 1; i < target_len + 1; i++) {
        table[0][i].score_ = i * gap;
        table[0][i].direction_ = kLeft;
    }

    CellComputer computer = CellComputer(query, query_len, target, target_len, table, match, mismatch, gap);

    for (int i = 1; i < query_len + 1; i++) {
        for (int j = 1; j < target_len + 1; j++) {
            computer.computeCell(i, j);
        }
    }

    if (cigar) {
        int i = query_len;
        int j = target_len;
        std::string cigar_tmp = "";
        while (i != 0 || j != 0) {
            switch (table[i][j].direction_) {
            case kDiagonal:
                if(table[i-1][j-1].score_ + mismatch == table[i][j].score_) {
                    cigar_tmp = "X" + cigar_tmp;
                } else {
                    cigar_tmp = "=" + cigar_tmp;
                }
                i--;
                j--;
                break;
            
            case kUp:
                cigar_tmp = "I" + cigar_tmp;
                i--;
                break;

            case kLeft:
                cigar_tmp = "D" + cigar_tmp;
                j--;
                break;

            case kNone:
                throw "cell doesn't have a computed direction";
                break;
            }
        }



        std::string cigar_result = "";
        char letter = cigar_tmp[0];
        int cnt = 1;
        for (int i = 1; i < cigar_tmp.size(); i++) {
            if (cigar_tmp[i] != letter) {
                cigar_result += std::to_string(cnt) + letter;
                letter = cigar_tmp[i];
                cnt = 1;
            } else {
                cnt++;
            }
        }
        cigar_result += std::to_string(cnt) + letter;
        *cigar = cigar_result;
    }

    if (target_begin) {
        *target_begin = 0;
    }
    
    // for (int i = 0; i < query_len + 1; i++) {
    //     for (int j = 0; j < target_len + 1; j++) {
    //         std::cout << " " << table[i][j].score_;
    //         if (table[i][j].direction_ == kDiagonal) {
    //             std::cout << "D";
    //         } else if (table[i][j].direction_ == kLeft) {
    //             std::cout << "L";
    //         } else if (table[i][j].direction_ == kUp) {
    //             std::cout << "T";
    //         }
    //     }
    //     std::cout << std::endl;
    // }

    return table[query_len][target_len].score_;

}

int Align(
    const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    AlignmentType type,
    int match,
    int mismatch,
    int gap,
    std::string* cigar = nullptr,
    unsigned int* target_begin = nullptr) {

    std::cout << "Target sequence:\n";
    for(int i = 0; i < target_len; i++) {
        std::cout << target[i];
    }
    std::cout << "\n";
    
    std::cout << "Query sequence:\n";
    for(int i = 0; i < query_len; i++) {
        std::cout << query[i];
    }
    std::cout << "\n";

    int result;
    switch (type) {
    case kLocal:
        std::cout << "Local" << std::endl;
        break;

    case kGlobal:
        std::cout << "Global ";
        result = needlemanWunsch(query, query_len, target, target_len, match, mismatch, gap, cigar, target_begin);
        break;

    case kSemiGlobal:
        std::cout << "SemiGlobal ";
        result = semiGlobal(query, query_len, target, target_len, match, mismatch, gap, cigar, target_begin);
        break;
    
    default:
        break;
    }

    return result;   
}

}
}

#endif 