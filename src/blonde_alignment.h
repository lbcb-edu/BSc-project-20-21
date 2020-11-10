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

private:
    void computeCell(int row, int col) {
        // Possible scores
        int diagonal_score = table[row - 1][col - 1].score_;
        diagonal_score += (query[row - 1] == target[col - 1]) ? match : mismatch;
        int top_score = table[row - 1][col].score_ + gap;
        int left_score = table[row][col - 1].score_ + gap;

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
        table[row][col].score_ = max_score;
        table[row][col].direction_ = dir;
    }

public:
    CellComputer(
        const char* query, unsigned int query_len,
        const char* target, unsigned int target_len,
        std::vector<std::vector<Cell>>& table,
        int match, int mismatch, int gap
    ) : query(query), query_len(query_len), target(target), target_len(target_len), table(table), match(match), mismatch(mismatch), gap(gap) {

    }

    void computeAllCells(AlignmentType type) {
        for (int i = 1; i < table.size(); i++) {
            for (int j = 1; j < table[0].size(); j++) {
                computeCell(i, j);
                // Ovo je dio koda za Lokalni alignment koji jos nije napisan, nisam sigurran je li < 0 ili <= 0
                // if (type == kLocal && table[i][j].score_ < 0) {
                //     table[i][j].score_ = 0;
                //     table[i][j].direction_ = kNone;
                // }
            }
        }
    }
};


void initAlignmentTable(std::vector<std::vector<Cell>>& table, int init_penalty) {
    int num_of_rows = table.size();
    int num_of_cols = table[0].size();
    
    for (int i = 1; i < num_of_rows; i++) {
    table[i][0].score_ = i * init_penalty;
    table[i][0].direction_ = kUp;
    }
    for (int i = 1; i < num_of_cols; i++) {
        table[0][i].score_ = i * init_penalty;
        table[0][i].direction_ = kLeft;
    }
}

void compressCigar(std::string& uncompressed_cigar, std::string& cigar_result) {
    cigar_result = "";
    char letter = uncompressed_cigar[0];
    int cnt = 1;
    for (int i = 1; i < uncompressed_cigar.size(); i++) {
        if (uncompressed_cigar[i] != letter) {
            cigar_result += std::to_string(cnt) + letter;
            letter = uncompressed_cigar[i];
            cnt = 1;
        } else {
            cnt++;
        }
    }
    cigar_result += std::to_string(cnt) + letter;
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

    int init_penalty = 0;
    if (type == kGlobal)
        init_penalty = gap;

    std::vector<std::vector<Cell>> table = std::vector<std::vector<Cell>> (query_len + 1, std::vector<Cell>(target_len + 1));
    initAlignmentTable(table, init_penalty);

    CellComputer computer = CellComputer(query, query_len, target, target_len, table, match, mismatch, gap);
    computer.computeAllCells(type);

    int align_score;
    int target_begin_result;
    std::string cigar_result = "";
    switch (type) {
    case kLocal: { // Treba napisati
        std::cout << "Local" << std::endl;
        align_score = 0;
        break;
    }

    case kGlobal: {
        std::cout << "Global ";
        
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
            compressCigar(cigar_tmp, cigar_result);
        }
        target_begin_result = 0;
        align_score = table[query_len][target_len].score_;
        break;
    }

    case kSemiGlobal: {
        std::cout << "SemiGlobal ";
        
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
            target_begin_result = j;
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
            compressCigar(cigar_tmp, cigar_result);
        }
        align_score = maximum;
        break;
    }
    
    default:
        break;
    }
    
    //Rezultati
    if (cigar) *cigar = cigar_result;
    if (target_begin) *target_begin = target_begin_result;
    return align_score;   
}

}
}

#endif 