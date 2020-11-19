#ifndef BLONDE_ALIGNMENT_H_
#define BLONDE_ALIGNMENT_H_

#include <iostream>
#include <string>
#include <getopt.h>
#include <vector>
#include <time.h>

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
                if (type == kLocal && table[i][j].score_ <= 0) {
                    table[i][j].score_ = 0;
                    table[i][j].direction_ = kNone;
                }
            }
        }
    }
};


void initAlignmentTable(std::vector<std::vector<Cell>>& table, int init_penalty, AlignmentType type) {
    int num_of_rows = table.size();
    int num_of_cols = table[0].size();
    if (type == kLocal || type == kSemiGlobal) init_penalty = 0;
    SrcDirection top_row_direction = kNone;
    SrcDirection first_col_direction = kNone;
    switch(type) {
    case kGlobal:
        top_row_direction = kLeft;
        first_col_direction = kUp;
        break;

    case kSemiGlobal:
        first_col_direction = kUp;
        break;

    case kLocal:
        break;
    }

    for (int i = 1; i < num_of_rows; i++) {
    table[i][0].score_ = i * init_penalty;
    table[i][0].direction_ = first_col_direction;
    }
    for (int i = 1; i < num_of_cols; i++) {
        table[0][i].score_ = i * init_penalty;
        table[0][i].direction_ = top_row_direction;
    }
}

void calcBacktrackPath(
    const std::vector<std::vector<Cell>> table,
    int mismatch,
    std::string& cigar_tmp,
    int& i, int& j) {

    while (table[i][j].direction_ != kNone) {
        switch (table[i][j].direction_) {
        case kDiagonal:
            if(table[i-1][j-1].score_ + mismatch == table[i][j].score_) {
                cigar_tmp += "X";
            } else {
                cigar_tmp += "=";
            }
            i--;
            j--;
            break;

        case kUp:
            cigar_tmp += "I";
            i--;
            break;

        case kLeft:
            cigar_tmp += "D";
            j--;
            break;

        case kNone:
            break;
        }
    }
}

void calcCigar(std::string& uncompressed_cigar, std::string& cigar_result) {
    std::reverse(uncompressed_cigar.begin(), uncompressed_cigar.end());
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

    std::vector<std::vector<Cell>> table = std::vector<std::vector<Cell>> (query_len + 1, std::vector<Cell>(target_len + 1));
    initAlignmentTable(table, gap, type);
    int row_cnt = table.size();
    int col_cnt = table[0].size();

    CellComputer computer = CellComputer(query, query_len, target, target_len, table, match, mismatch, gap);
    computer.computeAllCells(type);

    int align_score;
    int target_begin_result;
    std::string cigar_result = "";
    switch (type) {
    case kLocal: {
        //Find Maximum in whole table
        int maximum = table[0][0].score_;
        int max_indx_row = 0;
        int max_indx_col = 0;
        for (int i = 0; i < row_cnt; i++) {
            for (int j = 0; j < col_cnt; j++) {
                if (table[i][j].score_ >= maximum) {
                    maximum = table[i][j].score_;
                    max_indx_row = i;
                    max_indx_col = j;
                }
            }
        }
        if (cigar || target_begin) {
            int i = max_indx_row;
            int j = max_indx_col;
            std::string cigar_tmp = "";
            for(int k = i + 1; k < query_len + 1; k++) {
                cigar_tmp += "S";
            }
            calcBacktrackPath(table, mismatch, cigar_tmp, i, j);
            for(int k = i; k > 0; k--) {
                cigar_tmp += "S";
            }
            target_begin_result = j;
            calcCigar(cigar_tmp, cigar_result);

        }
        align_score = maximum;
        break;
    }

    case kGlobal: {
        if (cigar) {
            int i = query_len;
            int j = target_len;
            std::string cigar_tmp = "";
            calcBacktrackPath(table, mismatch, cigar_tmp, i, j);
            calcCigar(cigar_tmp, cigar_result);
        }
        target_begin_result = 0;
        align_score = table[query_len][target_len].score_;
        break;
    }

    case kSemiGlobal: {
        //Find Maximum in last row or column
        int maximum = table[0][target_len].score_;
        int max_indx_row = 0;
        int max_indx_col = target_len;
        for (int i = 0; i < row_cnt; i++) {
            if (table[i][target_len].score_ > maximum) {
                maximum = table[i][target_len].score_;
                max_indx_row = i;
                max_indx_col = target_len;
            }
        }
        for (int i = 0; i < col_cnt; i++) {
            if (table[query_len][i].score_ > maximum) {
                maximum = table[query_len][i].score_;
                max_indx_row = query_len;
                max_indx_col = i;
            }
        }

        if(cigar || target_begin) {
            int i = max_indx_row;
            int j = max_indx_col;
            std::string cigar_tmp = "";
            for(int k = i + 1; k < row_cnt; k++) {
                cigar_tmp += "I";
            }
            calcBacktrackPath(table, mismatch, cigar_tmp, i, j);
            target_begin_result = j;
            calcCigar(cigar_tmp, cigar_result);
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
