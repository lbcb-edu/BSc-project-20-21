#ifndef BLONDE_ALIGNMENT_H_
#define BLONDE_ALIGNMENT_H_

#include <iostream>
#include <string>
#include <getopt.h>
#include <vector>
#include <time.h>
#include <algorithm>

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
    void computeCell(int row, int col);

public:
    CellComputer(
        const char* query, unsigned int query_len,
        const char* target, unsigned int target_len,
        std::vector<std::vector<Cell>>& table,
        int match, int mismatch, int gap
    );

    void computeAllCells(AlignmentType type);
};

void initAlignmentTable(std::vector<std::vector<Cell>>& table, int init_penalty, AlignmentType type);

void calcBacktrackPath(
    const std::vector<std::vector<Cell>> table,
    int mismatch,
    std::string& cigar_tmp,
    int& i, int& j);
    
void calcCigar(std::string& uncompressed_cigar, std::string& cigar_result);

int Align(
    const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    AlignmentType type,
    int match,
    int mismatch,
    int gap,
    std::string* cigar,
    unsigned int* target_begin);

}
}

#endif
