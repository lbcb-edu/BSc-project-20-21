#include "alignment/blue_alignment.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace blue {
namespace alignment {

enum Operation { kMatch, kMismatch, kDelete, kInsert, kNone };

struct Cell {
  int value;
  Operation operation;
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

  Aligner(const char* query, unsigned int query_len, const char* target,
          unsigned int target_len, int match, int mismatch, int gap,
          std::string* cigar, unsigned int* target_begin)
      : query(query),
        query_len(query_len),
        target(target),
        target_len(target_len),
        match(match),
        mismatch(mismatch),
        gap(gap),
        cigar(cigar),
        target_begin(target_begin),
        matrix_(query_len + 1, std::vector<Cell>(target_len + 1)) {
    matrix_[0][0] = {0, Operation::kNone};
  }

  int NeedlemanWunsch() {
    for (int i = 1; i < query_len + 1; i++)
      matrix_[i][0] = {gap * i, Operation::kDelete};
    for (int j = 1; j < target_len + 1; j++)
      matrix_[0][j] = {gap * j, Operation::kInsert};

    for (int i = 1; i < query_len + 1; i++) {
      for (int j = 1; j < target_len + 1; j++) {
        Cell diagonal;
        if (query[i - 1] == target[j - 1])  // match
          diagonal = {matrix_[i - 1][j - 1].value + match, Operation::kMatch};
        else
          diagonal = {matrix_[i - 1][j - 1].value + mismatch,
                      Operation::kMismatch};

        Cell top = {matrix_[i - 1][j].value + gap, Operation::kDelete};
        Cell left = {matrix_[i][j - 1].value + gap, Operation::kInsert};

        matrix_[i][j] = std::max({diagonal, top, left}, [](Cell a, Cell b) {
          return a.value < b.value;
        });
      }
    }
    // TODO: cigar
    return matrix_[query_len][target_len].value;
  }

  int SmithWaterman() {
    for (int i = 1; i < query_len + 1; i++)
      matrix_[i][0] = {0, Operation::kNone};
    for (int j = 1; j < target_len + 1; j++)
      matrix_[0][j] = {0, Operation::kNone};

    Cell neutral = {0,Operation::kNone};
    int max_value = 0;

    for (int i = 1; i < query_len + 1; i++) {
      for (int j = 1; j < target_len + 1; j++) {
        Cell diagonal;
        if(query[i - 1] == target[j - 1]) // match
          diagonal = {matrix_[i - 1][j - 1].value + match, Operation::kMatch};
        else
          diagonal = {matrix_[i - 1][j - 1].value + mismatch, 
                      Operation::kMismatch};
        
        Cell top = {matrix_[i - 1][j].value + gap, Operation::kDelete};
        Cell left = {matrix_[i][j - 1].value + gap, Operation::kInsert};


        matrix_[i][j] = std::max({neutral, diagonal, top, left}, [](Cell a, Cell b) {
          return a.value < b.value;
        });

        if (matrix_[i][j].value > max_value) {
            max_value = matrix_[i][j].value;
        } 
      }
    }
    // TODO: cigar
    return max_value;
  }

  int SemiGlobal() { 
    for (int i = 1; i < query_len + 1; i++)
      matrix_[i][0] = {0, Operation::kNone};
    for (int j = 1; j < target_len + 1; j++)
      matrix_[0][j] = {0, Operation::kNone};

    
    for (int i = 1; i < query_len + 1; i++) {
      for (int j = 1; j < target_len + 1; j++) {
        Cell diagonal;
        if(query[i - 1] == target[j - 1]) // match
          diagonal = {matrix_[i - 1][j - 1].value + match, Operation::kMatch};
        else
          diagonal = {matrix_[i - 1][j - 1].value + mismatch, 
                      Operation::kMismatch};
        
        Cell top = {matrix_[i - 1][j].value + gap, Operation::kDelete};
        Cell left = {matrix_[i][j - 1].value + gap, Operation::kInsert};


        matrix_[i][j] = std::max({diagonal, top, left}, [](Cell a, Cell b) {
          return a.value < b.value;
        });
      }
    }

    int max_value_row, max_value_column = 0;
    int max_value;

    for (int i = 0; i < query_len + 1; i++) {
        if (matrix_[i][target_len].value > max_value_column) {
          max_value_column = matrix_[i][target_len].value;
        }
    }
    for (int j = 0; j < target_len + 1; j++) {
        if (matrix_[query_len][j].value > max_value_row) {
          max_value_row = matrix_[query_len][j].value;
        }
    }

    max_value = max_value_column > max_value_row ? max_value_column : max_value_row;

    // TODO: cigar
    return max_value;
  } 

 private:
  std::vector<std::vector<Cell>> matrix_;
};

}  // namespace alignment
}  // namespace blue

int blue::Align(const char* query, unsigned int query_len, const char* target,
                unsigned int target_len, AlignmentType type, int match,
                int mismatch, int gap, std::string* cigar,
                unsigned int* target_begin) {
  alignment::Aligner aligner(query, query_len, target, target_len, match,
                             mismatch, gap, cigar, target_begin);

  switch (type) {
    case AlignmentType::kGlobal: return aligner.NeedlemanWunsch();
    case AlignmentType::kLocal: return aligner.SmithWaterman();
    case AlignmentType::kSemiGlobal: return aligner.SemiGlobal();
    default:
      throw std::invalid_argument(
          "[blue::alignment] error: invalid alignment type");
  }
}
