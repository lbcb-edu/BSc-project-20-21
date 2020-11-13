#include "alignment/blue_alignment.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

namespace blue {
namespace alignment {

enum Operation { kMatch = 0, kMismatch, kDelete, kInsert, kNone };

struct Cell {
  int value;
  Operation operation;
};

struct MaxCell {
  int value;
  int row_index;
  int col_index;
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
        matrix_(query_len + 1,
                std::vector<Cell>(target_len + 1, {0, Operation::kNone})) {}

  int NeedlemanWunsch() {
    for (int i = 1; i < query_len + 1; i++)
      matrix_[i][0] = {gap * i, Operation::kDelete};
    for (int j = 1; j < target_len + 1; j++)
      matrix_[0][j] = {gap * j, Operation::kInsert};

    ComputeMatrix();

    if (cigar) {
      std::string initial_cigar = "";
      int i = query_len;
      int j = target_len;
      *cigar = CigarCreator(i, j, initial_cigar);
      *target_begin = 1;
    }

    return matrix_[query_len][target_len].value;
  }

  int SmithWaterman() {
    ComputeMatrix({0, Operation::kNone});

    if (cigar)
      ClippedCigarCreator(smith_waterman_max_.row_index,
                          smith_waterman_max_.col_index);

    return smith_waterman_max_.value;
  }

  int SemiGlobal() {
    ComputeMatrix();

    if (cigar)
      ClippedCigarCreator(semi_global_max_.row_index,
                          semi_global_max_.col_index);

    return semi_global_max_.value;
  }

 private:
  std::vector<std::vector<Cell>> matrix_;

  // matrix metadata
  MaxCell smith_waterman_max_ = {0, -1, -1};
  MaxCell semi_global_max_ = {std::numeric_limits<int>::min(), -1, -1};

  // Computes all cells in the matrix and sets matrix metadata
  void ComputeMatrix(Cell min_cell_limit = {std::numeric_limits<int>::min(),
                                            Operation::kNone}) {
    MaxCell last_row_max = {std::numeric_limits<int>::min(), -1, -1};
    MaxCell last_col_max = {std::numeric_limits<int>::min(), -1, -1};

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

        matrix_[i][j] =
            std::max({min_cell_limit, diagonal, top, left},
                     [](Cell a, Cell b) { return a.value < b.value; });

        if (matrix_[i][j].value > smith_waterman_max_.value)
          smith_waterman_max_ = {matrix_[i][j].value, i, j};
        if (i == query_len && matrix_[i][j].value > last_row_max.value)
          last_row_max = {matrix_[i][j].value, i, j};
        if (j == target_len && matrix_[i][j].value > last_col_max.value)
          last_col_max = {matrix_[i][j].value, i, j};
      }
    }

    semi_global_max_ =  // max(last_row_max, last_col_max)
        last_row_max.value > last_col_max.value ? last_row_max : last_col_max;
  }

  // Computes CIGAR for clipped alignment
  void ClippedCigarCreator(int row, int col) {
    std::string initial_cigar = "";
    int i = row;
    int j = col;
    for (int k = i + 1; k < query_len + 1; k++) {
      initial_cigar = "S" + initial_cigar;
    }
    std::string cigar_result = CigarCreator(i, j, initial_cigar);
    if (i != 0) {
      cigar_result = std::to_string(i) + "S" + cigar_result;
    }
    *target_begin = j + 1;
    *cigar = cigar_result;
  }

  std::string CigarCreator(int& i, int& j, std::string initial_cigar) {
    while (matrix_[i][j].operation != kNone) {
      switch (matrix_[i][j].operation) {
        case kMatch:
          initial_cigar = "=" + initial_cigar;
          i--;
          j--;
          break;
        case kMismatch:
          initial_cigar = "X" + initial_cigar;
          i--;
          j--;
          break;
        case kDelete:
          initial_cigar = "D" + initial_cigar;
          i--;
          break;
        case kInsert:
          initial_cigar = "I" + initial_cigar;
          j--;
          break;
        default: break;
      }
    }

    std::string cigar_result = "";
    char letter = initial_cigar[0];
    int count = 1;
    for (int i = 1; i < initial_cigar.length(); i++) {
      if (initial_cigar[i] == letter) {
        count++;
      } else {
        cigar_result += std::to_string(count) + letter;
        count = 1;
        letter = initial_cigar[i];
      }
    }
    cigar_result += std::to_string(count) + letter;

    return cigar_result;
  }
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
