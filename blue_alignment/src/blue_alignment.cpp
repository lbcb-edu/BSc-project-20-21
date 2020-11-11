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



  std::string CigarCreator(int& i, int& j, std::string initial_cigar) {
    while (matrix_[i][j].operation != kNone) {
      switch(matrix_[i][j].operation) {
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
    
    if (cigar) {
      std::string initial_cigar = "";
      int i = query_len;
      int j = target_len;
      *cigar = CigarCreator(i,j,initial_cigar);
      *target_begin = 1;
    }

    return matrix_[query_len][target_len].value;
  }



  int SmithWaterman() {
    for (int i = 1; i < query_len + 1; i++)
      matrix_[i][0] = {0, Operation::kNone};
    for (int j = 1; j < target_len + 1; j++)
      matrix_[0][j] = {0, Operation::kNone};

    Cell neutral = {0,Operation::kNone};
    int max_value = 0;
    int max_row_index = 0;
    int max_column_index = 0;

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


        matrix_[i][j] = std::max({diagonal, top, left, neutral}, [](Cell a, Cell b) {
          return a.value < b.value;
        });

        if (matrix_[i][j].value > max_value) {
            max_value = matrix_[i][j].value;
            max_row_index = i;
            max_column_index = j;
        } 
      }
    }

    if (cigar) {
      std::string initial_cigar = "";
      int i = max_row_index;
      int j = max_column_index;
      for (int k = i+1; k < query_len + 1; k++) {
        initial_cigar = "S" + initial_cigar;
      }
      std::string cigar_result = CigarCreator(i,j,initial_cigar);
      if (i != 0) {
        cigar_result = std::to_string(i) + "S" + cigar_result;
      }
      *target_begin = j+1;
      *cigar = cigar_result;
    }
    
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

    int max_row_index = 0;
    int max_column_index = 0;

    int max_value = matrix_[0][target_len].value;

    for (int i = 0; i < query_len + 1; i++) {
        if (matrix_[i][target_len].value > max_value) {
          max_value = matrix_[i][target_len].value;
          max_row_index = i;
          max_column_index = target_len;
        }
    }
    for (int j = 0; j < target_len + 1; j++) {
        if (matrix_[query_len][j].value > max_value) {
          max_value = matrix_[query_len][j].value;
          max_row_index = query_len;
          max_column_index = j;
        }
    }

    if (cigar) {
      std::string initial_cigar = "";
      int i = max_row_index;
      int j = max_column_index;
      for (int k = i+1; k < query_len + 1; k++) {
        initial_cigar = "S" + initial_cigar;
      }
      std::string cigar_result = CigarCreator(i,j,initial_cigar);
      if (i != 0) {
        cigar_result = std::to_string(i) + "S" + cigar_result;
      }
      *target_begin = j+1;
      *cigar = cigar_result;
    }

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
