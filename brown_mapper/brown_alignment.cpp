#include <iostream>
#include "brown_alignment.hpp"

namespace brown {

    //enum AlignmentType {GLOBAL, LOCAL, SEMIGLOBAL};

    int bzvz() {
        return 1;
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
            return 2;
    }
}