
#ifndef BLONDE_ALIGNMENT_H_
#define BLONDE_ALIGNMENT_H_

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
    std::int32_t score_;
    SrcDirection direction_;
};

int Align(
    const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    AlignmentType type,
    int match,
    int mismatch,
    int gap,
    std::string* cigar = nullptr,
    unsigned int* target_begin = nullptr) {

    switch (type) {
    case kLocal:
        std::cout << "Local" << std::endl;
        break;

    case kGlobal:
        std::cout << "Global" << std::endl;
        break;

    case kSemiGlobal:
        std::cout << "SemiGlobal" << std::endl;
        break;
    
    default:
        break;
    }

    return target_len + query_len;   
}

}
}

#endif 