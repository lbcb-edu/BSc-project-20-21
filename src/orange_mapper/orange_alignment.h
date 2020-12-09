#pragma once
#include <string>

//#define WIN32_LEAN_AND_MEAN

namespace orange
{
	class Alignment
	{
	public:

        enum Operation { sMatch, sMismatch, sDelete, sInsert, sNone };
        struct Cell {
            int value;
            Operation op;
        };
        enum AlignmentType { global, local, semiGlobal };
        static int Align(const char* query, unsigned int query_len,
            const char* target, unsigned int target_len,
            AlignmentType type,
            int match,
            int mismatch,
            int gap,
            std::string* cigar = nullptr,
            unsigned int* target_begin = nullptr);
	};
}