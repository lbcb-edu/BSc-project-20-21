#pragma once
#include <vector>
#include<string>
#include<iostream>
#include<algorithm>

using namespace std;
namespace orange
{
	enum AlignmentType { global, local, semiGlobal };
	enum Operation { sMatch, sMismatch, sDelete, sInsert, sNone };
        struct Cell {
            int value;
            Operation op;
        };
	class Alignment
	{
	public:
	const char* query;
        unsigned int query_len;
        const char* target;
        unsigned int target_len;
        AlignmentType type;
        int match;
        int mismatch;
        int gap;
        std::string* cigar = nullptr;
        unsigned int* target_begin = nullptr;
	vector<vector<Cell>> mainMatrix;
        Alignment(const char* query, unsigned int query_len,
            char* target, unsigned int target_len,
            orange::AlignmentType type,
            int match,
            int mismatch,
            int gap,
            string* cigar,
            unsigned int* target_begin):
            query(query), query_len(query_len), target(target), target_len(target_len), type(type), match(match), 
            mismatch(mismatch), gap(gap), cigar(cigar), target_begin(target_begin), 
            mainMatrix((query_len + 1), std::vector<Cell>(target_len + 1, { 0, Operation::sNone })){};

	string generateCigar(int maxRow, int maxColumn, string initialVal, const std::vector<vector<Cell>> mainMatrix);
	void calculateMax(int i, int j, const char* query,
            const char* target, AlignmentType type,
            int match,
            int mismatch,
            int gap, vector<vector<Cell>> mainMatrix);
        int Align(const char* query, unsigned int query_len,
            const char* target, unsigned int target_len,
            AlignmentType type,
            int match,
            int mismatch,
            int gap,
            std::string* cigar = nullptr,
            unsigned int* target_begin = nullptr);
	};
}
