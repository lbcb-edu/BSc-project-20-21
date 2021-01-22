#pragma once

#include <vector>
#include<string>
#include<iostream>

using namespace std;
namespace orange
{
	class Minimizer
	{
	public:
		const char* sequence;
		unsigned int sequence_len;
		unsigned int kmer_len;
		unsigned int window_len;
		Minimizer(const char* sequence, unsigned int sequence_len, unsigned int kmer_len, unsigned int window_len) :sequence(sequence),
			sequence_len(sequence_len), kmer_len(kmer_len), window_len(window_len) {};

		vector<tuple<unsigned int, unsigned int, bool>> Minimize(
			const char* sequence, unsigned int sequence_len, unsigned int kmer_len, unsigned int window_len);
	};
}