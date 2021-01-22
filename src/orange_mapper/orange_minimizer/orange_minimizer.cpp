#include "orange_minimizer.h"
#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <list>
using namespace std;

namespace orange {
    //format: [value, position, origin]
    using Kmer = std::tuple<unsigned int, unsigned int, bool>;

    unsigned int assignValue(char c) {
        switch (toupper(c)) {
        case 'C': return 0;
        case 'A': return 1;
        case 'T': return 2;
        case 'G': return 3;
        }
    }
    unsigned int complement(unsigned int n) {
        return 3 - n;
    }
    //compares values of two kmers
    bool compare(Kmer kmer1, Kmer kmer2) {
        if (std::get<0>(kmer1) < std::get<0>(kmer2)) return true;
        return false;
    }
    //returns current kmer value
    // if bool is true -> complement
    unsigned long int kmerMask(const char* sequence, bool isCompl, unsigned int k, int position) {
        unsigned long int kmer = 0;

        for (unsigned int i = 0; i < k; i++) {
            if (!isCompl) {
                kmer = kmer << 2 | assignValue(sequence[i + position]);
            }
            else {
                kmer = kmer << 2 | complement(assignValue(sequence[i + position]));
            }
        }
        return kmer;
    }
    //function that returns minimizers from a window
    //s = sequence length
    //k = kmer length
    //w = window length
    Kmer checkWindow(const char* sequence, unsigned int s, unsigned int k, unsigned int w, int position) {
        Kmer result;
        vector<Kmer> allKmers;
        Kmer kmer;
        Kmer kmerCompl;
        unsigned int i = position;
        unsigned int value;
        for (; i < position + w; i++) {
            value = 0;
            value = kmerMask(sequence, false, k, i);
            kmer = make_tuple(value, i, true);
            value = kmerMask(sequence, true, k, i);
            kmerCompl = make_tuple(value, i, false);
            result = compare(kmer, kmerCompl) ? kmer : kmerCompl;
            allKmers.push_back(result);
        }
        sort(allKmers.begin(), allKmers.end());
        if (!allKmers.empty())
            result = allKmers.front();
        return result;
    }

    //function that finds end minimizers
    vector<Kmer> endMinimizers(const char* sequence, unsigned int s, unsigned int k, unsigned int w) {
        vector<Kmer> results;
        Kmer temp;
        Kmer beginKmer;
        Kmer endKmer;
        unsigned long int value;
        unsigned long int complValue;
        //begin minimizer
        value = kmerMask(sequence, false, k, 0);
        complValue = kmerMask(sequence, true, k, 0);
        if (value < complValue) {
            beginKmer = make_tuple(value, 0, 1);
        }
        else {
            beginKmer = make_tuple(complValue, 0, 0);
        }
        //end minimizer
        value = kmerMask(sequence, false, k, s - k);
        complValue = kmerMask(sequence, true, k, s - k);
        if (value < complValue) {
            endKmer = make_tuple(value, s - k, 1);
        }
        else {
            endKmer = make_tuple(complValue, s - k, 0);
        }
        results.push_back(beginKmer);
        results.push_back(endKmer);
        return results;
    }

    //function that finds interior minimizers
    vector<Kmer> checkAll(const char* sequence, unsigned int s, unsigned int k, unsigned int w) {
        vector<Kmer> results;
        Kmer temp;
        for (unsigned int i = 0; i <= s - (w + k - 1) + 1; i++) {
            //adding minimizers from the window to result
            temp = checkWindow(sequence, s, k, w, i);
            results.push_back(temp);
        }
        return results;
    }

	vector<tuple<unsigned int, unsigned int, bool>> Minimizer::all_minimizers(
        const char* sequence, unsigned int sequence_len, unsigned int kmer_len, unsigned int window_len) {
        vector<Kmer> all;
        vector<Kmer> temp;

        //interior minimizers
        temp = checkAll(sequence, sequence_len, kmer_len, window_len);
        all.insert(end(all), begin(temp), end(temp));

        //end minimizers
        temp = endMinimizers(sequence, sequence_len, kmer_len, window_len);
        all.insert(end(all), begin(temp), end(temp));

        sort(all.begin(), all.end());

        return all;
	}

    vector<tuple<unsigned int, unsigned int, bool>> Minimizer::Minimize(
        const char* sequence, unsigned int sequence_len, unsigned int kmer_len, unsigned int window_len) {
        vector<Kmer> minimizers;
        vector<Kmer> temp;

        //interior minimizers
        temp = checkAll(sequence, sequence_len, kmer_len, window_len);
        minimizers.insert(end(minimizers), begin(temp), end(temp));

        //end minimizers
        temp = endMinimizers(sequence, sequence_len, kmer_len, window_len);
        minimizers.insert(end(minimizers), begin(temp), end(temp));

        sort(minimizers.begin(), minimizers.end());
        minimizers.erase(unique(minimizers.begin(), minimizers.end()), minimizers.end());
        return minimizers;
    };
}
