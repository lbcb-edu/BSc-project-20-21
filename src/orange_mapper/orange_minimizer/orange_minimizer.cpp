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
    //checks ih the last one is better than current kmer
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
    vector<Kmer> checkWindow(const char* sequence, unsigned int s, unsigned int k, unsigned int w, int position) {
        Kmer result;
        vector<Kmer> results;
        Kmer kmer;
        Kmer kmerCompl;
        unsigned int i = position;
        //while(i < s)???
        for (; i < w; i++) {
            std::get<0>(kmer) = kmerMask(sequence, false, k, i);
            std::get<1>(kmer) = i;
            std::get<2>(kmer) = true;
            std::get<0>(kmerCompl) = kmerMask(sequence, true, k, i);
            std::get<1>(kmerCompl) = i;
            std::get<2>(kmerCompl) = false;
            result = compare(kmer, kmerCompl) ? kmer : kmerCompl;
            while (!results.empty() && compare(result, results.back())) {
                results.pop_back();
            }
            results.push_back(result);

        }
        return results;
    }
    //function that finds end minimizers
    vector<Kmer> beginingMinimizers(const char* sequence, unsigned int s, unsigned int k, unsigned int w) {
        vector<Kmer> results;

        return results;
    }

    //function that finds end minimizers
    vector<Kmer> endMinimizers(const char* sequence, unsigned int s, unsigned int k, unsigned int w) {
        vector<Kmer> results;
        Kmer minKmer;
        Kmer kmer;
        unsigned int value;
        Kmer kmerCompl;
        for (unsigned int i = k; i < s; i++) {
            //first pass
            if (i == k) {
                value = kmerMask(sequence, false, k, i - k);
                kmer = make_tuple(value, i - k, false);
                value = kmerMask(sequence, true, k, i - k);
                kmerCompl = make_tuple(value, i - k, true);
                minKmer = compare(kmer, kmerCompl) ? kmer : kmerCompl;
                results.push_back(minKmer);
            }
            else {
                for (int j = i - k + 1; j > 0; j--) {
                    value = kmerMask(sequence, false, k, j);
                    kmer = make_tuple(value, j, false);
                    value = kmerMask(sequence, true, k, j);
                    kmerCompl = make_tuple(value, j, true);
                    kmer = compare(kmer, kmerCompl) ? kmer : kmerCompl;
                    if (compare(kmer, minKmer)) {
                        minKmer = kmer;
                        results.push_back(minKmer);
                    }
                }
            }
        }
        return results;
    }

    //function that finds interior minimizers
    vector<Kmer> checkAll(const char* sequence, unsigned int s, unsigned int k, unsigned int w) {
        vector<Kmer> results;
        vector<Kmer> temp;
        Kmer kmer;
        Kmer kmerCompl;
        for (unsigned int i = 0; i <= s - (w + k - 1); i++) {
            //adding minimizers from the window to result
            temp = checkWindow(sequence, s, k, w, i);
            results.insert(end(results), begin(temp), end(temp));
        }
        //sorting vector and removing duplicates
        sort(results.begin(), results.end());
        results.erase(unique(results.begin(), results.end()), results.end());
        return results;
    }


    vector<tuple<unsigned int, unsigned int, bool>> Minimizer::Minimize(
        const char* sequence, unsigned int sequence_len, unsigned int kmer_len, unsigned int window_len) {
        vector<Kmer> minimizers;
        vector<Kmer> temp;

        //TODO : begininig minimizers
        temp = beginingMinimizers(sequence, sequence_len, kmer_len, window_len);
        minimizers.insert(end(minimizers), begin(temp), end(temp));

        //interior minimizers
        temp = checkAll(sequence, sequence_len, kmer_len, window_len);
        minimizers.insert(end(minimizers), begin(temp), end(temp));

        //TODO: end minimizers
        temp = endMinimizers(sequence, sequence_len, kmer_len, window_len);
        minimizers.insert(end(minimizers), begin(temp), end(temp));

        sort(minimizers.begin(), minimizers.end());
        minimizers.erase(unique(minimizers.begin(), minimizers.end()), minimizers.end());
        return minimizers;
    };
}
