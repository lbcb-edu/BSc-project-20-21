#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
using namespace std;

namespace pink{

struct kmer{
    unsigned long int current_kmer = 0;
    unsigned long int last_kmer = 0;
    unsigned int pos = 0;
};

int nucleotide_value(char c, bool reverse){
    if(c == 'C') return reverse ? 3 : 0;
    else if(c == 'A') return reverse ? 2 : 1;
    else if(c == 'T') return reverse ? 1 : 2;
    else if(c == 'G') return reverse ? 0 : 3;
}

tuple<unsigned int, unsigned int, bool> findMinimizers(unsigned long int orig, unsigned long int revComp, unsigned int sequence_len, 
    unsigned int kmer_len, unsigned int window_len, int i, unsigned long int mask){

        tuple<unsigned int, unsigned int, bool> minimizer;
        kmer origKmer;
        kmer revKmer;
        
        for(int j = 0; j < window_len; j++){
            origKmer.current_kmer = (orig & (mask << (window_len - 1 - j) * 2))  >> (window_len - 1 - j) * 2;
            revKmer.current_kmer = (revComp & (mask << (window_len - 1 - j) * 2))  >> (window_len - 1 - j) * 2;
            if(j == 0){
                origKmer.last_kmer = origKmer.current_kmer;
                revKmer.last_kmer = revKmer.current_kmer;
                origKmer.pos = i + j;
                revKmer.pos = i + j;
            }else{ 
                if(origKmer.current_kmer < origKmer.last_kmer){
                    origKmer.last_kmer = origKmer.current_kmer;
                    origKmer.pos = i + j;
                }
                if(revKmer.current_kmer < revKmer.last_kmer){
                    revKmer.last_kmer = revKmer.current_kmer;
                    revKmer.pos = i + j;
                }
            }
        }
        if(origKmer.last_kmer <= revKmer.last_kmer){ 
            minimizer = make_tuple(origKmer.last_kmer, origKmer.pos, true);
        }else{
            minimizer = make_tuple(revKmer.last_kmer, revKmer.pos, false);
        }

        return minimizer;
}

tuple<unsigned int, unsigned int, bool> findEndMinimizers(unsigned long int orig, unsigned long int revComp, unsigned int sequence_len, 
    unsigned int kmer_len, unsigned int window_len, unsigned long int len, unsigned int mask, bool begin){

        tuple<unsigned int, unsigned int, bool> endMinimizer;
        kmer origKmer;
        kmer revKmer;

        if(begin){ //find minimizers at the beginning
            for(int i = 0; i < (len - kmer_len + 1); i++){
                origKmer.current_kmer = (orig & (mask << (len - kmer_len + i) * 2)) >> ((len - kmer_len + i) * 2);
                revKmer.current_kmer = (revComp & (mask << (len - kmer_len + i) * 2)) >> ((len - kmer_len + i) * 2);
                if(i == 0){
                    origKmer.last_kmer = origKmer.current_kmer;
                    revKmer.last_kmer = revKmer.current_kmer;
                    origKmer.pos = i;
                    revKmer.pos = i;

                }else{ 
                    if(origKmer.current_kmer < origKmer.last_kmer){
                        origKmer.last_kmer = origKmer.current_kmer;
                        origKmer.pos = i;
                    }
                    if(revKmer.current_kmer < revKmer.last_kmer){
                        revKmer.last_kmer = revKmer.current_kmer;
                        revKmer.pos = i;
                    }
                }
            }
        }else{ //find minimizers at the end
            unsigned long int maskEnd = 0;
            for(int i = 0; i < len; i++){
                maskEnd = maskEnd << 2 | 3;
            }
            orig = orig & maskEnd;
            revComp = revComp & maskEnd;
            for(int i = 0; i < (len - kmer_len + 1); i++){
                origKmer.current_kmer = (orig & (mask << (len - kmer_len + i) * 2)) >> ((len - kmer_len + i) * 2);
                revKmer.current_kmer = (revComp & (mask << (len - kmer_len + i) * 2)) >> ((len - kmer_len + i) * 2);
                if(i == 0){
                    origKmer.last_kmer = origKmer.current_kmer;
                    revKmer.last_kmer = revKmer.current_kmer;
                    origKmer.pos = sequence_len - len;
                    revKmer.pos = sequence_len - len;
                }else{ 
                    if(origKmer.current_kmer < origKmer.last_kmer){
                        origKmer.last_kmer = origKmer.current_kmer;
                        origKmer.pos = sequence_len - len + i;
                    }
                    if(revKmer.current_kmer < revKmer.last_kmer){
                        revKmer.last_kmer = revKmer.current_kmer;
                        revKmer.pos = sequence_len - len + i;
                    }
                }
            }
        }

        if(origKmer.last_kmer < revKmer.last_kmer){ 
            endMinimizer = make_tuple(origKmer.last_kmer, origKmer.pos, true);
        }else{
            endMinimizer = make_tuple(revKmer.last_kmer, revKmer.pos, false);
        }
    return endMinimizer;
}

vector<tuple<unsigned int, unsigned int, bool>> Minimize(
    const char* sequence, unsigned int sequence_len, unsigned int kmer_len, unsigned int window_len){

        unsigned long int mask = 0;
        unsigned long int orig = 0;
        unsigned long int revComp = 0;
        vector<tuple<unsigned int, unsigned int, bool>> minimizers;

        //mask for finding candidates for minimizers
        for(int i = 0; i < kmer_len ; i++){
            mask = mask << 2 | 3;
        }

        //creating kmers and finding minimizers
        for(int i = 0; i < (sequence_len - (window_len + kmer_len - 1) + 1); i++){
            for(int j = 0; j < (window_len + kmer_len - 1); j++){
                orig = (orig << 2) | nucleotide_value(sequence[i+j], false);
                revComp = (revComp << 2) | nucleotide_value(sequence[i+j], true);
            }
            minimizers.push_back(findMinimizers(orig, revComp, sequence_len, kmer_len, window_len, i, mask));
            orig = 0;
            revComp = 0;
        }

        unsigned int len = kmer_len;
        orig = 0;
        revComp = 0;
        while(len < (window_len + kmer_len - 1)){
            for(int i = 0; i < len; i++){
                orig = (orig << 2) | nucleotide_value(sequence[i], false);
                revComp = (revComp << 2) | nucleotide_value(sequence[i], true);
            }
            minimizers.push_back(findEndMinimizers(orig, revComp, sequence_len, kmer_len, window_len, len, mask, true));
            orig = 0;
            revComp = 0;
            len++;
        }

        len = kmer_len;
        while(len < (window_len + kmer_len - 1)){
            for(int i = 0; i < len; i++){
                orig = (orig << 2) | nucleotide_value(sequence[sequence_len - len + i], false);
                revComp = (revComp << 2) | nucleotide_value(sequence[sequence_len - len + i], true);
            }
            minimizers.push_back(findEndMinimizers(orig, revComp, sequence_len, kmer_len, window_len, len, mask, false));
            orig = 0;
            revComp = 0;
            len++;
        }

        sort(minimizers.begin(), minimizers.end());
        minimizers.erase(unique(minimizers.begin(), minimizers.end()), minimizers.end());
        return minimizers;
};
} //namespace pink
//int main(){
//    vector<tuple<unsigned int, unsigned int, bool>> mins;

    /* char sequence[] = {'T', 'G', 'A', 'C', 'G', 'T', 'A', 'C', 'A', 'T', 'G', 'G', 'A', 'C', 'A'};
    mins = Minimize(sequence, sizeof(sequence), 3, 3); */
    //ocekivani minimizeri: 2 10 0, 6 4 0, 6 7 1, 11 1 0, 11 11 0, 17 12 1, 18 0 0

    /* char sequence[] = {'T', 'A', 'G', 'G', 'A', 'A', 'C', 'T', 'G', 'A', 'T', 'T', 'G', 'C', 'C', 'T', 'T', 'A'};
    mins = Minimize(sequence, sizeof(sequence), 4, 3); */
    //ocekivani minizeri: 10 2 0, 10 13 1, 37 8 0, 41 14 1, 43 3 0, 45 6 1, 61 12 0, 79 11 0, 96 0 0

    /* char sequence[] = {'A', 'G', 'C', 'T', 'T', 'T', 'T', 'C', 'A', 'T', 'T', 'C', 'T', 'G', 'A', 'C', 'T', 'G', 'C', 'A', 'A', 'C', 'G', 'G', 'G', 'C', 'A', 'A', 'T', 'A'};
    mins = Minimize(sequence, sizeof(sequence), 5, 3); */
    //ocekivani ispis: 14 22 0, 58 23 0, 83 18 1, 89 25 1, 104 7 1, 170 2 1, 177 15 1, 180 11 1, 180 13 0, 213 1 0, 252 21 1, 343 3 0, 350 4 0, 372 9 0, 458 0 1

    /* char sequence[] = {'T', 'G', 'T', 'C', 'T', 'C', 'T', 'G', 'T', 'G', 'T', 'G', 'G', 'A', 'T', 'T', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'G', 'A', 'G', 'T', 'G', 'T', 'C'};
    mins = Minimize(sequence, sizeof(sequence), 3, 3); */

    /* char sequence[] = {'A','G','C','T','T','T','T','C','A','T','T','C','T','G','A','C','T','G','C','A','A','C','G','G','G','C','A','A','T','A','T','G','T','C','T','C','T','G','T','G','T','G','G','A','T','T','A','A','A','A','A','A','A','G','A','G','T','G','T','C','T','G','A','T','A','G','C','A','G','C','T','T','C','T','G','A','A','C','T','G'};
    mins = Minimize(sequence, sizeof(sequence), 5, 3); */

    /* char sequence[] = {'T', 'G', 'G', 'T', 'T', 'T', 'T', 'C', 'A', 'T', 'T', 'C', 'T', 'G', 'A', 'C', 'T', 'G', 'C', 'A', 'A', 'C', 'G', 'G', 'G', 'C', 'A', 'A', 'T', 'A'};
    mins = Minimize(sequence, sizeof(sequence), 5, 3); */

    /* char sequence[] = {'A','G','C','T','T','T','T','C','A','T','T','C','T','G','A','C','T','G','C','A','A','C','G','G','G','C','A','A','T','A','T','G','T','C','T','C','T','G','T','G','T','G','G','A','T','T','A','A','A','A','A','A','A','G','A','G','T','G','T','C','T','G','A','T','A','G','C','A','G','C','T','T','C','T','G','A','A','C','T','G'};
    mins = Minimize(sequence, sizeof(sequence), 15, 5); */

/*    for(auto min : mins){
        cout << get<0>(min) << " " << get<1>(min) << " " << get<2>(min) << endl;
    }
    return 0;
}*/