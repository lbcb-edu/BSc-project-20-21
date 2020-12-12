#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
using namespace std;

namespace pink {

struct kmer{
    unsigned int k_mer = 0;
    unsigned int current_kmer = 0;
    unsigned int last_kmer = 0;
    unsigned int pos = 0;
};

int nucleotide_value(char c, bool reverse){
    if(c == 'C') return reverse ? 3 : 0;
    else if(c == 'A') return reverse ? 2 : 1;
    else if(c == 'T') return reverse ? 1 : 2;
    else if(c == 'G') return reverse ? 0 : 3;
    else return -1;
}

vector<tuple<unsigned int, unsigned int, bool>> findMinimizers(unsigned long int orig, unsigned long int revComp, unsigned int sequence_len, 
    unsigned int kmer_len, unsigned int window_len, int mask, int mask2){

        vector<tuple<unsigned int, unsigned int, bool>> minimizers;
        kmer origKmer;
        kmer revKmer;
        unsigned int k = 0;
        
        for(int i = 0; i < sequence_len - (window_len + kmer_len - 1) + 1; i++){
            origKmer.k_mer = (orig & (mask << (sequence_len - (window_len + kmer_len + i) + 1) * 2)) >> ((sequence_len - (window_len + kmer_len + i) + 1) * 2);
            revKmer.k_mer = (revComp & (mask << (sequence_len - (window_len + kmer_len + i) + 1) * 2)) >> ((sequence_len - (window_len + kmer_len + i) + 1) * 2);
            for(int j = 0; j < window_len; j++){
                origKmer.current_kmer = (origKmer.k_mer & (mask2 << (window_len - 1 - j) * 2))  >> (window_len - 1 - j) * 2;
                revKmer.current_kmer = (revKmer.k_mer & (mask2 << (window_len - 1 - j) * 2))  >> (window_len - 1 - j) * 2;
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
                minimizers.push_back(make_tuple(origKmer.last_kmer, origKmer.pos, true));
            }else{
                minimizers.push_back(make_tuple(revKmer.last_kmer, revKmer.pos, false));
            }
        }
        return minimizers;
}

vector<tuple<unsigned int, unsigned int, bool>> findEndMinimizers (unsigned long int orig, unsigned long int revComp, unsigned int sequence_len, 
    unsigned int kmer_len, unsigned int window_len, int mask, int mask2){

        unsigned int len = kmer_len;
        vector<tuple<unsigned int, unsigned int, bool>> endMinimizers;
        kmer origKmer;
        kmer revKmer;

        //find minimizers at the beginning
        while(len < (window_len + kmer_len - 1)){
            origKmer.k_mer = orig >> ((sequence_len - len) * 2);
            revKmer.k_mer = revComp >> ((sequence_len - len) * 2);
            for(int i = 0; i < (len - kmer_len + 1); i++){
                origKmer.current_kmer = (origKmer.k_mer & (mask2 << (len - kmer_len - i) * 2)) >> ((len - kmer_len - i) * 2);
                revKmer.current_kmer = (revKmer.k_mer & (mask2 << (len - kmer_len - i) * 2)) >> ((len - kmer_len - i) * 2);
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
            if(origKmer.last_kmer < revKmer.last_kmer){ 
                endMinimizers.push_back(make_tuple(origKmer.last_kmer, origKmer.pos, true));
            }else{
                endMinimizers.push_back(make_tuple(revKmer.last_kmer, revKmer.pos, false));
            }
            len++;
        }

        //find minimizers at the end
        len = kmer_len;
        while(len < (window_len + kmer_len - 1)){
            int maskEnd = 0;
            for(int i = 0; i < len; i++){
                maskEnd = maskEnd << 2 | 3;
            }
            origKmer.k_mer = orig & maskEnd;
            revKmer.k_mer = revComp & maskEnd;
            for(int i = 0; i < (len - kmer_len + 1); i++){
                origKmer.current_kmer = (origKmer.k_mer & (mask2 << (len - kmer_len - i) * 2)) >> ((len - kmer_len - i) * 2);
                revKmer.current_kmer = (revKmer.k_mer & (mask2 << (len - kmer_len - i) * 2)) >> ((len - kmer_len - i) * 2);
                if(i == 0){
                    origKmer.last_kmer = origKmer.current_kmer;
                    revKmer.last_kmer = revKmer.current_kmer;
                    origKmer.pos = sequence_len - len - i;
                    revKmer.pos = sequence_len - len - i;
                }else{ 
                    if(origKmer.current_kmer < origKmer.last_kmer){
                        origKmer.last_kmer = origKmer.current_kmer;
                        origKmer.pos = sequence_len - len - i;
                    }
                    if(revKmer.current_kmer < revKmer.last_kmer){
                        revKmer.last_kmer = revKmer.current_kmer;
                        revKmer.pos = sequence_len - len - i;
                    }
                }
            }
            if(origKmer.last_kmer < revKmer.last_kmer){ 
                endMinimizers.push_back(make_tuple(origKmer.last_kmer, origKmer.pos, true));
            }else{
                endMinimizers.push_back(make_tuple(revKmer.last_kmer, revKmer.pos, false));
            }
            len++;
        }
    return endMinimizers;
}

vector<tuple<unsigned int, unsigned int, bool>> Minimize(
    const char* sequence, unsigned int sequence_len, unsigned int kmer_len, unsigned int window_len){

        unsigned int mask = 0;
        unsigned int mask2 = 0;
        unsigned long int orig = 0;
        unsigned long int revComp = 0;
        vector<tuple<unsigned int, unsigned int, bool>> minimizers;
        vector<tuple<unsigned int, unsigned int, bool>> endMinimizers;

        //mask for creating kmers
        for(int i = 0; i < (window_len + kmer_len - 1) ; i++){
            mask = mask << 2 | 3;
        }

        //mask for finding candidates for minimizers
        for(int i = 0; i < kmer_len ; i++){
            mask2 = mask2 << 2 | 3;
        }

        //assign nucleotide values for original strand and reverse complement
        for(int i = 0; i < sequence_len; i++){
            orig = (orig << 2) | nucleotide_value(sequence[i], false);
            revComp = (revComp << 2) | nucleotide_value(sequence[i], true);
        }
        
        endMinimizers = findEndMinimizers(orig, revComp, sequence_len, kmer_len, window_len, mask, mask2);
        endMinimizers.erase(unique(endMinimizers.begin(), endMinimizers.end()), endMinimizers.end());

        minimizers = findMinimizers(orig, revComp, sequence_len, kmer_len, window_len, mask, mask2);
        minimizers.erase(unique(minimizers.begin(), minimizers.end()), minimizers.end());

        for(auto e : endMinimizers){
            minimizers.push_back(e);
        } 

        sort(minimizers.begin(), minimizers.end());
        minimizers.erase(unique(minimizers.begin(), minimizers.end()), minimizers.end());
        return minimizers;
};
} //namespace pink

/*int main(){
    vector<tuple<unsigned int, unsigned int, bool>> mins;

    char sequence[] = {'T', 'G', 'A', 'C', 'G', 'T', 'A', 'C', 'A', 'T', 'G', 'G', 'A', 'C', 'A'};
    mins = Minimize(sequence, sizeof(sequence), 3, 3);
    //ocekivani minimizeri: 2 10 0, 6 4 0, 6 7 1, 11 1 0, 11 11 0, 17 12 1, 18 0 0

    /* char sequence[] = {'T', 'A', 'G', 'G', 'A', 'A', 'C', 'T', 'G', 'A', 'T', 'T', 'G', 'C', 'C', 'T', 'T', 'A'};
    mins = Minimize(sequence, sizeof(sequence), 4, 3); */
    //ocekivani minizeri: 10 2 0, 10 13 1, 37 8 0, 41 14 1, 43 3 0, 45 6 1, 61 12 0, 79 11 0, 96 0 0

    /*for(auto min : mins){
        cout << get<0>(min) << " " << get<1>(min) << " " << get<2>(min) << endl;
    }
    return 0;
}*/
