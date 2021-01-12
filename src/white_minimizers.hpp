#include <iostream>
#include <vector>
#include <string>
#include <tuple>

namespace white {

    const char* complementSeq(const char* sequence, unsigned int sequence_len) {
        //A -> T
        //T -> A
        //G -> C
        //C -> G

        const char* compSeq = "";
        for (int i = 0; i < sequence_len; i++) {
            if (sequence[i] == 'A') {
                compSeq += 'T';
            }
            if (sequence[i] == 'T') {
                compSeq += 'A';
            }
            if (sequence[i] == 'G') {
                compSeq += 'C';
            }
            if (sequence[i] == 'C') {
                compSeq += 'G';
            }
        }

        return compSeq;
    }

    std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
    const char* sequence, unsigned int sequence_len,
    unsigned int kmer_len, //u windowu gleda po k slova
    unsigned int window_len //duljina prozora u kojem gleda po k slova
    ) {

        std::string chToStr(sequence);
        std::string chToStrComp(complementSeq(sequence, sequence_len));
        //window_len = 7, kmer_len = 3
        int tupleCounter = 1;
        std::string min = "";
        std::string minComp = "";
        std::tuple<unsigned int, unsigned int, bool> minTuple;
        std::tuple<unsigned int, unsigned int, bool> minTupleComp;
        std::vector<std::tuple<unsigned int, unsigned int, bool>> returnVec;
            for (int i = 0; i < sequence_len / window_len; i++) {

                std::string subSeq = chToStr.substr(window_len*i, window_len*(i+1)); //next window
                for (int j = 0; j < window_len-(kmer_len-1); j++) {
                    std::string kmer = subSeq.substr(j, kmer_len-1 + j); //next kmer
                    if (min == "" || kmer.compare(min) < 0) {
                        min = kmer;
                        minTuple = {-1, window_len*i + j, true};
                    }
                }


                subSeq = chToStrComp.substr(window_len*i, window_len*(i+1));
                for (int j = 0; j < window_len-(kmer_len-1); j++) { //complement loop
                    std::string kmer = subSeq.substr(j, kmer_len-1 + j);
                    if (minComp == "" || kmer.compare(minComp) < 0) {
                        minComp = kmer;
                        minTupleComp = {-1, window_len*i + j, false};
                    }
                }

                //compare which kmer is alphabetically smaller
                if (min.compare(minComp) < 0) {
                    minTuple = {tupleCounter, std::get<1>(minTuple), std::get<2>(minTuple)};
                } 
                else {
                    minTuple = {tupleCounter, std::get<1>(minTupleComp), std::get<2>(minTupleComp)};
                }

                tupleCounter++;
                returnVec.push_back(minTuple);
            }

    }
  
}