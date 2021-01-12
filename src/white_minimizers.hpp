#include <iostream>
#include <vector>
#include <string>
#include <tuple>

namespace white {

    std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
    const char* sequence, unsigned int sequence_len,
    unsigned int kmer_len, //u windowu gleda po k slova
    unsigned int window_len //duljina prozora u kojem gleda po k slova
    )
    {
        //rezultat
         std::vector<std::tuple<unsigned int, unsigned int, bool>> returnVec;

        std::string chToStr(sequence);
        std::string chToStrComp(complementSeq(sequence, sequence_len));

        std::tuple<unsigned int, unsigned int, bool> minTuple;

        //ocitaj sve kmerove
        // ima ih sequence_len - kmer_len + 1
        std::vector<std::tuple<unsigned int, unsigned int, bool>> all_kmers (sequence_len - kmer_len + 1);

        //end-minimizers
        //TODO: smanjii kod
            //original
        std::string kmer_first = kmerNumbers(chToStr.substr(0, kmer_len-1)); 
        std::string kmer_last = kmerNumbers(chToStr.substr(sequence_len-kmer_len, sequence_len-1));
            //complement
        std::string kmer_first_c = kmerNumbers(chToStrComp.substr(0, kmer_len-1)); 
        std::string kmer_last_c = kmerNumbers(chToStrComp.substr(sequence_len-kmer_len, sequence_len-1));

        for (int i = 1; i < window_len - 1; i++) {

            std::string new_first = kmerNumbers(chToStr.substr(i, kmer_len + (i - 1)));
            std::string new_first_c = kmerNumbers(chToStrComp.substr(i, kmer_len + (i - 1)));

            if (std::stoul(new_first) < std::stoul(kmer_first)) {
                kmer_first = new_first;
            }
            if (std::stoul(new_first_c) < std::stoul(kmer_first_c)) {
                kmer_first_c = new_first_c;
            }

            std::string new_last = kmerNumbers(chToStr.substr(sequence_len-kmer_len-i, sequence_len-1-i));
            std::string new_last_c = kmerNumbers(chToStrComp.substr(sequence_len-kmer_len-i, sequence_len-1-i));

            if (std::stoul(new_last) < std::stoul(kmer_last)) {
                kmer_last = new_last;
            }
            if (std::stoul(new_last) < std::stoul(kmer_last)) {
                kmer_last_c = new_last_c;
            }

            if(std::stoul(kmer_first_c) < std::stoul(kmer_first)) {
                minTuple = {std::stoul(kmer_first_c), i , false};
            } else {
                minTuple = {std::stoul(kmer_first), i , true};
            }
            returnVec.push_back(minTuple);

            if(std::stoul(kmer_last_c) < std::stoul(kmer_last)) {
                minTuple = {std::stoul(kmer_last_c), i , false};
            } else {
                minTuple = {std::stoul(kmer_last), i , true};
            }
            returnVec.push_back(minTuple);
        }


        //interior-minimizers TODO!!

        //fraction is win_len + kmer_len - 1 "a set of w consecutive k-mers covers a string of exactly w +k −1 letters"
        // for i od 0 do seq_len - fraction
        //zadnja pozicija u rj > i onda od (seq tj pocetak) + i + wind_len - 1 substring

        int set_fraction = window_len + kmer_len - 1;
        
        std::string kmer = "";
        std::string kmer_Comp = "";

        for (int i = 0; i < sequence_len - set_fraction; i++) {

            int poz = i + window_len - 1;
            kmer = kmerNumbers(chToStr.substr(poz, kmer_len + poz)); //next kmer numeric
            kmer_Comp = kmerNumbers(chToStrComp.substr(poz, kmer_len + poz));

            if (std::stoul(kmer_Comp) < std::stoul(kmer)) {
                minTuple = {std::stoul(kmer_Comp), i , false};
            } else
            {
                minTuple = {std::stoul(kmer), i , true};
            }
            returnVec.push_back(minTuple);
            
        }

    }

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

    const std::string kmerNumbers(std::string kmer) {
        std::string rez = "";
        for (int a = 0; a < kmer.length(); a++) {
            switch(kmer[a]) {
            case 'C':
                rez += '0';
            case 'A':
                rez += '1';
            case 'T':
                rez += '2';
            case 'G':
                rez += '3';
            default:
                return "kriva baza"; //GREŠKA
            }
        }
        return rez;  
    }

    int letterMapping(char letter) {
        switch(letter) {
            case 'C':
                return 0;
            case 'A':
                return 1;
            case 'T':
                return 2;
            case 'G':
                return 3;
            default:
                return -1; //GREŠKA
        }
    }


    int compareLetter(std::string first, std::string second) {
        int counter = 0;
        for (int i = 0; i < first.length(); i++) {
            if (letterMapping(first[i]) < letterMapping(second[i])) {
                return -1;
            }
            if (letterMapping(first[i]) == letterMapping(second[i])) {
                continue;
            }
            else {
                return 1;
            }
        }
        return 0;
    }
}