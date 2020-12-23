#include <iostream>
#include "brown_minimizer.hpp"
#include <string.h>
#include <tuple>


int getKmerValue(char* kmer, bool origin) {
    std::string numbers = "";
    int i = 0;
    while(kmer[i] != '\0') {
        switch(kmer[i]) {
            case 'C':
                origin == true ? numbers += '1' : numbers += '4';
                break;
            case 'A':
                origin == true ? numbers += '2' : numbers += '3';
                break;
            case 'U':
                origin == true ? numbers += '3' : numbers += '2';
                break;
            case 'T':
                origin == true ? numbers += '3' : numbers += '2';
                break;
            case 'G':
                origin == true ? numbers += '4' : numbers += '1';
                break;
            default:
                std::cerr << "Wrong base in sequence! "<< std::endl;
                std::exit(EXIT_FAILURE);
            }
        i++;
    }
    return std::stoi(numbers);
}

    std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
        const char* sequence, unsigned int sequence_len,
        unsigned int kmer_len,
        unsigned int window_len) {  //pokusat nac ljepsi nacin za izvodenje
                    
        std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;
        int fraction_length;
        char *fraction;
        bool origin;

        char* kmer = (char*) malloc(kmer_len + 1);
        kmer[kmer_len] = '\0';
        unsigned int kmer_noreverse_value, kmer_reverse_value, kmer_value;

        //stavljanje pocetnih minimizera
        for (unsigned int i = 1; i < window_len; i++) {
            strncpy(kmer, sequence + i - 1, kmer_len);
            kmer_noreverse_value = getKmerValue(kmer, true);
            kmer_reverse_value = getKmerValue(kmer, false);
            if (kmer_noreverse_value < kmer_reverse_value) {
                kmer_value = kmer_noreverse_value;
                origin = true;
            } else {
                kmer_value = kmer_reverse_value;
                origin = false;
            }
            if (i == 1 || kmer_value <= std::get<0>(minimizers[minimizers.size() - 1])) {
                minimizers.push_back(std::make_tuple(kmer_value, i - 1, origin));
            }
        }

        //stavljanje unutarnjih minimizera
        fraction_length = window_len + kmer_len - 1;
        for (unsigned int i = 0; i < sequence_len - fraction_length; i++) {
            if (i == 0 || std::get<1>(minimizers[minimizers.size() - 1]) >= i) {
                strncpy(kmer, sequence + i + fraction_length - 1 - kmer_len, kmer_len);
                kmer_noreverse_value = getKmerValue(kmer, true);
                kmer_reverse_value = getKmerValue(kmer, false);
                if (kmer_noreverse_value < kmer_reverse_value) {
                    kmer_value = kmer_noreverse_value;
                    origin = true;
                } else {
                    kmer_value = kmer_reverse_value;
                    origin = false;
                }
                if (kmer_value <= std::get<0>(minimizers[minimizers.size() - 1])) {
                    minimizers.push_back(std::make_tuple(kmer_value, i + fraction_length - kmer_len - 1, origin));
                }
            } else {
                int minimizer_value, minimizer_position;
                bool minimizer_origin;
                for (int j = i; j + kmer_len <= fraction_length; j++) {                    
                    strncpy(kmer, sequence + j, kmer_len);
                    kmer_noreverse_value = getKmerValue(kmer, true);
                    kmer_reverse_value = getKmerValue(kmer, false);
                    if (kmer_noreverse_value < kmer_reverse_value) {
                        kmer_value = kmer_noreverse_value;
                        origin = true;
                    } else {
                        kmer_value = kmer_reverse_value;
                        origin = false;
                    }
                    if (j == i || kmer_value < minimizer_value) {
                        minimizer_value = kmer_value;
                        minimizer_position = j;
                        minimizer_origin = origin;
                    }
                }
                minimizers.push_back(std::make_tuple(minimizer_value, minimizer_position, minimizer_origin));
                
            }
        }

        //stavljanje krajnih minimizera, treba popravit
        for (unsigned int i = sequence_len - window_len - 1; i < sequence_len; i++) {
            if (std::get<1>(minimizers[minimizers.size() - 1]) < i) {
                int minimizer_value, minimizer_position;
                bool minimizer_origin;
                fraction_length = sequence_len - i - 1 + kmer_len - 1;
                for (int j = i; j + kmer_len <= fraction_length; j++) {                    
                    strncpy(kmer, sequence + j, kmer_len);
                    kmer_noreverse_value = getKmerValue(kmer, true);
                    kmer_reverse_value = getKmerValue(kmer, false);
                    if (kmer_noreverse_value < kmer_reverse_value) {
                        kmer_value = kmer_noreverse_value;
                        origin = true;
                    } else {
                        kmer_value = kmer_reverse_value;
                        origin = false;
                    }
                    if (j == i || kmer_value < minimizer_value) {
                        minimizer_value = kmer_value;
                        minimizer_position = j;
                        minimizer_origin = origin;
                    }
                }
                minimizers.push_back(std::make_tuple(minimizer_value, minimizer_position, minimizer_origin));
            }
            
        }

        return minimizers;
    }





                    

