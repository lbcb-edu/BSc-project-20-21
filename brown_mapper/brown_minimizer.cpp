#include <iostream>
#include "brown_minimizer.hpp"
#include <string.h>
#include <tuple>

unsigned int getKmerValue(char* kmer, bool origin) {
    std::string numbers = "";
    int i = 0;
    while(kmer[i] != '\0') {
        switch(kmer[i]) {
            case 'A':
                origin == true ? numbers += '1' : numbers += '4';
                break;
            case 'C':
                origin == true ? numbers += '2' : numbers += '3';
                break;
            case 'G':
                origin == true ? numbers += '3' : numbers += '2';
                break;
            case 'U':
                origin == true ? numbers += '4' : numbers += '1';
                break;
            case 'T':
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

namespace brown {

    std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
        const char* sequence, unsigned int sequence_len,
        unsigned int kmer_len,
        unsigned int window_len) {  
                    
        std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;
        unsigned int fraction_length;
        bool origin;

        char* kmer = (char*) malloc(kmer_len);
        unsigned int kmer_noreverse_value, kmer_reverse_value, kmer_value;

        //stavljanje pocetnih minimizera
        for (unsigned int i = 0; i < window_len - 1; i++) {
            strncpy(kmer, sequence + i, kmer_len);
            kmer_noreverse_value = getKmerValue(kmer, true);
            kmer_reverse_value = getKmerValue(kmer, false);
        
            if (kmer_noreverse_value < kmer_reverse_value) {
                kmer_value = kmer_noreverse_value;
                origin = true;
            } else {
                kmer_value = kmer_reverse_value;
                origin = false;
            }
            
            if (i == 0 || kmer_value <= std::get<0>(minimizers.back())) 
                minimizers.push_back(std::make_tuple(kmer_value, i ,origin));
        
        }

        //stavljanje unutarnjih minimizera
        fraction_length = window_len + kmer_len - 1;
        for (unsigned int i = 0; i <= sequence_len - fraction_length; i++) {
            if (std::get<1>(minimizers.back()) >= i) {
                strncpy(kmer, sequence + i + fraction_length - kmer_len, kmer_len);
                kmer_noreverse_value = getKmerValue(kmer, true);
                kmer_reverse_value = getKmerValue(kmer, false);
                
                if (kmer_noreverse_value < kmer_reverse_value) {
                    kmer_value = kmer_noreverse_value;
                    origin = true;
                } else {
                    kmer_value = kmer_reverse_value;
                    origin = false;
                }
                
                if (kmer_value <= std::get<0>(minimizers.back())) 
                    minimizers.push_back(std::make_tuple(kmer_value, i + fraction_length - kmer_len, origin));
            
            } else {
                unsigned int minimizer_value, minimizer_position;
                bool minimizer_origin;
                
                for (unsigned int j = 0; j <= fraction_length -  kmer_len; j++) {                    
                    strncpy(kmer, sequence + i + j, kmer_len);
                    kmer_noreverse_value = getKmerValue(kmer, true);
                    kmer_reverse_value = getKmerValue(kmer, false);
                    
                    if (kmer_noreverse_value < kmer_reverse_value) {
                        kmer_value = kmer_noreverse_value;
                        origin = true;
                    } else {
                        kmer_value = kmer_reverse_value;
                        origin = false;
                    }
                    
                    if (j == 0 || kmer_value < minimizer_value) {
                        minimizer_value = kmer_value;
                        minimizer_position = j;
                        minimizer_origin = origin;
                    }
                }
                minimizers.push_back(std::make_tuple(minimizer_value, minimizer_position, minimizer_origin));
            
            }
        }

        //stavljanje krajnih minimizera
        for (unsigned int i = sequence_len - window_len - kmer_len + 2; i <= sequence_len - kmer_len; i++) {
            if (std::get<1>(minimizers.back()) < i) {
                unsigned int minimizer_value, minimizer_position;
                bool minimizer_origin;
                
                for (unsigned int j = i; j <= sequence_len - kmer_len; j++) {                    
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

        free(kmer);

        return minimizers;
    }
}




        





                    

