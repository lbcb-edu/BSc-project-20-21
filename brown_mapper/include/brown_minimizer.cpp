#include <iostream>
#include "brown_minimizer.hpp"
#include <string.h>
#include <tuple>
#include <cmath>

unsigned int getKmerValue(char* kmer, unsigned int kmer_len) {
    unsigned int value = 0;
    unsigned int i = 0;
    while(i < kmer_len) {
        switch(kmer[i]) {
            case 'A':
                value = (value << 2) + 0;
                break;
            case 'C':
                value = (value << 2) + 1;
                break;
            case 'G':
                value = (value << 2) + 2;
                break;
            case 'U':
                value = (value << 2) + 3;
                break;
            case 'T':
                value = (value << 2) + 3;
                break;
            default:
                std::cerr << "Wrong base in sequence -> " << kmer[i] << " at " << i << std::endl;
                std::exit(EXIT_FAILURE);
            }
        i++;
    }
    return value;
}

namespace brown {

    unsigned int getReversedComplKmerValue(unsigned int value, unsigned int length) {
        unsigned int ans = 0;;
        for(int i = length * 2 - 2; i >= 0; i-= 2){
            ans |= (value & 3) << i;
            value >>= 2;
        }
        unsigned int n = ((1 << (length * 2)) - 1) ^ ans;
        return n;
    }

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
            kmer_noreverse_value = getKmerValue(kmer, kmer_len);
            //kmer_reverse_value = getReversedKmerValue(kmer_noreverse_value);
        
            //if (kmer_noreverse_value < kmer_reverse_value) {
                kmer_value = kmer_noreverse_value;
                origin = true;
           /*  } else {
                kmer_value = kmer_reverse_value;
                origin = false;
            } */
            
            if (i == 0 || kmer_value <= std::get<0>(minimizers.back())) 
                minimizers.push_back(std::make_tuple(kmer_value, i ,origin));
        
        }

        //std::cout << "gotovi pocetni minimizeri\n" ;

        //stavljanje unutarnjih minimizera
        fraction_length = window_len + kmer_len - 1;
        for (unsigned int i = 0; i <= sequence_len - fraction_length; i++) {
            if (std::get<1>(minimizers.back()) >= i) {
                strncpy(kmer, sequence + i + fraction_length - kmer_len, kmer_len);
                kmer_noreverse_value = getKmerValue(kmer, kmer_len);
                //kmer_reverse_value = getReversedKmerValue(kmer_noreverse_value);
                
                //if (kmer_noreverse_value < kmer_reverse_value) {
                    kmer_value = kmer_noreverse_value;
                    origin = true;
                /* } else {
                    kmer_value = kmer_reverse_value;
                    origin = false;
                } */
                
                if (kmer_value <= std::get<0>(minimizers.back())) 
                    minimizers.push_back(std::make_tuple(kmer_value, i + fraction_length - kmer_len, origin));
            
            } else {
                unsigned int minimizer_value, minimizer_position;
                bool minimizer_origin;
                
                for (unsigned int j = 0; j <= fraction_length -  kmer_len; j++) {                    
                    strncpy(kmer, sequence + i + j, kmer_len);
                    kmer_noreverse_value = getKmerValue(kmer, kmer_len);
                    //kmer_reverse_value = getReversedKmerValue(kmer_noreverse_value);
                    
                    //if (kmer_noreverse_value < kmer_reverse_value) {
                        kmer_value = kmer_noreverse_value;
                        origin = true;
                    /* } else {
                        kmer_value = kmer_reverse_value;
                        origin = false;
                    } */
                    
                    if (j == 0 || kmer_value < minimizer_value) {
                        minimizer_value = kmer_value;
                        minimizer_position = j;
                        minimizer_origin = origin;
                    }
                }
                minimizers.push_back(std::make_tuple(minimizer_value, minimizer_position + i, minimizer_origin));
            
            }
        }

        //stavljanje krajnih minimizera
        for (unsigned int i = sequence_len - window_len - kmer_len + 2; i <= sequence_len - kmer_len; i++) {
            if (std::get<1>(minimizers.back()) < i) {
                unsigned int minimizer_value, minimizer_position;
                bool minimizer_origin;
                
                for (unsigned int j = i; j <= sequence_len - kmer_len; j++) {                    
                    strncpy(kmer, sequence + j, kmer_len);
                    kmer_noreverse_value = getKmerValue(kmer, kmer_len);
                    //kmer_reverse_value = getReversedKmerValue(kmer_noreverse_value);
                    
                    //if (kmer_noreverse_value < kmer_reverse_value) {
                        kmer_value = kmer_noreverse_value;
                        origin = true;
                    /* } else {
                        kmer_value = kmer_reverse_value;
                        origin = false;
                    } */

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




        





                    

