#include <iostream>
#include "brown_minimizer.hpp"
#include <string.h>
#include <tuple>
#include <cmath>

unsigned int getKmerValue(char* kmer) {
    unsigned int value = 0;
    int i = 0;
    while(kmer[i + 1] != '\0') {
        switch(kmer[i]) {
            case 'A':
                value = (value + 0) << 2;
                break;
            case 'C':
                value = (value + 1) << 2;
                break;
            case 'G':
                value = (value + 2) << 2;
                break;
            case 'U':
                value = (value + 3) << 2;
                break;
            case 'T':
                value = (value + 3) << 2;
                break;
            default:
                std::cerr << "Wrong base in sequence! "<< std::endl;
                std::exit(EXIT_FAILURE);
            }
        i++;
    }
    switch(kmer[i]) {
        case 'A':
            value += 0;
            break;
        case 'C':
            value += 1;
            break;
        case 'G':
            value += 2;
            break;
        case 'U':
            value += 3;
            break;
        case 'T':
            value += 3;
            break;
        default:
            std::cerr << "Wrong base in sequence! "<< std::endl;
            std::exit(EXIT_FAILURE);
        }
    //std::cout << "ide radit stoi\n";
    return value;
}

namespace brown {

    unsigned int getReversedComplKmerValue(unsigned int value, unsigned int length) {
        int number_of_bits = floor(log2(value)) + 1; 

        unsigned int n = ((1 << number_of_bits) - 1) ^ value; 

        unsigned int ans = 0;
        for(int i = length * 2 - 2; i >= 0; i-= 2){
            ans |= (n & 3) << i;
            n >>= 2;
        }
        return ans;
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

        //std::cout << "krece radit minimizere\n";
        //stavljanje pocetnih minimizera
        for (unsigned int i = 0; i < window_len - 1; i++) {
            strncpy(kmer, sequence + i, kmer_len);
            kmer_noreverse_value = getKmerValue(kmer);
            //std::cout << "treba zavrsit stoi\n";
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
                kmer_noreverse_value = getKmerValue(kmer);
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
                    kmer_noreverse_value = getKmerValue(kmer);
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
                    kmer_noreverse_value = getKmerValue(kmer);
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




        





                    

