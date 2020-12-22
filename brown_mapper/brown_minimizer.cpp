#include <iostream>
#include "brown_minimizer.hpp"
#include <string.h>
#include <tuple>

namespace brown {

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
                    int minimizerValue, minimizerPosition;
                    bool origin;

                    //minimizeri s pocetka
                    for (int i = 1; i < window_len; i++) {
                        fraction_length = i + kmer_len - 1;
                        fraction = (char*) malloc(fraction_length + 1);
                        strncpy(fraction, sequence, fraction_length);
                        fraction[fraction_length] = '\0';

                        for (int j = 0; j + kmer_len <= fraction_length; j++) {
                            char *kmer = (char*) malloc(kmer_len + 1);
                            strncpy(kmer, fraction + j, kmer_len);
                            kmer[kmer_len] = '\0';
                            int kmerValueNoReverse = getKmerValue(kmer, true);
                            int kmerValueReverse = getKmerValue(kmer, false);
                            int kmerValue = kmerValueReverse < kmerValueNoReverse ? kmerValueReverse : kmerValueNoReverse;
                            if (j == 0 || kmerValue < minimizerValue) {
                                minimizerValue = kmerValue;
                                minimizerPosition = j;
                                origin = kmerValueReverse < kmerValueNoReverse ? false : true;
                            }
                            free(kmer);
                        }
                        std::tuple<unsigned int, unsigned int, bool> tuple(minimizerValue, minimizerPosition, origin); 
                        minimizers.push_back(tuple);
                        free(fraction);
                    }

                    //obicni minimizeri
                    fraction_length = kmer_len + window_len - 1;
                    for (int i = 0; i < sequence_len - fraction_length; i++) {
                        fraction = (char*) malloc(fraction_length + 1);
                        strncpy(fraction, sequence + i, fraction_length);
                        fraction[fraction_length] = '\0';

                        for (int j = 0; j + kmer_len < fraction_length; j++) {
                            char *kmer = (char*) malloc(kmer_len + 1);
                            strncpy(kmer, fraction + j, kmer_len);
                            kmer[kmer_len] = '\0';
                            int kmerValueNoReverse = getKmerValue(kmer, true);
                            int kmerValueReverse = getKmerValue(kmer, false);
                            int kmerValue = kmerValueReverse < kmerValueNoReverse ? kmerValueReverse : kmerValueNoReverse;
                            if (j == 0 || kmerValue < minimizerValue) {
                                minimizerValue = kmerValue;
                                minimizerPosition = i + j;
                                origin = kmerValueReverse < kmerValueNoReverse ? false : true;
                            }
                            free(kmer);
                        }

                        std::tuple<unsigned int, unsigned int, bool> tuple(minimizerValue, minimizerPosition, origin); 
                        minimizers.push_back(tuple);
                        free(fraction);
                    }

                    //minimizeri s kraja
                    for (int i = window_len - 1; i > 0; i--) {
                        fraction_length = i + kmer_len - 1;
                        fraction = (char*) malloc(fraction_length + 1);
                        strncpy(fraction, sequence + sequence_len - 1 - fraction_length, fraction_length);
                        fraction[fraction_length] = '\0';
                        
                        for (int j = 0; j + kmer_len <= fraction_length; j++) {
                            char *kmer = (char*) malloc(kmer_len + 1);
                            strncpy(kmer, fraction + j, kmer_len);
                            kmer[kmer_len] = '\0';
                            int kmerValueNoReverse = getKmerValue(kmer, true);
                            int kmerValueReverse = getKmerValue(kmer, false);
                            int kmerValue = kmerValueReverse < kmerValueNoReverse ? kmerValueReverse : kmerValueNoReverse;
                            if (j == 0 || kmerValue < minimizerValue) {
                                minimizerValue = kmerValue;
                                minimizerPosition = sequence_len - 1 - fraction_length + j;
                                origin = kmerValueReverse < kmerValueNoReverse ? false : true;
                            }
                            free(kmer);
                        }

                        std::tuple<unsigned int, unsigned int, bool> tuple(minimizerValue, minimizerPosition, origin); 
                        minimizers.push_back(tuple);
                        free(fraction);
                    }

                    return minimizers;

                }
    
}