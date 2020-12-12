#include <iostream>
#include "brown_minimizer.hpp"
#include <string.h>
#include <tuple>

namespace brown {

        int getKmerValueReverse(char* kmer) {
            std::string numbers = "";
            int i = 0;
            while(kmer[i] != '\0') {
                switch(kmer[i]) {
                    case 'G':
                        numbers += '1';
                        break;
                    case 'T':
                        numbers += '2';
                        break;
                    case 'U':
                        numbers += '2';
                        break;
                    case 'A':
                        numbers += '3';
                        break;
                    case 'C':
                        numbers += '4';
                        break;
                    default:
                        std::cerr << "Wrong base in sequence! "<< std::endl;
                        std::exit(EXIT_FAILURE);
                }
                i++;
            }
            return std::stoi(numbers);
        }

        int getKmerValue(char* kmer) {
            std::string numbers = "";
            int i = 0;
            while(kmer[i] != '\0') {
                switch(kmer[i]) {
                    case 'C':
                        numbers += '1';
                        break;
                    case 'A':
                        numbers += '2';
                        break;
                    case 'U':
                        numbers += '3';
                        break;
                    case 'T':
                        numbers += '3';
                        break;
                    case 'G':
                        numbers += '4';
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
                unsigned int window_len) {
                    
                    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;
                    int fraction_length = kmer_len + window_len - 1;
                    char *fraction = (char*) malloc(fraction_length + 1);
                    int minimizerValue, minimizerPosition;
                    bool origin;

                    int i;
                    for (i = 0; i < sequence_len; i += fraction_length) {
                        //const char* fraction(sequence + i, sequence + i + fraction_length);
                        //char *fraction = (char*) malloc(fraction_length + 1);
                        strncpy(fraction, sequence + i, fraction_length);
                        fraction[fraction_length] = '\0';
                        //int minimizerValue, minimizerPosition;
                        //bool origin;

                        for (int j = 0; j < fraction_length; j++) {
                            //std::string kmer(fraction + j, fraction + j + kmer_len)
                            char *kmer = (char*) malloc(kmer_len + 1);
                            strncpy(kmer, fraction + j, kmer_len);
                            kmer[kmer_len] = '\0';
                            int kmerValueNoReverse = getKmerValue(kmer);
                            int kmerValueReverse = getKmerValueReverse(kmer);
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

                    i -= fraction_length;
                    bool first_entry;
                    for(; i < sequence_len; i++) {
                        char *kmer = (char*) malloc(kmer_len + 1);
                        strncpy(kmer, sequence + i, kmer_len);
                        kmer[kmer_len] = '\0';
                        int kmerValueNoReverse = getKmerValue(kmer);
                        int kmerValueReverse = getKmerValueReverse(kmer);
                        int kmerValue = kmerValueReverse < kmerValueNoReverse ? kmerValueReverse : kmerValueNoReverse;
                        if (!first_entry || kmerValue < minimizerValue) {
                            first_entry = true;
                            minimizerValue = kmerValue;
                            minimizerPosition = i;
                            origin = kmerValueReverse < kmerValueNoReverse ? false : true;
                        }
                        free(kmer);
                    }
                    std::tuple<unsigned int, unsigned int, bool> tuple(minimizerValue, minimizerPosition, origin); 
                    minimizers.push_back(tuple);

                    return minimizers;

                }
    
}