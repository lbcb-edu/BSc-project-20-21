
#include <iostream>
#include <vector>
#include <string>
#include <tuple>

namespace white {

 /**
     * @brief komplementira danu sekvencu
     * 
     * @param sequence sekvenca koju zelimo komplementirati
     * @param sequence_len duljina predane sekvence
     * @return komplement predane sekvence
     */
    const char* complementSeq(const char* sequence, unsigned int sequence_len);

    
    unsigned int mapLetter(char letter);

    unsigned int Complement(unsigned int base);

    /**
     * @brief 
     * 
     * @param sequence sekvenca nad kojom provodimo trazenje minimizera
     * @param sequence_len duljina predane sekvence
     * @param kmer_len duljina trazenih K-merova
     * @param window_len velicina prozora za trazenje K-merova
     * @return vektor K-merova
     */
    std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
    const char* sequence, unsigned int sequence_len,
    unsigned int kmer_len, //u windowu gleda po k slova
    unsigned int window_len //duljina prozora u kojem gleda po k slova
    );
  
}