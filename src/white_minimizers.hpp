
#include <iostream>
#include <vector>
#include <string>
#include <tuple>

namespace white {
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