
#include <iostream>
#include <math.h>
#include <vector>
#include <limits>

namespace white {
    enum AlignmentType {GLOBAL, LOCAL, SEMIGLOBAL};
    enum Operation {kMatch = 0, kMismatch, kDelete, kInsert, kNone};

    /**
     * @brief model polja u matrici, sadrzi operaciju te dosadasnju vrijednost poravnanja
     * 
     */
    struct Cell {
        Operation operation;
        int value;
    };
    
    class Aligner {

        public:
            const char* query;
            unsigned int query_len;
            const char* target;
            unsigned int target_len;
            int match;
            int mismatch;
            int gap;
            std::string* cigar;
            unsigned int* target_begin;
            std::vector<std::vector<Cell>> Mat;

        /**
         * @brief Konstruira novi Aligner objekt
         * 
         * @param query prva sekvenca
         * @param query_len duljina prve sekvence
         * @param target druga sekvenca
         * @param target_len duljina druge sekvence
         * @param match vrijednost koja se zbraja na ukupnu sumu ukoliko se pojavi match
         * @param mismatch vrijednost koja se zbraja na ukupnu sumu ukoliko se pojavi mismatch
         * @param gap vrijednost koja se zbraja na ukupnu sumu ukoliko se pojavi gap
         * @param cigar string u koji se sprema rezultat poravnanja u CIGAR formatu
         * @param target_begin pokazuje od koje pozicije pocinjemo sa poravnanjem
         */
        Aligner(const char* query, unsigned int query_len, const char* target,
                unsigned int target_len, int match, int mismatch, int gap,
                std::string* cigar, unsigned int* target_begin);

        /**
         * @brief Needleman-Wunsch metoda poravnanja
         * 
         * @return ukupan trosak poravnanja, ovisno o parametrima zadanim u Aligneru
         */
        int NeedlemanWunsch();

        /**
         * @brief Smith-Waterman metoda poravnanja
         * 
         * @return ukupan trosak poravnanja, ovisno o parametrima zadanim u Aligneru
         */
        int SmithWaterman();

        /**
         * @brief Semi-globalna metoda poravnanja
         * 
         * @return ukupan trosak poravnanja, ovisno o parametrima zadanim u Aligneru 
         */
        int SemiGlobal();

        /**
         * @brief funkcija koja gradi string u CIGAR formatu
         * 
         * @param row broj retka celije s maksimalnom vrijednosti 
         * @param col broj stupca celije s maksimalnom vrijednosti
         * @param starting_cigar string nad kojim funkcija nastavlja gradnju
         * @return string u CIGAR formatu
         */
        std::string CigarBuilder(int& row, int& col, std::string starting_cigar);

        /**
         * @brief 
         * 
         * @param row 
         * @param col 
         */
        void ClippedCigarBuilder(int row, int col);

        /**
         * @brief funkcija koja odabire metodu poravnanja
         * 
         * @param type tip poravnanja koji koristimo
         * @return ukupan trosak poravnanja
         */
        int Align(AlignmentType type);
    };   
}
