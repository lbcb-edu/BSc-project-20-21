#include <iostream>
#include <algorithm>


using namespace std;

namespace pink {

enum AlignmentType {global, local, semiglobal};
enum ComingFromDirection {fromNone, fromLeft, fromTop, fromDiagonal};
//global => match:0, mismatch:+1, gap:+1, alignment score: vrijednost u polju zadnjeg retka i zadnjeg stupca
//local => match:+2, mismatch:-1, gap:-2, alignment score: maksimalna vrijednost matrice
//semiglobal => match:+1, mismatch:-1, gap:-1, alignment score: najveća vrijednost od zadnjeg retka i zadnjeg stupca

int MATCH = 0, INSERTION = 0, DELETION = 0;

struct Cell{
    int value;
    ComingFromDirection direction;
    bool operator<(const Cell &other) const {
        return value < other.value;
    };

    bool operator>(const Cell &other) const {
        return value > other.value;
    };

    bool operator==(const Cell &other) const {
        return value == other.value;
    };
};

int Align(const char* query, unsigned int query_len,
	const char* target, unsigned int target_len,
	AlignmentType type,
    	int match, int mismatch, int gap,
    	string* cigar = nullptr, unsigned int* target_begin = nullptr){
        
        Cell** V = new Cell * [query_len+1];
        for(unsigned int i = 0; i <= query_len; i++) V[i] = new Cell[target_len+1];

        //initialization
        if(type == global){
            V[0][0] = {0, fromNone}; 
            for(int i = 1; i <= query_len; i++){
                V[i][0].value = gap * i;
                V[i][0].direction = fromTop;
            }
            for(int j = 1; j <= target_len; j++){
                V[0][j].value = gap * j;
                V[0][j].direction = fromLeft;
            }
        }else{
            for(int i = 0; i <= query_len; i++){
                for(int j = 0; j <= target_len; j++){
                    V[i][j].value = 0;
                    V[i][j].direction = fromNone;
                } 
            }
        } 
        for(int i = 1; i <= query_len; i++){
            for(int j = 1; j <= target_len; j++){
                MATCH = V[i-1][j-1].value + ((query[i-1] == target[j-1]) ? match : mismatch);
                INSERTION = V[i][j-1].value + gap;
                DELETION = V[i-1][j].value + gap;
                if(type == local) V[i][j] = max({Cell{0, fromNone}, Cell{MATCH, fromDiagonal}, Cell{INSERTION, fromLeft}, Cell{DELETION, fromTop}});
                else if(type == semiglobal) V[i][j] = max({Cell{MATCH, fromDiagonal}, Cell{INSERTION, fromLeft}, Cell{DELETION, fromTop}});
                else V[i][j] = min({Cell{INSERTION, fromLeft}, Cell{MATCH, fromDiagonal}, Cell{DELETION, fromTop}}); 
            }
        }
        
        //find alignment score
        int alignmentScore = 0;
        unsigned int i_scoreIndex, j_scoreIndex;
        if(type == local){
            for(int i = 1; i <= query_len; i++){
                for(int j = 1; j <= target_len; j++){
                    if(V[i][j].value > alignmentScore){
                        alignmentScore = V[i-1][j-1].value;
                        i_scoreIndex = i - 1;
                        j_scoreIndex = j - 1;
                    }
                }
            }
        }else if(type == semiglobal){
            for(int i = 1; i < query_len; i++){
                if(V[i][target_len].value > alignmentScore){
                    alignmentScore = V[i][target_len].value;
                    i_scoreIndex = i;
                    j_scoreIndex = target_len;
                }
            }
            for(int j = 1; j < target_len; j++){
                if(V[query_len][j].value > alignmentScore){
                    alignmentScore = V[query_len][j].value;
                    i_scoreIndex = query_len;
                    j_scoreIndex = j;
                }
            }
        }else{
            alignmentScore = V[query_len][target_len].value;
            i_scoreIndex = query_len;
            j_scoreIndex = target_len;
        }

        //create cigar string
        string tmp = "";
        
        while(V[i_scoreIndex][j_scoreIndex].direction != fromNone){ 
            switch (V[i_scoreIndex][j_scoreIndex].direction)
            {
            case fromLeft:
                //Insertion
                tmp = tmp + 'I';
                j_scoreIndex--;
                break;
            case fromDiagonal:
                //Match and mismatch
                if(query[i_scoreIndex -1] == target[j_scoreIndex -1]){
                    tmp = tmp + '=';
                }
                else{
                    tmp = tmp + 'X';
                }
                j_scoreIndex--;
                i_scoreIndex--;
                break;
            case fromTop:
                //Deletion
                tmp = tmp + 'D';
                i_scoreIndex--;
                break;
            default:
                break;
            }
        }
        reverse(tmp.begin(), tmp.end());

        string result = "";
        char c = tmp.at(0);
        int n = 1;
        for(int i = 1; i < tmp.size(); i++){
            if(tmp.at(i) == c){
                n++;
            }else{
                result = result + to_string(n) + c;
                n = 1;
                c = tmp.at(i);
            }
        }
        result = result + to_string(n) + c;
        target_begin = &j_scoreIndex;
        cigar = &result;
        
        for (int i = 0; i <= query_len; i++) {
            delete[] V[i];
        }
        delete[] V;
        cout << *cigar << "\n";
        return alignmentScore;
};
} //namespace pink
/*int main(){
    //primjeri iz obrade informacije
    char query[] = {'G', 'T', 'A','C', 'C'};
    char target[] = {'T', 'T', 'C', 'A', 'C', 'G', 'T', 'T', 'A'};
    cout << Align(query, sizeof(query), target, sizeof(target), semiglobal, 1, -1, -1);
    cout << endl;*/


/*     char queryg[] = {'G', 'T', 'A','C', 'C'};
    char targetg[] = {'G', 'A', 'T', 'A', 'C', 'G', 'T', 'T', 'A'};
    cout << Align(queryg, sizeof(queryg), targetg, sizeof(targetg), global, 0, 1, 1);
    cout << endl; */


/*     char queryl[] = {'G', 'A', 'T','C', 'A', 'T', 'A', 'T', 'T'};
    char targetl[] = {'T', 'C', 'G', 'T', 'A', 'G', 'C', 'G'};
    cout << Align(queryl, sizeof(queryl), targetl, sizeof(targetl), local, 2, -1, -2);
    cout << endl; */ 

    /*return 0;
}*/
